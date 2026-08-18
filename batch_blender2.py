#!/usr/bin/env python3
"""Run blender2 over the BAM/guide pairs in a CSV manifest."""

from __future__ import annotations

import argparse
import csv
import multiprocessing as mp
import os
from dataclasses import dataclass
from pathlib import Path
import subprocess
import sys
from typing import Sequence


@dataclass(frozen=True)
class Job:
    bam: Path
    guide: str
    output: Path
    line_number: int


@dataclass(frozen=True)
class JobResult:
    job: Job
    returncode: int


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Batch run blender2 for BAM/guide pairs in a CSV manifest"
    )
    batch = parser.add_argument_group("batch options")
    batch.add_argument("--manifest", required=True, type=Path,
                       help="CSV with BAM and guide columns (an optional third site column is ignored)")
    batch.add_argument("--blender", required=True, type=Path,
                       help="path to the executable blender2.py")
    batch.add_argument("-j", "--ncpu", dest="ncpu", type=positive_int, default=8,
                       help="number of concurrent blender processes (default: 8)")
    batch.add_argument("--output-dir", type=Path, default=Path.cwd(),
                       help="existing output directory (default: current directory)")
    batch.add_argument("--filter-mode", choices=("strict", "weak", "none"),
                       default="strict", help="filtering applied to every run (default: strict)")

    blender = parser.add_argument_group("blender2 options")
    blender.add_argument("-n", "--nuclease", required=True,
                         help="nuclease to search for")
    blender.add_argument("-r", "--reference-fasta", required=True, type=Path,
                         help="indexed reference genome in FASTA format")
    blender.add_argument("-c", "--control-bam", type=Path,
                         help="control BAM used for every sample")
    blender.add_argument("-w", "--window-size", type=int,
                         help="override the nuclease-default score window size")
    blender.add_argument("-p", "--pams", nargs="+",
                         help="override default PAMs with a space-separated list")
    blender.add_argument("-t", "--threshold", type=int, default=3,
                         help="read threshold (default: 3)")
    blender.add_argument("-m", "--max-mismatches", type=int, default=8,
                         help="maximum guide mismatches (default: 8)")
    blender.add_argument("-s", "--score-min", type=int, default=3,
                         help="minimum hit score (default: 3)")
    blender.add_argument("-b", "--blacklist", type=Path,
                         help="BED blacklist used for every sample. VERY STRONGLY RECOMMENDED")
    blender.add_argument("-q", "--mapq", type=int, default=10,
                         help="minimum read MAPQ (default: 10)")
    blender.add_argument("--verbose", action="store_true")
    blender.add_argument("--debug", action="store_true")
    return parser


def _require_file(parser: argparse.ArgumentParser, path: Path, label: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        parser.error(f"{label} is not a file: {path}")
    return path


def validate_args(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    args.manifest = _require_file(parser, args.manifest, "manifest")
    args.blender = _require_file(parser, args.blender, "blender executable")
    if not os.access(args.blender, os.X_OK):
        parser.error(f"blender is not executable: {args.blender}")
    args.reference_fasta = _require_file(parser, args.reference_fasta, "reference FASTA")
    if args.control_bam is not None:
        args.control_bam = _require_file(parser, args.control_bam, "control BAM")
    if args.blacklist is not None:
        args.blacklist = _require_file(parser, args.blacklist, "blacklist")
    args.output_dir = args.output_dir.expanduser().resolve()
    if not args.output_dir.is_dir():
        parser.error(f"output directory is not a directory: {args.output_dir}")


def read_manifest(manifest: Path, output_dir: Path) -> list[Job]:
    jobs: list[Job] = []
    errors: list[str] = []
    with manifest.open(newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle), start=1):
            if not row or all(not value.strip() for value in row):
                continue
            if len(row) not in (2, 3):
                errors.append(f"line {line_number}: expected 2 or 3 columns, found {len(row)}")
                continue
            bam_text, guide = (value.strip() for value in row[:2])
            if not bam_text or not guide:
                errors.append(f"line {line_number}: BAM and guide must both be non-empty")
                continue
            bam = Path(bam_text).expanduser()
            if not bam.is_absolute():
                bam = manifest.parent / bam
            bam = bam.resolve()
            if not bam.is_file():
                errors.append(f"line {line_number}: BAM is not a file: {bam}")
                continue
            guide = guide.upper()
            output = output_dir / f"{guide}_{bam.stem}_blender.tsv"
            jobs.append(Job(bam=bam, guide=guide, output=output, line_number=line_number))

    if not jobs and not errors:
        errors.append("manifest contains no jobs")
    outputs: dict[Path, int] = {}
    for job in jobs:
        if job.output in outputs:
            errors.append(
                f"lines {outputs[job.output]} and {job.line_number}: duplicate output path {job.output}"
            )
        else:
            outputs[job.output] = job.line_number
    if errors:
        raise ValueError("Invalid manifest:\n  " + "\n  ".join(errors))
    return jobs


def common_blender_args(args: argparse.Namespace) -> list[str]:
    command = [str(args.blender), "-n", args.nuclease,
               "-r", str(args.reference_fasta)]
    if args.control_bam is not None:
        command.extend(("-c", str(args.control_bam)))
    if args.window_size is not None:
        command.extend(("-w", str(args.window_size)))
    if args.pams:
        command.append("-p")
        command.extend(pam.upper() for pam in args.pams)
    command.extend(("-t", str(args.threshold), "-m", str(args.max_mismatches),
                    "-s", str(args.score_min), "-q", str(args.mapq)))
    if args.blacklist is not None:
        command.extend(("-b", str(args.blacklist)))
    if args.filter_mode == "strict":
        command.append("--filter_strict")
    elif args.filter_mode == "weak":
        command.append("--filter")
    if args.verbose:
        command.append("--verbose")
    if args.debug:
        command.append("--debug")
    return command


def run_job(job: Job, common_args: Sequence[str]) -> JobResult:
    command = list(common_args)
    command.extend(("-f", str(job.bam), "-g", job.guide, "-o", str(job.output)))
    try:
        result = subprocess.run(command, shell=False, check=False)
        returncode = result.returncode
    except OSError as error:
        print(f"Could not start blender for manifest line {job.line_number}: {error}",
              file=sys.stderr)
        returncode = 127
    return JobResult(job=job, returncode=returncode)


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(parser, args)
    try:
        jobs = read_manifest(args.manifest, args.output_dir)
    except (OSError, csv.Error, ValueError) as error:
        parser.error(str(error))

    shared_args = common_blender_args(args)
    with mp.Pool(processes=args.ncpu) as pool:
        pending = [pool.apply_async(run_job, (job, shared_args)) for job in jobs]
        results = [result.get() for result in pending]

    failures = [result for result in results if result.returncode != 0]
    if failures:
        print(f"{len(failures)} of {len(results)} blender jobs failed:", file=sys.stderr)
        for result in failures:
            print(f"  line {result.job.line_number}: {result.job.bam} "
                  f"(exit {result.returncode})", file=sys.stderr)
        return 1
    print(f"Completed {len(results)} blender jobs successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
