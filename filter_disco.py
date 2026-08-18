#!/usr/bin/env python3
"""Filter BLENDER/Disco tab-separated results by score and mismatches."""

import argparse
import csv
import sys
from contextlib import nullcontext
from pathlib import Path
from typing import TextIO


REQUIRED_COLUMNS = (
    "Mismatches",
    "Discoscore"
)

DEFAULT_THRESHOLDS = ((7, 4), (5, 3), (3, 2))
STRICT_THRESHOLDS = ((7, 5), (5, 4), (3, 3))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter Disco results using score and mismatch thresholds."
    )
    parser.add_argument(
        "input",
        type=Path,
        help="tab-separated Disco results file, or - to read from standard input",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="output file (default: standard output)",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="use the stricter thresholds (originally from filter_pool.pl)",
    )
    return parser.parse_args()


def keep_row(discoscore: float, mismatches: int, strict: bool = False) -> bool:
    thresholds = STRICT_THRESHOLDS if strict else DEFAULT_THRESHOLDS
    return any(
        mismatches <= max_mismatches and discoscore >= min_score
        for max_mismatches, min_score in thresholds
    )


def filter_rows(input_file: TextIO, output_file: TextIO, strict: bool) -> None:
    reader = csv.DictReader(input_file, delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError("input file is empty or has no header")

    missing = [column for column in REQUIRED_COLUMNS if column not in reader.fieldnames]
    if missing:
        raise ValueError(f"input is missing required columns: {', '.join(missing)}")

    writer = csv.DictWriter(
        output_file,
        fieldnames=reader.fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()

    for line_number, row in enumerate(reader, start=2):
        try:
            discoscore = float(row["Discoscore"])
            mismatches = int(row["Mismatches"])
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"line {line_number} has an invalid Discoscore or Mismatches value"
            ) from error

        if keep_row(discoscore, mismatches, strict):
            writer.writerow(row)


def main() -> int:
    args = parse_args()
    try:
        if (
            str(args.input) != "-"
            and args.output is not None
            and args.input.resolve() == args.output.resolve()
        ):
            raise ValueError("input and output must be different files")

        input_context = (
            nullcontext(sys.stdin)
            if str(args.input) == "-"
            else args.input.open("r", encoding="utf-8", newline="")
        )
        output_context = (
            nullcontext(sys.stdout)
            if args.output is None
            else args.output.open("w", encoding="utf-8", newline="")
        )
        with input_context as input_file, output_context as output_file:
            filter_rows(input_file, output_file, args.strict)
    except (OSError, ValueError) as error:
        print(f"filter_disco.py: error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
