#!/usr/bin/env python
import pysam
import argparse
import bisect
import logging
import multiprocessing as mp
import sys
import pandas as pd


log = logging.getLogger('AutoBLENDER')


# yyX X X X X X X|X X X NGGzz
# yyX X X X X|X X X X X NCCzz
# protospacer_seq = XXXXXXXXXX
# _protospacer_len = len(XXXXXXXXXX) = 10
# pam_seqs = (GG, AG)
# _pam_len = 2
# pam_loc: 1
# cut_sites = (-3,-5)
# cut_separation = abs(self.cut_sites[0] - self.cut_sites[1]) = 2


class Nuclease():
    name: str
    protospacer_seq: str
    _protospacer_len: int
    pam_seqs: tuple  # all possible pams (should all be the same length)
    _pam_len: int
    pam_loc: int  # start of the PAM relative to the 3' end of the protospacer
    cut_sites: tuple  # relative to the 3' end of the protospacer:  (top,bottom)
    _cut_separation: int

    def __str__(self):
        return f"name={self.name}, protospacer={self.protospacer_seq}, pams={self.pam_seqs}, pam_location={self.pam_loc}, cut_sites={self.cut_sites}"

    def __init__(self, name: str, protospacer_seq: str):
        self.name = name
        self.protospacer_seq = protospacer_seq
        self._protospacer_len = len(self.protospacer_seq)

        if "Cas9" in self.name:
            if "Spy" in self.name:
                self.pam_seqs = ("AG", "GG")
                self.pam_loc = 1
                self._pam_len = 2
                self.cut_sites = (-3, -3)
            elif "Sa" in self.name:
                self.pam_seqs = ("GGG", "GAG", "GGA", "GAA")
                self.pam_loc = 2
                self._pam_len = 3
                self.cut_sites = (-3, -3)
            else:
                raise TypeError(f"Nuclease {self.name} not recognized!")
        elif "Cas12a" in self.name:
            if "As" in self.name:
                self.pam_seqs = ("TTT", "TCC", "CTT", "TCT", "TTC")
                self.pam_loc = -1*max([len(x) for x in self.pam_seqs]) -1 + -1*len(self.protospacer_seq) # constant PAM + degenerate base + protospacer
                self.cut_sites = (-5, -1)
            elif "Lb" in self.name:
                self.pam_seqs = ("TTT", "TCC", "CTT", "TCT", "TTC")
                self.pam_loc = -1*max([len(x) for x in self.pam_seqs]) -1 + -1*len(self.protospacer_seq) # constant PAM + degenerate base + protospacer
                self.cut_sites = (-6, 2)
            else:
                raise TypeError(f"Nuclease {self.name} not recognized!")
        elif "Cas12f" in self.name:
            if "Un1" in self.name:
                self.pam_seqs = ("TTTG", "TTTA")
                self.pam_loc = -1*max([len(x) for x in self.pam_seqs]) + -1*len(self.protospacer_seq) # constant PAM + degenerate base + protospacer
                self.cut_sites = (-2, 4)
            else:
                raise TypeError(f"Nuclease {self.name} not recognized!")
        elif "TnpB" in self.name:
            if "ISDra2" in self.name:
                # [[TGC][TCA][GT]A[TC]
                self.pam_seqs = (
                    "TTGAT", "TTGAC", "TTTAT", "TTTAC",
                    "TCGAT", "TCGAC", "TCTAT", "TCTAC",
                    "TAGAT", "TAGAC", "TATAT", "TATAC",
                    "GTGAT", "GTGAC", "GTTAT", "GTTAC",
                    "GCGAT", "GCGAC", "GCTAT", "GCTAC",
                    "GAGAT", "GAGAC", "GATAT", "GATAC",
                    "CTGAT", "CTGAC", "CTTAT", "CTTAC",
                    "CCGAT", "CCGAC", "CCTAT", "CCTAC",
                    "CAGAT", "CAGAC", "CATAT", "CATAC"
                )
                self.pam_loc = max([len(x) for x in self.pam_seqs])*-1 + len(self.protospacer_seq)*-1 # constant PAM + protospacer
                self.cut_sites = (-1, 5)
            else:
                raise TypeError(f"Nuclease {self.name} not recognized!")                
        else:
            raise TypeError("Nuclease not recognized!")
        self._cut_separation = abs(self.cut_sites[0] - self.cut_sites[1])
        self._pam_len = max([len(x) for x in self.pam_seqs])

def parse_arguments():
    parser = argparse.ArgumentParser(prog="autoBLENDER", description="find Cas on- and off-targets using AutoDisco")
    parser.add_argument('-n', '--nuclease', type=str, default=None, required=True,
                        help='Nuclease to search for. Possibilities are: SpyCas9, SaCas9, LbCas12a, AsCas12a, Cas12f, TnpB')
    parser.add_argument('-w', '--window_size', type=int,
                        help='Override nuclease-default window size for score summing. Cas9 default = 5, Cas12 default = 10')
    parser.add_argument('-f', '--experiment_bam', required=True, help='experimental BAM (required)')
    parser.add_argument('-o', '--output', required=True, help='output file (required)')
    parser.add_argument('-c', '--control_bam', help='control BAM (optional, but highly recommended)')
    parser.add_argument('--filter', action='store_true', help='Weak filter sites based on score and number of mismatches (default: False)')
    parser.add_argument('--filter_strict', action='store_true', help='Strict filter sites based on score and number of mismatches. Recommended! (default: False)')
    parser.add_argument('-g', '--guide', required=True,
                        help="On-target guide RNA sequence. Required. provided 5'-3' without the PAM sequence")
    parser.add_argument('-t', '--threshold', type=int, default=3,
                        help='Number of reads to consider for threshold (default 3)')
    parser.add_argument('-p', '--pams', metavar='PAM', nargs='*',
                        help='Override nuclease-default PAMs with a space-separated list e.g. "-p GG AG"')
    parser.add_argument('-r', '--reference_fasta', required=True,
                        help='Indexed (faidx) reference genome (fasta format). Index should be called <reference>.fai')
    parser.add_argument('-m', '--max_mismatches', type=int, default=8,
                        help="Maximum number of mismatches to allow to the guide sequence (default 8)")
    parser.add_argument('-s', '--score_min', type=int, default=3, help="Minimum score to consider a hit (default 3)")
    parser.add_argument('-b', '--blacklist', help='Blacklist to use for filtering hits, e.g. from ENCODE (BED3 format). STRONGLY RECOMMENDED')
    parser.add_argument('--verbose', action='store_true', default=False, help="verbose output")
    parser.add_argument('--debug', action='store_true', default=False, help='debug output')
    parser.add_argument('-q', '--mapq', type=int, default=10, help='MAPQ value of reads (default 10)')
    parser.add_argument('-j', '--workers', type=int, default=1,
                        help='Number of chromosomes to process in parallel (default 1)')
    args = parser.parse_args()
    args.pams = [x.upper() for x in args.pams] if args.pams else None  # convert to uppercase
    args.guide = args.guide.upper()  # convert to uppercase
    if args.filter_strict and args.filter:
        log.error(f"Cannot use both --filter and --filter_strict. Exiting.")
        sys.exit(1)
    if args.workers < 1:
        parser.error("--workers must be at least 1")

    return args

def check_read(read, min_MQ = None, blacklist = None):
    if min_MQ == None:
        min_MQ = 10
    if blacklist == None:
        blacklist = {}
    if read.is_unmapped:
        return False
    if read.mapping_quality <= min_MQ:
        return False
    if blacklist != {}:  # if there is a blacklist, check if the read is in it
        if read_in_blacklist(read, blacklist):
            return False
    return True

def _location_in_intervals(location: int, intervals) -> bool:
    """Return whether a 1-based location is in merged, sorted BED intervals."""
    if not intervals:
        return False
    starts, ends = intervals
    index = bisect.bisect_right(starts, location) - 1
    return index >= 0 and location <= ends[index]


def read_in_blacklist(read, blacklist: dict):
    return _location_in_intervals(read.reference_start + 1,
                                  blacklist.get(read.reference_name))

def location_in_blacklist(chromosome: str, start: int, blacklist: dict):
    return _location_in_intervals(start + 1, blacklist.get(chromosome))

def combine_starts(nuclease: Nuclease, for_starts: dict, rev_starts: dict):
    """Combine paired cut ends, evaluating each canonical end pair once."""
    both_starts = {}
    separation = nuclease._cut_separation + 1
    anchors = {start - separation for start in for_starts if start >= separation}
    anchors.update(start for start in rev_starts if start >= separation)
    for anchor in anchors:
        other = anchor + separation
        both = for_starts.get(other, 0) + rev_starts.get(anchor, 0)
        both_starts[anchor] = max(both, both_starts.get(anchor, 0))
        both_starts[other] = max(both, both_starts.get(other, 0))
    return both_starts

def get_pam(nuclease: Nuclease, chromosome: str, location: int, strand: str, fastaref):
    ref_pam = ""
    s = -1  # 0-based chromosomal coordinate relative to the cutsite
    e = -1
    if strand == "minus":  # Cas12 xxxBTTT, Cas9 CCNxxx
        s = location + nuclease.cut_sites[1] - nuclease.pam_loc - nuclease._pam_len
    elif strand == "plus":  # Cas12 TTTVxxx, Cas9 xxxNGG
        s = location + -1*nuclease.cut_sites[0] + nuclease.pam_loc
    e = s + nuclease._pam_len  # python closed end notation
    if s < 1:
        return ""
    ref_pam = fastaref.fetch(reference=chromosome, start=s+1, end=e+1).upper() # convert to 1-based indexing
    if strand == "minus":
        ref_pam = revcomp(ref_pam)
    return ref_pam

def get_sequence(chromosome: str, start: int, end: int, fastaref):
    return fastaref.fetch(reference=chromosome, start=start, end=end + 1).upper()


def sum_window(for_starts: dict, rev_starts: dict, start: int, window_size: int):
    x = 0
    for i in range((start - window_size - 1), (start + 1 + 1)):
        x += rev_starts.get(i, 0)
    for i in range(start, (start + window_size + 1 + 1)):
        x += for_starts.get(i, 0)
    return x

def n_mm(seq1: str, seq2: str):  # number of mismatches between two sequences
    if len(seq1) != len(seq2):
        log.warning("Length of sequences are not equal: " + seq1 + " " + str(len(seq1)) + " / " + seq2 + " " + str(len(seq2)))
    mm = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
    return mm

def revcomp(seq: str):
    rc = seq[::-1];  # reverse slicing
    table = str.maketrans("ABCDGHMNRSTUVWXYabcdghmnrstuvwxy", "TVGHCDKNYSAABWXRtvghcdknysaabwxr")
    rc = rc.translate(table)
    return rc

def get_blacklist(blacklist_fname):
    blacklist = {}
    if blacklist_fname == None:
        return {}
    f = open(blacklist_fname)
    for line in f:
        (chromosome, start, end) = line.rstrip().split()[:3]
        blacklist.setdefault(chromosome, []).append((int(start), int(end)))
    for chromosome, locations in blacklist.items():
        merged = []
        for start, end in sorted(locations):
            if merged and start <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        blacklist[chromosome] = (
            tuple(start for start, unused_end in merged),
            tuple(end for unused_start, end in merged),
        )
    log.debug(f"BLACKLIST: {blacklist}")
    f.close()
    return blacklist


def check_reverse_sequence(nuc: Nuclease, input_guide: str, chromosome: str, start: int, pam_minus: str, max_mismatches: int, reference_fasta: pysam.FastaFile) -> dict: # {'success': True/False, 'e': int, 's': int, 'guide': str, 'pam_fullseq': str, 'mm': int}
    if nuc.pam_loc < 0:
        pam_fullseq = pam_minus + abs(nuc.pam_loc + len(pam_plus) + nuc._protospacer_len ) * 'N'
    else:
        pam_fullseq = 'N' * nuc.pam_loc + pam_minus
    # s : Cas9 start-3,   Cas12 start-1
    s = start + nuc.cut_sites[1] + 1
    e = s + nuc._protospacer_len - 1
    if e < 0 or s < 0:  # too close to an end to be our sequence
        return {'success': False, 'e': e, 's': s, 'guide': False, 'pam_fullseq': pam_fullseq, 'mm': False}
    guide = revcomp(get_sequence(chromosome, s, e, reference_fasta))
    mm = n_mm(input_guide, guide)
    log.debug(f"MINUS {guide} {pam_minus} {input_guide} {mm} {len(guide)}")
    if mm > max_mismatches:
        return {'success': False, 'e': e, 's': s, 'guide': guide, 'pam_fullseq': pam_fullseq, 'mm': mm}
    return {'success': True, 'e': e, 's': s, 'guide': guide, 'pam_fullseq': pam_fullseq, 'mm': mm}

def check_forward_sequence(nuc: Nuclease, input_guide: str, chromosome: str, start: int, pam_plus: str, max_mismatches: int, reference_fasta: pysam.FastaFile) -> dict:
    if nuc.pam_loc < 0:
        pam_fullseq = pam_plus + abs(nuc.pam_loc + len(pam_plus) + nuc._protospacer_len ) * 'N'
    else:
        pam_fullseq = 'N' * nuc.pam_loc + pam_plus
    # s: Cas9 start -17, Cas12 start-19
    s = start - nuc.cut_sites[0] - nuc._protospacer_len + 1
    e = start - nuc.cut_sites[0]
    if e < 0 or s < 0:  # cannot be our sequence
        return {'success': False, 'e': e, 's': s, 'guide': False, 'pam_fullseq': pam_fullseq, 'mm': False}
    guide = get_sequence(chromosome, s, e, reference_fasta)
    mm = n_mm(input_guide, guide)
    log.debug(f"PLUS {guide} {pam_plus} {input_guide} {mm} {len(guide)}")
    if mm > max_mismatches:
        return {'success': False, 'e': e, 's': s, 'guide': guide, 'pam_fullseq': pam_fullseq, 'mm': mm}
    return {'success': True, 'e': e, 's': s, 'guide': guide, 'pam_fullseq': pam_fullseq, 'mm': mm}

OUTPUT_COLUMNS = (
    "Chr:Start-End", "Cutsite", "Discoscore", "Cutsite Ends", "Strand",
    "PAM", "Guide sequence", "Mismatches",
)


def process_chromosome(chromosome, args, nuc, window_size, blacklist,
                       edited_bamfile, reference_fasta, control_bamfile=None):
    """Process one chromosome and return (output rows, good-read count)."""
    nsites = 0
    n_goodreads = 0
    rows = []
    log.info(f"Working on {chromosome}")
    for_starts = {}
    rev_starts = {}
    chromosome_blacklist = blacklist.get(chromosome)

    for read in edited_bamfile.fetch(contig=chromosome):
        if read.is_unmapped or read.mapping_quality <= args.mapq:
            continue
        if chromosome_blacklist and _location_in_intervals(
                read.reference_start + 1, chromosome_blacklist):
            continue
        n_goodreads += 1
        start = read.reference_start
        end = read.reference_end
        if read.mate_is_unmapped:
            if not read.is_reverse:
                for_starts[start] = for_starts.get(start, 0) + 1
            else:
                rev_starts[end - 1] = rev_starts.get(end - 1, 0) + 1
        elif read.template_length > 0:
            for_starts[start] = for_starts.get(start, 0) + 1
        elif read.template_length < 0:
            rev_starts[end - 1] = rev_starts.get(end - 1, 0) + 1

    both_starts = combine_starts(nuc, for_starts, rev_starts)
    if args.debug:
        for key, value in for_starts.items():
            log.debug(f"{chromosome} FORWARD START {key} {value}")
        for key, value in rev_starts.items():
            log.debug(f"{chromosome} REVERSE START {key} {value}")
        for key, value in both_starts.items():
            log.debug(f"{chromosome} BOTH START {key} {value}")
        log.debug(f"{chromosome} testing {len(both_starts)} starts")

    for site, n_ends in both_starts.items():
            if blacklist != {}:
                if location_in_blacklist(chromosome, site, blacklist):
                    if args.verbose:
                        log.info(f"{chromosome}:{site}-{site}\tFILTERED:blacklisted")
                    continue

            score = sum_window(for_starts, rev_starts, site, window_size=window_size)
            if n_ends < args.threshold:
                continue

            # CALCULATE DISCO SCORE
            # filter if the control bam has many reads at this site
            if control_bamfile is not None:
                control_count = control_bamfile.count(
                    contig=chromosome,
                    start=site - nuc._cut_separation - 1,
                    stop=site + nuc._cut_separation + 1,
                    read_callback=check_read,
                )
                if control_count > 10:
                    if args.verbose:
                        log.info(f"CONTROL skipping {chromosome}:\t{site + 1}\t{control_count}")
                    continue

            if score < args.score_min:  # doesn't pass discover-score cutoff
                if args.verbose:
                    log.info(f"{chromosome}:x-x\t{site + 1}\t{site}\t{n_ends}\tFILTERED:fails disco score {args.score_min}")
                continue

            # START TESTS AGAINST REFERENCE SEQUENCE
            pam_plus = get_pam(nuc, chromosome, site, "plus", reference_fasta)
            pam_minus = get_pam(nuc, chromosome, site, "minus", reference_fasta)
            log.debug(f"{chromosome}, {site}, {n_ends}, {pam_minus}, {pam_plus}, {nuc.pam_seqs}")

            # MINUS STRAND SITES
            strand = "" # should always get set. leave empty for debugging
            if pam_minus in nuc.pam_seqs:
                strand = "antisense"
                result = check_reverse_sequence(nuc, args.guide, chromosome, site, pam_minus, args.max_mismatches, reference_fasta)
                site1 = site + 1  # convert to 1-based indexing for output (pysam uses 0-based indexing)
                s1 = result['s'] + 1 # convert to 1-based indexing for output
                e1 = result['e'] + 1 # convert to 1-based indexing for output
                pam_fullseq = result['pam_fullseq']
                guide = result['guide']
                mm = result['mm']

                if result['success']:
                    nsites += 1 # we found a site
                    sitedict = {"Chr:Start-End":f"{chromosome}:{s1}-{e1}", "Cutsite":site1, "Cutsite Ends":n_ends, "Strand":strand, "PAM":pam_fullseq, "Guide sequence":guide, "Mismatches":mm, "Discoscore":score} 
                    rows.append(sitedict)
                    if args.verbose:
                        log.info(sitedict)
                else:
                    if args.verbose:
                        if result['mm'] and result['mm'] >= args.max_mismatches:
                            log.info(f"{chromosome}:{s1}-{e1}\t{site1}\t{score}\t{n_ends}\t{strand}\t{pam_fullseq}\t{guide} FILTERED: {mm} mismatches")
                        if result['s'] < 0 or result['e'] < 0:
                            log.info(f"{chromosome}:{s1}-{e1}\t{site1}\t{score}\t{n_ends}\t{strand}\t{pam_fullseq}\t{guide} FILTERED: too close to start of chromosome")

            # PLUS STRAND SITES
            if pam_plus in nuc.pam_seqs:
                strand = "sense"
                result = check_forward_sequence(nuc, args.guide, chromosome, site, pam_plus, args.max_mismatches, reference_fasta)
                site1 = site + 1  # convert to 1-based indexing for output (pysam uses 0-based indexing)
                s1 = result['s'] + 1 # convert to 1-based indexing for output
                e1 = result['e'] + 1 # convert to 1-based indexing for output
                pam_fullseq = result['pam_fullseq']
                guide = result['guide']
                mm = result['mm']

                if result['success']:
                    nsites += 1 # we found a site
                    sitedict = {"Chr:Start-End":f"{chromosome}:{s1}-{e1}", "Cutsite":site1, "Cutsite Ends":n_ends, "Strand":strand, "PAM":pam_fullseq, "Guide sequence":guide, "Mismatches":mm, "Discoscore":score} 
                    rows.append(sitedict)
                    if args.verbose:
                        log.info(sitedict)
                else:
                    if args.verbose:
                        if result['mm'] and result['mm'] >= args.max_mismatches:
                            log.info(f"{chromosome}:{s1}-{e1}\t{site1}\t{score}\t{n_ends}\t{strand}\t{pam_fullseq}\t{guide} FILTERED: {mm} mismatches")
                        if result['s'] < 0 or result['e'] < 0:
                            log.info(f"{chromosome}:{s1}-{e1}\t{site1}\t{score}\t{n_ends}\t{strand}\t{pam_fullseq}\t{guide} FILTERED: too close to start of chromosome")
    log.info(f"Found {nsites} unfiltered sites on {chromosome}")
    return rows, n_goodreads


def _process_chromosome_worker(payload):
    chromosome, args, nuc, window_size, blacklist = payload
    if args.debug:
        log.setLevel(logging.DEBUG)
    with pysam.AlignmentFile(args.experiment_bam, "rb") as edited_bamfile, \
            pysam.FastaFile(filename=args.reference_fasta) as reference_fasta:
        if args.control_bam:
            with pysam.AlignmentFile(args.control_bam, "rb") as control_bamfile:
                return process_chromosome(
                    chromosome, args, nuc, window_size, blacklist,
                    edited_bamfile, reference_fasta, control_bamfile,
                )
        return process_chromosome(
            chromosome, args, nuc, window_size, blacklist,
            edited_bamfile, reference_fasta,
        )


def main():
    logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                        format='%(asctime)s %(levelname)s : %(message)s')
    log.info(f"Running AutoBLENDER \'{' '.join(sys.argv)}\'")
    args = parse_arguments()
    if args.debug:
        log.setLevel(logging.DEBUG)

    nuc = Nuclease(args.nuclease, args.guide)
    if args.pams:
        nuc.pam_seqs = tuple(args.pams)
    if args.window_size:
        window_size = args.window_size
    elif "Cas12a" in nuc.name or "Cas12f" in nuc.name or "TnpB" in nuc.name:
        window_size = 10
    else:
        window_size = 5
    if args.verbose:
        log.info(f"nuclease parameters: {nuc}")

    blacklist = get_blacklist(args.blacklist)
    with pysam.AlignmentFile(args.experiment_bam, "rb") as bam_header, \
            pysam.FastaFile(filename=args.reference_fasta) as fasta_header:
        chromosomes = []
        fasta_references = set(fasta_header.references)
        for chromosome in bam_header.references:
            if chromosome in fasta_references:
                chromosomes.append(chromosome)
            else:
                log.warning(
                    f"{chromosome} not in reference fasta {args.reference_fasta}. Skipping."
                )

    if args.workers == 1:
        with pysam.AlignmentFile(args.experiment_bam, "rb") as edited_bamfile, \
                pysam.FastaFile(filename=args.reference_fasta) as reference_fasta:
            if args.control_bam:
                with pysam.AlignmentFile(args.control_bam, "rb") as control_bamfile:
                    results = [process_chromosome(
                        chromosome, args, nuc, window_size, blacklist,
                        edited_bamfile, reference_fasta, control_bamfile,
                    ) for chromosome in chromosomes]
            else:
                results = [process_chromosome(
                    chromosome, args, nuc, window_size, blacklist,
                    edited_bamfile, reference_fasta,
                ) for chromosome in chromosomes]
    else:
        payloads = [
            (chromosome, args, nuc, window_size, blacklist)
            for chromosome in chromosomes
        ]
        with mp.Pool(processes=args.workers) as pool:
            results = pool.map(_process_chromosome_worker, payloads)

    rows = [row for chromosome_rows, unused_count in results for row in chromosome_rows]
    n_goodreads = sum(count for unused_rows, count in results)
    df = pd.DataFrame.from_records(rows, columns=OUTPUT_COLUMNS)

    if args.filter:
        log.info(f"Normal filtering sites based on score and number of mismatches")
        filter1 = df[(df['Mismatches'] <= 7) & (df['Discoscore'] >= 4)]
        filter2 = df[(df['Mismatches'] <= 5) & (df['Discoscore'] >= 3)]
        filter3 = df[(df['Mismatches'] <= 3) & (df['Discoscore'] >= 2)]
        df = pd.merge(filter1, filter2, how = "outer")
        df = pd.merge(df, filter3, how = "outer")
    elif args.filter_strict:
        log.info(f"Strict filtering sites based on score and number of mismatches")
        filter1 = df[(df['Mismatches'] <= 7) & (df['Discoscore'] >= 5)]
        filter2 = df[(df['Mismatches'] <= 5) & (df['Discoscore'] >= 4)]
        filter3 = df[(df['Mismatches'] <= 3) & (df['Discoscore'] >= 3)]
        df = pd.merge(filter1, filter2, how = "outer")
        df = pd.merge(df, filter3, how = "outer")

    # NORMALIZE TO FRACTION OF MAX
    df['DPM'] = df['Discoscore'] / (n_goodreads / 1e6)
    df['normDPM'] = df['DPM'] / df['DPM'].max()
    
    # SORT AND OUTPUT
    df.sort_values(by=['normDPM'], ascending=False, inplace=True)
    df.to_csv(args.output, index=False, sep='\t', float_format='%.6f')
    log.info(f"{len(df.index)} site(s) written to {args.output}")


if __name__ == '__main__':
    main()
