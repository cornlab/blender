#!/usr/bin/env python
import pysam
import re
import argparse
import warnings

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
    pam_seqs: tuple # all possible pams (should all be the same length)
    _pam_len: int
    pam_loc: int # relative to the 3' end of the protospacer
    cut_sites: tuple # relative to the 3' end of the protospacer:  (top,bottom)
    _cut_separation: int
    
    def __str__(self):
        return f"name={self.name}, protospacer={self.protospacer_seq}, pams={self.pam_seqs}, pam_location={self.pam_loc}, cut_sites={self.cut_sites}"
    
    def __init__(self, name: str, protospacer_seq: str, pam_seqs: tuple, pam_loc: int, cut_sites: tuple):
        self.name = name
        self.protospacer_seq = protospacer_seq
        self._protospacer_len: len(self.protospacer_seq)
        self.pam_seqs = pam_seqs
        self.pam_loc = pam_loc
        self._pam_len = max([len(x) for x in self.pam_seqs])
        self.cut_sites = cut_sites
        self._cut_separation = abs(self.cut_sites[0] - self.cut_sites[1])
        
        if "N" in self.pam_seqs:
            m = re.match(r'(N*)([ABCDGHMNRSTUVWXYabcdghmnrstuvwxy]+)(N*)')
            if m.group(1):
                assert(self.pam_loc < 0)
                self.pam_loc = -1*len(m.group(1))
            if m.group(2):
                self.pam_seqs = m.group(2)
                self._pam_len = len(m.group(2))
            if m.group(3):
                assert(self.pam_loc > 0)
                self.pam_loc = len(m.group(3))

    def __init__(self, name: str, protospacer_seq: str):
        self.name = name
        self.protospacer_seq = protospacer_seq
        self._protospacer_len = len(self.protospacer_seq)
        
        if "Cas9" in self.name:
            if "Spy" in self.name:
                self.pam_seqs = ("AG", "GG")
                self.pam_loc = 1
                self._pam_len = 2
                self.cut_sites = (-3,-3)
            elif "Sa" in self.name:
                self.pam_seqs = ("GGG", "GAG", "GGA", "GAA")
                self.pam_loc = 2
                self._pam_len = 3
                self.cut_sites = (-3,-3)
            else: 
                raise TypeError(f"Nuclease {self.name} not recognized!")
        elif "Cas12a" in self.name:
            if "As" in self.name:
                self.pam_seqs = ("TTT", "TCC", "CTT", "TCT", "TTC")
                self.pam_loc = (max([len(x) for x in self.pam_seqs]) + 1 + len(self.protospacer_seq))*-1
                self.cut_sites = (-5,-1)
            elif "Lb" in self.name:
                self.pam_seqs = ("TTT", "TCC", "CTT", "TCT", "TTC")
                self.pam_loc = (max([len(x) for x in self.pam_seqs]) + 1 + len(self.protospacer_seq))*-1
                self.cut_sites = (-5,-1)
            else: 
                raise TypeError(f"Nuclease {self.name} not recognized!")
        elif "Cas12f" in self.name:
            pass
        elif "TnpB" in self.name:
            pass
        else:
            raise TypeError("Nuclease not recognized!")
        self._cut_separation = abs(self.cut_sites[0] - self.cut_sites[1])
        self._pam_len = max([len(x) for x in self.pam_seqs])

def parse_arguments():
    parser = argparse.ArgumentParser(prog="autoBLENDER", description = "find Cas on- and off-targets using AutoDisco")
    parser.add_argument('-n', '--nuclease', type=str, default=None, required=True, help='Nuclease to search for. Possibilities are: SpyCas9, SaCas9, LbCas12a, AsCas12a, Cas12f, TnpB')
    parser.add_argument('-w', '--window_size', type=int, help='Override nuclease-default window size for score summing. Cas9 default = 5, Cas12 default = 10')
    parser.add_argument('-f', '--file', required=True, help='experimental BAM (required)')
    parser.add_argument('-c', '--control', help='control BAM (optional, but highly recommended)')
    parser.add_argument('-g', '--guide', required=True, help="On-target guide RNA sequence. Required. provided 5'-3' without the PAM sequence")
    parser.add_argument('-t', '--threshold', type=int, default=3, help='Number of reads to consider for threshold (dfault 3)')
    parser.add_argument('-p', '--pams', metavar='PAM', nargs='*', help='Override nuclease-default PAMs with a space-separated list e.g. "-p GG AG"')
    parser.add_argument('-r', '--reference', required=True, help='Indexed (faidx) reference genome (fasta format). Index should be called <reference>.fai')
    parser.add_argument('-m', '--max_mismatches', type=int, default=8, help="Maximum number of mismatches to allow to the guide sequence (default 8)")
    parser.add_argument('-s', '--score_min', type=int, default=3, help="Minimum score to consider a hit (default 3)")
    parser.add_argument('-b', '--blacklist', help='Blacklist to use for filtering hits, e.g. from ENCODE (BED3 format)')
    parser.add_argument('--verbose', action='store_true', default=False, help="verbose output")
    parser.add_argument('--debug', action='store_true', default=False, help='debug output')
    args = parser.parse_args()
    return args

def check_read(read, min_MQ: int = 25):
    if read.is_unmapped:
        return False
    if read.mapping_quality <= min_MQ:
        return False
    return True

def read_in_blacklist( read, blacklist: dict ): # backlist format = {str(chr): (start,end)}
    try:
        b_locations = blacklist[read.reference_name]
    except KeyError:
        return False
    else:
        for location in b_locations:
            (b_start, b_end) = location
            if (read.reference_start+1 >= b_start) and (read.reference_start+1 <= b_end): # add 1 to match python 0-indexing to bedfile 1-indexing
                return True
        return False

def location_in_blacklist( chromosome: str, start: int, blacklist: dict ): # backlist format = {str(chr): (start,end)}
    try:
        b_locations = blacklist[chromosome]
    except KeyError:
        return False
    else:
        for location in b_locations:
            (b_start, b_end) = location
            if (start+1 >= b_start) and (start+1 <= b_end): # add 1 to match python 0-indexing to bedfile 1-indexing
                return True
        return False

def combine_starts(nuclease: Nuclease, for_starts: dict, rev_starts: dict, threshold: int):
    both_starts = {}
    for start in sorted(for_starts.keys()):
        if start - nuclease._cut_separation - 1 < 0:
            continue
        both = for_starts.get(start, 0) + rev_starts.get(start - nuclease._cut_separation - 1, 0)
        if both >= threshold:
            both_starts[start] = both
        if both > for_starts[start]: # blunt ends exist on the reverse strand
            both_starts[start - nuc._cut_separation - 1] = both
    for start in sorted(rev_starts.keys()):
        if start not in both_starts.keys():
            both = for_starts.get(start + nuclease._cut_separation + 1, 0) + rev_starts.get(start,0)
            if both >= threshold:
                both_starts[start] = both
            if both > rev_starts[start]: # blunt ends exist on the forward strand
                both_starts[start + nuc._cut_separation + 1] = both
    return both_starts

def get_pam(nuclease: Nuclease, chromosome: str, location: int, strand: str, fastaref):
    ref_pam = ""
    s = None # relative to the cutsite
    e = None # relative to the cutsite
    if strand == "minus": #+23 Cas12, -5 Cas9 CCNxxx   543210
        s = location + nuclease.cut_sites[1] - nuclease.pam_loc - nuclease._pam_len + 1
    elif strand == "plus": #-19 Cas12, +5 Cas9
        s = location + -1*nuclease.cut_sites[0] + nuclease.pam_loc + nuclease._pam_len - 1
    e = s+nuclease._pam_len # python closed end notation
    if s < 1:
        return
    ref_pam = fastaref.fetch(reference=chromosome, start=s, end=e).upper()
    if strand == "minus":
        ref_pam = revcomp(ref_pam)
    return ref_pam

def get_sequence(chromosome: str, start: int, end: int, fastaref):
    return fastaref.fetch(reference=chromosome, start=start, end=end+1).upper()

def sum_window(for_starts: dict, rev_starts: dict, start: int, window_size: int):
    x = 0
    for i in range( (start-window_size-1), (start+1+1) ):
        x += rev_starts.get(i, 0)
    for i in range( start, (start + window_size+1+1)):
        x += for_starts.get(i, 0)
    return x

def n_mm(seq1: str, seq2: str): # number of mismatches between two sequences
    if len(seq1) != len(seq2):
        if verbose:
            warnings.warn("Length of sequences are not equal: " + seq1 + " " + str(len(seq1)) + " / " + seq2 + " " + str(len(seq2)), RuntimeWarning)
    mm = sum(c1!=c2 for c1,c2 in zip(seq1,seq2))
    return mm

def revcomp(seq: str):
    rc = seq[::-1]; # reverse slicing
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
        if chromosome not in blacklist:
            blacklist[chromosome] = [(int(start),int(end))]
        else:
            blacklist[chromosome].append((int(start),int(end)))
    if debug:
        print(blacklist)
    f.close()
    return blacklist

if __name__ == '__main__':
    args = parse_arguments()
    edited_fname = args.file
    control_fname = args.control
    input_guide = args.guide
    reference_fname = args.reference
    threshold = args.threshold
    verbose = args.verbose
    debug = args.debug
    max_mismatches = args.max_mismatches
    score_min = args.score_min
    blacklist_fname = args.blacklist
    
    nuc = Nuclease(args.nuclease, args.guide)
    window_size = 5 # default
    if args.pams:
        nuc.pam_seqs = tuple(args.pams)
    
    if args.window_size:
        window_size = args.window_size
    elif "Cas9" in nuc.name:
        window_size = 5
    elif "Cas12a" in nuc.name:
        window_size = 10
    if verbose:
        print(f"nuclease parameters: {nuc}")

    if (verbose):
        print (f"arguments {args}")
    
    blacklist = get_blacklist(blacklist_fname)

    # Note that all output needs to be 1-offset to change from python 0-indexing to 1-indexing!    
    print("Chr:Start-End\tCutsite\tDiscoscore\tCutsite Ends\tStrand\tPAM\tGuide sequence\tMismatches")

    edited_bamfile = pysam.AlignmentFile(edited_fname, "rb")
    for chromosome in edited_bamfile.references:
        if chromosome[:4] == "chrUn" or chromosome[:4] == "chrM":
            continue
        if debug:
            print("Working on " + chromosome)
        for_starts = {}
        rev_starts = {}
        both_starts = {}
        for read in edited_bamfile.fetch(contig=chromosome, multiple_iterators=True):
            goodread = check_read(read)
            if goodread:
                start = read.reference_start
                end = read.reference_end
                if read.template_length == 0: # mate is unmapped
                    if not (read.flag and read.mate_is_unmapped): # this combo only happens if the read is not the 2nd in the pair
                        for_starts[start] = for_starts.get(start, 0) + 1
                    else:
                        rev_starts[end + nuc._cut_separation - 1] = rev_starts.get(end + nuc._cut_separation - 1, 0) + 1
                elif read.template_length > 0: # first in pair
                    for_starts[start] = for_starts.get(start, 0) + 1
                elif read.template_length < 0: # second in pair
                    rev_starts[end-1] = rev_starts.get(end-1, 0) + 1
        both_starts = combine_starts(nuc, for_starts, rev_starts, threshold)
        if debug:
            print(for_starts)
            print(rev_starts)
            print(both_starts)
        if debug:
            print(chromosome + " testing " + str(len(both_starts.keys())) + " starts")

        edited_count = {}
        if control_fname:
            ctrl_count = {}
            control_bamfile = pysam.AlignmentFile(control_fname, "rb")
        for start in both_starts.keys():
            edited_count[start] = edited_bamfile.count(contig=chromosome, start=start, stop=start)
            if control_fname:
                ctrl_count[start] = control_bamfile.count(contig=chromosome, start=start, stop=start)
        if control_fname:
            control_bamfile.close()

        reference_fasta = pysam.Fastafile(filename=reference_fname)
        
        output = {}
        for start in both_starts.keys():
            if blacklist != {}:
                if location_in_blacklist(chromosome, start, blacklist ):
                    if verbose: 
                        print (read.reference_name + ":" + str(read.reference_start) + "-" + str(read.reference_end) + "\t read FILTERED:blacklisted")
                    continue
            if control_fname:
                if ctrl_count[start] > 10:
                    if verbose:
                        print("CONTROL skipping " + chromosome + ":" + str(start+1) + " " + str(ctrl_count[start]))
                    continue

            if edited_count[start] > 0:
                if(both_starts[start]/edited_count[start] < 0.25):
                    if verbose: 
                        print(chromosome + ":x-x\t" +
                                str(start+1) + "\t" +
                                str(score) + "\t" + 
                                str(both_starts[start]) + "\t" +
                                "\tFILTERED: deep area")
                    continue

            score = sum_window(for_starts, rev_starts, start, window_size=window_size)
            if score < score_min: # doesn't pass discover-score cutoff
                if verbose: 
                        print(chromosome + ":" + "x" + "-" + "x" + "\t" +
                            str(start+1) + "\t" +
                            str(score) + "\t" + 
                            str(both_starts[start]) + "\t" + 
                            "\tFILTERED:fails disco score " + str(score_min))
                continue

            start1 = start+1 # convert to 1-based indexing (pysam fetch uses genome reference)
            pam_plus = get_pam(nuc, chromosome, start, "plus", reference_fasta)
            pam_minus = get_pam(nuc, chromosome, start, "minus", reference_fasta)
            if debug:
                print(chromosome, start, both_starts[start], pam_minus, pam_plus, nuc.pam_seqs)

            if pam_minus in nuc.pam_seqs:
                pam_fullseq = ""
                if nuc.pam_loc < 0:
                    pam_fullseq = pam_minus + abs(nuc.pam_loc+nuc._protospacer_len+nuc._pam_len)*'N'
                else:
                    pam_fullseq = 'N'*nuc.pam_loc+pam_minus


                #s : Cas9 start-3,   Cas12 start-1
                s = start + nuc.cut_sites[1] +1
                #e =  Cas9 start + 17
                #e = start + nuc._protospacer_len + nuc.cut_sites[1]
                e = s + nuc._protospacer_len - 1
                
                s1 = s+1 
                e1 = e+1

                if e1 < 1 or s1 < 1: # too close to and end to be our sequence
                    continue
                guide = get_sequence(chromosome, s, e, reference_fasta)
                guide = revcomp(guide)
                mm = n_mm(input_guide, guide)
                if debug:
                    print(f"MINUS {guide} {pam_minus} {input_guide} {mm} {len(guide)}")
                if mm > max_mismatches:
                    if verbose: 
                        print(chromosome + ":" + str(s1) + "-" + str(e1) + "\t" +
                            str(start1) + "\t" +
                            str(score) + "\t" + 
                            str(both_starts[start]) + "\t" + 
                            "antisense\t" +
                            pam_fullseq+"\t" +
                            guide + 
                            " FILTERED: " + str(mm) +  " mismatches")
                else:
                    outstr = chromosome + ":" + str(s1) + "-" + str(e1) + "\t" + str(start1) + "\t" + str(score) + "\t" +  str(both_starts[start]) + "\t" + "antisense\t" + pam_fullseq+"\t" + guide + "\t" + str(mm)
                    output[chromosome+str(start)+guide] = outstr

            if pam_plus in nuc.pam_seqs:
                pam_fullseq = ""
                if nuc.pam_loc < 0:
                    pam_fullseq = pam_plus + abs(nuc.pam_loc+nuc._protospacer_len+nuc._pam_len)*'N'
                else:
                    pam_fullseq = 'N'*nuc.pam_loc+pam_plus
                # s: Cas9 start -17, Cas12 start-19
                s = start - nuc.cut_sites[0] - nuc._protospacer_len + 1
                e = start - nuc.cut_sites[0]
                #e = s + nuc._protospacer_len - 1

                s1 = s+1 
                e1 = e+1

                if e1 < 1 or s1 < 1: # cannot be our sequence
                    continue
                guide = get_sequence(chromosome, s, e, reference_fasta)
                mm = n_mm(input_guide, guide)
                if debug:
                    print(f"PLUS {guide} {pam_plus} {input_guide} {mm} {len(guide)}")
                if mm > max_mismatches:
                    if verbose: 
                        print(chromosome + ":" + str(s1) + "-" + str(e1) + "\t" +
                        str(start1) + "\t" +
                        str(score) + "\t" + 
                        str(both_starts[start]) + "\t" + 
                        "sense\t" +
                        pam_fullseq+"\t" +
                        guide + 
                        " FILTERED: " + str(mm) +  " mismatches")
                else:
                    outstr = chromosome + ":" + str(s1) + "-" + str(e1) + "\t" + str(start1) + "\t" + str(score) + "\t" + str(both_starts[start]) + "\t" + "sense\t" + pam_fullseq+"\t" + guide + "\t" + str(mm)
                    output[chromosome+str(start)+guide] = outstr
        for site in sorted(output.keys()):
            print(output[site])
