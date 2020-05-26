#!/usr/bin/python2

import argparse
from collections import defaultdict
import gzip
import sys


def parse_exons(reg2transcript_file):
    """ Parse through the exons of nirvana 2.0.3 
    
    Returns dict of dict of dict:
    {
        refseq: {
            "position": {
                chrom: [
                    (start1, end1, exon_nb1),
                    (start2, end2, exon_nb2)
                ]
            }
        }
    }
    """

    exons = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    with open(reg2transcript_file) as f:
        for line in f:
            chrom, start, end, gene_symbol, refseq, exon_nb = line.strip().split()
            exons[refseq]["position"][chrom].append((int(start)+1, int(end), int(exon_nb)))

    return exons


def parse_coverage_file(coverage_file):
    """ Parse coverage file

    Returns dict:
    {
        chrom: [
            (start, end),
            (start, end),
        ]
    }
    """

    cov_data = {}

    with gzip.open(coverage_file, "rb") as f:
        for line in f:
            line = line.decode().strip("\n")

            if not line.startswith("#"):
                (chrom, start, end) = line.split("\t")[0:3]
                cov_data.setdefault(chrom, []).append((int(start), int(end)))

    return cov_data        


def main(transcripts, cov_file, exon_file):
    """ Print the exons of given transcripts in a perl hash format

    Returns exons_covered
    """

    exons_covered = []

    exons_data = parse_exons(exon_file)
    coverage_data = parse_coverage_file(cov_file)

    # loop through transcripts
    for transcript in transcripts:
        # check if the transcript is in the exon file
        if transcript in exons_data:
            # loop through the chrom and its exons
            for chrom, exons in exons_data[transcript]["position"].items():
                # loop through the exons
                for exon in exons:
                    exon_start, exon_end, exon_nb = exon
                    
                    # loop through the regions of that chromosome
                    for region in coverage_data[chrom]:
                        region_start, region_end = region

                        # if the exon == the region
                        if (exon_start == region_start and exon_end == region_end):
                            data = {
                                "refseq": transcript,
                                "region": {
                                   "chrom": chrom,
                                   "start": exon_start,
                                   "end": exon_end
                                },
                                "exon_nr": exon_nb
                            }

                            if data not in exons_covered:
                                exons_covered.append(data)
                    
    print("(")
    for exon in exons_covered:
        print("\t{")

        for field, val in exon.items():
            if field == "region":
                print("\tregion=>{")

                for pos_field, pos_data in val.items():
                    print("\t\t\t{}=>\"{}\",".format(pos_field, pos_data))

                print("\t\t}")
            else:
                print("\t{}=>\"{}\",".format(field, val))

        print("\t},")

    print(")")

    return exons_covered



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Get the region coverage data from ref bed and transcripts")

    parser.add_argument("-c", "--coverage", help="Coverage file from region_coverage.py")
    parser.add_argument("-t", "--transcripts", nargs="+", help="Transcript(s) to define regions")
    parser.add_argument("-e", "--exon_file", help="Dump of nirvana GFF with all exons")

    args = parser.parse_args()

    transcripts = args.transcripts
    coverage_file = args.coverage
    exon_file = args.exon_file

    main(transcripts, coverage_file, exon_file)
