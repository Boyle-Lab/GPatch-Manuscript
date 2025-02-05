#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser
from Bio import SeqIO
from Bio.Seq import Seq

__version__ = "0.0.1"

"""
Filter out unpatched contig records from ragtag.patch.fasta
and replace scaffold names with chromosome names, based on
data in the agp files.
"""

def parse_agp(agp_f, c=0):
    """
    Parse an AGP file into a dict keyed on the given columm, c.
    """
    ret = {}
    with open(agp_f, "r") as agp:
        for line in agp:
            if re.match("#", line):
                # Comment/header line. 
                continue
            rec = line.strip().split()
            ret[rec[c]] = rec
    return ret


def parse_patched_agp(agp_f):
    """
    Parse ragtag.patch.agp into a dict, keyed on scaffold name
    where each key contains a list of records on the given scaffold.
    """
    ret = {}
    with open(agp_f, "r") as agp:
        for line in agp:
            rec = line.strip().split()
            if not re.match("scf", rec[0]):
                # Skip all "seq" records in the patched agp as these are
                # not used in the patched sequence records.
                continue
            elif rec[0] in ret:
                ret[rec[0]].append(rec)
            else:
                ret[rec[0]] = [rec]
    return ret


def get_chrom_name(reference_agp, scf_agp, ctg_agp):
    """
    Given the set of records in scf_agp, find the first
    patch record and retrieve the correspinding chromosome
    name from the reference assembly agp.
    """
    qseq = {}
    for rec in scf_agp:
        if re.match("qseq", rec[5]):
            # We have a patch record!
            #return reference_agp[rec[5]][5]
            if reference_agp[rec[5]][5] not in qseq:
                qseq[reference_agp[rec[5]][5]] = reference_agp[rec[5]][5]
    if len(qseq.keys()) == 1:
        # Single match
        return list(qseq.keys())[0]
    elif len(qseq.keys()) > 1:
        sys.stderr.write("Warning: %s matches multiple reference chromosomes: %s. Using the first matching record.\n" % (scf_agp[0][0], qseq))
        return list(qseq.keys())[0]
    # If we've reached this point, there are no patch records. We
    # will make our best guess on the correct chromosome based on
    # the source chrom for the first contig record.
    #sys.stderr.write("%s\t%s\n" % (ctg_agp[scf_agp[0][5]][0], scf_agp[0][0]))
    for rec in scf_agp:
        ctg = ctg_agp[rec[5]]
        chrom = ctg[0].split("_")[0]
        if chrom not in qseq:
            qseq[chrom] = chrom
    if len(qseq.keys()) == 1:
        # Single match
        return list(qseq.keys())[0]
    elif len(qseq.keys()) > 1:
        sys.stderr.write("Warning: %s matches multiple reference chromosomes: %s. Using the first matching record.\n" % (scf_agp[0][0], qseq))
        return list(qseq.keys())[0]
    else:
        return None


def main():
    parser = ArgumentParser(description='Filter out unpatched contig records from ragtag.patch.fasta and replace scaffold names with chromosome names, based on data in the agp files.')
    parser.add_argument('-f', '--fasta', metavar='FASTA', type=str,
                        required=True, help='Path to patched genome assembly fasta from RagTag patch.')
    parser.add_argument('-r', '--reference_agp', metavar='AGP', type=str,
                        required=True, help='AGP file relating original to renamed chromosome names in the reference genome (usually ragtag.patch.rename.agp).')
    parser.add_argument('-p', '--patched_agp', metavar='AGP', type=str,
                        required=True, help='AGP file containing final order and orientation of contigs and patches in the patched genome. (usually ragtag.patch.agp).')
    parser.add_argument('-c', '--contigs_agp', metavar='AGP', type=str,
                        required=True, help='AGP file relating original to renamed contig names (usually ragtag.patch.ctg.agp).')
    
    args = parser.parse_args()

    reference_agp = parse_agp(args.reference_agp, c=0)
    patched_agp = parse_patched_agp(args.patched_agp)
    contigs_agp = parse_agp(args.contigs_agp, c=5)

    for fasta in SeqIO.parse(args.fasta, "fasta"):
        if re.match("seq", fasta.name):
            # Unplaced contig
            continue
        else:
            chrom = get_chrom_name(reference_agp, patched_agp[fasta.name], contigs_agp)
            fasta.id = chrom
            sys.stdout.write("%s\n" % fasta.format("fasta"))
        
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
