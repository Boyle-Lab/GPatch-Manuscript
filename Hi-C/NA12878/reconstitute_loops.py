#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser

__version__ = "0.0.1"

"""
Given a BED file with loop anchor annotations, reconstitute
all loop records for which both anchors are available and print
results in BEDPE format.
"""

def parse_anchors_into_loop_dict(bed):
    """
    Given a set of BED records for loop anchors, generate a
    BEDPE within a dict using all available pairs of loop anchors.
    """
    ret = {}
    with open(bed, "r") as bed_f:
        for line in bed_f:
            rec = line.strip().split()
            # Get the root loop ID and anchor #
            ID_rec = rec[3].split('_')
            loop_id = ID_rec[0]
            anchor_id = ID_rec[1]
            if loop_id in ret:
                ret[loop_id].append(rec)
            else:
                ret[loop_id] = [rec]
    return ret


def reconstitute_bedpe(loop_id, rec):
    """
    Given a record from the loop_records dict, reconstitute a BEDPE
    line. Prints '.' as placeholders for missing anchor records.
    """
    if len(rec) < 2:
        # Singleton. See if we have the 5' or 3' anchor and construct
        # a BEDPE record with blanks for the opposite anchor.
        anchor_id = rec[0][3].split('_')[1]
        if anchor_id == 1:
            sys.stdout.write("%s\t%d\t%d\t.\t.\t.\t%s\n" % (rec[0][0], int(rec[0][1]), int(rec[0][2]), loop_id))
        else:
            sys.stdout.write(".\t.\t.\t%s\t%d\t%d\t%s\n" % (rec[0][0], int(rec[0][1]), int(rec[0][2]), loop_id))
    else:
        sys.stdout.write("%s\t%d\t%d\t%s\t%d\t%d\t%s\n" % (rec[0][0], int(rec[0][1]), int(rec[0][2]), rec[1][0], int(rec[1][1]), int(rec[1][2]), loop_id))


def main():
    parser = ArgumentParser(description='Given a BED file with loop anchor annotations, reconstitute all loop records for which both anchors are available and print output in BEDPE format.')
    parser.add_argument('-b', '--bed', metavar='BED', type=str,
                        required=True, help='Path to BED file containing loop anchor records.')
    
    
    args = parser.parse_args()

    loop_records = parse_anchors_into_loop_dict(args.bed)
    for loop_id in loop_records.keys():
        reconstitute_bedpe(loop_id, loop_records[loop_id])
    

    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
