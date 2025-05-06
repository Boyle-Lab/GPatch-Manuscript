#!/usr/bin/env python3

import sys, os, re
from argparse import ArgumentParser
import coolbox
from coolbox.api import *

__version__ = "0.0.1"

"""
Given a Hi-C matrix in .hic format, produce a set of matrix
heatmaps for all chromosomes.
"""

def parse_chromsizes(chrom_sizes_f):
    """
    Parse a chrom.sizes file into a dict.
    """
    ret = {}
    with open(chrom_sizes_f, "r") as infile:
        for line in infile:
            rec = line.strip().split()
            ret[rec[0]] = rec
    return ret

STDERR = ''

def new_stderr(old):
    """
    This is a hack to detect when the chromosome range in chrom.sizes does not
    match what CoolBox expects. Since no exception is thrown, we need to check
    for stderr output. This captures writes to sys.stderr.write in the STDERR 
    variable. See https://stackoverflow.com/questions/31816190/how-to-redirect-stderr-to-variable-in-python-2-7
    """
    def new(*args):
        # put your code here, you will intercept writes to stderr
        #print('Intercepted: ' + repr(args))
        global STDERR   # add new write to STDERR 
        STDERR += args[0]
        # Uncomment below to still print STDERR content to terminal.
        #old(*args)
    return new


def find_end_coord(brw, chrom, cs_end_coord):
    """
    Find the end coordinate for a chromosome as it exists in CoolBox.
    Reads global STDERR, which captures output redirected from STDERR,
    to determine if brw.goto was successful or threw an error, since
    there is no excpetion raised. Returns the first coordinate that
    does not throw an error as the end coordinate.
    """
    global STDERR
    ec = cs_end_coord
    # This loop setup assumes we enter the loop after being called
    # due to an error with brw.goto, which should leave the "is not valid"
    # string in stderr content.
    while re.search("is not valid.", STDERR):
        # We must reset STDERR each time or we will never break the loop!
        STDERR = ''
        # We assume here that the end coordinate given was the
        # coordinate that raised the original error, and decrement first.
        ec -= 1
        #sys.stdout.write("ec: %d\n" % (ec))
        coords = chrom + ":1-" + str(ec)
        brw.goto(coords)
    # Reset STDERR before we go, to be safe!
    STDERR = ''
    return ec


def main():
    parser = ArgumentParser(description='Given a Hi-C matrix in .hic format, produce a set of matrix heatmaps for all chromosomes.')
    parser.add_argument('-m', '--matrix', metavar='FILE', type=str,
                        required=True, help='Path to Hi-C matrix, in .hic format.')
    parser.add_argument('-s', '--chrom_sizes', metavar='FILE', type=str,
                        required=True, help='Path to chrom.sizes file for the reference genome used in Hi-C analysis.')
    parser.add_argument('-b', '--contigs_bed', metavar='BED', type=str, default="",
                        required=False, help='Path to BED file with contig coordinates.')
    parser.add_argument('-k', '--skip', metavar='FILE', type=str,
                        required=False, default="", help='Comma-delimited list of chromosomes to skip plotting.')
    parser.add_argument('-t', '--plot', metavar='FILE', type=str,
                        required=False, default="", help='Comma-delimited list of chromosomes to plot. (All others are skipped.)')
    parser.add_argument('-p', '--prefix', metavar='STR', type=str,
                        required=False, default="", help='Prefix for outfile names.')
    parser.add_argument('-f', '--format', metavar='FILE', type=str,
                        required=False, default="png", help='Output file format. Default="png"')
    parser.add_argument('-r', '--resolution', metavar='INT', type=int,
                        required=False, default=80000,  help='Resolution to plot. Must be present in the Hi-C matrix or will throw an error. Default=80000kb')
    parser.add_argument('-c', '--cmap', metavar='INT', type=str,
                        required=False, default="JuiceBoxLike",  help='Color map to use. See https://gangcaolab.github.io/CoolBox/api.html#coolbox.core.track.HicMatBase for available options. Default="JuiceBoxLike"')
    parser.add_argument('-n', '--normalize', metavar='NORM', type=str,
                        required=False, default="False", help='Type of normalization to apply. See https://gangcaolab.github.io/CoolBox/api.html#coolbox.core.track.HicMatBase for options. Default=False (no normalization)')
    parser.add_argument('-l', '--balance', action='store_true',
                        required=False, default=False, help='Use the balanced contact matrix.')
    
    args = parser.parse_args()

    chrom_sizes = parse_chromsizes(args.chrom_sizes)

    to_skip = []
    if args.skip != "":
        to_skip = args.skip.split(',')

    normalize = False
    if args.normalize != "False":
        normalize = args.normalize
    hic_matrix = HiCMat(args.matrix, resolution=args.resolution, cmap=args.cmap, style='triangular', normalize=normalize, balance=args.balance)
    if args.contigs_bed != "":
        contigs_bed = BED(args.contigs_bed, display='interlaced', labels=False)
    
    # Loop over chromosomes, producing a plot for each.
    chroms = chrom_sizes.keys()
    if args.plot != "":
        chroms = args.plot.split(',')
    for chrom in chroms:
        if chrom in to_skip:
            continue
        sys.stderr.write("Plotting %s\n" % (chrom))

        # Get the length of the chromosome and set up the viewer frame.
        frame = hic_matrix + Title(chrom) + XAxis()
        if args.contigs_bed != "":
            frame = frame + contigs_bed + Title("Contigs")
        #frame.plot(coords_str)
        # Create the browser with the current frame.
        brw = Browser(frame)
        global STDERR
        old_stderr = sys.stderr.write
        sys.stderr.write = new_stderr(sys.stderr.write)
        coords = chrom + ":1-" + chrom_sizes[chrom][1]
        brw.goto(coords)
        if re.search("is not valid.", STDERR):
            # The chrom.sizes end coordinate does not match what CoolBox expects.
            ec = find_end_coord(brw, chrom, int(chrom_sizes[chrom][1]))
            coords = chrom + ":1-" + str(ec)
            #sys.stdout.write("%s\n" % (coords))
            brw.goto(coords)
        sys.stderr.write = old_stderr
        # Finally, save the image with an appropriate filenmame.
        fname = chrom + '.' + args.format
        if args.prefix != "":
            fname = args.prefix + '.' + fname
        brw.save(fname)
    
    sys.stderr.write("Done!\n")
        
if __name__ == '__main__':
    main()
    exit(0)
