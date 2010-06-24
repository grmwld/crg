#!/usr/bin/env python2.6

from __future__ import with_statement
import os
import sys
import optparse
from AGBio.io.Fasta import *

def main():

    parser = optparse.OptionParser()

    parser.add_option('-i', '--in',
                      dest='inputfile',
                      help='Input file.',
                      metavar='FILE')

    parser.add_option('-o', '--out',
                      dest='outputfile',
                      help='Output file.',
                      metavar='FILE')

    (opt, args) = parser.parse_args()

    with open(opt.inputfile) as iff:
        seqs = loadSequences(iff)

    checked = set()

    with open(opt.outputfile, 'w') as off:
        for seq in seqs:
            if seq.sequence not in checked:
                checked.add(seq.sequence)
                seq.prints(off)

if __name__ == '__main__':
    main()
