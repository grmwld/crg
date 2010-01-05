#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import sys
import os
import re
import optparse
import AGBio.IO.Fasta as Fasta


def main():

    parser = optparse.OptionParser()

    parser.add_option( '-i', '--inputfile',
                       dest='inputfilename',
                       help='file containing the alignments that will be used to build the PSSM using prepare_alignment_selenoprofiles.py.',
                       metavar='FILE' )

    (options, args) = parser.parse_args()
    
    infile = open(options.inputfilename, 'r')

    sequences = Fasta.loadSequences( infile ).uniques(method='headers')

    sequences.save(sys.stdout)

    infile.close()

if __name__ == '__main__':

    main()
