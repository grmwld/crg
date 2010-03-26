#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import sys
import os
import optparse
from string import maketrans
from AGBio.io.Fasta import *

def main():

    parser = optparse.OptionParser()

    parser.add_option('-i', '--infile',
                      dest='infile',
                      help='input fasta file',
                      metavar='FILE')

    parser.add_option('-t', '--trans_table',
                      dest='trans_table',
                      help='characters to be replaced in the format abcd:efgh')

    parser.add_option('-d', '--deltete_chars',
                      dest='delete_chars',
                      help='characters to be deleted')

    parser.add_option('-c', '--case_insensitive',
                      action='store_true', dest='insensitive', default=False,
                      help='case insensitive match')

    parser.set_defaults(delete_chars = None)

    (opt, args) = parser.parse_args()

    if opt.trans_table:
        before, after = opt.trans_table.split(':')
        if len(before) != len(after):
            parser.error('number of characters in before and after differs.')
        print before, type(before), after, type(after)
        transtable = maketrans(before, after)
    else:
        transtable = None

    with open(opt.infile) as iff:
        sequences = loadSequences(iff)

    for s in sequences:
        t = Sequence(s.header,
                     s.sequence.translate(transtable))
        t.prints()


if __name__ == '__main__':
    main()
