#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import sys
import os
import optparse
sys.path.append('/users/rg/agrimaldi/Code/crg/python/libs')
from AGBio.io.Fasta import *
from AGBio.selenoprofiles_tools.files_analysers.ali import *

def parse_file(filename):
    pass

def main():

    parser = optparse.OptionParser()

    parser.add_option('-i', '--input',
                      dest='infilenames',
                      help='input filenames as a stream',
                      metavar='FILES')

    parser.add_option('-s', '--symbols',
                      dest='symbols',
                      help='symbols to keep track of during the parsing')

    parser.add_option('-p', '--display',
                      dest='display',
                      help='display the statistics relative to this symbol')

    parser.set_defaults(symbols = 'UC-',
                        display = 'U')

    (opts, args) = parser.parse_args()

    ali = AliData(opts.infilenames)

    ali.find_x_positions(aas=opts.symbols)

    hitUcount = 0
    for seq in ali.ali:
        if 'hit' in seq.header and 'U' in seq.sequence:
            hitUcount += 1

    print hitUcount
    ali.u_redundant.prints(opts.display)

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(1)
