#!/usr/bin/env python2.6
# -*- coding:utf-8 -*-

import sys
import os
import optparse
from AGBio.ncbi.BlastWrappers import *


if __name__ == '__main__':

    parser = optparse.OptionParser()
    
    parser.add_option( '-i', '--gi',
                       dest='ginum',
                       help='GI accession number(s) to use. If several values \
                       are provided, must be coma separated.',
                       metavar='GI' )

    parser.add_option( '-o', '--out',
                       dest='outfilename',
                       help='Name of the output filename.',
                       metavar='NAME' )

    parser.add_option( '-d', '--db',
                       dest='db',
                       help='Blast database to use to retrieve the sequence.',
                       metavar='DATABASE' )
    
    parser.add_option( '-w', '--web',
                       dest='web',
                       action='store_true',
                       help='',
                       default=False )

    parser.add_option( '-v', '--verbose',
                       dest='verbose',
                       help='Be a bit verbose',
                       action='store_true',
                       default=False)
    
    (options, args) = parser.parse_args()

    if 'ginum' not in options.__dict__:
        parser.error('You must provide a GI number, or "all" to get everythong.')

    if options.web and options.db:
        parser.error('-w [--web] and -d [--db] are mutually exclusive.\n' \
                     +'You must choose how to fetch the sequence, either ' \
                     +'online or using a local database.')

    fetcher = BlastDbCmdWrapper(entry=options.ginum.split(','),
                                db=options.db,
                                outfile=options.outfilename)

    if options.verbose:
        print fetcher.cline
        
    fetcher.run()
    
