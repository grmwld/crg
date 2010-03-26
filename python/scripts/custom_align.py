#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import optparse
from AGBio.Utilities import *
from AGBio.io.Fasta import *

def main():

    parser = optparse.OptionParser()

    parser.add_option('-i', '--in',
                      dest='infile',
                      help='input file')

    parser.add_option('-o', '--out',
                      dest='outfile',
                      help='output file')

    parser.add_option('-m', '--method',
                      dest='method',
                      help='alignment method. tcoffee|mafft|auto.' \
                      'If set to auto, will use mafft if tcoffee fails')

    parser.add_option('-a', '--n_core',
                      dest='ncore', type='int',
                      help='If method is tcoffee, how many threads.' \
                      'If it is mafft, how many parallel jobs to launch.')

    parser.add_option('-u', '--replace',
                      dest='replace',
                      help='which symbol to replace. Must be formated as o:r')

    parser.add_option('-t', '--temp',
                      dest='temp',
                      help='temporary root directory')

    parser.set_defaults(temp = '/tmp',
                        method = 'auto',
                        ncore = 1)

    (opts, args) = parser.parse_args()

    if opts.method not in ['auto', 'tcoffee', 'mafft']:
        parser.error('Wrong method.')

    
