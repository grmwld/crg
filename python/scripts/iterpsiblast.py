#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import subprocess
import optparse

GENCHK = 'blastpgp -C chk.asn'
BLAST = 'blastpgp -F F -b 5000 -I'

def main():

    parser = optparse.OptionParser()

    parser.add_option( '-i', '--inputfile',
                       dest='inputfilename',
                       help='file containing the alignments that will be used to build the PSSM using prepare_alignment_selenoprofiles.py.',
                       metavar='FILE' )

    parser.add_option( '-o', '--outputfile',
                       dest='outputfilename',
                       help='base name used for outputs',
                       metavar='NAME' )

    parser.add_option( '-d', '--database',
                       dest='database',
                       help='database used by blast.',
                       metavar='File' )

    parser.add_option( '-k', '--checkpoint',
                       dest='checkpoint',
                       help='base name of the generated checkpoints.',
                       metavar='NAME' )

    parser.add_option( 'j', '--n_iteration',
                       dest='niteration',
                       help='number of successive iterations to perform.',
                       metavar='INTEGER' )

    parser.add_option( '-c', '--n_core',
                       dest='ncore',
                       help='number of cores to use during blast',
                       metavar='INTEGER' )

    parser.set_defaults( ncore = '1',
                         niteration = '1' )

    (options, args) = parser.parse_args()

    inputfile = options.inputfilename
    outputfile = options.outputfilename
    database = options.database
    ncore = options.ncore
    niteration = options.niteration

    for i in range( niteration ):

        
