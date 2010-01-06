#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import sys
import os
import optparse
sys.path.append('/users/rg/agrimaldi/Code/crg/python/libs')
import AGBio.IO.Fasta as FastaLib


def findSelenoproteins( sequences ):

    return sequences.findPattern('U', mode='full')


def main():

    parser = optparse.OptionParser()

    parser.add_option( '-i', '--inputfile',
                       dest='inputfilename',
                       help='fasta file in which selenoproteins should be looked for.',
                       metavar='FILE' )

    parser.add_option( '-o', '--outputfile',
                       dest='outputfilename',
                       help='fasta file containing the selenoproteins',
                       metavar='FILE' )

    parser.add_option( '-v', '--verbose',
                       dest='verbosity',
                       help='verbosity level : 0=none ; 1=standard ; 2=detailed ; 3=full',
                       metavar='INTEGER' )

    parser.set_defaults( verbosity = '1' )

    (options, args) = parser.parse_args()

    stdoutflag = False

    verbosity = int( options.verbosity )

    if options.inputfilename:

        inputfilenames = options.inputfilename.split(',')
        infiles = []
        for i in inputfilenames:
            infiles.append( open( i, 'r' ) )

    else: sys.exit( 'You must provide an input filename.')

    if options.outputfilename:

        outfile = open( options.outputfilename, 'w' )
        stdoutflag = True

    else: outfile = sys.stdout

    for f in infiles:
        if verbosity >= 1:
            print
            print '>>> Searching for selenoproteins in file ' + f.name
            print

        if verbosity >= 2:
            print '>>> Loading sequences ...'

        sequences = FastaLib.loadSequences( f )

        if verbosity >= 2:
            print '>>> ... Done.'
            print

        if verbosity >= 2:
            print '>>> Searching for U containing sequences ...'

        selenoproteins = findSelenoproteins( sequences )

        if verbosity >= 2:
            print '>>> ... Done.'
            print

        FastaLib.saveSequences(selenoproteins, outfile)

        for selP in selenoproteins:

            if verbosity >= 3 and stdoutflag:
                print selP.header.strip()
                print selP.sequence.strip()

        if verbosity >= 1:
            print
            print 'Found ' + str( len( selenoproteins ) ) + ' selenoproteins'
            print

    for i in infiles:
        i.close()
    outfile.close()


if __name__ == '__main__':

    main()
