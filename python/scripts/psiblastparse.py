#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import sys
import os
import optparse
import AGBio.IO.Fasta as Fasta
import AGBio.IO.NCBI as NCBI


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

    parser.add_option( '-e', '--evalue',
                       dest='evalue',
                       help='e-value threshold.',
                       metavar='FLOAT' )

    parser.add_option( '-a', '--header_pattern',
                       dest='hpattern',
                       help='pattern that should match in the headers of the results.',
                       metavar='PATTERN' )

    parser.add_option( '-b', '--sequence_pattern',
                       dest='spattern',
                       help='pattern that should match in the sequences of the results.',
                       metavar='PATTERN' )

    parser.add_option( '-v', '--verbose',
                       dest='verbosity',
                       help='verbosity level : 0=none ; 1=standard ; 2=detailed ; 3=full',
                       metavar='INTEGER' )

    parser.set_defaults( verbosity = '1',
                         evalue = '10',
                         hpattern = '',
                         spattern = '' )
    
    (options, args) = parser.parse_args()

    verbosity = int( options.verbosity )
    evalue = float( options.evalue )
    hpattern = options.hpattern
    spattern = options.spattern

    infile = open(options.inputfilename, 'r')

    if options.outputfilename:
        outfile = open(options.outputfilename + '.fasta', 'w')
    else:
        outfile = sys.stdout
    
    if verbosity >= 1:
        sys.stderr.write('>>> Parsing the blast output file.\n\n')

    psi_parser = NCBI.PsiBlastParser()
    psiblast_record = psi_parser.parse(infile)

    if verbosity >= 1:
        sys.stderr.write('>>> Extracting required data.\n\n')

    sequences = psi_parser.extractData( psiblast_record,
                                        evalue=evalue,
                                        spattern=spattern)
                                        #hpattern=hpattern )#.uniques(method='headers')

    if verbosity >= 1:
        sys.stderr.write('>>> Saving results.\n\n')

    sequences.save(outfile)

    infile.close()
    if outfile == sys.stdout:
        outfile.close()

    if verbosity >= 1:
        sys.stderr.write('>>> done.\n')


if __name__ == '__main__':
    
    main()
