#!/usr/bin/env python2.6

import os
import sys
import optparse
import AGBio.io.Fasta as Fasta

def main():

    parser = optparse.OptionParser()

    parser.add_option( '-i', '--inputfile',
                       dest='inputfilename',
                       help='file with incomplete headers.',
                       metavar='FILE' )

    parser.add_option( '-o', '--outputfile',
                       dest='outputfilename',
                       help='outputfile.',
                       metavar='FILE' )

    parser.add_option( '-p', '--pattern',
                       dest='patternfilename',
                       help='pattern file containing the complete headers.',
                       metavar='FILE' )

    parser.set_defaults( outputfilename = None )

    (options, args) = parser.parse_args()

    if not (options.inputfilename and options.patternfilename):
        parser.error('You have to provide two files, check help.')

    with open(options.inputfilename, 'r') as iff:
        inlines = Fasta.loadSequences(iff)
    with open(options.patternfilename, 'r') as pff:
        patlines = [line for line in pff.readlines() \
                    if line.startswith('>gi|') \
                    or line.startswith('gi|')]

    if not options.outputfilename:
        outfile = sys.stdout
    else:
        outfile = open(options.outputfilename, 'w')

    for iseq in inlines:
        nofound = True
        for phead in patlines:
            try:
                if iseq.header.split('|')[1].split()[0] == phead.split('|')[1].split()[0]:
                    tmpseq = Fasta.Sequence(phead, iseq.sequence)
                    tmpseq.prints(outfile)
                    nofound = False
                    break
            except IndexError as e:
                sys.stderr.write( '\nError while processing the files:\n' )
                sys.stderr.write( pline + '\n' )
                sys.stderr.write( line + '\n' )
                break
        if nofound:
            sys.stderr.write('\n' + iseq.header + '\n')

    outfile.close()


if __name__ == '__main__':
    main()
            

