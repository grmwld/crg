#!/usr/bin/env python2.6

import os
import sys
import optparse
import re
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

    parser.add_option( '-l', '--inplace',
                       action='store_true',
                       dest='inplace' )

    parser.add_option( '-g', '--gi',
                       action='store_false',
                       dest='inplace' )

    parser.set_defaults( outputfilename = None )

    (options, args) = parser.parse_args()

    if not (options.inputfilename and options.patternfilename):
        parser.error('You have to provide two files, check help.')

    with open(options.inputfilename, 'r') as iff:
        inlines = Fasta.loadSequences(iff)
    with open(options.patternfilename, 'r') as pff:
        patlines = [line for line in pff.readlines() \
                    if line.startswith('>')]

    if not options.outputfilename:
        outfile = sys.stdout
    else:
        outfile = open(options.outputfilename, 'w')

    if not options.inplace:
        GI_REGEX = re.compile(r'gi\|(\d+)\|')

        for iseq in inlines:
            nofound = True
            for phead in patlines:
                try:
                    giq = GI_REGEX.search(iseq.header).group(1)
                    gis = GI_REGEX.search(phead).group(1)
                    if giq == gis:
                        tmpseq = Fasta.Sequence(phead, iseq.sequence)
                        tmpseq.prints(outfile)
                        nofound = False
                        break
                except AttributeError as e:
                    sys.stderr.write(iseq.header + ' ' + phead)
                    sys.exit(-1)
                except IndexError as e:
                    sys.stderr.write( '\nError while processing the files:\n' )
                    sys.stderr.write( pline + '\n' )
                    sys.stderr.write( line + '\n' )
                    break
            if nofound:
                sys.stderr.write('\n' + iseq.header + '\n')
    else:
        if len(inlines) != len(patlines):
            raise Exception, 'Different number of sequences'
        for seq, pat in zip(inlines, patlines):
            Fasta.Sequence(pat, seq.sequence).prints(outfile, 60)

    outfile.close()


if __name__ == '__main__':
    main()
            

