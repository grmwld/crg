#!/usr/bin/env python2.6

import os
import sys
import optparse
import AGBio.io.Fasta as Fasta
from AGBio.io.NCBI import PsiBlastXMLParser
from AGBio.UtilityWrappers import *
from AGBio.ncbi.BlastWrappers import *

TEMP = '/tmp/'

def uniq(infile):
    '''Use remove duplicate headers, keeping the one with the lowest evalue.
    '''
    headict={}
    with open(infile) as ifile:
        for line in ifile:
            sline = line.split()
            if sline[0] in headict:
                if float(sline[1]) < float(headict[sline[0]]):
                    headict[sline[0]] = sline[1]
            else:
                headict[sline[0]] = sline[1]
    with open(infile, 'w') as ofile:
        for header, evalue in headict.items():
            ofile.write(header + ' ' + evalue + '\n')


def main():

    parser = optparse.OptionParser()

    parser.add_option( '-i', '--inputfile',
                       dest='inputfilename',
                       help='fasta file in which selenoproteins should be looked for.',
                       metavar='FILE' )

    parser.add_option( '-o', '--outputfile',
                       dest='outputfilename',
                       help='base output filename',
                       metavar='FILE' )

    parser.add_option( '-e', '--evalue',
                       dest='evalue',
                       help='e-value threshold.',
                       metavar='FLOAT' )

    parser.add_option( '-p', '--pattern',
                       dest='pattern',
                       help='pattern to look for in sequences.',
                       metavar='REGEX')

    parser.add_option( '-b', '--blast_version',
                       dest='blastversion',
                       help='set the blast version to use, either `legacy` or `plus`.',
                       metavar='VERSION' )

    parser.add_option( '-v', '--verbose',
                       dest='verbosity',
                       help='verbosity level : 0=none ; 1=standard ; 2=detailed ; 3=full',
                       metavar='INTEGER' )

    parser.set_defaults( verbosity = '1',
                         evalue = '10',
                         pattern = None,
                         blastversion = 'legacy' )

    (options, args) = parser.parse_args()

    verbosity = int(options.verbosity)
    evalue = options.evalue
    pattern = options.pattern

    blastindexfile = ''.join(( options.outputfilename, '.index.0' ))
    blastfastafile = ''.join(( options.outputfilename, '.fasta.0' ))

    os.system(' '.join(( 'touch', blastindexfile )))
    os.system(' '.join(( 'touch', blastfastafile )))

    if options.blastversion == 'legacy':
        fetcher = FastaCmdWrapper( entry=[],
                                   db='/seq/databases/nr_uncompressed/nr',
                                   outfile=blastfastafile )
    else:
        fetcher = BlastDbCmdWrapper( entry=[],
                                     db='nr',
                                     outfile=blastfastafile )


    ## Parse the blast output file.
    if verbosity >= 1:
        sys.stderr.write( '\n' )
        sys.stderr.write( '>>> Parsing blast output.\n' )
    with open(options.inputfilename, 'r') as infile:
        blastparser = PsiBlastXMLParser(infile)
        blastparser.parse()
        if verbosity >= 1:
            sys.stderr.write('>>> Extracting required data.\n\n')
        sequences = blastparser.extractData( evalue=evalue,
                                             fmt='id,evalue',
                                             outfile=blastindexfile)

    ## Only keep one copy of a header, the one with the best evalue.
    if verbosity >= 1:
        sys.stderr.write( '\n' )
        sys.stderr.write( '>>> Keeping only best evalues.\n' )
    niq(blastindexfile)

    ## Gather all GIs in list
    entries = []
    with open(blastindexfile, 'r') as bif:
        for line in bif:
            entries.append(line.split('|')[1])
    fetcher.entry = entries

    ## Fetch the sequences from the local databases
    ## TODO : Fetch failed from the web.
    if verbosity >= 1:
        sys.stderr.write( '\n' )
        sys.stderr.write( '>>> Building fasta.0 file.\n' )
    fetcher.run()

if __name__ == '__main__':
    main()
