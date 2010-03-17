#!/usr/bin/env python2.6

import os
import sys
import optparse
import AGBio.io.Fasta as Fasta
from AGBio.io.NCBI import PsiBlastXMLParser
from AGBio.UtilityWrappers import *
from AGBio.ncbi.BlastWrappers import *

def uniq(infile):
    '''Use to remove duplicate headers, keeping the one with the lowest evalue.
    '''
    headict={}
    with open(infile) as ifile:
        for line in ifile:
            sline = (' '.join(line.split()[:-1]), line.split()[-1])
            if sline[0] in headict:
                if float(sline[1]) < float(headict[sline[0]]):
                    headict[sline[0]] = sline[1]
            else:
                headict[sline[0]] = sline[1]
    with open(infile, 'w') as ofile:
        for header, evalue in headict.items():
            ofile.write(header + ' ' + evalue + '\n')

def getTopSeqs(seqs, maxnumseqs=400, startevalue=10, pattern=None,
               verbose=False ):
    '''Function used to find the evalue yielding a number of sequences
    inferior to the specified threshold.

    Arguments:
    - `seqs`: the sequences.
    - `maxnumseqs`: the maximum number of sequences that is needed.
    - `startevalue`: iterative evalue decrementing should start with this power.
                     (e.g. startevalue=3 ==> evalue = 1e3)
    - `keepU`: should any U containing sequence be kept ?
    - `verbose`: should the function be verbose ?

    Returns a tuple (sequences, evalue) containing the sequences found and the
    corresponding evalue.
    '''
    topseqs = Fasta.SequenceList()
    evalue = startevalue + 1
    count = 0
    while (len(topseqs) > maxnumseqs or evalue == startevalue + 1) and evalue >= -400:
        prev_top_seqs = topseqs
        topseqs = Fasta.SequenceList()
        evalue -= 1
        for seq in seqs:
            if (pattern and pattern in seq.sequence) or\
                   float(seq.header.split()[-1]) <= float('1e' + str(evalue)):
                topseqs.append(seq)
        if verbose:
            sys.stderr.write( '        >>> evalue : 1e' + str(evalue) + '.\n' )
            sys.stderr.write( '            ' + str(len(topseqs)) + ' matching sequences\n' )
    return (topseqs, evalue)

def fmtOptPat(param, args):
    '''format the patterns.
    '''
    if args:
        return tuple([tuple((param, arg)) for arg in args.split(',')])
    return None

def main():

    parser = optparse.OptionParser()

    parser.add_option( '-i', '--inputfile',
                       dest='inputfilename',
                       help='blast output file, in xml format.',
                       metavar='FILE.xml' )

    parser.add_option( '-o', '--outputfile',
                       dest='outputfilename',
                       help='base output filename',
                       metavar='FILE' )

    parser.add_option( '-d', '--db',
                       dest='database',
                       help='database from which the sequences should be fetched.',
                       metavar='FILE' )

    parser.add_option( '-e', '--evalue',
                       dest='evalue',
                       type='float',
                       help='e-value threshold.',
                       metavar='FLOAT' )

    parser.add_option( '-E', '--start_expo_evalue',
                       dest='startexpoeval',
                       type='int',
                       help='exponent of the evalue threshold used when refiltering.',
                       metavar='INT' )
    
    parser.add_option( '-b', '--blast_version',
                       dest='blastversion',
                       help='set the blast version to use, either `legacy` or `plus`.',
                       metavar='VERSION' )
    
    parser.add_option( '-f', '--filter',
                       action='store_true', dest='dofilter', default=False,
                       help='do the filter step.')

    parser.add_option( '-p', '--keep_patterns_iff',
                       dest='keeppatiff',
                       help='Keep only if patterns match exactly. The patterns should be coma seperated.',
                       metavar='pat1,pat2,pat3' )

    parser.add_option( '-q', '--keep_patterns',
                       dest='keeppat',
                       help='Keep patterns that match exactly, no matter what. The patterns should be coma seperated.',
                       metavar='pat1,pat2,pat3' )

    parser.add_option( '-F', '--format',
                       dest='formatop',
                       help='format of the output. default is `header,evalue`',
                       metavar='INTEGER' )
                       
    parser.add_option( '-M', '--max_num_start_seq',
                       dest='maxnumstartseq',
                       type='int',
                       help='maximum number of sequences in the first alignement to be' +\
                       'processed. If set, a new input file with the top sequences ordered' +\
                       'by evalue is created and used.',
                       metavar='INTEGER' )
    
    parser.add_option( '-k', '--keep_U',
                       action='store_true', dest='keepu', default=False,
                       help='Should U containing sequences be kept regardless of their evalues ?.'+\
                       'Use in conjunction of -M')

    parser.add_option( '-T', '--temp',
                       dest='temp',
                       help='set the temp folder to use.',
                       metavar='FOLDER' )

    parser.add_option( '-S', '--simple',
                       dest='do_simple', action='store_true', default=False,
                       help='do not do extra fancy steps. Just parse the file and return the disired output in a file.' )
    
    parser.add_option( '-v', '--verbose',
                       dest='verbosity',
                       type='int',
                       help='verbosity level : 0=none ; 1=standard ; 2=detailed ; 3=full',
                       metavar='INTEGER' )

    parser.set_defaults( verbosity = 1,
                         database = 'nr',
                         evalue = 10,
                         startexpoeval = -10,
                         keeppat = None,
                         blastversion = 'legacy',
                         temp = '/tmp/',
                         maxnumstartseq = None,
                         formatop = 'header,evalue')

    (options, args) = parser.parse_args()

    verbosity = options.verbosity
    database = options.database
    evalue = options.evalue
    temp = options.temp
    maxnumstartseq = options.maxnumstartseq

    blastindexfile = ''.join(( options.outputfilename, '.index.0' ))
    blastfastafile = ''.join(( options.outputfilename, '.fasta.0' ))

    os.system(' '.join(( 'touch', blastindexfile )))
    os.system(' '.join(( 'touch', blastfastafile )))

    if options.blastversion == 'legacy':
        fetcher = FastaCmdWrapper( entry=[],
                                   db=database,
                                   outfile=blastfastafile )
    else:
        fetcher = BlastDbCmdWrapper( entry=[],
                                     db=database,
                                     outfile=blastfastafile )

    ## Parse the blast output file.
    if verbosity >= 1:
        sys.stderr.write( '\n' )
        sys.stderr.write( '>>> Parsing blast output : ' +\
                          options.inputfilename + '\n' )
    with open(options.inputfilename, 'r') as infile:
        blastparser = PsiBlastXMLParser(infile)
        blastparser.parse()
        if verbosity >= 2:
            sys.stderr.write('    >>> Extracting required data.\n')
        if options.dofilter:
            sequences = blastparser.extractData( evalue=evalue,
                                                 fmt=options.formatop,
                                                 outfile=blastindexfile,
                                                 includepatternsiff=fmtOptPat('title', options.keeppatiff),
                                                 includepatterns=fmtOptPat('title', options.keeppat),
                                                 excludepatterns=(('title', 'hypothetical'),
                                                                  ('title', 'predicted'),
                                                                  ('title', 'PREDICTED')))
        else:
            sequences = blastparser.extractData( evalue=evalue,
                                                 fmt=options.formatop,
                                                 outfile=blastindexfile )
        if options.do_simple:
            sys.exit(0)
        
    ## Only keep one copy of a header, the one with the best evalue.
    if verbosity >= 1:
        sys.stderr.write( '\n' )
        sys.stderr.write( '>>> Keeping only best evalues.\n' )
    uniq(blastindexfile)
    
    ## Gather all GIs in list
    if verbosity >= 2:
        sys.stderr.write( '\n' )
        sys.stderr.write( '>>> Gathering all Gis.\n' )
    entries = []
    with open(blastindexfile, 'r') as bif:
        for line in bif:
            entries.append(line.split('|')[1])
    fetcher.entry = entries

    ## Fetch the sequences from the local databases.
    ## TODO : Fetch failed from the web.
    if verbosity >= 1:
        sys.stderr.write( '\n' )
        sys.stderr.write( '>>> Building fasta.0 file by fetching sequences from local database.\n' )
    fetcher.run()

    ## Apply final filters : keep only top evalues and U containing until a threshold is reached
    if maxnumstartseq:
        if verbosity >= 1:
            sys.stderr.write( '\n' )
            sys.stderr.write( '>>> Applying final filters on ' + \
                              blastfastafile + '.\n' )
        if verbosity >= 2:
            sys.stderr.write( '    >>> Adding evalue to headers.\n' )
        ### TODO : use .fasta.fh in tmp dir.
        tmpfullheadfasta = blastfastafile + '.fh'
        addheaders = AddFullHeadersWrapper2(blastfastafile,
                                           tmpfullheadfasta,
                                           blastindexfile)
        addheaders.run()
        if verbosity >= 3:
            sys.stderr.write( '        >>> Loading sequences.\n' )
        with open(tmpfullheadfasta, 'r') as ff:
            allseqs = Fasta.loadSequences(ff)
        if verbosity >= 2:
            sys.stderr.write( '    >>> Keeping valid sequences.\n' )
        tmppat = None
        if options.keepu:
            tmppat = 'U'
        validseqs = getTopSeqs(seqs=allseqs,
                               maxnumseqs=maxnumstartseq,
                               startevalue=options.startexpoeval,
                               pattern=tmppat,
                               verbose=verbosity>=4 )
        keptseqs = '.'.join(( options.outputfilename,
                              str(validseqs[1]),
                              str(len(validseqs[0])),
                              'fasta' ))
        if verbosity >= 2:
            sys.stderr.write( '    >>> Found ' + str(len(validseqs[0])) + \
                              ' sequences with evalue <= 1e' + \
                              str(validseqs[1]) + '\n' )
        with open(keptseqs, 'w') as ff:
            validseqs[0].save(ff)

    sys.stderr.write( '\n' )


if __name__ == '__main__':
    main()
