#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import sys
import os
import time
import shutil
import subprocess
import optparse
from AGBio.Utilities import *
import AGBio.io.Fasta as Fasta
import AGBio.UtilityWrappers as UtilityWrappers


def spDiff( fileA, fileB ):

    sfileA = open( fileA, 'r' )
    sfileB = open( fileB, 'r' )

    seqsA = Fasta.loadSequences( sfileA )
    seqsB = Fasta.loadSequences( sfileB )

    sfileA.close()
    sfileB.close()

    spA = seqsA.findPattern('U', mode='full')
    spB = seqsB.findPattern('U', mode='full')
    
    return  spA.symetric_difference(spB, method='raw')


def removeGaps( sequence ):

        return sequence.replacePattern('-', '')
    

def main():

    parser = optparse.OptionParser()

    parser.add_option( '-i', '--inputfile',
                       dest='inputfilename',
                       help='file containing the alignments that will be used to build the PSSM using prepare_alignment_selenoprofiles.py.',
                       metavar='FILE' )

    parser.add_option( '-r', '--datadir',
                       dest='datadir',
                       help='directory containing, for each familly FAM, a directory FAM.blast and a directory FAM.selenoprofiles.prep',
                       metavar='DIR' )

    parser.add_option( '-o', '--outputfile',
                       dest='outputfilename',
                       help='base name used for outputs',
                       metavar='NAME' )

    parser.add_option( '-a', '--n_core',
                       dest='ncore',
                       type='int',
                       help='number of cores to use during the various operations.',
                       metavar='INTEGER' )

    parser.add_option( '-M', '--mafft',
                       action='store_true', dest='domafft', default=False,
                       help='do the mafft step.')
    
    parser.add_option( '-T', '--trimal',
                       action='store_true', dest='dotrimal', default=False,
                       help='do the trimal step.')

    parser.add_option( '-C', '--tcoffee',
                       action='store_true', dest='dotcoffee', default=False,
                       help='do the t_coffee step.')

    parser.add_option( '-B', '--headers',
                       action='store_true', dest='doheaders', default=False,
                       help='do the addheaders step.')

    parser.add_option( '-p', '--patternfile',
                       dest='patternfile',
                       help='pattern file to use if the -D option is used.',
                       metavar='FILE' )

    parser.add_option( '-F', '--filter',
                       action='store_true', dest='dofilter', default=False,
                       help='do the filter step.')

    parser.add_option( '-P', '--prepare',
                       action='store_true', dest='doprepal', default=False,
                       help='do the prepare_alignment_selenoprofiles step.')

    parser.add_option( '-g', '--tag_threshold',
                       dest='tagthreshold',
                       type='float',
                       help='tag threshold to use if the -P or --prepare is used.',
                       metavar='FLOAT' )

    parser.add_option( '-A', '--all',
                       action='store_true', dest='doall', default=False,
                       help='do all steps.')

    parser.add_option( '-Y', '--dry',
                       action='store_true', dest='dryrun', default=False,
                       help="Prints the commands without executing them.")
    
    parser.add_option( '-D', '--debug',
                       action='store_true', dest='debug', default=False,
                       help="Debug mode. Nothing is cleaned.")

    parser.add_option( '-t', '--temp',
                       dest='temp',
                       help='set the temp folder to use.',
                       metavar='FOLDER' )

    parser.add_option( '-v', '--verbose',
                       dest='verbosity',
                       type='int',
                       help='verbosity level : 0=none ; 1=standard ; 2=detailed ; 3=full',
                       metavar='INTEGER' )



    parser.set_defaults( verbosity = 1,
                         ncore = 1,
                         tagthreshold = 0.5,
                         temp = '/tmp/',
                         patternfile = 'None' )

    (options, args) = parser.parse_args()

    if options.doall:
        options.doheaders = True
        options.dofilter = True
        options.domafft = True
        options.dotrimal = True
        options.dotcoffee = True
        
    infile = options.inputfilename
    tmpinitfilename = genTempfilename(options.temp, 'ungapped_')
    with open(infile, 'r') as iff:
        tmpseqs = Fasta.loadSequences(iff)
    with open(tmpinitfilename, 'w') as ugf:
        for seq in tmpseqs:
            removeGaps(seq).prints(ugf)
    tmpinfile = tmpinitfilename

    mafftoutfile = ''.join((options.outputfilename, '_mafft.fasta'))
    trimaloutfile1 = ''.join((options.outputfilename, '_trimmed_native.fasta'))
    trimaloutfile2 = ''.join((options.outputfilename, '_trimmed_spadded.fasta'))
    trimaloutfile = trimaloutfile1
    tcoffeeoutfile = ''.join((options.outputfilename, '_tcoffee.fasta'))
    fullheadoutfile = ''.join((options.outputfilename, '.det.fasta'))
#    patternfile = ''.join(('.'.join(options.inputfilename.split('.')[:2]), '.index.0'))
    patternfile = options.patternfile
    filteroutfile = ''.join((options.outputfilename, '.filt.fasta'))

    ncore = options.ncore
    verbosity = options.verbosity
    temp = options.temp

    addheaders = UtilityWrappers.AddFullHeadersWrapper(tmpinfile,
                                                       fullheadoutfile,
                                                       patternfile)

    filterseqs = UtilityWrappers.FilterWrapper(tmpinfile,
                                               filteroutfile,
                                               inverse=True,
                                               titlematch=('PREDICTED', 'predicted', 'hypothetical'))

    mafft = UtilityWrappers.MafftWrapper(tmpinfile,
                                         mafftoutfile,
                                         auto=True)

    trimal = UtilityWrappers.TrimalWrapper(tmpinfile,
                                           trimaloutfile1,
                                           clusters=100)

    tcoffee = UtilityWrappers.TcoffeeWrapper(trimaloutfile2,
                                             tcoffeeoutfile,
                                             ncore=ncore)

    prepsp = UtilityWrappers.SelenoprofilesPreWrapper(tmpinfile,
                                                      options.outputfilename,
                                                      all=True,
                                                      tagthreshold=options.tagthreshold,
                                                      temp=temp)

    try:

        if options.dryrun:
            print('\nThis is a dry run. Relaunch the command without the option -Y to do the actual stuff.\n')



        ## Add full headers
    ##     if options.doheaders:
    ##         addheader.infile = tmpinfile
    ##         tmpinfile = fullheadoutfile
    ##         if options.dryrun:
    ##             print addheaders.cline
    ##         else:
    ##             if verbosity >= 1:
    ##                 sys.stderr.write('\n    >>> Adding headers\n\n')
    ##             addheaders.run()

        ## Filter out the 'fake' proteins
        if options.dofilter:
            time.sleep(0.5)
            filterseqs.infile = tmpinfile
            tmpinfile = filteroutfile
            if options.dryrun:
                print filterseqs.cline
            else:
                if verbosity >= 1:
                    sys.stderr.write('\n    >>> Filtering out\n\n')
                filterseqs.run()

        ## run mafft
        numseqinmafftoutput = 0 
        if options.domafft:
            time.sleep(0.5)
            mafft.infile = tmpinfile
            tmpinfile = mafftoutfile
            if options.dryrun:
                print mafft.cline
            else:
                if verbosity >= 1:
                    sys.stderr.write('\n    >>> Running Mafft\n\n')
                mafft.run()
                with open(mafftoutfile, 'r') as mfo:
                    seqs = Fasta.loadSequences(mfo)
                    numseqinmafftoutput = len(seqs)

        ## run trimal
        if options.dotrimal and numseqinmafftoutput > 200:
            time.sleep(0.5)
            trimal.infile = tmpinfile
            tmpinfile = trimaloutfile1
            if options.dryrun:
                print trimal.cline
            else:
                if verbosity >= 1:
                    sys.stderr.write('\n    >>> Running Trimal\n\n')
                trimal.run()

        if not options.dryrun and options.dotcoffee and options.dotrimal:
            if verbosity >= 1:
                sys.stderr.write('\n    >>> Removing gaps\n\n')

            ti = open(tmpinfile, 'r')
            tmpinfile = trimaloutfile2
            to = open(tmpinfile, 'w')

            si = Fasta.loadSequences(ti)
            ti.close()
            refs = Fasta.SequenceList()

            ## saves the sequences with no gaps
            for s in si:
                refs.append(removeGaps(s))
            Fasta.saveSequences(refs, to)

            if options.dotrimal and numseqinmafftoutput > 200:
                if verbosity >= 1:
                    sys.stderr.write('\n    >>> Adding ommited selenoproteins\n')
                ## Gather the non intersecting proteins from the 2 files
                diffSelenoproteins = spDiff( mafftoutfile,
                                             trimaloutfile1 )

                spDiffr = Fasta.SequenceList()
                ## remove gaps from selenoproteins
                for s in diffSelenoproteins:
                    spDiffr.append(removeGaps(s))
                ## append to the file the selenoproteins that were not present
                Fasta.saveSequences(spDiffr, to)
            to.close()

        ## run t_coffee
        if options.dotcoffee:
            time.sleep(0.5)
            tcoffee.infile = tmpinfile
            tmpinfile = tcoffeeoutfile
            if options.dryrun:
                print tcoffee.cline
            else:
                if verbosity >= 1:
                    sys.stderr.write('\n    >>> Running T_coffee\n\n')
                tcoffee.run()

        ## Add full headers
        if options.doheaders:
            time.sleep(0.5)
            addheaders.infile = tmpinfile
            tmpinfile = fullheadoutfile
            if options.dryrun:
                print addheaders.cline
            else:
                if verbosity >= 1:
                    sys.stderr.write('\n    >>> Adding headers\n\n')
                addheaders.run()

        ## prepare alignments for selenoprofiles
        if options.doprepal:
            time.sleep(0.5)
            prepsp.infile = tmpinfile
            if options.dryrun:
                print prepsp.cline
            else:
                if verbosity >= 1:
                    sys.stderr.write('\n    >>> preparing for selenoprofiles\n\n')
                prepsp.run()

    except KeyboardInterrupt:
        sys.exit('manual exit.')
    finally:
        if not options.debug:
            if verbosity >= 2:
                sys.stderr.write('\n    >>> Removing temporary file ' + tmpinitfilename +'\n\n')
            os.remove(tmpinitfilename)

if __name__ == '__main__':
        main()
