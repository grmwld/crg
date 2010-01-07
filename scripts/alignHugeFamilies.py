#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import sys
import os
import shutil
import subprocess
import optparse
#sys.path.append('/users/rg/agrimaldi/Code/crg/python/libs')
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

def getTopSeqs( filename, threshold, minevalue=10 ):
    num_seq = int(subprocess.Popen(["wc",
                                    "-l",
                                    changeFileExtension(filename, 'index.0', 2)],
                                   stdout=subprocess.PIPE).communicate()[0].split()[0])
    evalue = 
    #gawkprc = subprocess.Popen(["gawk", "{if ($2 < 1e"+str(evalue)+"){print $1}}", changeFileExtension(filename, 'index.0', 2)], stdout=subprocess.PIPE)
    #wcprc = subprocess.Popen(["wc", "-l"], stdin=gawkprc.stdout)

    while threshold < num_seq:
        evalue -= 1
        gawkprc = subprocess.Popen(["gawk", "{if ($2 < 1e"+str(evalue)+"){print $1}}", changeFileExtension(filename, 'index.0', 2)], stdout=subprocess.PIPE)
        num_seq = int(subprocess.Popen(["wc", "-l"], stdin=gawkprc.stdout, stdout=subprocess.PIPE).communicate()[0].split()[0])

    return ([gi.split('|')[1] for gi in subprocess.Popen(["gawk", "{if ($2 < 1e"+str(evalue)+"){print $1}}", changeFileExtension(filename, 'index.0', 2)], stdout=subprocess.PIPE).communicate()[0].split()])

def getTopSeqs2(indexfile, threshold, fastafile)

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

    parser.add_option( '-a', '--n_core',
                       dest='ncore',
                       help='number of cores to use during the various operations.',
                       metavar='INTEGER' )

    parser.add_option( '-m', '--mafft',
                       action='store_true', dest='domafft', default=False,
                       help='do the mafft step.')
    
    parser.add_option( '-t', '--trimal',
                       action='store_true', dest='dotrimal', default=False,
                       help='do the trimal step.')

    parser.add_option( '-c', '--tcoffee',
                       action='store_true', dest='dotcoffee', default=False,
                       help='do the t_coffee step.')

    parser.add_option( '-d', '--headers',
                       action='store_true', dest='doheaders', default=False,
                       help='do the addheaders step.')

    parser.add_option( '-f', '--filter',
                       action='store_true', dest='dofilter', default=False,
                       help='do the filter step.')

    parser.add_option( '-p', '--prepare',
                       action='store_true', dest='doprepal', default=False,
                       help='do the prepare_alignment_selenoprofiles step.')

    parser.add_option( '-M', '--max_num_start_seq',
                       dest='maxnumstartseq',
                       help='maximum number of sequences in the first alignement to be processed. If set, a new input file with the top sequences ordered by evalue is created and used.',
                       metavar='NAME' )

    parser.add_option( '-A', '--all',
                       action='store_true', dest='doall', default=False,
                       help='do all steps.')

    parser.add_option( '-Y', '--dry',
                       action='store_true', dest='dryrun', default=False,
                       help="Prints the commands without executing them.")
    
    parser.add_option( '-T', '--temp',
                       dest='temp',
                       help='set the temp folder to use.',
                       metavar='FOLDER' )

    parser.add_option( '-v', '--verbose',
                       dest='verbosity',
                       help='verbosity level : 0=none ; 1=standard ; 2=detailed ; 3=full',
                       metavar='INTEGER' )



    parser.set_defaults( verbosity = '1',
                         ncore = '1',
                         temp = '/tmp/')

    (options, args) = parser.parse_args()

    if options.doall:
        options.doheaders = True
        options.dofilter = True
        options.domafft = True
        options.dotrimal = True
        options.dotcoffee = True
        
    
    infile = options.inputfilename
    tmpinfile = infile

    mafftoutfile = ''.join((options.outputfilename, '_mafft.fasta'))
    trimaloutfile1 = ''.join((options.outputfilename, '_trimmed_native.fasta'))
    trimaloutfile2 = ''.join((options.outputfilename, '_trimmed_spadded.fasta'))
    trimaloutfile = trimaloutfile1
    tcoffeeoutfile = ''.join((options.outputfilename, '_tcoffee.fasta'))
    fullheadoutfile = ''.join((options.outputfilename, '.det.fasta'))
    patternfile = ''.join(('.'.join(options.inputfilename.split('.')[:2]), '.index.0'))
    filteroutfile = ''.join((options.outputfilename, '.filt.fasta'))

    ncore = options.ncore
    verbosity = options.verbosity

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
                                                      tagthreshold=0.8,
                                                      temp=temp)

    if options.dryrun:
        print('\nThis is a dry run. Relaunch the command without the option -Y to do the actual stuff.\n')



    ## Add full headers
    if options.doheaders:
        addheader.infile = tmpinfile
        tmpinfile = fullheadoutfile
        if options.dryrun:
            print addheaders.cline
        else:
            if verbosity >= 1:
                sys.stderr.write('\n    >>> Adding headers\n\n')
            addheaders.run()

    ## Filter out the 'fake' proteins
    if options.dofilter:
        filterseqs.infile = tmpinfile
        tmpinfile = filteroutfile
        if options.dryrun:
            print filterseqs.cline
        else:
            if verbosity >= 1:
                sys.stderr.write('\n    >>> Filtering out\n\n')
            filterseqs.run()
    
    ## Keep the N best sequences based on their evalues
    ## if options.maxnumstartseq:
##         seqs = getTopSeqs(infile, int(options.maxnumstartseq), keepU=True)
##         refetch()
    ##with open()
        
    ## run mafft
    if options.domafft:
        mafft.infile = tmpinfile
        tmpinfile = mafftoutfile
        if options.dryrun:
            print mafft.cline
        else:
            if verbosity >= 1:
                sys.stderr.write('\n    >>> Running Mafft\n\n')
            mafft.run()

    ## run trimal
    if options.dotrimal:
        trimal.infile = tmpinfile
        tmpinfile = trimaloutfile2
        if options.dryrun:
            print trimal.cline
        else:
            if verbosity >= 1:
                sys.stderr.write('\n    >>> Running Trimal\n\n')
            trimal.run()

    if not options.dryrun:
        if verbosity >= 1:
            sys.stderr.write('\n    >>> Removing gaps\n\n')

        ti = open(trimaloutfile1, 'r')
        to = open(trimaloutfile2, 'w')

        si = Fasta.loadSequences(ti)
        ti.close()
        refs = Fasta.SequenceList()

        ## saves the sequences with no gaps
        for s in si:
            refs.append(removeGaps(s))

        Fasta.saveSequences(refs, to)

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
        tcoffee.infile = tmpinfile
        tmpinfile = tcoffeeoutfile
        if options.dryrun:
            print tcoffee.cline
        else:
            if verbosity >= 1:
                sys.stderr.write('\n    >>> Running T_coffee\n\n')
            tcoffee.run()
    
    ## prepare alignments for selenoprofiles
    if options.doprepal:
        prepsp.infile = tmpinfile
        if options.dryrun:
            print prepsp.cline
        else:
            if verbosity >= 1:
                sys.stderr.write('\n    >>> preparing for selenoprofiles\n\n')
            prepsp.run()

    ## Add full headers
    if options.doheaders:
        if verbosity >= 1:
            sys.stderr.write('\n    >>> Adding headers\n\n')
#        addheader.
#        addheaders.run()

if __name__ == '__main__':
    main()
