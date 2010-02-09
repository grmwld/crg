#!/usr/bin/env python2.6

import os
import sys
import shutil
import optparse
from AGBio.Utilities import *
from AGBio.UtilityWrappers import *
from AGBio.io.Fasta import *

def incList(llist, d1, d2):
    return [llist[0]+d1, llist[1]+d2]

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

    parser.add_option( '-n', '--num_seqs',
                       dest='numseqs',
                       help='Number of sequences to achieve.',
                       type='int',
                       metavar='INT')

    parser.add_option( '-r', '--res_overlap',
                       dest='resoverlap',
                       help='Use fixed value of the resoverlap.',
                       type='float',
                       metavar='FLOAT')

    parser.add_option( '-p', '--percent_deviation',
                       dest='perdev',
                       help='percentage of devaition allowed.',
                       type='int',
                       metavar='INT')

    parser.add_option( '-c', '--max_identical_cycles',
                       dest='maxidcycles',
                       help='maximum number of cycles to perform in case the number of sequences does not varies anymore.',
                       type='int',
                       metavar='INT')

    parser.add_option( '-C', '--max_total_cycles',
                       dest='maxtotalcycles',
                       help='maximum number of total cycles to perform.',
                       type='int',
                       metavar='INT')

    parser.add_option( '-t', '--temp',
                       dest='temp',
                       help='temporary folder.',
                       metavar='FOLDER')

    parser.add_option( '-D', '--debug',
                       action='store_true', dest='debug', default=False,
                       help='Debug.')

    parser.set_defaults(temp='/home/agrimaldi/temp/',
                        perdev=15,
                        maxidcycles=5,
                        resoverlap=None,
                        maxtotalcycles=20)

    (options, args) = parser.parse_args()

    num_seqs = 99999
    range_seqs = range(options.numseqs*(100-options.perdev),
                       options.numseqs*(100+options.perdev))
    idcycles = 0
    prevnumcycles = num_seqs
    bestfile = options.inputfilename
    bestnum = num_seqs
    ncycle = 0

    tmpfname = options.inputfilename
    if not options.resoverlap:
        overlapvals = [0.5, 50]
        minvals = [0, 0]
        maxvals = [1, 100]
    else:
        overlapvals = [options.resoverlap, 50]
        minvals = [options.resoverlap, 0]
        maxvals = [options.resoverlap, 100]
        
    trimmer = TrimalWrapper(infile=options.inputfilename,
                            outfile='boo')

    tmpdir = genTempfilename(options.temp, 'CT_')
    os.mkdir(tmpdir)

    try:
        while idcycles < options.maxidcycles \
                  or num_seqs in range_seqs \
                  or options.maxtotalcycles >= ncycle :
            verb = str(overlapvals) + str(minvals) + str(maxvals) + str(ncycle) + ''
            sys.stdout.write(verb)
            tmpfname = genTempfilename(tmpdir)
            trimmer.outfile = tmpfname
            trimmer.scoreoverlap = overlapvals
            trimmer.run()
            with open(tmpfname, 'r') as tf:
                num_seqs = len(loadSequences(tf))
            if num_seqs > options.numseqs:
                minvals = overlapvals[:]
                overlapvals = incList(overlapvals,
                                      (maxvals[0]-overlapvals[0])/2,
                                      round((maxvals[1]-overlapvals[1])/2))
            elif num_seqs < options.numseqs:
                maxvals = overlapvals[:]
                overlapvals = incList(overlapvals,
                                      (minvals[0]-overlapvals[0])/2,
                                      round((minvals[1]-overlapvals[1])/2))
            else:
                break
            if num_seqs == prevnumcycles:
                idcycles += 1
            else:
                idcycles = 0
            prevnumcycles = num_seqs
            if num_seqs > 0 and abs(options.numseqs - num_seqs) < bestnum:
                bestnum = num_seqs
                bestfile = tmpfname
            ncycle += 1
            print idcycles, bestnum, num_seqs
    except IOError, KeyboardInterrupt:
        pass
    finally:
        print
        print 'copying', bestfile, 'to', options.outputfilename
        shutil.copy(bestfile, options.outputfilename)
        if not options.debug:
            shutil.rmtree(tmpdir)

if __name__ == '__main__':
    main()
