#!/usr/bin/env python2.6

import os
import sys
import shutil
import multiprocessing
import subprocess
import optparse
from AGBio.Utilities import *
from AGBio.UtilityWrappers import *
from AGBio.io.Fasta import *

def incList(llist, d1, d2):
    return [llist[0]+d1, llist[1]+d2]

def sliceInterval(minv, maxv, ncuts):
    ssize = (float(maxv) - float(minv)) / (ncuts + 1)
    cut_sites = []
    tmin = float(minv)
    for i in range(ncuts):
        cut_sites.append(tmin + ssize)
        tmin += ssize
    return cut_sites

def inWhichSliceIs(slices, num, reverse=False):
    if num in slices:
        if num == slices[0]:
            return (slices[:2], (0, 1))
        elif num == slices[-1]:
            return slices[-2:], (len(slices)-2, len(slices)-1)
        else:
            return ((slices[slices.index(num)-1],
                     slices[slices.index(num)+1]),
                    (slices.index(num)-1,
                     slices.index(num)+1))
    result = []
    dmin = abs(slices[0] - slices[-1])
    for i, s in enumerate(slices[:-1]):
        d = abs(num - s) + abs(num - slices[i+1])
        print d
        if d < dmin:
            dmin = d
            reslut = ((s, slices[i+1]), (i, i+1))

def inWhichSliceIs2(trimmers, num, reverse=False):
    t_num_seqs = [t.num_seqs for t in trimmers]
    if not reverse:
        if num < trimmers[0].num_seq:
            return ((0, trimmers[0]), (0, 0))
        elif num > trimmers[-1].num_seqs:
            return ((trimmers[-1], 0), (len(trimmers)-1, 0))
    else:
        if num > trimmers[0].num_seqs:
            return ((0, trimmers[0]), (0, 0))
        elif num < trimmers[-1].num_seqs:
            return ((trimmers[-1], 0), (len(trimmers)-1, 0))

    if num in [t_num_seqs]:
        if num == trimmers[0].num_seqs:
            return (trimmers[:2], (0, 1))
        elif num == trimmers[-1].num_seqs:
            return trimmers[-2:], (len(trimmers)-2, len(trimmers)-1)
        else:
            return ((trimmers[trimmers.index(num)-1],
                     trimmers[trimmers.index(num)+1]),
                    (trimmers.index(num)-1,
                     trimmers.index(num)+1))
    for i, s in enumerate(trimmers[:-1]):
        if not reverse:
            if num > s.num_seqs and num < trimmers[i+1].num_seqs:
                return ((s, trimmers[i+1]), (i, i+1))
        else:
            if num < s.num_seqs and num > trimmers[i+1].num_seqs:
                return ((s, trimmers[i+1]), (i, i+1))

def worker(trimmer):
    """worker function"""
    ## subprocess.call(['trimal_dev',
##                      '-in', 'nifh_frxc_mafft.fasta',
##                      '-out', tmpout,
##                      '-resoverlap', str(scores[0]),
##                      '-seqoverlap', str(scores[1])])
    print trimmer.cline
    trimmer.run()
    return


class Trimmer(TrimalWrapper):
    def __init__(self, infile, outfile='boo', scoreoverlap=(0.5, 50)):
        TrimalWrapper.__init__(self, infile, outfile, scoreoverlap=scoreoverlap)
        self.num_seqs = 0

    def update_num_seq(self):
        with open(self['outfile'][1], 'r') as iff:
            self.num_seqs = len(loadSequences(iff))


def main():

    parser = optparse.OptionParser()

    parser.add_option( '-i', '--inputfile',
                       dest='inputfilename',
                       help='input alignment.',
                       metavar='FILE' )

    parser.add_option( '-o', '--outputfile',
                       dest='outputfilename',
                       help='base output filename',
                       metavar='FILE' )

    parser.add_option( '-a', '--n_core',
                       dest='ncore',
                       help='Number of core to use. `0` means find the max for this machine.',
                       type='int',
                       metavar='NCORE' )

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
                        maxtotalcycles=20,
                        ncore=1)

    (options, args) = parser.parse_args()

    with open(options.inputfilename, 'r') as ifff:
        t_num_seqs = len(loadSequences(ifff))
    num_seqs = t_num_seqs
    range_seqs = range(options.numseqs*(100-options.perdev),
                       options.numseqs*(100+options.perdev))
    idcycles = 0
    prevnumcycles = num_seqs
    bestfile = options.inputfilename
    bestnum = num_seqs
    ncycle = 0

    tmpfname = options.inputfilename
    if not options.resoverlap:
        overlapvals = zip(*[sliceInterval(0, 1, options.ncore),
                       sliceInterval(0, 100, options.ncore)])
        minvals = [0, 0]
        maxvals = [1, 100]
    else:
        overlapvals = zip(*[sliceInterval(options.resoverlap, 1, options.ncore),
                       sliceInterval(options.resoverlap, 100, options.ncore)])
        minvals = [options.resoverlap, 0]
        maxvals = [options.resoverlap, 100]

    trimmers = []
    for i in range(options.ncore):
        trimmer = Trimmer(infile=options.inputfilename,
                                outfile='boo')
        trimmers.append(trimmer)

    tmpdir = genTempfilename(options.temp, 'CT_')
    os.mkdir(tmpdir)

    try:
        while idcycles < options.maxidcycles \
                  or num_seqs in range_seqs \
                  or options.maxtotalcycles >= ncycle :
            jobs = []
            tmpofiles = []
            for i, trimmer in enumerate(trimmers):
                verb = str(overlapvals) + str(minvals) + str(maxvals)
                tmpfname = genTempfilename(tmpdir)
                tmpofiles.append(tmpfname)
                trimmer.outfile = tmpfname
                trimmer.scoreoverlap = overlapvals[i]
                p = multiprocessing.Process(target=worker, args=(trimmer,))
                jobs.append(p)
            for j in jobs:
                j.start()
            for j in jobs:
                j.join()
            for trimmer in trimmers:
                trimmer.update_num_seq()
            extremTrimmers = inWhichSliceIs2(trimmers, options.numseqs, True)
            ###
            print extremTrimmers
            for tt in extremTrimmers[0]:
                try:
                    print tt.scoreoverlap
                    print tt.num_seqs
                except Exception as e:
                    print e
            ###
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
            print verb, idcycles, num_seqs
    except (IOError, KeyboardInterrupt) as err:
        print err
    finally:
        print
        print 'copying', bestfile, 'to', options.outputfilename
        shutil.copy(bestfile, options.outputfilename)
        if not options.debug:
            shutil.rmtree(tmpdir)

if __name__ == '__main__':
    main()
