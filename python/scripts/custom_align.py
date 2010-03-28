#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import optparse
import string
from AGBio.Utilities import *
from AGBio.io.Fasta import *
from AGBio.UtilityWrappers import *

def find_positions(sstr, symbol):
    opos = []
    for i, c in enumerate(sstr):
        if c == symbol:
            opos.append(i)
    return opos

def backplace(sstr, positions, o, r):
    output = ''
    ngaps = 0
    for i, c in enumerate(sstr):
        if c == '-':
            ngaps += 1
        if (i - ngaps) in positions:
            output += r
        else:
            output += c
    return output


def main():

    parser = optparse.OptionParser()

    parser.add_option('-i', '--in',
                      dest='infile',
                      help='input file')

    parser.add_option('-o', '--out',
                      dest='outfile',
                      help='output file')

    parser.add_option('-m', '--method',
                      dest='method',
                      help='alignment method. tcoffee|mafft|auto.' \
                      'If set to auto, will use mafft if tcoffee fails')

    parser.add_option('-a', '--n_core',
                      dest='ncore', type='int',
                      help='If method is tcoffee, how many threads.' \
                      'If it is mafft, how many parallel jobs to launch.')

    parser.add_option('-u', '--replace',
                      dest='replace',
                      help='which symbol to replace. Must be formated as o:r')

    parser.add_option('-t', '--temp',
                      dest='temp',
                      help='temporary root directory')

    parser.set_defaults(temp = '/tmp',
                        method = 'auto',
                        ncore = 1)

    (opts, args) = parser.parse_args()

    if opts.method not in ['auto', 'tcoffee', 'mafft']:
        parser.error('Wrong method.')

    data_dict = {}
    before, after = opts.replace.split(':')

    with open(opts.infile, 'r') as iff:
        sequences = loadSequences(iff)

    ## replaces the symbol in each sequence and save an index of headers
    ori_headers = genTempfilename(opts.temp)
    replaced_nonal = genTempfilename(opts.temp)
    with open(ori_headers, 'w') as hof:
        with open(replaced_nonal, 'w') as rof:
            for sequence in sequences:
                tpos = find_positions(sequence.sequence, before)
                sequence.replacePattern(before, after)
                data_dict[sequence.header] = tpos
                hof.write(sequence.header+'\n')
                sequence.prints(rof)

    ## write to a temp file the modified sequences
    with open(replaced_nonal, 'w') as off:
        sequences.prints(off)

    ## align the sequences and write to a temp file
    tmp_output_al = genTempfilename(opts.temp)
    if opts.method in ['tcoffee', 'auto']:
        tcoffee = TcoffeeWrapper(infile=replaced_nonal,
                                 outfile=tmp_output_al,
                                 ncore=opts.ncore)
        print tcoffee.cline
        tcoffee.run()
    elif opts.method in ['mafft', 'auto']:
        mafft = MafftWrapper(infile=replaced_nonal,
                             outfile=tmp_output_al)
        print mafft.cline
        mafft.run()

    ## replace the headers with those saved in the index
    tmp_output_filled = genTempfilename(opts.temp)
    header_filler = AddFullHeadersWrapper2(infile=tmp_output_al,
                                           outfile=tmp_output_filled,
                                           patternfile=ori_headers,
                                           method='inplace')
    print header_filler.cline
    header_filler.run()

    ## replace back the initial symbol in each sequence
    opts.outfile = os.path.abspath(os.path.expanduser(opts.outfile))
    with open(tmp_output_filled, 'r') as iff:
        sequences = loadSequences(iff)
        with open(opts.outfile, 'w') as off:
            for s in sequences:
                backplaced_seq = backplace(s.sequence,
                                           data_dict[s.header],
                                           after, before)
                Sequence(s.header, backplaced_seq).prints(off)

   
if __name__ == '__main__':
    main()
    
