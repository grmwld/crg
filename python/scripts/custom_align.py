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
        if c == o and (i - ngaps) in positions:
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
                      'If set to auto, mafft will be used if the number ' \
                      'of sequences exceeds 1000, otherwise, tcoffee will' \
                      'be used')

    parser.add_option('-r', '--replace_method',
                      dest='replace_method',
                      help='method used to replace headers.')
    
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

    parser.add_option('-d', '--debug',
                      dest='debug', action='store_true',
                      default=False, help='do not delete the temp files')

    parser.set_defaults(temp = '/tmp',
                        method = 'auto',
                        ncore = 1,
                        replace_method = 'gi')

    (opts, args) = parser.parse_args()

    if opts.method not in ['auto', 'tcoffee', 'mafft']:
        parser.error('Wrong method.')

    try:
        
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
        if opts.method != 'auto':
            al_method = opts.method
        else:
            if len(sequences) > 1000:
                al_method = 'mafft'
            else:
                al_method = 'tcoffee'
        tmp_output_al = genTempfilename(opts.temp)
        if al_method == 'tcoffee':
            tcoffee = TcoffeeWrapper(infile=replaced_nonal,
                                     outfile=tmp_output_al,
                                     ncore=opts.ncore)
            print tcoffee.cline
            tcoffee.run()
        elif al_method == 'mafft':
            mafft = MafftWrapper(infile=replaced_nonal,
                                 outfile=tmp_output_al)
            print mafft.cline
            mafft.run()

        ## replace the headers with those saved in the index
        tmp_output_filled = genTempfilename(opts.temp)
        header_filler = AddFullHeadersWrapper2(infile=tmp_output_al,
                                               outfile=tmp_output_filled,
                                               patternfile=ori_headers,
                                               method=opts.replace_method)
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

    except Exception, (e):
        raise e
    finally:
        if not opts.debug:
            os.remove(tmp_output_filled)
            os.remove(tmp_output_al)
            os.remove(ori_headers)
            os.remove(replaced_nonal)

        
if __name__ == '__main__':
    main()
    
