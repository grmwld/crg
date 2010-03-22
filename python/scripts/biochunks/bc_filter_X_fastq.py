#!/usr/bin/env python2.6

import os
import sys
import optparse

def filter_X(chunks, char='N', max_N=0):
    '''Filters chunks based on the composition of their sequences
    TODO : Make it more modular so that only the result of a function
    would be required to filter out the chunks
    '''
    for chunk in chunks:
        if chunk[1].count(char) <= max_N and chunk[5].count(char) <= max_N:
            yield chunk

def chop_reads(lines):
    for line in lines:
        if line.startswith('@'):
            chunk = [line]
            for i in xrange(7):
                chunk.append(lines.next())
            yield chunk

def main():

    parser = optparse.OptionParser()

    parser.add_option('-c', '--char',
                      dest='forbiden_char',
                      help='Which forbiden char should exclude reads ?',
                      metavar='CHAR')

    parser.add_option('-n', '--max_occurence',
                      dest='max_occurence',
                      help='Maximum number of occurence of the forbidden char',
                      type='int',
                      metavar='INT')

    (opt, args) = parser.parse_args()

    chunks = filter_X(chop_reads(sys.stdin.xreadlines()),
                      opt.forbiden_char, opt.max_occurence)

    for chunk in chunks:
        for line in chunk:
            sys.stdout.write(line)


if __name__ == '__main__':
    main()
