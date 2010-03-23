#!/usr/bin/env python2.6

import os
import sys
import optparse

def chop_reads(lines):
    for line in lines:
        if line[0] == '>':
            chunk = [line]
            chunk.append(lines.next())
            chunk.append(lines.next())
            chunk.append(lines.next())
            yield chunk

def trim(chunks, num):
    for chunk in chunks:
        yield [chunk[0], chunk[1][num:],
               chunk[2], chunk[3][num:]]

def main():

    parser = optparse.OptionParser()

    parser.add_option('-n', '--num',
                      dest='num', type='int',
                      help='Number of chars to be trimmed at the left end of reads',
                      metavar='INT')

    (opt, args) = parser.parse_args()

    chunks = trim(chop_reads(sys.stdin.xreadlines()),
                  opt.num)

    for chunk in chunks:
        map(sys.stdout.write, chunk)


if __name__ == '__main__':
    main()
