#!/usr/bin/env python2.6

import os
import sys
import optparse

def filter_X(lines, char='N', max_N=0):
    '''Filters chunks based on the composition of their sequences
    '''
    for line in lines:
        if line[0] == '@':
            chunk = [line]
            chunk.append(lines.next())        
            chunk.append(lines.next())        
            chunk.append(lines.next())        
            chunk.append(lines.next())        
            chunk.append(lines.next())        
            chunk.append(lines.next())        
            chunk.append(lines.next())        
            if chunk[1].count(char) <= max_N and chunk[5].count(char) <= max_N:
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

    chunks = filter_X(sys.stdin.xreadlines(),
                      opt.forbiden_char, opt.max_occurence)

    for chunk in chunks:
        map(sys.stdout.write, chunk)
 

if __name__ == '__main__':
    main()
