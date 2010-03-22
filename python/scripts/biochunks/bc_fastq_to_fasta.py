#!/usr/bin/env python2.6

import os
import sys

def chop_reads(lines):
    for line in lines:
        if line.startswith('@'):
            chunk = ['>'+line[1:]]
            chunk.append(lines.next())
            lines.next()
            lines.next()
            chunk.append('>'+lines.next()[1:])
            chunk.append(lines.next())
            yield chunk

def main():

    chunks = chop_reads(sys.stdin.xreadlines())

    for chunk in chunks:
        for line in chunk:
            sys.stdout.write(line)


if __name__ == '__main__':
    main()
