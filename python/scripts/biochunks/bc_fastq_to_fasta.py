#!/usr/bin/env python2.6

import os
import sys

def main():

    lines = sys.stdin.xreadlines()

    for line in lines:
        if line[0] == '@':
            sys.stdout.write('>'+line[1:])
            sys.stdout.write(lines.next())
            lines.next()
            lines.next()
            sys.stdout.write('>'+lines.next()[1:])
            sys.stdout.write(lines.next())


if __name__ == '__main__':
    main()
