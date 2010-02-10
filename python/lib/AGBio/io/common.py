#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

from __future__ import with_statement

def writeline(outf, line=''):
    outf.write(line + '\n')
    
def formatText(text, width=80, leftpad=0):
    result = ' ' * leftpad
    intext = ''.join((text.split()))
    for i, c in enumerate(intext):
        if i != 0 and i%width == 0:
            result += '\n'
            result += ' ' * leftpad
        result += c
    return result

def isFileEmpty(filename):
    with open(filename, 'r') as f:
        cont = f.read().strip()
        if cont:
            return False
    return True