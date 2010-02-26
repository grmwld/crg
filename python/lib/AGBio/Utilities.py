#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import os
import sys
import random


def getBaseName(filename, numextension=1):
    return '.'.join(filename.split('.')[:-numextension])

def changeFileExtension(filename, newextension, numextension=1):
    return getBaseName(filename, numextension) + '.' + newextension

class _GetchUnix:
    def __init__(self):
        import tty, sys

    def __call__(self, msg=None):
        import sys, tty, termios
        sys.stdout.write(msg+'\n')
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(sys.stdin.fileno())
            ch = sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch
getch = _GetchUnix()

def genTempfilename(dirr='./', prefix=''):
    '''Generate temporary filename
    '''
    hashn = str(hex(random.getrandbits(128))[2:-1])
    return os.path.join(dirr, prefix + hashn)

def flatten(x):
    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result


def prop(fcn):
    return property(**fcn())


class BiDict(dict):
    def __init__(self, dictionary={}):
        self.update(dictionary)
    @prop
    def reverse():
        def fget(self):
            output = {}
            for k, v in self.items():
                if v not in output:
                    output[v] = list(k)
                else:
                    output[v].append(k)
            return output
        return locals()
