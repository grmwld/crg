#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import os
import sys

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
