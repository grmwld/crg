#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import sys
import os
import re

filein = open(sys.argv[1], 'r')
file2 = open (sys.argv[2], 'r')

f2 = file2.readlines()
ff = ''.join(f2)

for i in filein:
    try:
        if not re.search(i.split('|')[1], ff):
            print i
    except IndexError:
        sys.stderr.write(i)

filein.close()
file2.close()
