#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import os

def getBaseName(filename, numextension=1):
    return '.'.join(filename.split('.')[:-numextension])

def changeFileExtension(filename, newextension, numextension=1):
    return getBaseName(filename, numextension) + '.' + newextension
