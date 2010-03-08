#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

from __future__ import with_statement
import os
import sys
import shutil
from AGBio.Utilities import *

class BaseFileWrapper(object):
    '''Wrapper for a file.
    '''
    def __init__(self, filename):
        self.name = filename
        self.abspath = os.path.abspath(filename)


class Openable(object):
    def openf(cls, mode):
        try:
            cls.fobj = cls.name.open(cls.name, mode)
        except Exception, (e):
            print e


class Closable(object):
    '''Interface for object that can be closed (ie : files)
    '''
    def closef(cls):
        try:
            cls.obj.close()
        except Exception, (e):
            print e


class Container(object):
    def sub_dirs(cls):
        return self._sub(os.path.isdir, FolderWrapper)

    def sub_files(cls):
        return self._sub(os.path.isfile, FileWrapper)

    def get_sub_links(cls):
        return self._sub(os.path.islink)

    def _sub(cls, evaluator, classtype):
        llist = []
        for dir in [d for d in os.listdir(cls.name) if evaluator(d)]:
            llist.append(classtype(dir))
        return llist
    

class FileWrapper(BaseFileWrapper, Openable, Closable):
    def __init__(self, name):
        BaseFileWrapper.__init__(self, name)

    def create(self):
        with open(self.abspath, 'w'):
            pass

    def delete(self):
        os.remove(self.abspath)


class FolderWrapper(BaseFileWrapper, Container):
    def __init__(self, name):
        BaseFileWrapper.__init__(self, name)

    def create(self):
        os.mkdir(self.abspath)

    def delete(self, recursive=False):
        if recursive:
            shutil.rmtree(self.abspath)
        else:
            os.rmdir(self.abspath)

    def sub_file(self, name=None):
        if name:
            return FileWrapper(os.path.join(self.abspath, name))
        else:
            return FileWrapper(genTempfilename(self.abspath))

    def sub_folder(self, name):
        return FolderWrapper(os.path.join(self.abspath, name))


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
