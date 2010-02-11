#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

from __future__ import with_statement
import os
import sys
import shutil

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


class FolderWrapper(BaseFileWrapper, Container):
    def __init__(self, name):
        BaseFileWrapper.__init__(self, name)

    def create_file(self, name):
        with open(name, 'w') as mkf:
            pass
        return FileWrapper(name)

    def create_folder(self, name):
        os.mkdir(name)
        return FolderWrapper(name)


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
