#!/usr/bin/env python2.6
# -*- coding:utf-8 -*-

from __future__ import with_statement
import os
import sys
import subprocess
import shutil
sys.path.append('/users/rg/mmariotti/libraries')
from MMlib import bbash
from AGBio.io.common import *
from AGBio.selenoprofiles_tools.files.p2g import *
from AGBio.selenoprofiles_tools.files.b_secisearch import *
            

class GenomeFolderParser(object):
    '''Parser of a genome folder (from selenoprofiles)

    Parses a set of folder inside the organisme folder.
    Typically, the folder to investigate is `output`,
    but could be `discarded`.

    Usage:

    parser = GenomeFolderParser('/path/to/folder')
    parser.parse(doall=True)
    print parser.cys
    print parser.sec
    ...
    '''
    def __init__(self, path, dirs=['output'], everywhere=False):
        self.rootdir = os.path.abspath(path)
        if everywhere:
            self.dirs = [os.path.join(self.rootdir, 'output'),
                         os.path.join(self.rootdir, 'discarded')]
        else:
            self.dirs = [os.path.join(self.rootdir, d) for d in dirs]
        self.notempty = []
        self.sec = {}
        self.cys = {}
        self.thr = {}
        self.arg = {}
        self.uga = {}
        self.ual = {}
        self.secis_std = {}
        self.secis_nonstd = {}
        self.secis_twil = {}
        self.secis_b = {}
        self.kw2dict = {'cysteine' : self.cys,
                        'selenocysteine' : self.sec,
                        'threonine' : self.thr,
                        'arginine' : self.arg,
                        'uga_containing' : self.uga,
                        'unaligned' : self.ual,
                        'std' : self.secis_std,
                        'non_std' : self.secis_nonstd,
                        'twil' : self.secis_twil,
                        'bsecis' : self.secis_b}

    def parse(self, doall=False,
              sec=False, cys=False, thr=False, arg=False,
              uga=False, ual=False, bsecis=False,
              stdsecis=False, nonstdsecis=False, twilsecis='False'):
        '''Method to parse the folder.

        The wanted informations are given as arguments. If `doall` is set to
        True, everything is parsed.

        Arguments:
        - `doall`: check for any output file.
        - `sec`: check for selenocysteine containing candidates.
        - `cys`: check for cysteine containing candidates.
        - `thr`: check for threonine containing candidates.
        - `arg`: check for arginine containing candidates.
        - `uga`: check for UGA containing candidates.
        - `ual`: check for unaligned candidates.
        - `stdsecis`: check for standard SECIS containing candidates.
        - `nonstdsecis`: check for non-std SECIS containing candidates.
        - `twilsecis`: check for twilight SECIS containing candidates.
        - `bsecis`: check for bSECIS containing candidates
        '''
        lukw = []
        if doall:
            lukw = self.kw2dict.keys()
        else:
            if sec: lukw.append('selenocysteine')
            if cys: lukw.append('cysteine')
            if thr: lukw.append('threonine')
            if arg: lukw.append('arginine')
            if uga: lukw.append('uga_containing')
            if ual: lukw.append('unaligned')
            if stdsecis : lukw.append('std')
            if nonstdsecis : lukw.append('non_std')
            if twilsecis : lukw.append('twil')
            if bsecis : lukw.append('bsecis')
        for d in self.dirs:
            resfiles = [f for f in os.listdir(d) if self._keep(f)]
            for ff in resfiles:
                for kw in lukw:
                    if kw in ff.split('.'):
                        self._update_dict(kw, ff)
        self.notempty = [pp for pp in [self.sec,
                                       self.cys,
                                       self.thr,
                                       self.arg,
                                       self.uga,
                                       self.ual,
                                       self.bsecis] if pp]

    def parseResultFiles(self, p2g=False, bsecisearch=False, force=False):
        '''Calls a parsing function on each file of a given result.

        bsecisearch is used on selenoproteins candidates.
        p2g files are parsed with an appropriate parser.
        '''
        if bsecisearch:
            for proteink, proteinv in self.sec.items():
                b_secisearcher = BSeciSearchWrapper('dull', 'dul', 1, 1)
                for hitk, hitv in proteinv.items():
                    outdir = '.'.join([proteink,
                                       str(hitk),
                                       'selenoproteine'
                                       'bsecis'])
                    os.mkdir(outdir)
                    cdsfile = [f for f in hitv if f.endswith('.cds')][0]
                    b_secisearcher.infile = cdsfile
                    b_secisearcher.outdir = outdir
                    b_secisearcher.run()
                    bs_output = os.path.join(outdir, 'bsecis_containing_sequences')
                    if not isFileEmpty(bs_output):
                        self._update_dict('bsecis', outdir)
                        self._update_dict('selenocysteine', outdir)
        if p2g:
            for case in self.notempty:
                for proteink, proteinv in case.items():
                    for hitk, hitv in protein.items():
                        p2gfile = [f for f in hitv if f.endswith('.p2g')][0]
                        p2g_parser = P2G_Parser(p2gfile)
                        p2g_parser.parse()

    def _keep(self, filename):
        trash = ['.ali', '.hit']
        for i in trash:
            if filename.endswith(i):
                return False
        return True

    def _get_hit_num(self, filename, case):
        spe_cases = ('std', 'non_std', 'twil')
        if case not in self.kw2dict.keys():
            raise 'Unknown case : ' + case
        ff = filename.split('.')
        if case not in spe_cases:
            index = ff.index(case) - 1
        else:
            index = ff.index(case) - 3
        return int(ff[index])

    def _get_prot_name(self, filename, case):
        spe_cases = ('std', 'non_std', 'twil')
        parts = filename.split('.')
        name = parts[0]
        if case not in spe_cases:
            index = parts.index(case) - 1
        else:
            index = parts.index(case) - 3
        return '.'.join(parts[:index])

    def _update_dict(self, keyword, filename):
        t_dict = {}
        t_dict.update(self.kw2dict[keyword])
        protname = self._get_prot_name(filename, keyword)
        hit_num = self._get_hit_num(filename, keyword)
        if protname in t_dict:
            if hit_num in t_dict[protname]:
                t_dict[protname][hit_num].append(filename)
            else:
                t_dict[protname][hit_num] = [filename]
        else:
            t_dict[protname] = {hit_num : [filename]}
        self.kw2dict[keyword].update(t_dict)


class ResultFolderParser(object):
    pass


def main():

    testpath = '/users/rg/agrimaldi/Data/gos/selenoprofiles_U_runs/2010_01_28_bact_U_failed/run/Zymomonas_mobilis_subsp._mobilis_ZM4_'

    parser = GenomeFolderParser(testpath)

    parser.parse(doall=True)
    print parser.cys
    pp = os.path.join(testpath, 'output', 'seld.1.unaligned.p2g')
    pparser = P2G_Parser(pp)
    pparser.parse()
    print pparser.result
    


if __name__ == '__main__':
    main()
