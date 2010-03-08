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
from AGBio.selenoprofiles_tools.files_analysers.p2g import *
from AGBio.selenoprofiles_tools.files_analysers.b_secisearch import *
from AGBio.selenoprofiles_tools.files_analysers.ali import *


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
        self.excluded = []
        self.sec = {}
        self.cys = {}
        self.thr = {}
        self.arg = {}
        self.uga = {}
        self.ual = {}
        self.hom = {}
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
                        'bsecis' : self.secis_b,
                        'homologue' : self.hom}
        self.p2g = []
        self.ali = []

    def parse(self, doall=False,
              sec=False, cys=False, thr=False, arg=False,
              uga=False, ual=False, hom=False, bsecis=False,
              stdsecis=False, nonstdsecis=False, twilsecis=False):
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
        - `hom`: check for homologues candidates.
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
            if hom: lukw.append('homologue')
            if ual: lukw.append('unaligned')
            if stdsecis : lukw.append('std')
            if nonstdsecis : lukw.append('non_std')
            if twilsecis : lukw.append('twil')
            if bsecis : lukw.append('bsecis')
        for d in self.dirs:
            resfiles = [f for f in os.listdir(d) \
                        if self._keep(os.path.join(self.rootdir,
                                                   'output', f))]
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
                                       self.hom,
                                       self.secis_b] if pp]


    def parseResultFiles(self, p2g=False, bsecisearch=False, force=False):
        '''Calls a parsing function on each file of a given result.

        bsecisearch is used on selenoproteins candidates.
        p2g files are parsed with an appropriate parser.
        '''
        if bsecisearch:
            for proteink, proteinv in self.sec.items():
                b_secisearcher = BSeciSearchWrapper('dull', 'dul', 1, 1)
                for hitk, hitv in proteinv.items():
                    outdir = os.path.join(self.rootdir,
                                          'output',
                                          '.'.join([proteink,
                                                    str(hitk),
                                                    'selenocysteine',
                                                    'bsecis']))
                    try:
                        os.mkdir(outdir)
                    except OSError, (e):
                        if e.errno == 17:
                            pass
                    cdsfile = os.path.join(self.rootdir,
                                           'output',
                                           [f for f in hitv if f.endswith('.cds')][0])
                    b_secisearcher.infile = cdsfile
                    b_secisearcher.outdir = outdir
                    print b_secisearcher.cline
                    b_secisearcher.run()
                    bs_output = os.path.join(outdir, 'bsecis_containing_sequences')
                    if not isFileEmpty(bs_output):
                        self._update_dict('bsecis', outdir)
                        self._update_dict('selenocysteine', outdir)
                    else:
                        print 'Removing', outdir
                        shutil.rmtree(outdir)

    def parse_p2gs(self):
        for case in self.notempty:
            for proteink, proteinv in case.items():
                for hitk, hitv in proteinv.items():
                    p2gfile = os.path.join(self.rootdir,
                                           'output',
                                           [f for f in hitv \
                                            if f.endswith('.p2g')][0])
                    p2g_parser = P2G_Parser(p2gfile)
                    self.p2g.append(p2g_parser)
                    self.p2g[-1].parse()

    def ali_stats(self):
        for case in self.notempty:
            for protname in case.keys():
                alifile = os.path.join(self.rootdir,
                                       'output',
                                       '.'.join([protname, 'ali']))
                alignment = AliData(alifile)
                if alignment.filename not in [s.filename for s in self.ali]:
                    self.ali.append(alignment)
                    self.ali[-1].find_x_positions()

    def isexcluded(self, case):
        ccase = getattr(self, case)
        if self.excluded:
            for prot, hit in self.excluded:
                if not (prot in ccase and hit in ccase[prot]):
                    return False
            return True
        else:
            return False

    def _keep(self, filename):
        trash = ['.hit']
        for i in trash:
            if filename.endswith(i) \
                   or (filename.endswith('.ali') \
                       and os.path.getsize(filename) == 0):
                return False
        return True

    def _get_hit_info(self, filename, case):
        spe_case1 = ('std', 'non_std', 'twil')
        spe_case2 = ('bsecis')
        if case not in self.kw2dict.keys():
            raise 'Unknown case : ' + case
        parts = filename.split('.')
        index = None
        if case in spe_case1:
            index = parts.index(case) - 3
        elif case in spe_case2:
            index = parts.index(case) - 2
        else:
            index = parts.index(case) - 1
        return ['.'.join(parts[:index]),
                int(parts[index])]

    def _update_dict(self, keyword, filename):
        t_dict = {}
        t_dict.update(self.kw2dict[keyword])
        protname, hit_num = self._get_hit_info(filename, keyword)
        exclusionfile = filename.split('.')[-1] == 'exclude_from_tree'
        if protname in t_dict:
            if hit_num in t_dict[protname]:
                t_dict[protname][hit_num].append(filename)
            else:
                t_dict[protname][hit_num] = [filename]
        else:
            t_dict[protname] = {hit_num : [filename]}
        self.kw2dict[keyword].update(t_dict)
        if exclusionfile:
            self.excluded.append([protname, hit_num])


class ResultFolderParser(object):
    '''Parser for whole result folder, containing several genomes.
    '''
    def __init__(self, folder, lookat='output'):
        rootfolder = os.path.abspath(folder)
        self._genomes = [d for d in os.listdir(folder) if os.path.isdir(d)]


def main():

    testpath = '/users/rg/agrimaldi/Data/gos/selenoprofiles_U_runs/2010_01_28_bact_U_failed/run/Zymomonas_mobilis_subsp._mobilis_ZM4_'

    parser = GenomeFolderParser(testpath)

    parser.parse(doall=True)
    print parser.cys
    pp = os.path.join(testpath, 'output', 'seld.1.unaligned.p2g')
    pparser = P2G_Parser(pp)
    pparser.parse()
    print pparser.result
    pparser.result.target.fasta().prints()
    


if __name__ == '__main__':
    main()
