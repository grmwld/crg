#!/usr/bin/env python2.6
# -*- coding:utf-8 -*-

import os
import sys
import subprocess
sys.path.append('/users/rg/mmariotti/libraries')
from MMlib import bbash
from AGBio.UtilityWrappers import *
from AGBio.io.common import *


class P2G_Parser(object):
    '''Parse a p2g output file, using appropriate parsing method
    according to where it comes from (genewise or exonerate).
    '''
    def __init__(self, filename):
        self._filename = filename
        self.parse = self._source_prog()
        self.result = ParserResult()

    def _parse_genewise(self):
        #proc = subprocess.Popen(['parse_genewise.py', '-i', self._filename],
        #                        stdout=subprocess.PIPE,
        #                        shell=True)
        op = bbash(' '.join(['parse_genewise.py', '-i', self._filename]), True)
        #raw_input()
        #os.system(' '.join(['parse_genewise.py', '-i', self._filename]))
        output = [c.strip().split() for c in op.split(';')]
        self.result.query.parse(output[0], output[5][1])
        self.result.target.parse(output[1], output[5][2])
        self.result.frameshifts = output[2][1]
        self.result.orig_pos = output[3]
        self.result.score = output[4][1]
        if len(output[6]) == 1:
            self.result.stop_codons = 0
        else: self.result.stop_codons = output[6][1]

    def _parse_exonarate(self):
        pass

    def _source_prog(self):
        with open(self._filename, 'r') as iff:
            for line in iff:
                if line.strip().startswith('genewise output'):
                     return self._parse_genewise
        return self._parse_exonarate


class ParserResult(object):
    '''Class to hold results from parsing.
    '''
    def __init__(self):
        self.score = 0
        self.introns = 0
        self.frameshifts = 0
        self.query = self._result_seq()
        self.target = self._result_seq()
        self.orig_pos = []
        self.consensus = ''

    def __str__(self):
        return '\n'.join((str(self.query), str(self.target), self.score,
                          self.frameshifts, str(self.orig_pos)))

    class _result_seq(object):
        def __init__(self):
            self.name = ''
            self.start = 0
            self.end = 0
            self.sequence = ''
        def __str__(self):
            return '\n'.join((self.name,
                              self.start + ' --- ' + self.end,
                              formatText(self.sequence, 50, 4)))
        def parse(self, info, seq):
            self.name = info[0]
            self.start = info[1]
            self.end = info[2]
            self.sequence = seq
            

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
        self.sec = {}
        self.cys = {}
        self.thr = {}
        self.arg = {}
        self.uga = {}
        self.ual = {}
        self.secis_std = {}
        self.secis_nonstd = {}
        self.secis_twil = {}
        self.kw2dict = {'cysteine' : self.cys,
                        'selenocysteine' : self.sec,
                        'threonine' : self.thr,
                        'arginine' : self.arg,
                        'uga_containing' : self.uga,
                        'unaligned' : self.ual,
                        'std' : self.secis_std,
                        'non_std' : self.secis_nonstd,
                        'twil' : self.secis_twil}

    def parse(self, doall=False,
              sec=False, cys=False, thr=False, arg=False,
              uga=False, ual=False,
              stdsecis=False, nonstdsecis=False, twilsecis='False'):
        '''Method to parse the folder.

        The wanted informations are given as arguments. If `doall` is set to
        True, everything is parsed.

        Arguments:
        - `self`:
        - `doall`:
        - `sec`:
        - `cys`:
        - `thr`:
        - `arg`:
        - `uga`:
        - `ual`:
        - `stdsecis`:
        - `nonstdsecis`:
        - `twilsecis`:
        '''
        lukw = []
        if doall:
            lukw.append('selenocysteine')
            lukw.append('cysteine')
            lukw.append('threonine')
            lukw.append('arginine')
            lukw.append('uga_containing')
            lukw.append('unaligned')
            lukw.append('std')
            lukw.append('non_std')
            lukw.append('twil')
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
        for d in self.dirs:
            resfiles = [f for f in os.listdir(d) if self._keep(f)]
            for ff in resfiles:
                for kw in lukw:
                    if kw in ff.split('.'):
                        self._update_dict(kw, ff)

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

    def _update_dict(self, keyword, filename):
        t_dict = {}
        t_dict.update(self.kw2dict[keyword])
        protname = filename.split('.')[0]
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
