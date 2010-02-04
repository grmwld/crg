#!/usr/bin/env python2.6
# -*- coding:utf-8 -*-

import os
import sys


class GenomeFolderParser(object):

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
    print parser.secis_std


if __name__ == '__main__':
    main()
