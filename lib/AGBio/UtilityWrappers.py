#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import os

TEMP = '/home/agrimaldi/temp/'
SBIN = '/soft/bin/'
HBIN = '/users/rg/agrimaldi/usr/bin/'
MBIN = '/users/rg/mmariotti/bin/'
MSCRIPT = '/users/rg/mmariotti/scripts/'
BLAST = SBIN + 'blastall'
MAFFT = SBIN + 'mafft'
TRIMAL = HBIN + 'trimal_dev'
TCOFFEE = HBIN + 't_coffee'
ADDFULLHEADERS = 'add_detail_to_titles2.py'
FILTER = MBIN + 'fetch_seq.g'
PRESELENOPROFILES = HBIN + 'prepare_alignment_selenoprofiles.py'
SELENOPROFILES = HBIN + 'selenoprofiles.py'
BLASTPARSER = 'blaster_parser.g'
DOWNLOADER = HBIN + 'get_sequence_from_gi.pl'

NOOPT = ['', '']


class BaseUtilityWrapper(dict):
    '''Base class for wrappers.

    Wrappers ease the use of common bioinformatic external tools.
    '''
    def __init__(self, utilitypath):
        '''Constructor.

        @param utilitypath : The absolute path of the tool.
        '''
        self.utilitypath = utilitypath
        self.cline = ''

    def getClineOpt(self, option):
        '''Getter of commandline options.

        @param option : The option to get
        @return : The relevent commandline option.
        '''
        return ''.join(self[option])

    def updateCline(self):
        '''Updates the commandline based on the options.
        '''
        self.cline = self.utilitypath
        for opt in self.keys():
            self.cline += (' ' + self.getClineOpt(opt))

    def run(self):
        '''Runs the commandline built with the various options.
        '''
        os.system(self.cline)


class UtilityWrapper(dict):
    '''Base class for wrappers.

    Wrappers ease the use of common bioinformatic external tools.
    '''
    def __init__(self, utilitypath, infile, outfile):
        '''Constructor.

        @param utilitypath : The absolute path of the tool.
        @param infile : The input file.
        @param outfile : The output file.
        '''
        self.utilitypath = utilitypath
        self.update({
            'infile' : ['', infile],
            'outfile' : ['', outfile]
            })
        self.cline = ''
##        self.updateCline()

    @property
    def infile(self):
        return self.getClineOpt('infile')
    @infile.setter
    def infile(self, infile):
        self['infile'][1] = infile
        self.updateCline()

    @property
    def outfile(self):
        return self.getClineOpt('outfile')    
    @outfile.setter    
    def outfile(self, outfile):
        self['outfile'][1] = outfile
        self.updateCline()
        
    def getClineOpt(self, option):
        '''Getter of commandline options.

        @param option : The option to get
        @return : The relevent commandline option.
        '''
        return ''.join(self[option])

    def updateCline(self):
        '''Updates the commandline based on the options.
        '''
        self.cline = self.utilitypath
        for opt in self.keys():
            self.cline += (' ' + self.getClineOpt(opt))

    def run(self):
        '''Runs the commandline built with the various options.
        '''
        os.system(self.cline)

class HasOption(dict):
    '''Base interface for providing option specification
    '''
    def __init__(self, keyarg='', carg='', value=None):
        self._keyarg = keyarg
        self._carg = carg
        self._value = value
        

    def _get_value(self):
        return self.getClineOpt(getattr(self._keyarg, self))
    def _set_value(self):
        if value:
            self[self._keyarg] = [self._carg, value]
        elif self._keyarg in self.keys():
            del self[self._keyarg]
        self.updateCline()


class HasInFile(dict):
    '''Interface providing infile input.
    '''
    def __init__(self, carg='', infile=None):
        '''Interface providing infile specification
        '''
        self._icarg = carg
        self.infile = infile

    @property
    def infile(self):
        return self.getClineOpt('infile')
    @infile.setter
    def infile(self, value):
        if value:
            self['infile'] = [self._icarg, value]
        elif 'infile' in self.keys():
            del self['infile']
        self.updateCline()


class HasOutFile(dict):
    '''Interface providing outfile specification.
    '''
    def __init__(self, carg='', outfile=None):
        '''
        '''
        self._ocarg = carg
        self.outfile = outfile
        
    @property
    def outfile(self):
        return self.getClineOpt('outfile')    
    @outfile.setter    
    def outfile(self, value):
        if value:
            self['outfile'] = [self._ocarg, value]
        elif 'outfile' in self.keys():
            del self['outfile']
        self.updateCline()


class MafftWrapper(UtilityWrapper):
    '''Wrapper for thr tool Mafft.
    '''
    def __init__(self, infile, outfile, auto=True):
        UtilityWrapper.__init__(self, MAFFT, infile, outfile)
        self['outfile'] = ['> ', outfile]
        if auto:
            self['auto'] = '--auto'
        self.updateCline()
    
    @property
    def auto(self):
        return self.getClineOpt('auto')
    @auto.setter
    def auto(self, auto):
        if auto:
            self['auto'] = '--auto'
        elif 'auto' in self.keys():
            del self['auto']
        self.updateCline()
        
    def updateCline(self):
        '''Updates the commandline based on the options.

        Overloaded in order to keep the order of the various parameters
        consistent : mafft options infile > outile
        '''
        self.cline = self.utilitypath
        for opt in self.keys():
            if opt != 'infile' and opt != 'outfile':
                self.cline += (' ' + self.getClineOpt(opt))
        self.cline += ' ' + (' '.join((self.infile, self.outfile)))


class TcoffeeWrapper(UtilityWrapper):
    '''Wrapper for the tool t_coffee.
    '''
    def __init__(self, infile, outfile=None, output='fasta', ncore=1):
        UtilityWrapper.__init__(self, TCOFFEE, infile, outfile)
        if outfile:
            self['outfile'] = ['-outfile=', outfile]
        self['output'] = ['-output=', output]
        self['ncore'] = ['-n_core=', str(ncore)]
        self.updateCline()

    @property
    def outfile(self):
        return UtilityWrapper.infile.fget(self)    
    @outfile.setter    
    def outfile(self, outfile):
        if outfile:
            self['outfile'] = ['-outfile=', outfile]
        elif 'outfile' in self.keys():
            del self['outfile']
        self.updateCline()

    @property
    def ncore(self):
        return self.getClineOpt('ncore')
    @ncore.setter
    def ncore(self, ncore):
        self['ncore'][1] = str(ncore)
        self.updateCline()

    def updateCline(self):
        '''Updates the commandline based on the options.
        Makes sure that the infile option comes first.
        '''
        self.cline = ' '.join((self.utilitypath, self.infile))
        for opt in self.keys():
            if opt != 'infile':
                self.cline += (' ' + self.getClineOpt(opt))
                

class TrimalWrapper(UtilityWrapper):
    '''Wrapper for the tool trimal.
    '''
    def __init__(self, infile, outfile, clusters=None, gapthreshold=None, scoreoverlap=None):
        UtilityWrapper.__init__(self, TRIMAL, infile, outfile) 
        self['infile'] = ['-in ', infile]
        self['outfile'] = ['-out ', outfile]
        if clusters:
            self['clusters'] = ['-clusters ', str(clusters)]
        if gapthreshold:
            self['gapthreshold'] = ['-gt ', str(gapthreshold)]
        if scoreoverlap:
            self['scoreoverlap'] = ['-resoverlap ', str(scoreoverlap[0]),
                                    ' -seqoverlap ', str(scoreoverlap[1])]
        self.updateCline()

    @property
    def clusters(self):
        return self.getClineOpt('clusters')
    @clusters.setter
    def clusters(self, clusters):
        if clusters:
            self['cluster'] = ['-clusters ', str(clusters)]
        elif 'clusters' in self.keys():
            del self['clusters']
        self.updateCline()

    @property
    def gapthreshold(self):
        return self.getClineOpt('gapthreshold')                
    @gapthreshold.setter
    def gapthreshold(self, gapthreshold):
        if gapthreshold:
            self['gapthreshold'] = ['-gapthreshold ', str(gapthreshold)]
        elif 'gapthreshold' in self.keys():
            del self['gapthreshold']
        self.updateCline()

    @property
    def scoreoverlap(self):
        return self.getClineOpt('scoreoverlap')                
    @gapthreshold.setter
    def scoreoverlap(self, value):
        if value:
            self['scoreoverlap'] = ['-resoverlap ', str(value[0]),
                                    ' -seqoverlap ', str(value[1])]
        elif 'scoreoverlap' in self.keys():
            del self['scoreoverlap']
        self.updateCline()


class AddFullHeadersWrapper(UtilityWrapper):
    '''Wrapper to use the add_detail_titles2.py utility,
    with fetches the full header of a sequence based on the gi
    '''
    def __init__(self, infile, outfile, patternfile):
        UtilityWrapper.__init__(self, ADDFULLHEADERS, infile, outfile)
        self['infile'] = ['-i ', infile]
        self['outfile'] = ['-o ', outfile]
        self['patternfile'] = ['-p ', patternfile]
        self.updateCline()

    @property
    def patternfile(self):
        return self.getClineOpt('patternfile')
    @patternfile.setter
    def patternfile(self, patternfile):
        self['patternfile'] = ['-p ', patternfile]
        self.updateCline()

    def updateCline(self):
        '''Updates the commandline based on the options.
        Makes sure that the outfile option comes last.
        '''
        self.cline = self.utilitypath
        for opt in self.keys():
            if opt != 'outfile':
                self.cline += (' ' + self.getClineOpt(opt))
        self.cline += ' ' + self.outfile


class AddFullHeadersWrapper2(BaseUtilityWrapper, HasInFile, HasOutFile):
    '''Wrapper to use the add_detail_titles3.py utility,
    with fetches the full header of a sequence based on the gi
    '''
    def __init__(self, infile, outfile, patternfile):
        BaseUtilityWrapper.__init__(self, 'add_detail_to_titles3.py')
        HasInFile.__init__(self, '-i ', infile)
        HasOutFile.__init__(self, '-o ', outfile)
        self['patternfile'] = ['-p ', patternfile]
        self.updateCline()

    @property
    def patternfile(self):
        return self.getClineOpt('patternfile')
    @patternfile.setter
    def patternfile(self, patternfile):
        self['patternfile'] = ['-p ', patternfile]
        self.updateCline()
        

class FilterWrapper(UtilityWrapper):
    '''Class to use the fetch_seq.g script
    '''
    def __init__(self, infile, outfile, titlematch, inverse=False):
        UtilityWrapper.__init__(self, FILTER, infile, outfile)
        self['outfile'] = ['> ', outfile]
        self['titlematch'] = ['-v TITLE=', '\"' + ';;;'.join(titlematch) + '\"']
        if inverse:
            self['match'] = ['-v V=1']
        else:
            self['match'] = ['-v ALL=1']
        self.updateCline()

    @property
    def titlematch(self):
        return self.getClineOpt('titlematch')
    @titlematch.setter
    def titlematch(self, titlematch):
        if titlematch:
            self['titlematch'] = ['-v TITLE=', '\"' + ';;;'.join(titlematch) + '\"']
        elif 'titlematch' in self.keys():
            del self['titlematch']
        self.updateCline()

    @property
    def inverse(self):
        return self.getClineOpt('inverse')
    @inverse.setter
    def titlematch(self, inverse):
        if inverse:
            self['match'] = ['-v V=1']
        else:
            self['match'] = ['-v ALL=1']
        self.updateCline()

    def updateCline(self):
        '''Updates the commandline based on the options.
        Makes sure that the infile and outfile options are
        the last.
        '''
        self.cline = self.utilitypath
        for opt in self.keys():
            if opt != 'infile' and opt != 'outfile':
                self.cline += (' ' + self.getClineOpt(opt))
        self.cline += ' ' + (' '.join((self.infile, self.outfile)))


class SelenoprofilesPreWrapper(UtilityWrapper):
    '''Wrapper for the tool prepare_alignment_selenoprofiles.py
    '''
    def __init__(self, infile, outfile, temp=None, all=True, tagthreshold=None):
        UtilityWrapper.__init__(self, PRESELENOPROFILES, infile, outfile)
        self['outfile'] = ['-output ', outfile]
        if temp:
            self['temp'] = ['-temp ', temp]
        if all:
            self['all'] = ['-all']
        if tagthreshold:
            self['tagthreshold'] = ['-tag_threshold ', str(tagthreshold)]
        self.updateCline()

    @property
    def temp(self):
        return self.getClineOpt('temp')
    @temp.setter
    def temp(self, temp):
        if temp:
            self['temp'] = ['-temp ', temp]
        elif 'temp' in self.keys():
            del self['temp']
        self.updateCline()

    @property
    def all(self):
        return self.getClineOpt('all')
    @all.setter
    def all(self, all):
        if all:
            self['all'] = '-all'
        elif 'all' in self.keys():
            del self['all']
        self.updateCline()

    @property
    def tagthreshold(self):
        return self.getClineOpt('tagthreshold')
    @tagthreshold.setter
    def tagthreshold(self, value):
        if value:
            self['tagthreshold'] = ['-tagthreshold ', str(value)]
        elif 'tagthreshold' in self.keys():
            del self['tagthreshold']
        self.updateCline()
        
    def updateCline(self):
        '''Updates the commandline based on the options.
        Makes sure that the infile option comes first.
        '''
        self.cline = ' '.join((self.utilitypath, self.infile, '-f'))
        for opt in self.keys():
            if opt != 'infile':
                self.cline += (' ' + self.getClineOpt(opt))


class BlasterParserWrapper(UtilityWrapper):
    '''Wrapper to use blaster_parser.g
    '''
    def __init__(self, infile, outfile, pgi=True, pevalue=True, evalue=10, pattern=None, sort=False):
        UtilityWrapper.__init__(self, BLASTPARSER, infile, outfile)
        self['outfile'] = ['> ', outfile]
        if evalue and pattern:
            self['gfilter'] = ['| ', "gawk '{if( $(NF-6)~/"+pattern+"/ || $(NF-4) <= "+evalue+")"]
        elif evalue and not pattern:
            self['gfilter'] = ['| ', "gawk '{if( $(NF-4) <= "+evalue+")"]
        if pgi and pevalue:
            self['output'] = ['', "{print \">\" $1, $(NF-4)}}'"]
        elif pgi:
            self['output'] = ['', "{print \">\" $1}}'"]
        elif pevalue:
            self['output'] = ['', "{print $(NF-4)}}'"]
        if sort:
            self['sort'] = ['| ', "sort"]
        else:
            self['sort'] = NOOPT
        self.updateCline()

    @property
    def gfilter(self):
        return self.getClineOpt('gfilter')

    @property
    def output(self):
        return self.getClineOpt('output')
    
    @property
    def sort(self):
        return self.getClineOpt('sort')
    @sort.setter
    def sort(self, sort):
        if sort:
            self['sort'] = ['| ', "sort"]
        else:
            self['sort'] = NOOPT

    def updateCline(self):
        '''Updates the command line making sure that
        option is in the correct order
        '''
        self.cline = ' '.join((self.utilitypath, '-v ALL=1',
                               self.infile,
                               self.gfilter,
                               self.output,
                               self.sort,
                               self.outfile))


class DownloaderWrapper(UtilityWrapper):
    '''Wrapper to use get_sequence_from_gi.pl
    '''
    def __init__(self, infile, outfile):
        UtilityWrapper.__init__(self, DOWNLOADER, infile, outfile)
        self['outfile'] = ['>> ', outfile]

    def updateCline(self):
        self.cline = ' '.join((self.utilitypath, self.infile, self.outfile))
        
    
class SelenoprofilesWrapper(UtilityWrapper):
    pass


class Blastall(UtilityWrapper):
    pass
