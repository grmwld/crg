#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import optparse
import shutil
import traceback
#sys.path.append('/soft/general/python-2.5.2/lib/python2.5/site-packages/MySQL_python-1.2.2-py2.5-linux-i686.egg')
#sys.path.append('/soft/general/python-2.5.2/lib/python2.5/site-packages/PyQt4')
#import PyQt4
#sys.path.append('/soft/general/python-2.5.2/lib/python2.5/site-packages')
sys.path.append('/users/rg/agrimaldi/usr/lib/python2.5/site-packages')
from ete2 import Tree, faces
from AGBio.selenoprofiles_tools.results_analyser import *
from AGBio.Utilities import *
from PyQt4 import QtGui


MAX_SN_LEN = 10
MAX_LN_LEN = 80

def load_equivalent_names(filename):
    output = {}
    with open(filename, 'r') as iff:
        for line in iff:
            sline = line.split()
            output[sline[0]] = sline[1]
    return output

def long2short(longname):
    ln = longname.split()
    return ln[0][0] + ln[1][:3]

def add_to_species_dict(longname, genomes):
    sln = sanitize(longname)
    if sln in genomes and sln not in species:
        species[sln] = sln

def sanitize(longname):
    t = longname
    for c in '\\\'\"/':
        t = t.replace(c, '')
    for c in '-_':
        t = t.rstrip(c)
    t = t.replace(' ', '_')
    return t + '_'

def parse_result_folders_file(ffile):
    info = []
    with open(ffile, 'r') as iff:
        for line in iff:
            ppath, names = line.split()
            pppath = os.path.abspath(os.path.expanduser(ppath))
            info.append({pppath:names.split(',')})
    return info

def getProts(ff):
    return ff[ff.keys()[0]]

def getFolder(ff):
    return ff.keys()[0]

def layout(node):
    if node.is_leaf():
        try:
            has_sec = False
            add_to_species_dict(node.name, g_genomes)
            node.img_style['size'] = 10
            species_separator = faces.TextFace(' ', ftype='courier', fsize=12)
            shortNameFace = faces.TextFace(long2short(node.name).ljust(MAX_SN_LEN),
                                           ftype='courier',
                                           fsize=12)
            pathNameFace = faces.TextFace(species[sanitize(node.name)],
                                          ftype='courier',
                                          fsize=12)
            longNameFace = faces.TextFace(node.name.ljust(MAX_LN_LEN),
                                          ftype='courier',
                                          fsize=12)

            for index, info_folder in enumerate(info):
                folder = getFolder(info_folder)
                protlist = getProts(info_folder)
                before = 0
                if index > 0:
                    before = len(getProts(info[index - 1]))
            
                fp = os.path.join(folder, species[sanitize(node.name)])
                sp_parser = GenomeFolderParser(fp)
                sp_parser.parse(sec=True, cys=True, thr=True, arg=True, bsecis=True)
                if bsecisearchoption:
                    sp_parser.parseResultFiles(p2g=False, bsecisearch=True)
                    sp_parser = GenomeFolderParser(fp)
                    sp_parser.parse(sec=True, cys=True, thr=True, arg=True, bsecis=True)

                protnames = set()

                for keyword in sp_parser.notempty:
                    protnames.update(keyword.keys())
                protnames = list(protnames)

                for col, protname in enumerate(protlist):
                    prot_count = 0
                    if protname in sp_parser.cys.keys() \
                           and not sp_parser.isexcluded('cys'):
                        prot_count += 1
                        faces.add_face_to_node(facecys, node,
                                               col + 2 + before,
                                               aligned=True)
                    if protname in sp_parser.sec.keys() \
                           and not sp_parser.isexcluded('sec'):
                        has_sec = True
                        prot_count += 1
                        if protname in sp_parser.secis_b.keys():
                            faces.add_face_to_node(facesec_b, node,
                                                   col + 2 + before,
                                                   aligned=True)
                        else:
                            faces.add_face_to_node(facesec, node,
                                                   col + 2 + before,
                                                   aligned=True)
                    if protname in sp_parser.thr.keys() \
                           and not sp_parser.isexcluded('thr'):
                        prot_count += 1
                        faces.add_face_to_node(facethr, node,
                                               col + 2 + before,
                                               aligned=True)
                    if protname in sp_parser.arg.keys() \
                           and not sp_parser.isexcluded('arg'):
                        prot_count += 1
                        faces.add_face_to_node(facearg, node,
                                               col + 2 + before,
                                               aligned=True)
                    if prot_count == 0:
                        faces.add_face_to_node(facenan, node,
                                               col + 2 + before,
                                               aligned=True)

                if has_sec:
                    node.img_style['fgcolor'] = '#75af51'
                    shortNameFace.fgcolor = QtGui.QColor('#479042')
    #                longNameFace.bgcolor = QtGui.QColor('#479042')
                else:
    #                node.img_style['fgcolor'] = '#af5b5b'
                    shortNameFace.bgcolor = QtGui.QColor('#9c3939')
                    longNameFace.bgcolor = QtGui.QColor('#9c3939')

            faces.add_face_to_node(shortNameFace, node, column=0, aligned=True)
            ##faces.add_face_to_node(pathNameFace, node, column=1, aligned=True)
            faces.add_face_to_node(longNameFace, node, column=1, aligned=True)
            ##faces.add_face_to_node(species_separator, node, column=1, aligned=True)
            
        except KeyError:
            print traceback.print_exc()
            node.delete()
        except OSError, e:
            if e.errno == 2:
                pass
        except Exception, e:
            print e
            print traceback.print_exc()

    else:
        node.img_style['size'] = 0
        

def main():

    parser = optparse.OptionParser()

    parser.add_option('-i', '--input_file',
                      dest='infilename',
                      help='input filename, newick tree.',
                      metavar='FILE')
    
    parser.add_option('-f', '--result_folders',
                      dest='result_folders',
                      help='file used to specify result folders and corresponding proteins. Each line should contain the path to a result folder and a coma separated list of proteins. The path and the list have to be separated by a space.',
                      metavar='FILE')

    parser.add_option('-g', '--gui',
                      action='store_true', dest='gui', default=False,
                      help='Use GUI for printing.')

    parser.add_option('-r', '--render',
                      action='store_true', dest='render', default=False,
                      help='render the tree in .png format.')

    parser.add_option('-b', '--bsecisearch',
                      action='store_true', dest='bsecisearch', default=False,
                      help='search for bSECIS elements before building the tree.')

    parser.add_option('-O', '--organism_names',
                      dest='org_names',
                      help='in case some organism names in the result folder do'+\
                      'not match the names in the newick tree, provide equivalences'+\
                      'in a file.',
                      metavar='FILE')
    
    parser.add_option('-R', '--resources',
                      dest='res_folder',
                      help='resources folder. (images)',
                      metavar='DIR')

    parser.set_defaults(res_folder = '.')

    (options, args) = parser.parse_args()

    global domain_of_life
    global g_genomes
    global species
    global resultfolders
    global info
    global bsecisearchoption
    global facecys, facesec, facesec_b, facearg, facethr, facenan
    facecys = faces.ImgFace(options.res_folder + 'cys.png')
    facesec = faces.ImgFace(options.res_folder + 'sec.png')
    facesec_b = faces.ImgFace(options.res_folder + 'sec-secis.png')
    facearg = faces.ImgFace(options.res_folder + 'arg.png')
    facethr = faces.ImgFace(options.res_folder + 'thr.png')
    facenan = faces.ImgFace(options.res_folder + 'nan.png')

    species = {}
    if options.org_names:
        species = load_equivalent_names(options.org_names)
    info = parse_result_folders_file(options.result_folders)
    print info
    resultfolders = flatten([i.keys() for i in info])
    print resultfolders
    g_genomes = list(set(flatten([os.listdir(folder) for folder in resultfolders])))
    print g_genomes

    bsecisearchoption = options.bsecisearch
    
    t = Tree(options.infilename)

    if options.gui:
        t.show(layout)
    else:
        print t

    if options.render:
        t.render('tree.png', layout)


if __name__ == '__main__':
    try:
        main()
    except Exception, (e):
        print e
        
