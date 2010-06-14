#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import optparse
import shutil
import traceback
import Image, ImageDraw, ImageFont
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

def parse_nrfile(nrfile):
    output = {}
    with open(nrfile) as iff:
        for line in iff:
            i = line.strip()
            org = i.split('/')[-3]
            prot, num, typ = i.split('/')[-1].split('.')
            if org in output:
                if prot in output[org]:
                    if typ in output[org][prot]:
                        output[org][prot][typ].append(i)
                    else: output[org][prot][typ] = [i]
                else: output[org][prot] = {typ:[i]}
            else: output[org] = {prot:{typ:[i]}}
    return output

def gen_png(base_im, num_candidates):
    im = Image.open(base_im)
    draw = ImageDraw.Draw(im)
    draw.text((8, 4), str(num_candidates))
    tmp_png = genTempfilename('/home/agrimaldi/temp')
    im.save(tmp_png, 'PNG')
    temp_bin.append(tmp_png)
    return tmp_png

def load_equivalent_names(filename):
    output = {}
    with open(filename, 'r') as iff:
        for line in iff:
            sline = line.strip().split(':::')
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

def get_num_nr(keyword, protname, orgname):
    if protname in nr_p2g[orgname] and keyword in nr_p2g[orgname][protname]:
        return len(nr_p2g[orgname][protname][keyword])
    return 0

def layout(node):
    if node.is_leaf():
        try:
            sane_node_name = sanitize(node.name)
            has_sec = False
            add_to_species_dict(node.name, g_genomes)
            node.img_style['size'] = 10
            species_separator = faces.TextFace(' ', ftype='courier', fsize=12)
            ## shortNameFace = faces.TextFace(long2short(node.name).ljust(MAX_SN_LEN),
##                                            ftype='courier',
##                                            fsize=12)
##             pathNameFace = faces.TextFace(species[sane_node_name],
##                                           ftype='courier',
##                                           fsize=12)
            longNameFace = faces.TextFace(node.name.ljust(MAX_LN_LEN),
                                          ftype='courier',
                                          fsize=12)

            for index, info_folder in enumerate(info):
                folder = getFolder(info_folder)
                protlist = getProts(info_folder)
                before = 0
                if index > 0:
                    before = len(getProts(info[index - 1]))
            
                fp = os.path.join(folder, species[sane_node_name])
                sp_parser = GenomeFolderParser(fp.strip())
                sp_parser.parse(sec=True, cys=True, thr=True, arg=True,
                                uga=True, ual=True, oth=True, bsecis=True)
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
                           and not sp_parser.isexcluded('cys') \
                           and (not nr_p2g \
                                or nr_p2g and get_num_nr('cysteine',
                                                      protname,
                                                      species[sane_node_name])):
                        prot_count = 1
                        n_cand = len(sp_parser.cys[protname])
                        if nr_p2g:
                            n_cand = get_num_nr('cysteine',
                                                protname,
                                                species[sane_node_name])
                        tmp_png = gen_png(imcys, n_cand)
                        faces.add_face_to_node(faces.ImgFace(tmp_png),
                                               node, col + 2 + before,
                                               aligned=True)
                    if protname in sp_parser.sec.keys() \
                           and not sp_parser.isexcluded('sec') \
                           and (not nr_p2g \
                                or nr_p2g and get_num_nr('selenocysteine',
                                                      protname,
                                                      species[sane_node_name])):
                        has_sec = True
                        prot_count = 1
                        n_cand = len(sp_parser.sec[protname])
                        if nr_p2g:
                            n_cand = get_num_nr('selenocysteine',
                                                protname,
                                                species[sane_node_name])
                        if protname in sp_parser.secis_b.keys():
                            faces.add_face_to_node(facesec_b, node,
                                                   col + 2 + before,
                                                   aligned=True)
                        else:
                            tmp_png = gen_png(imsec, n_cand)
                            faces.add_face_to_node(faces.ImgFace(tmp_png), node,
                                                   col + 2 + before,
                                                   aligned=True)
                    if protname in sp_parser.thr.keys() \
                           and not sp_parser.isexcluded('thr') \
                           and (not nr_p2g \
                                or nr_p2g and get_num_nr('threonine',
                                                      protname,
                                                      species[sane_node_name])):
                        prot_count = 1
                        n_cand = len(sp_parser.thr[protname])
                        if nr_p2g:
                            n_cand = get_num_nr('threonine',
                                                protname,
                                                species[sane_node_name])
                        tmp_png = gen_png(imthr, n_cand)
                        faces.add_face_to_node(faces.ImgFace(tmp_png), node,
                                               col + 2 + before,
                                               aligned=True)
                    if protname in sp_parser.arg.keys() \
                           and not sp_parser.isexcluded('arg') \
                           and (not nr_p2g \
                                or nr_p2g and get_num_nr('arginine',
                                                      protname,
                                                      species[sane_node_name])):
                        prot_count = 1
                        n_cand = len(sp_parser.arg[protname])
                        if nr_p2g:
                            n_cand = get_num_nr('arginine',
                                                protname,
                                                species[sane_node_name])
                        tmp_png = gen_png(imarg, n_cand)
                        faces.add_face_to_node(faces.ImgFace(tmp_png), node,
                                               col + 2 + before,
                                               aligned=True)
                    n_cand = 0
                    for kw in [('ual', 'unaligned'), ('uga', 'uga_containing'),
                               ('oth', 'other')]:
                        if protname in getattr(sp_parser, kw[0]).keys() \
                           and not sp_parser.isexcluded(kw[0]) \
                           and (not nr_p2g \
                                or nr_p2g and get_num_nr(kw[1],
                                                      protname,
                                                      species[sane_node_name])):
                            prot_count = 1
                            n_cand += len(getattr(sp_parser, kw[0])[protname])
                            if nr_p2g:
                                n_cand += get_num_nr(kw[1],
                                                     protname,
                                                     species[sane_node_name])
                    if n_cand:
                        tmp_png = gen_png(imoth, n_cand)
                        faces.add_face_to_node(faces.ImgFace(tmp_png), node,
                                               col + 2 + before,
                                               aligned=True)
                    if prot_count == 0:
                        faces.add_face_to_node(facenan, node,
                                               col + 2 + before,
                                               aligned=True)

                if has_sec:
                    node.img_style['fgcolor'] = '#75af51'
    #                shortNameFace.fgcolor = QtGui.QColor('#479042')
                    longNameFace.bgcolor = QtGui.QColor('#479042')
                else:
    #                node.img_style['fgcolor'] = '#af5b5b'
    #                shortNameFace.bgcolor = QtGui.QColor('#9c3939')
                    longNameFace.bgcolor = QtGui.QColor('#9c3939')

            ##faces.add_face_to_node(shortNameFace, node, column=0, aligned=True)
            ##faces.add_face_to_node(pathNameFace, node, column=1, aligned=True)
            faces.add_face_to_node(longNameFace, node, column=1, aligned=True)
            ##faces.add_face_to_node(species_separator, node, column=1, aligned=True)
            
        except KeyError:
            print traceback.print_exc()
            node.delete()
        except OSError, e:
            if e.errno == 2:
                print traceback.print_exc()
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

    parser.add_option('-c', '--nr_file',
                       dest='nr_file',
                       help='file containing all the non redundant p2g files. The file should be a list of p2g files, one per line')

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

    global temp_bin
    global nr_p2g
    global resource_folder
    global domain_of_life
    global g_genomes
    global species
    global resultfolders
    global info
    global bsecisearchoption
    global imcys, imsec, imsec_b, imarg, imthr, imoth
    global facecys, facesec, facesec_b, facearg, facethr, facenan
    resource_folder = options.res_folder
    temp_bin = []
    imcys = resource_folder + 'cys.png'
    imsec = resource_folder + 'sec.png'
    imsec_b = resource_folder + 'sec-secis.png'
    imarg = resource_folder + 'arg.png'
    imthr = resource_folder + 'thr.png'
    imoth = resource_folder + 'oth.png'
    facesec_b = faces.ImgFace(resource_folder + 'sec-secis.png')
    facenan = faces.ImgFace(resource_folder + 'nan.png')

    species = {}
    nr_p2g = None
    if options.nr_file:
        nr_p2g = parse_nrfile(options.nr_file)
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
        t.render('tree.png', layout, w=2500, h=40*len(g_genomes))

    for t in temp_bin:
        os.remove(t)

if __name__ == '__main__':
    try:
        main()
    except Exception, (e):
        print e
        
