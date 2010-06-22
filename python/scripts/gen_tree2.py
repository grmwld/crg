#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import os
import optparse
import traceback
import Image, ImageDraw, ImageFont, ImageOps
from ete2 import Tree, faces
sys.path.append('/users/rg/agrimaldi/Code/python/python/lib/')
print sys.path
from AGBio.Utilities import *


MAX_SN_LEN = 10
MAX_LN_LEN = 35

def prop(fcn):
    return property(**fcn())


class GlobalInfo(object):
    def __init__(self):
        self.__resource_folder = ''
        self.tmp_dir = ''
        self.families = []
        self.species = {}
        self.result_folders = []
        self.results_data = []
        self.run_info = []
        self.img = {'cysteine':[None, None], 'selenocysteine':[None, None],
               'arginine':[None, None], 'threonine':[None, None], 'unknown':[None, None],
               'nan':[None, None], 'uga_containing':[None, None], 'fam_name':[None, None]}

    @prop
    def resource_folder():
        def fget(self):
            return self.__resource_folder    
        def fset(self, value):
            self.__resource_folder = value
            for key in self.img.keys():
                self.img[key][0] = os.path.join(value, key + '.png')
                self.img[key][1] = {}
            self.img['nan'][1][0] = self.img['nan'][0]
        return locals()

    def parse_p2g_list(self, p2g_list):
        output = {}
        with open(p2g_list) as iff:
            for line in iff:
                i = line.strip()
                org = self.species[i.split('/')[-3]]
                prot, num, label = i.split('/')[-1].split('.')
                if org in output:
                    if prot in output[org]:
                        if label in output[org][prot]:
                            output[org][prot][label].append(i)
                        else: output[org][prot][label] = [i]
                    else: output[org][prot] = {label:[i]}
                else: output[org] = {prot:{label:[i]}}
        self.results_data = output

    def gen_png(self, kw, num, angle=0, save=True):
        if num not in self.img[kw][1].keys():
            #font = ImageFont.load("arial.pil", 12)
            im = Image.open(self.img[kw][0])
            draw = ImageDraw.Draw(im)
            draw.text((2, 4), str(num))
            tmp_png = genTempfilename(self.tmp_dir)
            im = im.rotate(angle)
            im.save(tmp_png, 'PNG')
            if save:
                self.img[kw][1][num] = tmp_png
        else:
            tmp_png = self.img[kw][1][num]
        return tmp_png
    
    def parse_families(self, ffile):
        with open(ffile, 'r') as iff:
            self.families = [f.strip() for f in iff.readlines()]

    def deduce_families(self):
        for org in self.results_data.keys():
            for fam in org.keys():
                if fam not in self.families:
                    self.families.add(fam)

def load_equivalent_names(filename):
    output = {}
    with open(filename, 'r') as iff:
        for line in iff:
            sline = line.strip().split(':::')
            output[sline[1]] = sline[0]
    return output

def sanitize(longname):
    t = longname
    for c in '\\\'\"/':
        t = t.replace(c, '')
    for c in '-_':
        t = t.rstrip(c)
    t = t.replace(' ', '_')
    return t

def layout(node):
    if node.is_leaf():
        try:
            sane_node_name = sanitize(node.name)
            node.img_style['size'] = 10
            longNameFace = faces.TextFace(node.name.ljust(MAX_LN_LEN),
                                          ftype='courier',
                                          fsize=12)
            if node.name == 'Family Name':
                for family in globalInfo.families:
                    tmp_png = globalInfo.gen_png('fam_name', family, angle=90, save=False)
                    faces.add_face_to_node(faces.ImgFace(tmp_png),
                                           node, globalInfo.families.index(family)+2,
                                           aligned = True)
            else:
                organism = globalInfo.results_data[sane_node_name]
                print organism
                for family in globalInfo.families:
                    if family in organism.keys():
                        num_other = 0
                        for label in organism[family]:
                            print 'candidate :', organism[family][label], globalInfo.families.index(family), len(organism[family][label])
                            if label in globalInfo.img.keys():
                                tmp_png = globalInfo.gen_png(label,
                                                             len(organism[family][label]))
                                faces.add_face_to_node(faces.ImgFace(tmp_png),
                                                       node, globalInfo.families.index(family)+2,
                                                       aligned = True)

                            else:
                                num_other += 1
                        if num_other:
                            tmp_png = globalInfo.gen_png('unknown', num_other)
                            faces.add_face_to_node(faces.ImgFace(tmp_png),
                                           node, globalInfo.families.index(family)+2,
                                           aligned = True)
                    else:
                        tmp_png = globalInfo.gen_png('nan', '')
                        faces.add_face_to_node(faces.ImgFace(tmp_png),
                                               node, globalInfo.families.index(family)+2,
                                               aligned = True)

            faces.add_face_to_node(longNameFace, node, column=1, aligned=True)
                
        except Exception, e:
            print e
            print traceback.print_exc()


def main():

    parser = optparse.OptionParser()

    parser.add_option('-i', '--input_file',
                      dest='infilename',
                      help='input filename, newick tree.',
                      metavar='FILE')
    
    parser.add_option('-f', '--families',
                      dest='families',
                      help='file used to specify the families.',
                      metavar='FILE')

    parser.add_option('-p', '--p2g_list',
                       dest='p2g_list',
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
                      dest='resource_folder',
                      help='resources folder. (images)',
                      metavar='DIR')

    parser.add_option('-t', '--temp',
                      dest='tempfolder',
                      help='temporary folder')

    parser.set_defaults(resource_folder = '.',
                        tempfolder = '/tmp' )


    (options, args) = parser.parse_args()

    global globalInfo

    globalInfo = GlobalInfo()

    globalInfo.resource_folder = options.resource_folder
    globalInfo.tmp_dir = options.tempfolder
    
    if options.org_names:
        globalInfo.species = load_equivalent_names(options.org_names)

    print globalInfo.species

    print 'Loading data ...'
    globalInfo.parse_p2g_list(options.p2g_list)
    print '... Done'


    if options.families:
        globalInfo.parse_families(options.families)
    else:
        globalInfo.deduce_families()

    print globalInfo.families

    t = Tree(options.infilename)


    if options.gui:
        t.show(layout)
    else:
        print t

    if options.render:
        t.render('tree.png', layout, w=2500, h=40*len(globalInfo.species))

    for k, v in globalInfo.img.items():
        if k != 'nan':
            for vv in v[1].values():
                os.remove(vv)


if __name__ == '__main__':
    try:
        main()
    except Exception, e:
        print e
        print traceback.print_exc()
