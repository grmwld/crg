#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement
import sys
import optparse

def load_nodes(nodes):
    data = {}
    c = None
    p = None
    with open(nodes, 'r') as iff:
        for line in iff:
            c, p = [i.strip() for i in line.split('|')[:2]]
            if p in data:
                data[p].append(c)
            else:
                data[p] = [c]
    return data

## def get_tree(cp_couples):
##     pc_dict = {}
##     for couple in cp_couples:
##         if couple[1] in pc_dict:
##             pc_dict[couple[1]].append(couple[0])
##         else:
##             pc_dict[couple[1]] = [couple[0]]
##     return pc_dict

def get_sub_taxids(tree, taxid, subs=[], only_leaves=False):
    lsubs = subs
    try:
        for c in tree[taxid]:
            lsubs.append(get_sub_taxids(tree, c, lsubs, only_leaves))
        if not only_leaves:
            lsubs.append(taxid)
    except KeyError, (e):
        return taxid
    for i in lsubs:
        if type(i) != type(''):
            lsubs.remove(i)
    return lsub


def main():

    parser = optparse.OptionParser()

    parser.add_option('-i', '--taxid',
                      dest='taxid',
                      help='taxid from wich all children nodes or parents nodes should be returned',
                      metavar='ID',
                      type='int')

    parser.add_option('-l', '--only_leaves',
                      dest='only_leaves', action='store_true', default=False,
                      help='only return those children nodes that are leaves')

    parser.add_option('-n', '--nodes_file',
                      dest='nodes_file',
                      help='file containing relations between nodes from NCBI',
                      metavar='FILE')

    parser.set_defaults(nodes_file = '/users/rg/mmariotti/libraries/nodes.dmp')

    (opt, args) = parser.parse_args()

    if not opt.taxid:
        parser.error('You must at least provide a taxid as an input')

    tree = load_nodes(opt.nodes_file)

    sub_nodes = get_sub_taxids(tree, str(opt.taxid),
                               only_leaves=opt.only_leaves)
    
    for i in sub_nodes:
        sys.stdout.write(i+'\n')


if __name__ == '__main__':
    main()
