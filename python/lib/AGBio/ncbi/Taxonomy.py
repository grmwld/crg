#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement


def load_taxid_gi_map(filename):
    taxid_gi_dict = {}
    with open(filename, 'r') as iff:
        for line in iff:
            gi, taxid = [l.strip() for l in line.split()]
            if taxid in taxid_gi_dict:
                taxid_gi_dict[taxid].append(gi)
            else:
                taxid_gi_dict[taxid] = [gi]
    return taxid_gi_dict


class Node(object):
    def __init__(self, taxid=None, parent=None, children=None):
        self.taxid = taxid
        self.children = children
        self.parent = parent

    def get_sub_children(self, tree, taxid):
        lsubs = subs
        try:
            for c in tree[taxid]:
                lsubs.append(get_sub_taxids(tree, c, lsubs))
            lsubs.append(taxid)
        except KeyError, (e):
            return taxid
        for i in lsubs:
            if type(i) != type(''):
                lsubs.remove(i)
        return lsubs

    def add_child(self, node):
        self.children.append(node)
        

class Tree(dict):
    def __init__(self, nodesfile=None):
        self.filename = nodesfile
        self.root = Node()

    def load_nodes(self):
        c = None
        p = None
        t_dict = {}
        with open(self.filename, 'r') as iff:
            for line in iff:
                c, p = [i.strip() for i in line.split('|')[:2]]
                if p in self:
                    t_dict[p].append(c)
                else:
                    t_dict[p] = [c]
        self.root = self.__populate_tree(t_dict, Node('1', None, None))

    def __populate_tree(self, t_dict, node):
        curnode = node
        print 'curnode', node.taxid
        try:
            print t_dict[curnode.taxid]
            for c in t_dict[curnode.taxid]:
                newnode = Node(c, curnode, None)
                print 'newnode', newnode.taxid
                curnode.add_child(self.__populate_tree(t_dict,
                                                       node=newnode))
            self[curnode.taxid] = curnode
        except KeyError:
            return curnode
        return curnode
            


def get_sub_taxids(tree, taxid, subs=[]):
    lsubs = subs
    try:
        for c in tree[taxid]:
            lsubs.append(get_sub_taxids(tree, c, lsubs))
        lsubs.append(taxid)
    except KeyError, (e):
        return taxid
    for i in lsubs:
        if type(i) != type(''):
            lsubs.remove(i)
    return lsubs


def main():

    nodesdmp = '/users/rg/mmariotti/libraries/nodes.dmp'

    tt = Tree(nodesdmp)
    tt.load_nodes()

    print tt.keys()



if __name__ == '__main__':
    main()
