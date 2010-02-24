#!/usr/bin/env python

from __future__ import with_statement
import os
import sys
sys.path.append('/users/rg/agrimaldi/Code/python/python/lib')
sys.path.append('/users/rg/mmariotti/libraries')
import shutil
import subprocess
import optparse
import time
from AGBio.Utilities import *
from profiles_classes import parse_blast_output

def preformat(filename, pdir='/home/agrimaldi/temp/'):
    tmpfolder = genTempfilename(pdir)
    with open(tmpfolder, 'w') as tof:
        tof.write(subprocess.Popen(['alignthingie.pl',
                                    '-k',
                                    '-A',
                                    filename],
                                   stdout=subprocess.PIPE).communicate()[0])
    return tmpfolder

def gotoline(opened_file, line):
    c = 0
    while c < line:
        opened_file.readline()
        c += 1
    return opened_file

def get_full_hit(hit):
    full_hit = []
    with open(hit.ref.file, 'r') as iff:
        gotoline(iff, hit.ref.start_line)
        for i in range(hit.ref.end_line - hit.ref.start_line):
            full_hit.append(iff.readline())
    return full_hit

def stop_containing(hit):
    full_hit = get_full_hit(hit)
    for i in full_hit:
        if i.strip().startswith('Sbjct:') and '*' in i.split()[2]:
            return True
    return False

def main():

    try:
        time.sleep(0.1)
        infile = sys.argv[1]
        with open(sys.argv[1], 'r') as testif:
            pass
        prefmtfile = preformat(infile)
        blasthits = parse_blast_output(prefmtfile)
        filteredhits = []

        for hit in blasthits:
            if stop_containing(hit):
                filteredhits.append(hit)

        for hit in filteredhits:
            sys.stdout.write(''.join(get_full_hit(hit)))
            
    except IOError:
        sys.stderr.write('/nFile ' + infile + ' : not found./n/n')
    except Exception (e):
        print e
    finally:
        try:
            os.remove(prefmtfile)
        except Exception:
            pass

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.exit('man exit')
