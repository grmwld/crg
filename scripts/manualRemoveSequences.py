#!/usr/bin/env python2.6

import os
import sys
import optparse
import AGBio.io.Fasta as Fasta
from AGBio.Utilities import getch

def main():

    parser = optparse.OptionParser()

    parser.add_option( '-i', '--inputfile',
                       dest='inputfilename',
                       help='fasta file in which selenoproteins should be looked for.',
                       metavar='FILE' )

    parser.add_option( '-a', '--alignmentfile',
                       dest='alfilename',
                       help='alignment file used when details are requested.',
                       metavar='FILE' )

    parser.add_option( '-o', '--outputfile',
                       dest='outputfilename',
                       help='base output filename',
                       metavar='FILE' )

    parser.add_option( '-f', '--keep_prefilter',
                       dest='keepprefilter',
                       help='prefilters all sequences that have the given pattern in their name and keep them.',
                       metavar='PATTERN' )
    
    parser.add_option( '-F', '--throw_prefilter',
                       dest='throwprefilter',
                       help='prefilters all sequences that have the given pattern in their name and throw them.',
                       metavar='PATTERN' )

    parser.add_option( '-b', '--autothrow_abscents',
                       action='store_true', dest='atabscent', default=False,
                       help='Throw all sequences not present in the alignment provided.')

    parser.set_defaults( keepprefilter = False,
                         throwprefilter = False,
                         alfilename = False )

    (options, args) = parser.parse_args()

    with open(options.inputfilename, 'r') as inf:
        sequences = Fasta.loadSequences(inf)

    if options.alfilename:
        with open(options.alfilename, 'r') as alf:
            alignment = Fasta.Alignment(Fasta.loadSequences(alf))
        nrdetail = alignment.findPositions(('U','C','-'), False)
        rdetail = alignment.findPositions(('U','C','-'), True)
        
    if options.keepprefilter:
        kpatterns = options.keepprefilter.split(',')
    if options.throwprefilter:
        tpatterns = options.throwprefilter.split(',')

    kept_seq = Fasta.SequenceList()
    thrown_seq = Fasta.SequenceList()
    man_check_list = Fasta.SequenceList()

    for seq in sequences:
        kept = False
        thrown = False
        if options.keepprefilter:
            for pattern in kpatterns:
                if pattern in seq.header:
                    kept_seq.append(seq)
                    kept = True
        if options.throwprefilter:
            for pattern in tpatterns:
                if pattern in seq.header:
                    thrown_seq.append(seq)
                    thrown = True
        if not kept and not thrown:
            man_check_list.append(seq)

    idx = 0
    while idx < len(man_check_list):
        seq = man_check_list[idx]
        gi = seq.header.split('|')[1]
        choice = 'r'
        decided = False
        print seq.header
        while not decided:
            print len(kept_seq), len(thrown_seq)
            choice = getch('# '+str(idx+1)+' / '+str(len(man_check_list))+' -- Keep ? [Y/n]')
            if choice == 'b':
                if idx > 0:
                    idx -= 1
                    seq = man_check_list[idx]
                    gi = seq.header.split('|')[1]
                    print seq.header
                    try:
                        thrown_seq.remove(seq)
                    except:
                        pass
                    try:
                        kept_seq.remove(seq)
                    except:
                        pass
            elif choice in ('y', '\n'):
                kept_seq.append(seq)
                decided = True
                idx += 1
            elif choice == 'n':
                thrown_seq.append(seq)
                decided = True
                idx += 1
            elif choice == 's':
                os.system('fetch_seq.g -v TITLE="'+gi+'" -v ALL=1 '+options.inputfilename )
            elif choice == 'd' and options.alfilename:
                print
                print 'General Detail :'
                for pos in rdetail['U']:
                    sys.stdout.write('    '+str(pos) + ' ')
                    for xpos in rdetail:
                        try:
                            sys.stdout.write(str(xpos)+': ')
                            sys.stdout.write(str(len(rdetail[xpos][pos])) + ' ; ')
                        except KeyError:
                            sys.stdout.write('0 ; ')
                    sys.stdout.write('\n')
                print
                tmpseq = None
                for seqal in alignment:
                    if seqal.header == seq.header:
                        tmpseq = seqal
                if tmpseq:
                    tmppos = [i for i, x in enumerate(tmpseq.sequence) if x == 'U']
                    print 'In the sequence provided :'
                    print '    U :', tmppos
                    print '    U in those positions :', [len(rdetail['U'][(l,)]) for l in tmppos]
                    print '    C in those positions:', [len(rdetail['C'][(l,)]) for l in tmppos]
                    print '    - in those positions:', [len(rdetail['-'][tuple((l,))]) for l in tmppos]
                    print
                    print '    Symbols present at the positions of each U :'
                    for pos in [p for p in rdetail['U'] if p != ()]:
                        spos = str(pos[0])
                        print '        Position :', spos, '---', tmpseq.sequence[int(spos)]
                else:
                    print 'Not present in the alignment provided'
                print
            elif choice == 'q':
                cc = 'r'
                while cc not in ('y', 'n'):
                    cc = raw_input('Manual quit. Would you like to save your changes ? [y/N]')
                    if cc in 'y':
                        pass
                    if cc in 'n':
                        sys.exit('Quiting without saving.')
            else:
                print 'Wrong command'

    with open(options.outputfilename, 'w') as of:
        kept_seq.prints(of, 80)


if __name__ == '__main__':

    main()
