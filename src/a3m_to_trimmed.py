#!/usr/bin/env python
import sys

def convert(infile):
    with open(infile) as aln:
        counter = 0
        for l in aln:
            if '>' in l and not counter == 0:
                yield '\n>sequence{0:07d}/1-100\n'.format(counter)
                counter += 1
            elif '>' not in l:
                l = l.strip()
                upperseq = ''.join([c for c in l if not c.islower()])
                upperseq = upperseq.replace('X', '-')
                yield upperseq
            elif '>' in l and counter == 0:
                yield '>target/1-100\n'
                counter += 1
        yield '\n'


if __name__ == '__main__':

    infile = sys.argv[1]
    outfile_generator = convert(infile)
    for l in outfile_generator:
        sys.stdout.write(l)
