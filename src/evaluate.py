#!/usr/bin/env python

import sys
import os
import argparse
import subprocess
import numpy as np

from itertools import islice

import ppv
import plot_contact_map
import a3m_to_trimmed
import parse_contacts


WORKDIR = os.path.dirname(os.path.realpath(__file__))


def extract_seq(aln_file):
    with open(aln_file) as f:
        seq = list(islice(f, N))
    return '\n'.join(seq)


def reformat_alignment(aln_file, aln_file_a3m, reformat_method=''):
    """ Optional STEP 0: Convert alignment to a3m """
    aln_type = aln_file.split('.')[-1]
    if aln_type == 'fa' or aln_type == 'fasta':
        aln_type = 'fas'
    if not reformat_method:
        cmd = ['%s/reformat.pl' % WORKDIR]
    else:
        cmd = [reformat_method]
    cmd += [aln_type, 'a3m', aln_file, aln_file_a3m]
    subprocess.call(cmd, stdout=open(os.devnull, 'wb'))


def trim_alignment(aln_file, trimmed_aln_file):
    """ STEP 1: convert a3m to trimmed format """
    outfile_generator = a3m_to_trimmed.convert(aln_file)
    with open(trimmed_aln_file, 'w') as outf:
        for l in outfile_generator:
            outf.write(l)


def predict_contacts(trimmed_aln_file, cm_file, cm_method=''):
    """ STEP 2: run given contact prediction method """
    # default:
    if not cm_method:
        cmd = ['%s/run_gdca.sh' % WORKDIR]
    else:
        cmd = [cm_method]
    cmd += [trimmed_aln_file, cm_file]
    subprocess.call(cmd)


def get_ppv(seq_file, cm_file, native_file):
    """ STEP 3a: compare contact map to native """
    result = ppv.get_ppv(seq_file, cm_file, native_file)
    #plot_contact_map.plot_map(seq_file, cm_file, pdb_filename=native_file)
    return result[1]


def get_numc(cm_file, th=0.):
    """ STEP 3b: get (normalized) number of contacts above score
        REMARK: in case of GaussDCA the score is a log-likelihood,
        => score above 0.0 = more likely to be correct than false
        Threshold needs to be changed depending on contact predictor.
    """
    contacts = parse_contacts.parse(open(cm_file))
    contacts_np = parse_contacts.get_numpy_cmap(contacts)
    numc = len(np.where(contacts_np > th)[0])
    numc_norm = numc / float(pow(contacts_np.shape[0], 2))
    return numc, numc_norm


def get_maxc(cm_file):
    """ STEP 3c: get maximal contact score """
    contacts = parse_contacts.parse(open(cm_file))
    contacts_np = parse_contacts.get_numpy_cmap(contacts)
    return np.max(contacts_np)

    
def evaluate(aln_file, seq_file='', native_file='', cm_method='', reformat_method='', th=0.):
    """ Run evaluation pipeline on given alignment """
    if not seq_file:
        seq_file = '%s/%s.fa' % (os.path.dirname(aln_file), os.path.basename(aln_file)[:5])
    if not native_file:
        native_file = '%s/native.pdb' % os.path.dirname(aln_file)
    if not os.path.isfile(seq_file):
        sys.exit('Please provide an existing sequence file.')
    
    # STEP 0
    if not aln_file.endswith('.a3m'):
        aln_file_a3m = '.'.join(aln_file.split('.')[:-1]) + '.a3m'
        if not os.path.isfile(aln_file_a3m) or os.stat(aln_file_a3m).st_size == 0:
            reformat_alignment(aln_file, aln_file_a3m, reformat_method)
        aln_file = aln_file_a3m

    # STEP 1
    trimmed_aln_file = '.'.join(aln_file.split('.')[:-1]) + '.trimmed'
    if not os.path.isfile(trimmed_aln_file):
        trim_alignment(aln_file, trimmed_aln_file)

    # STEP 2
    cm_file = '.'.join(trimmed_aln_file.split('.')[:-1]) + '.cm'
    if not os.path.isfile(cm_file) or os.stat(cm_file).st_size == 0:
        predict_contacts(trimmed_aln_file, cm_file, cm_method=cm_method)
    
    # STEP 3
    ppv = get_ppv(seq_file, cm_file, native_file)
    numc, numc_norm = get_numc(cm_file, th=th)
    maxc = get_maxc(cm_file)

    return [ppv, numc, numc_norm, maxc]


if __name__ == '__main__':
    
    p = argparse.ArgumentParser(description='Run alignment quality\
            evaluation workflow.\nFor given alignment it outputs PPV,\
            (normalized) number of contacts above score threshold, and\
            maximum contact score.')
    p.add_argument('alignment', help='Input aligment file')
    p.add_argument('-s', '--seqfile', default='', help='Sequence file')
    p.add_argument('-n', '--native', default='', help='Reference pdb file to compare with')
    p.add_argument('-c', '--contact', default='', help='Path to contact predictor executable')
    p.add_argument('-t', '--threshold', default=0., type=float, help='Contact score threshold')
    p.add_argument('-r', '--reformat', default='', help='Path to reformat.pl script from HHsuite')
    p.add_argument('-o', '--output', default='', help='Save output in csv format')

    args = vars(p.parse_args(sys.argv[1:]))
    stats = evaluate(args['alignment'], seq_file=args['seqfile'],\
            native_file=args['native'], cm_method=args['contact'],\
            reformat_method=args['reformat'])

    out_file = args['output']
    if out_file:
        with open(out_file, 'w') as outf:
            outf.write('alignment_file,PPV,numc,maxc')
            outf.write('%s,%s\n' % (args['alignment'], ','.join(map(str, stats))))
    else:
        print '%s,%s' % (args['alignment'], ','.join(map(str, stats)))

