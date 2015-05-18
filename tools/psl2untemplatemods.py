#!/usr/bin/env python3
#
# Copyright (c) 2014 Hyeshik Chang
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# - Hyeshik Chang <hyeshik@snu.ac.kr>
#

import pandas as pd
import gzip
from itertools import dropwhile
from Bio import SeqIO
import re

PSL_HEADER = [
    'match', 'mismatch', 'rep. match', 'Ns', 'Q gap count', 'Q gap bases', 'T gap count',
    'T gap bases', 'strand', 'Q name', 'Q size', 'Q start', 'Q end', 'T name',
    'T size', 'T start', 'T end', 'block count', 'blockSizes', 'qStarts', 'tStarts',
    'Q seq', 'T seq'
]

MULTIMAP_SORT_COLUMNS = [
    'match', 'mismatch', 'Q gap count', 'T gap count', 'Q gap bases', 'T gap bases'
]
MULTIMAP_SORT_ORDER = [ # ascending?
    False,   True,       True,          True,          True,          True,
]


def iterseq_sequential(seqin, format='fasta'):
    nextid = yield
    for seq in SeqIO.parse(seqin, format=format):
        while seq.name == nextid:
            nextid = yield seq

def rescue_untemplated_mods(row, pat_tailtrimmer=re.compile('[at]+$')):
    qseqend = row['Q seq'].rstrip(',').split(',')[-1]
    tseqend = row['T seq'].rstrip(',').split(',')[-1]

    if qseqend == tseqend:
        return row

    tailmatch = pat_tailtrimmer.search(qseqend)
    if tailmatch is None:
        return row

    scanstart = tailmatch.start()
    qseqend_toscan = qseqend[scanstart:]
    tseqend_toscan = tseqend[scanstart:]
    if qseqend_toscan == tseqend_toscan:
        return row

    firstmismatch = next(dropwhile(lambda x: x[1][0] == x[1][1],
                            enumerate(zip(qseqend_toscan, tseqend_toscan))))[0]
    additional_clipping = len(qseqend_toscan) - firstmismatch
    if additional_clipping <= 0:
        return row

    row['T end'] -= additional_clipping
    row['addedSeq'] = qseqend_toscan[-additional_clipping:].upper() + row['addedSeq']

    return row

def run(options):
    t_rescue = options.t_rescue
    refhairpins = SeqIO.to_dict(SeqIO.parse(options.reffile, 'fasta'))

    tbl = pd.read_table(options.alnfile, compression='gzip', names=PSL_HEADER)
    sortedmaps = tbl[tbl['strand'] == '+'].sort(['Q name'] + MULTIMAP_SORT_COLUMNS,
                                         ascending=[True] + MULTIMAP_SORT_ORDER)
    uniquemaps = sortedmaps.drop_duplicates('Q name').copy()
    uniquemaps['qSeqNumber'] = uniquemaps['Q name'].apply(lambda x: int(x.split('-')[0]))
    uniquemaps['readCount'] = uniquemaps['Q name'].apply(lambda x: int(x.split('-')[1]))
    uniquemaps = uniquemaps.sort('qSeqNumber', ascending=True)

    seqgetter = iterseq_sequential(gzip.open(options.readfile, 'rt'))
    next(seqgetter)
    uniquemaps['readSeq'] = uniquemaps.apply(lambda x:
                                             str(seqgetter.send(x['Q name']).seq), axis=1)
    uniquemaps['addedSeq'] = uniquemaps.apply(lambda x: x['readSeq'][x['Q end']:], axis=1)

    if t_rescue:
        uniquemaps = uniquemaps.apply(rescue_untemplated_mods, axis=1)

    finaltbl = uniquemaps[['T name', 'qSeqNumber', 'readCount', 'readSeq', 'T end', 'addedSeq']]
    finaltbl.columns = ['hairpin', 'seqNo', 'readCount', 'readSeq', 'refEndPos', 'addedSeq']

    finaltbl.to_csv(options.output if options.output != '-' else '/dev/stdout', sep='\t',
                    index=False)

def parse_arguments():
    import argparse

    parser = argparse.ArgumentParser(description='Parse BLAT alignments to produce a '
                                                 'simplified table with untemplated additions')
    parser.add_argument('--read', dest='readfile', metavar='FILE', type=str,
                        help='Sequence read file in gzipped-FASTA format', required=True)
    parser.add_argument('--reference', dest='reffile', metavar='FILE', type=str,
                        help='Reference hairpin sequence file in FASTA format', required=True)
    parser.add_argument('--alignments', dest='alnfile', metavar='FILE', type=str,
                        help='Sequence alignments in gzipped-PSLx format', required=True)
    parser.add_argument('--untmpl-t-rescue', dest='t_rescue', action='store_true',
                        default=False, help='Rescue aggressively untemplate T-stretches '
                        'with partial match to reference sequence.')
    parser.add_argument('--output', dest='output', metavar='FILE', type=str,
                        help='Output file for resulting table', required=True)
    return parser.parse_args()

if __name__ == '__main__':
    options = parse_arguments()
    run(options)
