#!/usr/bin/env python3
"""
Filters and aligns fasta records based on cmsearch covariance model and
presence of invalid sequence characters.
"""
import argparse
import sys

from Bio import SeqIO, Alphabet


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument(
        'cmsearch',
        help='cmsearch alignments')
    p.add_argument(
        'fasta',
        help='cmsearch fasta file')
    p.add_argument(
        '--unknowns',
        required=True,  # mirrors ``taxit update_taxids``
        metavar='fasta',
        help=('fasta format output of sequences not '
              'aligned or with invalid sequence characters'))
    p.add_argument(
        '--out',
        default=sys.stdout,
        metavar='fasta',
        help='fasta output of sequences in forward orientation')

    args = p.parse_args()

    cmsearch = (row for row in open(args.cmsearch) if not row.startswith('#'))
    cmsearch = (row.split() for row in cmsearch)
    cmsearch = {row[0]: row[9] for row in cmsearch}

    with open(args.out, 'w') as out, open(args.unknowns, 'w') as unknowns:
        seqs = SeqIO.parse(args.fasta, 'fasta', Alphabet.IUPAC.ambiguous_dna)
        for s in seqs:
            if s.id in cmsearch and Alphabet._verify_alphabet(s.seq):
                if cmsearch[s.id] == '-':
                    s.seq = s.seq.reverse_complement()
                SeqIO.write(s, out, 'fasta')
            else:
                SeqIO.write(s, unknowns, 'fasta')


if __name__ == '__main__':
    main()
