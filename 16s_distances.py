#!/usr/bin/env python

import os
import argparse
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--silva-fasta',
                        type=str,
                        required=True,
                        help='Path to a FASTA file downloaded from SILVA. File should have been downloaded with the '
                             '"no gaps" option.')
    parser.add_argument('-l', '--taxonomy_level',
                        choices=['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'],
                        default='Genus',
                        help='Your desired taxonomy level.')
    parser.add_argument('-g1', '--group_one',
                        type=str,
                        required=True,
                        help='Name of first taxa you want to compare (as it appears in the SILVA Fasta).')
    parser.add_argument('-g2', '--group_two',
                        type=str,
                        required=True,
                        help='Name of second taxa you want to compare (as it appears in the SILVA Fasta).')
    parser.add_argument('-o', '--output_dir',
                        type=str,
                        required=True,
                        help='Name of desired output directory. Created if it does not exist.')
    args = parser.parse_args()
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    extract_fastas(taxonomy_level=args.taxonomy_level,
                   group=args.group_one,
                   silva_fasta=args.silva_fasta,
                   output_dir=args.output_dir)
    extract_fastas(taxonomy_level=args.taxonomy_level,
                   group=args.group_two,
                   silva_fasta=args.silva_fasta,
                   output_dir=args.output_dir)


def extract_fastas(taxonomy_level, group, silva_fasta, output_dir):
    sequences_to_write = list()
    for sequence in SeqIO.parse(silva_fasta, 'fasta'):
        tax_dict = header_to_taxonomy(sequence.description)
        if tax_dict:
            if tax_dict[taxonomy_level] == group:
                sequences_to_write.append(sequence)
    SeqIO.write(sequences_to_write, os.path.join(output_dir, group + '.fasta'), 'fasta')


def header_to_taxonomy(fasta_header):
    taxonomy = ''.join(fasta_header.split()[1:])
    tax_dict = dict()
    taxonomy = taxonomy.split(';')
    if len(taxonomy) != 7:
        return None
    tax_dict['Kingdom'] = taxonomy[0]
    tax_dict['Phylum'] = taxonomy[1]
    tax_dict['Class'] = taxonomy[2]
    tax_dict['Order'] = taxonomy[3]
    tax_dict['Family'] = taxonomy[4]
    tax_dict['Genus'] = taxonomy[5]
    tax_dict['Species'] = taxonomy[6]
    return tax_dict


if __name__ == '__main__':
    main()