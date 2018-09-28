#!/usr/bin/env python

import os
import shutil
import argparse


def find_primer_regions(forward_primer, reverse_primer, tmp_dir, output_fasta, input_fasta):
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)
    cmd = 'msa.sh in={input_fasta} out={forward_sam} literal={forward_primer}'.format(input_fasta=input_fasta,
                                                                                      forward_sam=os.path.join(tmp_dir, 'forward_primer.sam'),
                                                                                      forward_primer=forward_primer)
    os.system(cmd)
    cmd = 'msa.sh in={input_fasta} out={reverse_sam} literal={reverse_primer}'.format(input_fasta=input_fasta,
                                                                                      reverse_sam=os.path.join(tmp_dir, 'reverse_primer.sam'),
                                                                                      reverse_primer=reverse_primer)
    os.system(cmd)
    cmd = 'cutprimers.sh in={input_fasta} out={output_fasta} sam1={forward_sam} sam2={reverse_sam}'.format(input_fasta=input_fasta,
                                                                                                           output_fasta=output_fasta,
                                                                                                           forward_sam=os.path.join(tmp_dir, 'forward_primer.sam'),
                                                                                                           reverse_sam=os.path.join(tmp_dir, 'reverse_primer.sam'))
    os.system(cmd)
    shutil.rmtree(tmp_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--forward_primer',
                        type=str,
                        required=True,
                        help='DNA Sequence of forward primer.')
    parser.add_argument('-r', '--reverse_primer',
                        type=str,
                        required=True,
                        help='DNA Sequence of reverse primer.')
    parser.add_argument('-i', '--input_fasta',
                        type=str,
                        required=True,
                        help='Path to input fasta file where sequences need to be extracted from.')
    parser.add_argument('-o', '--output_fasta',
                        type=str,
                        required=True,
                        help='Path to output fasta file where amplified regions will be stored.')
    parser.add_argument('-t', '--tmp_dir',
                        type=str,
                        default='tmp',
                        help='Name of temporary directory used to store intermediate files.')
    args = parser.parse_args()

    find_primer_regions(forward_primer=args.forward_primer,
                        reverse_primer=args.reverse_primer,
                        input_fasta=args.input_fasta,
                        output_fasta=args.output_fasta,
                        tmp_dir=args.tmp_dir)
