#!/usr/bin/env python

import os
import argparse
import multiprocessing
from Bio import SeqIO, pairwise2
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from io import StringIO


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
    parser.add_argument('-t', '--threads',
                        type=int,
                        default=multiprocessing.cpu_count(),
                        help='Number of threads to run. Defaults to all cores on machine.')
    args = parser.parse_args()
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)

    # Extract the fasta sequences for each of the two groups of interest.
    extract_fastas(taxonomy_level=args.taxonomy_level,
                   group=args.group_one,
                   silva_fasta=args.silva_fasta,
                   output_dir=args.output_dir)
    extract_fastas(taxonomy_level=args.taxonomy_level,
                   group=args.group_two,
                   silva_fasta=args.silva_fasta,
                   output_dir=args.output_dir)

    make_blast_db(os.path.join(args.output_dir, args.group_two + '.fasta'))
    sequence_one_seqs = extract_sequences(os.path.join(args.output_dir, args.group_one + '.fasta'))
    p = multiprocessing.Pool(processes=args.threads)
    database_list = [os.path.join(args.output_dir, args.group_two + '.fasta')] * len(sequence_one_seqs)
    blast_results = p.starmap(blast_sequence, zip(sequence_one_seqs, database_list))
    p.close()
    p.join()
    percent_id_all = list()
    with open(os.path.join(args.output_dir, args.group_one + '_' + args.group_two + '_detail.csv'), 'w') as f:
        f.write('QuerySequence,SubjectSequence,PercentID\n')
        for blast_result_list in blast_results:
            for blast_result in blast_result_list:
                percent_id_all.append(blast_result.percent_id)
                f.write('{},{},{}\n'.format(blast_result.query_name, blast_result.subject_name, blast_result.percent_id))
    print(sum(percent_id_all)/len(percent_id_all))
    with open(os.path.join(args.output_dir, 'distances.csv'), 'a+') as f:
        f.write('{},{},{}\n'.format(args.group_one, args.group_two, str(sum(percent_id_all)/len(percent_id_all))))


def make_blast_db(fasta_file):
    cmd = 'makeblastdb -in {} -dbtype nucl'.format(fasta_file)
    os.system(cmd)


def extract_sequences(fasta_file):
    sequences = list()
    for sequence in SeqIO.parse(fasta_file, 'fasta'):
        sequences.append(sequence)
    return sequences


def find_percent_identities(sequences):
    # Modified from https://www.biostars.org/p/208540/
    # As it turns out, doing pairwise alignments in python is painfully slow, so this would end up taking
    # months to complete for an all-vs-all run of stuff. Will attempt some blast
    alignment = pairwise2.align.globalxx(sequences[0], sequences[1], score_only=True)
    seq_length = max(len(sequences[0]), len(sequences[1]))
    percent_id = (alignment/seq_length) * 100
    return percent_id


class BlastResult:
    def __init__(self, percent_id, query_name, subject_name):
        self.percent_id = percent_id
        self.query_name = query_name
        self.subject_name = subject_name


def blast_sequence(sequence, database):
    # BLAST a (16S) sequence against a database of 16S sequence from another genus.
    blastn = NcbiblastnCommandline(db=database, outfmt=5, max_target_seqs=100000)  # By default, only get top 500 matches. Change to ridiculous level
    stdout, stderr = blastn(stdin=str(sequence.seq))
    percent_identities = list()
    if stdout:
        for record in NCBIXML.parse(StringIO(stdout)):
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.align_length >= 0.9 * len(sequence.seq):
                        percent_id = 100 * hsp.identities/len(sequence.seq)
                        subject_name = alignment.title
                        query_name = sequence.description
                        blast_result = BlastResult(percent_id=percent_id,
                                                   subject_name=subject_name,
                                                   query_name=query_name)
                        # percent_identities.append(100 * hsp.identities/len(sequence)
                        percent_identities.append(blast_result)
    return percent_identities


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
