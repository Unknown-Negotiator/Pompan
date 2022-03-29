#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#@created: 07.05.2021
#@author: SWW 2021 team
#@contact: support@noname_tool.com

import argparse
import re
import os
import numpy as np
import fastaparser
import pprint
from Bio import Seq, SeqUtils, SeqIO


def count_GGG_triplets(sequence):
    '''Count GGG triplets in given sequence.
    
    Parameters:
    sequence (string): DNA sequence

    Returns:
    int: A number of GGG triplets
    '''
    return sequence.count('GGG')


def fasta_reader(fasta_file):
    """
    Args:
        fasta_file: path to genome in fasta format

    Returns:
         dictionary with key == fasta header, value == sequence
    """
    with open(fasta_file) as fh:
        parser = fastaparser.Reader(fh, parse_method='quick')
        seq_dict = {}
        for seq in parser:
            seq_dict[seq.header] = seq.sequence
    return seq_dict


def genome_N50(fasta_file):
    """
    Args:
        fasta_dictionary: dictionary with key == fasta header, value == sequence (fasta_reader function)

    Returns:
        N50 characteristics of genome
    """
    SeqLen = []
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            SeqLen.append(len(record.seq))
    ReverseLen = sorted(SeqLen, reverse=True)
    csum = np.cumsum(ReverseLen)
    MedianCsum = csum[-1] // 2
    csum2 = min(csum[csum >= MedianCsum])
    for index, item in enumerate(csum):
        if item == csum2:
            indexN50 = index
            N50 = ReverseLen[indexN50]
    return N50


def genome_length(fasta_file):
    """
    Args:
        fasta_dictionary: dictionary with key == fasta header, value == sequence (fasta_reader function)

    Returns:
        genome_length
    """
    SeqLen = []
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            SeqLen.append(len(record.seq))
    return int(sum(SeqLen))


def genome_gc(fasta_dictionary):
    """
    Args:
        fasta_dictionary: dictionary with key == fasta header, value == sequence (fasta_reader function)

    Returns:
        genome gc content (%)
    """
    total_CG = 0
    total_len = 0
    for fasta in fasta_dictionary.values():
        total_CG += fasta.count('C') + fasta.count('G')
        total_len += len(fasta)
    return round(total_CG / total_len * 100, 1)


def nucleotide_frequency(fasta_dictionary):
    """
    Args:
        fasta_dictionary: dictionary with key == fasta header, value == sequence (fasta_reader function)

    Returns:
        dictionary with frequency of each nucleotide, key == k-mer (sequence), value == frequency
    """
    A, T, C, G, N = 0, 0, 0, 0, 0
    for fasta in fasta_dictionary.values():
        A += fasta.count('A')
        T += fasta.count('T')
        C += fasta.count('C')
        G += fasta.count('G')
        N += fasta.count('N')
    return {'A': A, 'T': T, 'C': C, 'G': G, 'N': N}


def kmer_frequency(fasta_dictionary, k):
    """
    Args:
        fasta_dictionary: dictionary with key == fasta header, value == sequence (fasta_reader function)
        k: length of k-mer sequence

    Returns:
        dictionary with k-mer frequency, key == k-mer (sequence), value == frequency
    """
    k_mers = {}
    for heading in fasta_dictionary.keys():
        sequence = fasta_dictionary[heading]
        stop = int(k)
        for start in range(len(sequence)):
            if stop <= len(sequence):
                k_mer = sequence[start:stop]
                stop += 1
                if k_mer not in k_mers.keys():
                    k_mers[k_mer] = 1
                else:
                    k_mers[k_mer] += 1

    return k_mers


def unique_k_mers(k_mers_dict):
    uniq_count = 0
    for k_mer in k_mers_dict.keys():
        if k_mers_dict[k_mer] == 1:
            uniq_count += 1
    return uniq_count


def repetitive_k_mers(k_mers_dict):
    rep_count = 0
    for k_mer in k_mers_dict.keys():
        if k_mers_dict[k_mer] != 1:
            rep_count += 1
    return rep_count


def orf_finder(fasta_dictionary):
    """
    Args:
        fasta_dictionary: dictionary with key == fasta header, value == sequence (fasta_reader function)

    Returns:
        dictionary with key == number of ORF, value == array with positions of start and stop codon
    """
    pass


def codon_frequency(fasta_dictionary):
    """
    Args:
        fasta_dictionary: dictionary with key == fasta header, value == sequence (fasta_reader function)

    Returns:
        dictionary with key == codon (e.g. "ATG"), value == frequency of codone
    """
    dict_codons = {}
    for header, seq in fasta_dictionary.items():
        for i in range(0, len(seq) - 3, 3):
            if seq[i:i + 3] not in dict_codons.keys():
                dict_codons[seq[i:i + 3]] = 1
            else:
                dict_codons[seq[i:i + 3]] += 1
        return dict_codons


def aminoacid_frequency(fasta_dictionary):
    """
    Args:
        fasta_dictionary: dictionary with key == fasta header, value == sequence (fasta_reader function)

    Returns:
        dictionary with key == aminoacid (e.g. "M"), value == frequency of aminoacid
    """
    protein = {"TTT": "F", "CTT": "L", "ATT": "I", "GTT": "V",
               "TTC": "F", "CTC": "L", "ATC": "I", "GTC": "V",
               "TTA": "L", "CTA": "L", "ATA": "I", "GTA": "V",
               "TTG": "L", "CTG": "L", "ATG": "M", "GTG": "V",
               "TCT": "S", "CCT": "P", "ACT": "T", "GCT": "A",
               "TCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
               "TCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
               "TCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",
               "TAT": "Y", "CAT": "H", "AAT": "N", "GAT": "D",
               "TAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
               "TAA": "STOP", "CAA": "Q", "AAA": "K", "GAA": "E",
               "TAG": "STOP", "CAG": "Q", "AAG": "K", "GAG": "E",
               "TGT": "C", "CGT": "R", "AGT": "S", "GGT": "G",
               "TGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
               "TGA": "STOP", "CGA": "R", "AGA": "R", "GGA": "G",
               "TGG": "W", "CGG": "R", "AGG": "R", "GGG": "G"
               }

    def translate(seq, protein):
        trans_seq = ''
        for i in range(0, len(seq) - 3, 3):
            trans_seq += protein[seq[i: i + 3]]
        return trans_seq

    def count_aminas(trans_seq, list_aminas):
        aminas = {}
        for amina in list_aminas:
            aminas[amina] = trans_seq.count(amina)
        return aminas

    aminas_answer = {}
    list_aminas = np.unique(list(protein.values()))
    for header, seq in fasta_dictionary.items():
        trans_seq = translate(seq, protein)
        aminas = count_aminas(trans_seq, list_aminas)
        if aminas_answer == {}:
            aminas_answer = aminas
        else:
            for k, v in aminas_answer.items():
                aminas_answer[k] += aminas[k]
    return aminas_answer


def main(input_files, output_file, k):
    with open(output_file, 'w') as fw:
        fw.write("%s\n" % "\t".join(["Genome", "N50", "Genome length", "Genome GC (%)",
                                     f"Unique k-mers_{k}", f"Repetitive k-mers {k}"]))
        for file in input_files:
            # считаем все статистики по файлам
            # результат записываем в словарь, например results["N50"] = genome_N50
            genome_name = os.path.basename(file)
            fasta_dictionary = fasta_reader(file)
            N50 = genome_N50(file)
            length = genome_length(file)
            gc = genome_gc(fasta_dictionary)
            k_mers = kmer_frequency(fasta_dictionary, k)
            k_mers_unique= unique_k_mers(k_mers)
            k_mers_rep = repetitive_k_mers(k_mers)

            results = [genome_name, str(N50), str(length), str(gc), str(k_mers_unique), str(k_mers_rep)]
            fw.write("%s\n" % "\t".join(results))


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Super noname tool!')
    parser.add_argument('-f','--fasta', help='path to input file(s) in fasta format', nargs="+", required=True)
    parser.add_argument('-k', '--kmer', help='length of k for k-mers statistics', required=True)
    parser.add_argument('-o','--output', help='path to output file', required=True)
    args = vars(parser.parse_args())

    input_files = args["fasta"]
    output_file = args["output"]
    k = args["kmer"]

    main(input_files, output_file, k)
