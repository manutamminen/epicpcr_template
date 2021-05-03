#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import gzip
from itertools import groupby


bact_stds = [[">Mock1_16S;tax=d:Mock1_16S,p:Mock1_16S,c:Mock1_16S,o:Mock1_16S,f:Mock1_16S,g:Mock1_16S,s:Mock1_16S\n",
"""GGCTACCTTGTTACGACTTAGGTGTTAGCAGCTTACCTCAGGATCCCAACTGATCCTCGTACGATTTACCAAATCGGTCTT
AGGGGTTTGAGCGTCCGTGGCATCCCTGTTGGGCCGGTAGACTACCAGGGTATCTAATCTTCAGACGTGATTTATTGTACA
CGGAGCTGGTACCTTGGGCGACGATTTTGTCGGGAGCGTGCTTCTCATACGCTCCAAGCTACTTCATCCGTACTATAACGG
GTTCAGACGGTTTTGCGGCCTTATTCGGCCTTCCAGGAAATGCTAAGTATTCCCACTACTATAAACATCCTTATAGCCTTC
GCTCGTAACTATAACGGTCCTAAGGTAGCGAAGCTTAGTCCGGAATTACCGCGGCGGCTG""".replace("\n", "")],
[">Mock2_16S;tax=d:Mock2_16S,p:Mock2_16S,c:Mock2_16S,o:Mock2_16S,f:Mock2_16S,g:Mock2_16S,s:Mock2_16S\n",
"""GGCTACCTTGTTACGACTTAAATCGGAGGCTACAAGCATCTTGAACCTAAGGCTAGGAGACGTGCTAGTGGTCCTGTGTGA
AAACGTCTCCCAGATCTTGGGCCGACGTTAAGGTTGAGAGACTACCAGGGTATCTAATAGGAATTTTCGATGGCCGGGCCG
ATCTAAGCTAGTCCCTCCGGTAGCTTATGAGACGCAAGTAGACGCTTCTCAACTAACGATCCTCTCGGAACAGTGAACTAA
GTAATAGCACGGCGACGCTGTTGCTAAGTGATGGAGCAATACGATCACTGGCCCTTTAATCCCGTTGTAGCTGATGTCCAG
CGATCTAACTATAACGGTCCTAAGGTAGCGAAATGATCGAGTGAATTACCGCGGCGGCTG""".replace("\n", "")],
[">Mock3_16S;tax=d:Mock3_16S,p:Mock3_16S,c:Mock3_16S,o:Mock3_16S,f:Mock3_16S,g:Mock3_16S,s:Mock3_16S\n",
"""GGCTACCTTGTTACGACTTTAAGGGGTAATAGTCTGGCCTGCAGGCTAGAGACACACTGCCCACGGATATTGCTCCTGCTA
GGACAGTATCCCGTGTCACCTGGTCTCGTTCTATCACAAGACTACCAGGGTATCTAATTTGTGGTGACACCTCACTTATCT
CGCTCCAATTAATTGAGTCTACTGGTGCTTCGACAACCCTCTAAATTATTGAAGCAAGTCTAACATCCGCACGGGTGTTGG
TCGGAGACCCAACTACCTCACTTGGACTGTGGCCCAAGCAACAGACCCAGAACAAATGCTTACCAGCAACTTCGGTTACCG
TGGCTTAACTATAACGGTCCTAAGGTAGCGAAATCCGATAGTGAATTACCGCGGCGGCTG""".replace("\n", "")],
[">Mock4_18S_16S;tax=d:Mock4_18S_16S,p:Mock4_18S_16S,c:Mock4_18S_16S,o:Mock4_18S_16S,f:Mock4_18S_16S,g:Mock4_18S_16S,s:Mock4_16S\n",
"""GGCTACCTTGTTACGACTTTAAGCTAGCTAAGCGGTCGCAGCTGGAAGCCGAAAGGCCCATCAATAATTAAGCGGGGATCC
GGGCCGTCTTTGTGCAGAGTTATCCTTTGCTAATGCTAGACTACCAGGGTATCTAATGCGGTCTTAATAGATCCCGGAGTC
GCTACAGTGCAACACCCTGGACGATTTCGTCATGAGTAATCCGACGGCGAGCGAAACAGAGGGCGGACATACCCTTGTTGT
GAGATCGGACCTCGGGGGCAAAACGGTGGCAACCTGCAACGCGGTCGCTGACACGTCCCCGTGACGCAAGAATCTATGCAA
AGATAGAAGGACCGTTGTGAGTATCAGGGTTAGTGTATCTCGAATTACCGCGGCGGCTGGCGACGGGCGGTGTGTACAAGA
GCTAGCATTGTTGGAAGTCCTCAAGTGCAGGACTGATTACCTCGGTGTGTATCGTGTTGAATGTGTACAGATCGAATGTTC
GGGTGAATGATATCGTCGAACTCAGTCAGTCTGTAAGACTTATTTCTGTCTCTACCTTCCTGCAATCCAGCTGGGCCGAGG
CTGGATGAATAAGGCCGCGGTTAAGTCAGGTTATTTCTCATTATATCGGTGGATATTTAGGGGCATCACAGACCTGTTATT
CAGGCCGCCTCTATGCGTGCTTGTAGTGACGTATCCGCCCTTTATATGCTCAAAACTGCCCAATAGTGAACTCCAGGAACG
TCTGGATTCATAATGAACAAGAACGGCCATGCACCACCA""".replace("\n", "")],
[">Mock5_18S_16S;tax=d:Mock5_18S_16S,p:Mock5_18S_16S,c:Mock5_18S_16S,o:Mock5_18S_16S,f:Mock5_18S_16S,g:Mock5_18S_16S,s:Mock5_16S\n",
"""GGCTACCTTGTTACGACTTGCCAGCCTACTGGGACGTGTATGTATCGCGTTCTATATAGATACAACATACGGTCTACGCCT
AGTCTGTGAATGCCCGGCCGACTAGGTTATGTAGCACTAGACTACCAGGGTATCTAATGACGTCGAGGAGCCTTTGAGCAG
CAAACGTTTTAGGGGCCTGCCCGTGTGTCAGCAATTACCGTGGAGGGCTATATTTTTCCACCTATCCGCGCTGCTCAGGAC
GTCCGGCTACTCCTGCATGGGAAACCCACGCCCTTTGGAAATCATCCATTGCGTGTTTCCGAAAAATAGTTTGTACGTGCG
CGGTTGAATTCTTATGCCTTACTGTAGAGCTCCTCTTCGAAAGAATTACCGCGGCGGCTGGCGACGGGCGGTGTGTACCAT
GCCTCATTGAGGTCTTAGGTTCGCATAGTTGATGTATGTAGTCAGCCCACCGGATTACAGAGGGAGTAGGAAGTACACTCG
TACGAGTAATGCAGAAATCGCAGTGTTCGTAAGATCATGCGGGAACCACACTACCCACCGAATCTTTGCTAACTTGGCATC
TTTAACGTTGTCAATAGCCCGCCACGGTGCCGGAACTAGAGGCAGTCACTGCGAGGAGCCATGGCATCACAGACCTGTTAT
CTTTACAGACAGAGCGGCCTCAGTAATACTAGGGCCACGTTGTTGGAAAGTACCGAAAGGTCGGCACGTGAAGATCCCGGA
GGTGCACTCAATGCAACTGAAGAACGGCCATGCACCACCA""".replace("\n", "")],
[">Mock6_18S_16S;tax=d:Mock6_18S_16S,p:Mock6_18S_16S,c:Mock6_18S_16S,o:Mock6_18S_16S,f:Mock6_18S_16S,g:Mock6_18S_16S,s:Mock6_16S\n",
"""GGCTACCTTGTTACGACTTGCCAGGAGCCGGGCTGTAGCATCTGTTCAACACAAATCCGCTCTGGTAAATACTAGTAATAG
GTGTTGCTATTGGCTGTTGAATGGATAGAATGGCTCGGAGACTACCAGGGTATCTAATTCGCCCGTAAAACTGATAGGCCT
TCAGTACGGCGTTGGGGTCTTTACCACGGGCCCTTACCACCTTAGCACCTGCCTCTTCTCGGGCACCTAACAAGCAAGTCA
TTGTGTGTGGCACGCGTAGATAACCCGCAGGATTCACGCCGAATACGTCGACACCCATTCACATGGCCCGCGCCTTAAACA
AGGGATCTCAGGTAGGAAATCCTATCACTAACAAGTAAGCAAGAATTACCGCGGCGGCTGGCGACGGGCGGTGTGTACATC
TCGAACTACTTCCGACAATTATCGCGTCATGTCCCAAAAATCTACGGATGCCAGAGAACATGAAGAGCGTGAAGTGCGACG
TCTGAGTCGGCGCCATGAGCATGGAACTGGATGATGTTCATCGCATGGTTTTTGTGAACGCAGCAGTACGCTATTTCCGGA
TTATGGCATCTCTGGTTACGTAAAACTTACTCTTAACCATCTAAGAAGACTAAGCCATATGAGGCATCACAGACCTGTTAT
CCCCTAGCAATCCACGTATGAAATTCGAGTCATGTGCACCTTAACGGAACGCTGAATCCTCAGTCACCTCCCGCTCGGGTT
TGTGTCTGATATGACATCTAAGAACGGCCATGCACCACCA""".replace("\n", "")]]

euk_stds = [[">Mock1_18S;tax=d:Mock1_18S,p:Mock1_18S,c:Mock1_18S,o:Mock1_18S,f:Mock1_18S,g:Mock1_18S,s:Mock1_16S\n",
"""GCGACGGGCGGTGTGTACCAATGGGTATTAACTATAACGGTCCTAAGGTAGCGAAGCTGTAGTGAATGCCAAGTGTTTTAC
TTGAGGCGGACGGATTGCTGTACGCGGTAGGACAAATTTCCAGCTAGTACGTTTGCGTTGGCTCGATAAATCCATAGGAGC
TGCAAACTCAGTCCGCATCCGGAACATTCGGAGCTCTAATAGAACTATCACAGCGACACGGCTTTGTGCAGTCCATATCGG
GGGGCATCACAGACCTGTTATGCGCACCGCAGCATTGGCTATGTTCTAAAACTCTGACTAATCCCTTACTTGCGTGCATAC
GTCGTGCACCTCGAGGCCGGAGGATAAAAGACTGGGGCATAAGAACGGCCATGCACCACCA""".replace("\n", "")],
[">Mock2_18S;tax=d:Mock2_18S,p:Mock2_18S,c:Mock2_18S,o:Mock2_18S,f:Mock2_18S,g:Mock2_18S,s:Mock2_16S\n",
"""GCGACGGGCGGTGTGTACTAAAACTCGTTAACTATAACGGTCCTAAGGTAGCGAAGACATATGTGAAAAACTTGGTCTAAT
GGGCTCGAGTCGTTTTTTATATACATGCTTATGTGTGGAACTCCAGCTTTGCCGGTGTTCGTTAAGCGTGGTGTCACGACC
CTGAAAAAATCCCGCACCTTTACCGTCAGAAGAAGTTCAATTTAAACGCGGGACATGCCACAAACTCAGCTCCGTTAAGAT
GCGGCATCACAGACCTGTTATCTCCGTCGACAATTTCATCGGTCGCAGAAACCAAAAGATTTATTGCTTCTCAAGTGTGCA
GACCCTGTGCTTCTGTTAAAATCTATAGACTGGCATGCATAAGAACGGCCATGCACCACCA""".replace("\n", "")],
[">Mock3_18S;tax=d:Mock3_18S,p:Mock3_18S,c:Mock3_18S,o:Mock3_18S,f:Mock3_18S,g:Mock3_18S,s:Mock3_16S\n",
"""GCGACGGGCGGTGTGTACTCTTCGAGAGTAACTATAACGGTCCTAAGGTAGCGAAAGCCTCGGGCATTCTCTGCTTCACAT
CATCGGATGTTTGCGAACTAGGACTCTACTACATCCGCTATTGAAGTTTCGTAGCTCATATTGGGCTGTGATTTTGCCGAT
TGTCGCTCGACAAGGGATACGCAGTAGCCTAGGGAGGCAAAGGACAGATGATAGGCGGTCGACCGCATATCAATTGCGTCT
CCGGCATCACAGACCTGTTATCTACTAGTGATTGCCCTAATCTAATAGCTTCCCCATCGCAGGGAGGGTAAGATATGTGTG
CTCTTAGCTCACATGCCCTTCAGTAAGGTAGATGAGGACCAAGAACGGCCATGCACCACCA""".replace("\n", "")],
[">Mock4_18S_16S;tax=d:Mock4_18S_16S,p:Mock4_18S_16S,c:Mock4_18S_16S,o:Mock4_18S_16S,f:Mock4_18S_16S,g:Mock4_18S_16S,s:Mock4_16S\n",
"""GGCTACCTTGTTACGACTTTAAGCTAGCTAAGCGGTCGCAGCTGGAAGCCGAAAGGCCCATCAATAATTAAGCGGGGATCC
GGGCCGTCTTTGTGCAGAGTTATCCTTTGCTAATGCTAGACTACCAGGGTATCTAATGCGGTCTTAATAGATCCCGGAGTC
GCTACAGTGCAACACCCTGGACGATTTCGTCATGAGTAATCCGACGGCGAGCGAAACAGAGGGCGGACATACCCTTGTTGT
GAGATCGGACCTCGGGGGCAAAACGGTGGCAACCTGCAACGCGGTCGCTGACACGTCCCCGTGACGCAAGAATCTATGCAA
AGATAGAAGGACCGTTGTGAGTATCAGGGTTAGTGTATCTCGAATTACCGCGGCGGCTGGCGACGGGCGGTGTGTACAAGA
GCTAGCATTGTTGGAAGTCCTCAAGTGCAGGACTGATTACCTCGGTGTGTATCGTGTTGAATGTGTACAGATCGAATGTTC
GGGTGAATGATATCGTCGAACTCAGTCAGTCTGTAAGACTTATTTCTGTCTCTACCTTCCTGCAATCCAGCTGGGCCGAGG
CTGGATGAATAAGGCCGCGGTTAAGTCAGGTTATTTCTCATTATATCGGTGGATATTTAGGGGCATCACAGACCTGTTATT
CAGGCCGCCTCTATGCGTGCTTGTAGTGACGTATCCGCCCTTTATATGCTCAAAACTGCCCAATAGTGAACTCCAGGAACG
TCTGGATTCATAATGAACAAGAACGGCCATGCACCACCA""".replace("\n", "")],
[">Mock5_18S_16S;tax=d:Mock5_18S_16S,p:Mock5_18S_16S,c:Mock5_18S_16S,o:Mock5_18S_16S,f:Mock5_18S_16S,g:Mock5_18S_16S,s:Mock5_16S\n",
"""GGCTACCTTGTTACGACTTGCCAGCCTACTGGGACGTGTATGTATCGCGTTCTATATAGATACAACATACGGTCTACGCCT
AGTCTGTGAATGCCCGGCCGACTAGGTTATGTAGCACTAGACTACCAGGGTATCTAATGACGTCGAGGAGCCTTTGAGCAG
CAAACGTTTTAGGGGCCTGCCCGTGTGTCAGCAATTACCGTGGAGGGCTATATTTTTCCACCTATCCGCGCTGCTCAGGAC
GTCCGGCTACTCCTGCATGGGAAACCCACGCCCTTTGGAAATCATCCATTGCGTGTTTCCGAAAAATAGTTTGTACGTGCG
CGGTTGAATTCTTATGCCTTACTGTAGAGCTCCTCTTCGAAAGAATTACCGCGGCGGCTGGCGACGGGCGGTGTGTACCAT
GCCTCATTGAGGTCTTAGGTTCGCATAGTTGATGTATGTAGTCAGCCCACCGGATTACAGAGGGAGTAGGAAGTACACTCG
TACGAGTAATGCAGAAATCGCAGTGTTCGTAAGATCATGCGGGAACCACACTACCCACCGAATCTTTGCTAACTTGGCATC
TTTAACGTTGTCAATAGCCCGCCACGGTGCCGGAACTAGAGGCAGTCACTGCGAGGAGCCATGGCATCACAGACCTGTTAT
CTTTACAGACAGAGCGGCCTCAGTAATACTAGGGCCACGTTGTTGGAAAGTACCGAAAGGTCGGCACGTGAAGATCCCGGA
GGTGCACTCAATGCAACTGAAGAACGGCCATGCACCACCA""".replace("\n", "")],
[">Mock6_18S_16S;tax=d:Mock6_18S_16S,p:Mock6_18S_16S,c:Mock6_18S_16S,o:Mock6_18S_16S,f:Mock6_18S_16S,g:Mock6_18S_16S,s:Mock6_16S\n",
"""GGCTACCTTGTTACGACTTGCCAGGAGCCGGGCTGTAGCATCTGTTCAACACAAATCCGCTCTGGTAAATACTAGTAATAG
GTGTTGCTATTGGCTGTTGAATGGATAGAATGGCTCGGAGACTACCAGGGTATCTAATTCGCCCGTAAAACTGATAGGCCT
TCAGTACGGCGTTGGGGTCTTTACCACGGGCCCTTACCACCTTAGCACCTGCCTCTTCTCGGGCACCTAACAAGCAAGTCA
TTGTGTGTGGCACGCGTAGATAACCCGCAGGATTCACGCCGAATACGTCGACACCCATTCACATGGCCCGCGCCTTAAACA
AGGGATCTCAGGTAGGAAATCCTATCACTAACAAGTAAGCAAGAATTACCGCGGCGGCTGGCGACGGGCGGTGTGTACATC
TCGAACTACTTCCGACAATTATCGCGTCATGTCCCAAAAATCTACGGATGCCAGAGAACATGAAGAGCGTGAAGTGCGACG
TCTGAGTCGGCGCCATGAGCATGGAACTGGATGATGTTCATCGCATGGTTTTTGTGAACGCAGCAGTACGCTATTTCCGGA
TTATGGCATCTCTGGTTACGTAAAACTTACTCTTAACCATCTAAGAAGACTAAGCCATATGAGGCATCACAGACCTGTTAT
CCCCTAGCAATCCACGTATGAAATTCGAGTCATGTGCACCTTAACGGAACGCTGAATCCTCAGTCACCTCCCGCTCGGGTT
TGTGTCTGATATGACATCTAAGAACGGCCATGCACCACCA""".replace("\n", "")]]

def read_fasta(fasta, remove_commas=False):
    """ Read a fasta file and yield pairs of fasta ids and DNA sequences.
    """
    with gzip.open(fasta, 'rt') as f:
        grouped = groupby(f, lambda x: x[0] == ">")
        for cond, entry in grouped:
            if cond:
                fasta_id = next(entry)
                if remove_commas:
                    fasta_id = fasta_id.replace(",", "")
                _, seq_iter = next(grouped)
                seq = ''.join([line.strip() for line in seq_iter]).upper()
                yield ([fasta_id, seq])


def prepare_database(input_file, output_bact, output_euks):
    with open(output_bact, "w") as bact, open(output_euks, "w") as euks:
        tax_ids = ['tax=k:', 'p:', 'c:', 'o:', 'f:', 'g:', 's:']
        for id, seq in bact_stds:
            bact.write(id)
            bact.write(seq + "\n")
        for id, seq in euk_stds:
            euks.write(id)
            euks.write(seq + "\n")
        for fasta_id, seq in read_fasta(input_file):
            fasta_id, tax = fasta_id.split(" ", 1)
            tax = tax.strip().split(";")
            tax = list(map(lambda x: "".join(x), zip(tax_ids, tax)))
            tax = ",".join(tax)
            fid = fasta_id + ";" + tax
            seq = seq.replace("U", "T")
            if "Eukaryota" in fid:
                euks.write(fid + "\n" + seq + "\n")
            elif "Bacteria" in fid:
                bact.write(fid + "\n" + seq + "\n")


if __name__ == '__main__':
    log_fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_fmt)
    logger = logging.getLogger(__name__)
    logger.info('Setting up the 16S and 18S rRNA databases...')
    prepare_database(snakemake.input[1],
            snakemake.output[0],
            snakemake.output[1])