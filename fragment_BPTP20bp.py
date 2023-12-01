# -*- encoding: utf-8 -*-


'''Get fragements' insert size and motif.'''


import argparse
import multiprocessing
import os
import pysam

import pandas as pd

from pandarallel import pandarallel


def do_fragment(bam, fasta, output, min_mapq=30, processes=20):
    '''Get fragements' insert size and motif.'''
    chrom_list = [f'chr{i}' for i in range(1, 23)]
    args_iter = [(bam, fasta, min_mapq, chrom) for chrom in chrom_list]
    pool = multiprocessing.Pool(processes=processes)
    fragment_list = []
    for df in pool.map(fragment, args_iter):
        fragment_list.append(df)
    fragment_df = pd.concat(fragment_list, ignore_index=True)
    fragment_df.to_csv(output, sep='\t', index=False, header=None, 
        float_format='%.2f')


def fragment(args):
    '''Wrapper for parallel.'''
    return fragment_generator(*args)


def fragment_generator(bam, fasta, min_mapq, chrom):
    '''BPM generator.'''
    def filter_read(read):
        '''True if the given read should be process.'''
        return not (read.is_duplicate or read.is_secondary or read.is_unmapped 
            or read.is_qcfail or read.mapq < min_mapq or 'S' in read.cigarstring
            or 'H' in read.cigarstring or read.is_supplementary 
            or not read.is_proper_pair)
    fa = pysam.FastaFile(fasta)
    read1_names, read1_chrom, read1_start, read1_end = [], [], [], []
    read1_sizes, read1_edm, read1_bpm = [], [], []
    read2_names, read2_chrom, read2_start, read2_end = [], [], [], []
    read2_sizes, read2_edm, read2_bpm = [], [], []
    frag_gc = []
    with pysam.Samfile(bam, 'rb') as bamfile:
        for read in bamfile.fetch(reference=chrom):
            if filter_read(read):
                if read.flag in [99, 163]:
                    ref_seq = fa.fetch(read.reference_name, 
                        read.reference_start, 
                        read.reference_start + abs(int(read.template_length)
                        )).upper()
                    edm_seq = read.get_reference_sequence().upper()[:4]
                    try:
                        bpm_seq = fa.fetch(read.reference_name, 
                            read.reference_start - 20, 
                            read.reference_start + 20).upper()
                    except:
                        bpm_seq = 'N'
                    if 'N' not in edm_seq and 'N' not in bpm_seq:
                        read1_names.append(read.query_name)
                        read1_chrom.append(read.reference_name)
                        read1_start.append(read.reference_start)
                        read1_end.append(read.reference_end)
                        read1_sizes.append(abs(int(read.template_length)))
                        read1_edm.append(edm_seq)
                        read1_bpm.append(bpm_seq)
                        frag_gc.append((ref_seq.count('G') + 
                            ref_seq.count('C')) / len(ref_seq))
                elif read.flag in [83, 147]:
                    edm_seq = read.get_reference_sequence().upper()[-4:]
                    try:
                        bpm_seq = fa.fetch(read.reference_name, 
                            read.reference_end - 20, 
                            read.reference_end + 20).upper()
                    except:
                        bpm_seq = 'N'
                    if 'N' not in edm_seq and 'N' not in bpm_seq:
                        read2_names.append(read.query_name)
                        read2_chrom.append(read.reference_name)
                        read2_start.append(read.reference_start)
                        read2_end.append(read.reference_end)
                        read2_sizes.append(abs(int(read.template_length)))
                        read2_edm.append(reverse_complement(edm_seq))
                        read2_bpm.append(reverse_complement(bpm_seq))
    read1_df = pd.DataFrame({'name': read1_names, 'chrom': read1_chrom, 
        'read1_start': read1_start, 'read1_end': read1_end, 
        'insert_size': read1_sizes, 'read1_edm': read1_edm, 
        'read1_bpm': read1_bpm, 'gc': frag_gc})
    read2_df = pd.DataFrame({'name': read2_names, 'chrom': read2_chrom, 
        'read2_start': read2_start, 'read2_end': read2_end, 
        'insert_size': read2_sizes, 'read2_edm': read2_edm, 
        'read2_bpm': read2_bpm})
    df = pd.merge(read1_df, read2_df, on=['name', 'chrom', 'insert_size'])
    df.drop(columns=['name'], inplace=True)
    fa.close()
    return df


def reverse_complement(seq):
    '''Reverse complement sequence.'''
    base_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    seq_list = list(seq)
    seq_list.reverse()
    rc_seq = ''
    for base in seq_list:
        rc_seq += base_dict[base]
    return rc_seq


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=do_fragment.__doc__)
    parser.add_argument('-b', '--bam', help='Dedup bam file.', required=True)
    parser.add_argument('-f', '--fasta', help='UCSC hg19 fasta.', required=True)
    parser.add_argument('-o', '--output', help='Output file.', required=True)
    args = parser.parse_args()

    dirname = os.path.dirname(args.output)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    do_fragment(args.bam, args.fasta, args.output)
