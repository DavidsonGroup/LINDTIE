'''
Module      : LINDTIE_get_novel_contigs
Description : Gets novel transcripts from transcript count matrix
Copyright   : (c) Jia Wei Tan, Dec 2025
License     : MIT
Maintainer  : https://github.com/jiawei-tan
Portability : POSIX
Reads the transcript count matrix, identifies novel transcripts,
filters by CPM, and writes the novel transcripts to a TSV file.

Adapted from the MINTIE pipeline (https://github.com/Oshlack/MINTIE)
'''

import pandas as pd
import numpy as np
import os
import sys
from Bio import SeqIO
from argparse import ArgumentParser

def parse_args(args):
    '''
    Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Get novel transcripts'
    parser = ArgumentParser(description = description)
    parser.add_argument(dest='tcm_file',
                        metavar='TCM_FILE',
                        type=str,
                        help='''transcript_counts_matrix.tsv file.''')
    parser.add_argument(dest='ref_tx_fasta',
                        metavar='REF_FASTA',
                        type=str,
                        help='''Transcriptome reference fasta.''')
    parser.add_argument(dest='denovo_fasta',
                        metavar='DENOVO_FASTA',
                        type=str,
                        help='''Sample de novo filtered fasta assembly file.''')
    # NEW ARGUMENT HERE
    parser.add_argument('--min_cpm',
                        dest='min_cpm',
                        type=float,
                        default=1.0,
                        help='''Minimum CPM threshold for filtering (default: 1.0).''')
    return parser.parse_args(args)

def get_ref_txs(ref_tx_fasta):
    '''
    Get all reference transcript IDs
    from reference transcriptome fasta
    '''
    handle = open(ref_tx_fasta, 'r')
    ref_txs = []
    for record in SeqIO.parse(handle, 'fasta'):
        ref_txs.append(record.id)
    handle.close()
    return ref_txs

def process_transcripts(tcm, ref_txs, outdir, min_cpm):
    '''
    Process transcript count matrix to identify novel transcripts
    and prepare output table.
    '''
    print('Processing transcript count matrix...')
    
    # Create working copy
    full_tx_table = tcm.copy()
    
    # Identify the sample column (assuming it's the one that isn't 'transcript_id')
    sample_columns = [col for col in tcm.columns if col != 'transcript_id']
    if sample_columns:
        case_sample = sample_columns[0]
        
        # --- CALCULATE CPM ---
        # Calculate library size (sum of all reads in the sample)
        library_size = full_tx_table[case_sample].sum()
        
        if library_size > 0:
            # Formula: (Counts / Library Size) * 1 Million
            full_tx_table['CPM'] = (full_tx_table[case_sample] / library_size) * 1e6
        else:
            full_tx_table['CPM'] = 0.0

        # Calculate logCPM (log2) for output, adding 1 to avoid log(0)
        full_tx_table['logCPM'] = np.log2(full_tx_table['CPM'] + 1)

        # Rename the sample column to 'num_read_case'
        full_tx_table = full_tx_table.rename(columns={case_sample: 'num_read_case'})
    
    # Identify novel transcripts (not in reference)
    full_tx_table['is_novel'] = ~full_tx_table['transcript_id'].isin(ref_txs)

    # Filter to keep ONLY novel transcripts (excludes known references)
    print('Filtering for novel transcripts...')
    full_tx_table = full_tx_table[full_tx_table['is_novel']]

    # --- FILTER BY CPM ---
    print(f'Filtering for transcripts with CPM >= {min_cpm}...')
    full_tx_table = full_tx_table[full_tx_table['CPM'] >= min_cpm]

    if len(full_tx_table) == 0:
        print('Warning: No novel transcripts found!')

    # Add dummy DE columns so LINDTIE_post_process.py does not crash.
    # use logFC=10 and FDR=0 to ensure these are treated as "significant" downstream.
    full_tx_table['logFC'] = 10.0
    full_tx_table['FDR'] = 0.0
    full_tx_table['PValue'] = 0.0
    # Note: 'logCPM' is now calculated above, so we don't hardcode it to 1.0 anymore
    
    # Write the main significant transcript file
    print(f'Writing output to {outdir}/DE_transcript_significant.txt...')
    full_tx_table.to_csv('%s/DE_transcript_significant.txt' % outdir, sep='\t', index=False)

def main():
    args = parse_args(sys.argv[1:])
    try:
        print('Reading in transcript count matrix file...')
        tcm = pd.read_csv(args.tcm_file, sep = '\t')
        print('Fetching reference transcripts...')
        ref_txs = get_ref_txs(args.ref_tx_fasta)
    except IOError as message:
        print("{} ERROR: {}, exiting".format("get_novel_contigs", message), file=sys.stderr)
        sys.exit(1)

    # Process and write transcripts
    outdir = os.path.dirname(args.tcm_file)
    # Pass the min_cpm argument to the processing function
    process_transcripts(tcm, ref_txs, outdir, args.min_cpm)
    
if __name__ == '__main__':
    main()