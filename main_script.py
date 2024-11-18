# main_script.py

import os
import numpy as np
import pandas as pd
from config import get_config
from sequence_utils import complement, dna_to_mrna
from match_analysis import try_all_matches, score_match
from thermo_analysis import analyse_thermodynamics

# Get configuration settings
baseBuff, scoreLim, hpLEN, dataDIR, results_dir = get_config()

# Load MMP9 RNA sequence and convert to DNA
mmp9_rna = open(dataDIR + 'rna_mmp9.txt', 'r').read().replace('\n', '')
mmp9_dna = dna_to_mrna(mmp9_rna).replace('U', 'T')  # Convert RNA to DNA

# Define DNA sequences as a DataFrame
dna_seq = pd.DataFrame({
    'Name': ['sq1', 'sq2', 'sq3', 'sq4', 'sq5', 'sq6', 'sq7'],
    'Sequence': [
        'GGGCAGAGGTGTCT', 'TCGCTGTCAAAGTTCGAGGTG', 'TTCAGGGCGAGGACCATAGAG',
        'AATTCAGGGCGAGGACCATAG', 'GCCTCAGAGAATCGCCAGTACT', 
        'TTGAGGTCGCCCTCAAAGGTT', 'CCTCAAGTGGGACCATCATAA'
    ],
    'Complement': [
        complement('GGGCAGAGGTGTCT', 'DNA'), complement('TCGCTGTCAAAGTTCGAGGTG', 'DNA'),
        complement('TTCAGGGCGAGGACCATAGAG', 'DNA'), complement('AATTCAGGGCGAGGACCATAG', 'DNA'),
        complement('GCCTCAGAGAATCGCCAGTACT', 'DNA'), complement('TTGAGGTCGCCCTCAAAGGTT', 'DNA'),
        complement('CCTCAAGTGGGACCATCATAA', 'DNA')
    ],
    'Length': [14, 21, 21, 21, 22, 21, 21],
    'MW': [4359.9, 6477.2, 6520.3, 6504.3, 6704.4, 6437.2, 6399.2],
    'Tm(°C)': [43.2, 54.4, 56.3, 54.4, 56.7, 54.4, 52.4]
})

# Add the MMP9 DNA sequence column for thermodynamic analysis
dna_seq['MMP9_Seq'] = mmp9_dna

# Initialise lists to collect data for each sequence and its matches
dataMATCH = []

# Process each sequence in `dna_seq['Complement']`
for seq in dna_seq['Complement']:
    # Get the corresponding row from dna_seq for this `seq`
    dna_row = dna_seq[dna_seq['Complement'] == seq].iloc[0]
    
    # Call try_all_matches with a score limit of at least half length of hairpin
    matches = try_all_matches(seq, mmp9_dna, round(len(seq) * scoreLim))
    
    # Check if any matches were found
    if matches[0][0]:  # If there's at least one valid match
        hyb_col = True
        
        # Process each match
        for match in matches:
            _, query_start, query_end, subject_len, score = match
            
            # Calculate start and end positions for extracted slice with baseBuff buffer
            if query_start <= baseBuff:
                q_st = 0
                q_ed = query_end + baseBuff
            elif query_end >= len(mmp9_dna) - baseBuff:
                q_st = query_start - baseBuff
                q_ed = query_end + abs(len(mmp9_dna) - query_end)
            else:
                q_st = query_start - baseBuff
                q_ed = query_end + baseBuff
            
            # Extract the matched sequence from mmp9_dna with buffer
            matched_segment = mmp9_dna[q_st:q_ed]
            
            # Append data for each match, including all columns from `dna_seq`
            dataMATCH.append({
                'hyb_col': hyb_col,
                'mmp9_seqs': matched_segment,
                'search_seqs': seq,
                'query_start': query_start,
                'query_end': query_end,
                'subject_len': subject_len,
                'score': score,
                'Name': dna_row['Name'],
                'Sequence': dna_row['Sequence'],
                'Complement': dna_row['Complement'],
                'Length': dna_row['Length'],
                'MW': dna_row['MW'],
                'Tm(°C)': dna_row['Tm(°C)']
            })
    else:
        # If no matches met the score limit, add NaN entry for unmatched sequence, including dna_seq data
        dataMATCH.append({
            'hyb_col': False,
            'mmp9_seqs': np.nan,
            'search_seqs': seq,
            'query_start': None,
            'query_end': None,
            'subject_len': None,
            'score': None,
            'Name': dna_row['Name'],
            'Sequence': dna_row['Sequence'],
            'Complement': dna_row['Complement'],
            'Length': dna_row['Length'],
            'MW': dna_row['MW'],
            'Tm(°C)': dna_row['Tm(°C)']
        })

# Convert collected data to DataFrame
dfMATCH = pd.DataFrame(dataMATCH)

# Thermodynamic analysis on the hairpin sequences
hairpin_results = analyse_thermodynamics(dna_seq)

# Save the DataFrame to an Excel file
hairpin_results.to_excel(os.path.join(results_dir, 'hairpin_thermo_results.xlsx'), index=True)
dfMATCH.to_excel(os.path.join(results_dir, 'rn_seqDATA.xlsx'), index=True)