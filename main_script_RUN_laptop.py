#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 18:56:43 2024

@author: michaelbh23
"""

import platform
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as sc
import seaborn as sns
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SeqUtils import molecular_weight, MeltingTemp as mt
from nupack import *
import itertools
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math
from config import get_config
import matplotlib.cm as cm
from tqdm import tqdm  # Import the tqdm library
from difflib import SequenceMatcher  # For quick similarity scoring
import random


# Get the system/OS name
os_name = platform.system()
doTEST = 1 #this is for tesing the code so we dont have to loop through everything
WrongHP = 1 #this is for seeing what happens when we make bad hairpin
# Set the flag to determine loop type
doFullGene = 1  # Set to 1 for processing both hairpin and gene, 0 for just hairpin
# Get configuration settings
scoreLim, hpLEN, dataDIR, results_dir, doPLOT, numTOP, lowest_probability = get_config()
basePAD = hpLEN # if 0  then no pad. Was default of 15
plt.close('all')

# Define physical model
temperature = 37 
mmp9_model = Model(material='dna04', celsius=temperature, sodium=0.5, magnesium=0.2)
#mmp9_model = Model(material='dna04')
# Define the energy gap (kcal/mol) for calculating suboptimal structures
gap = 0.5 # 1.1 orginal

def _complement(sequence,verified):
    '''This function returns a reverse complement 
    of a DNA or RNA strand'''
    #verified = verify_na(sequence)
    if verified == "DNA":
        # complement strand
        seq = sequence.replace("A", "t").replace(
            "C", "g").replace("T", "a").replace("G", "c")
        seq = seq.upper()
        # reverse strand
        seq = seq[::-1]
        return seq
 
    elif verified == "RNA":
       
        # complement strand
        seq = sequence.replace("A", "u").replace(
            "C", "g").replace("U", "a").replace("G", "c")
        seq = seq.upper()
         
        # reverse strand
        seq = seq[::-1]
        return seq
    else:
        return False
    
    
class dna_sequence:
    def __init__(self, name, sequence, length, mw, tm):
        self.name = name
        self.sequence = sequence
        self.complement = _complement(self.sequence,'DNA')
        self.length = length
        self.mw = mw
        self.tm = tm
        
 # Best two sequences for each of the siRNA
 # Add first sequence
 # Control siRNA which does not hybridise with any mRNA with the same hairpin sticky ends: TTCTCCGAACGTGTCACGT
 # Check with control if it forms a complex with the mMMP9 and probabilities
 

       
sq6 = dna_sequence('sq6', 'TTGAGGTCGCCCTCAAAGGTT', 21, 6437.2, 54.4)
sq2 = dna_sequence('sq2', 'TCGCTGTCAAAGTTCGAGGTG', 21, 6477.2, 54.4)
sq3 = dna_sequence('sq3', 'TTCAGGGCGAGGACCATAGAG', 21, 6520.3, 56.3)
sq4 = dna_sequence('sq4', 'AATTCAGGGCGAGGACCATAG', 21, 6504.3, 54.4)
sq5 = dna_sequence('sq5', 'GCCTCAGAGAATCGCCAGTACT', 22, 6704.4, 56.7)
sq1 = dna_sequence('sq1', 'GGGCAGAGGTGTCT', 14, 4359.9, 43.2)
sq7 = dna_sequence('sq7', 'CCTCAAGTGGGACCATCATAA', 21, 6399.2, 52.4)
sq8 = dna_sequence('sq8', 'CGGGGGGGGGGGGA', 14, 6399.2, 52.4)
sq9 = dna_sequence('sq9', 'GACTCGACGGTGATGGGGGGC', 21, 6399.2, 52.4)
sq10= dna_sequence('sq10', 'GGGGGGGGGGGG', 14, 6399.2, 52.4)
if doTEST:
    sequences = [sq8,sq2,sq9]
    
else:   
    # List of sequence objects
    sequences = [sq1, sq2, sq3, sq4, sq5, sq6, sq7]    


def wrong_complement(sequence,verified):
    '''This function returns a reverse complement 
    of a DNA or RNA strand'''
    #verified = verify_na(sequence)
    if verified == "DNA":
        # complement strand
        # seq = sequence.replace("A", "t").replace(
        #     "C", "g").replace("T", "a").replace("G", "c")
        seq = sequence.upper()
        # reverse strand
        seq = seq[::-1]
        return seq
 
    elif verified == "RNA":
       
        # # complement strand
        # seq = sequence.replace("A", "u").replace(
        #     "C", "g").replace("U", "a").replace("G", "c")
        seq = sequence.upper()
         
        # reverse strand
        seq = seq[::-1]
        return seq
    else:
        return False


def verify_na(sequence):
    '''This code verifies if a sequence is a DNA or RNA'''
    # set the input sequence
    seq = set(sequence)
     
    # confirm if its elements is equal to the
    # set of valid DNA bases
    # Use a union method to ensure the sequence is
    # verified if does not contain all the bases
    if seq == {"A", "T", "C", "G"}.union(seq):
        return "DNA"
    elif seq == {"A", "U", "C", "G"}.union(seq):
        return "RNA"
    elif seq == {"A", "C", "G"}.union(seq):
        return "RNA"
    elif seq == {"C", "G"}.union(seq):
        return "RNA"
    else:
        return False

def sum_diagonals(matrix):
    matrix = np.array(matrix)
    rows, cols = matrix.shape
    diagonal_data = {}
    max_sum = float('-inf')
    max_mean = float('-inf')
    max_sum_index = None
    # diagVAL = np.fliplr(check).diagonal()

    # Combine upper and lower diagonals in a single loop
    for i in range(-rows + 1, cols):
        diagonal = np.diagonal(matrix, offset=i)
        diagonal_sum = np.nansum(diagonal)
        diagonal_mean = np.nanmean(diagonal)
        non_nan_count = np.count_nonzero(~np.isnan(diagonal))
        diagonal_data[f"Diagonal {i}"] = (diagonal_sum, non_nan_count)

        # Check if this is the highest sum so far
        if diagonal_sum > max_sum:
            max_sum = diagonal_sum
            max_sum_index = i

        # Check if this is the highest mean so far
        if diagonal_mean > max_mean:
            max_mean = diagonal_sum
            max_mean_index = i

    # return diagonal_data, max_sum_index, max_sum, max_mean, (main_diagonal_sum, main_diagonal_mean), (middle_diagonal_sum, middle_diagonal_mean)
    return diagonal_data, max_sum_index, max_sum, max_mean

   
def dna_to_mrna(inputDNA):
    '''This function reverts the DNA sequences to RNA transcripts'''
    MRNA = ''
    for code in inputDNA:
        if code == 'T':
            MRNA = MRNA + 'U'
        if code == 'C':
            MRNA = MRNA + 'C'
        if code == 'A':
            MRNA = MRNA + 'A'
        if code == 'G':
            MRNA = MRNA + 'G'

    return MRNA


def score_match(subject, query, subject_start, query_start, length):
    score = 0
    # for each base in the match
    for i in range(0,length):
        # first figure out the matching base from both sequences
        subject_base = subject[subject_start + i]
        query_base = query[query_start + i]
        # then adjust the score up or down depending on 
        # whether or not they are the same
        if subject_base == query_base:
            score = score + 1
        else:
            score = score - 1
    return score

   
def try_all_matches(subject, query, score_limit):
    matches = []  # List to store all matches meeting the score limit
    subject_len = len(subject)
    max_start = len(query) - subject_len + 1  # Maximum starting index for a full-length subject match within query
    
    # Slide `subject` across `query` in single increments
    for query_start in range(max_start):
        score = score_match(subject, query, 0, query_start, subject_len)  # Full-length match of `subject`
        
        
        # Check if the score meets or exceeds the threshold
        if score >= score_limit:
            query_end = query_start + subject_len
            matches.append((True, query_start, query_end, subject_len, score))
    
    # Return the list of matches if any were found, otherwise an empty list
    return matches if matches else [(False, None, None, None, None)]

   
def pretty_print_match(subject, query, subject_start, query_start, length):

    # first print the start/stop positions for the subject sequence
    print(str(subject_start) + (' ' * length) + str(subject_start+length))

    # then print the bit of the subject that matches
    print(' ' + subject[subject_start:subject_start+length])

    # then print the bit of the query that matches
    print(' ' + query[query_start:query_start+length])

    # finally print the start/stop positions for the query
    print(str(query_start) + (' ' * length) + str(query_start+length))

    print('\n--------------------\n')


def seqIsolate(matrix, hpLEN):
    # Step 1: Get rid of first and last hairpin rows and columns 
    matrix_step1 = matrix.copy()
    # matrix_step1[:hpLEN, :] = np.nan
    matrix_step1[:, :hpLEN] = np.nan
    # matrix_step1[-hpLEN:, :] = np.nan
    matrix_step1[:, -hpLEN:] = np.nan
    
    # Step 2: Do the opposite - set everything except first and last 7 rows and columns to NaN
    matrix_step2 = matrix.copy()
    matrix_step2[hpLEN:-hpLEN, hpLEN:-hpLEN] = np.nan
    
    return matrix_step1, matrix_step2

def process_results(my_result):
    results = {}
    for my_complex, complex_result in my_result.complexes.items():
        P = complex_result.pairs.to_array()
        s = 'Expected number of unpaired nucleotides at equilibrium in complex %s = %.2f'
        print(s % (my_complex.name, np.diagonal(P).sum()))
        results[my_complex.name] = P
    return results

# Function to find the best matching sequence in MMP9 for a given hairpin
def find_best_match(sequence, target_dna):
    window_size = len(sequence)
    best_match_score = -1
    best_match_segment = None

    # Slide a window over the target DNA to find the best match
    for start in range(len(target_dna) - window_size + 1):
        snippet = target_dna[start:start + window_size]
        match_score = SequenceMatcher(None, sequence, snippet).ratio()
        
        if match_score > best_match_score:
            best_match_score = match_score
            best_match_segment = snippet

    return best_match_segment, best_match_score



# Create the DataFrame dynamically using a dictionary comprehension
dna_seq = pd.DataFrame({
    'Name': [seq.name for seq in sequences],
    'Sequence': [seq.sequence for seq in sequences],
    'Complement': [seq.complement for seq in sequences],
    'Length': [seq.length for seq in sequences],
    'MW': [seq.mw for seq in sequences],
    'Tm(°C)': [seq.tm for seq in sequences]
})

    
print(dna_seq)


# To start with, we would need to double-check if the above sequences have a comlementariety towards the MMP9 RNA sequence
# First we will upload the sequences of the targetted genes.
# Once, we have double checked that the given sequences do indeed have a complement sequence towards the MMP9 gene
# We will extract and save the correponding sequences +,- 15bp to then perform the hybridization analysis
# The Nupack library used for the complex analysis of hybridization works with sequences up to a certain range, hence the need to split the MMP9 into the sequences which perform the binding

mmp9_rna = open(dataDIR+'rna_mmp9.txt', 'r').read()
mmp9_rna = mmp9_rna.replace('\n', '')
mmp9_dna = mmp9_rna.replace('U', 'T')
#need to convert mmp9_rna to dna


###############################################################################
# MIKES CHANGES
# Initialise lists to collect data for each sequence and its matches
dataMATCH = []

# Iterate through each sequence in `dna_seq['Complement']`
for seq in dna_seq['Complement']:
    # Get the corresponding row from dna_seq for this `seq`
    dna_row = dna_seq[dna_seq['Complement'] == seq].iloc[0]
    
    # Call try_all_matches with a score limit of at least half the length of hairpin
    matches = try_all_matches(seq, mmp9_dna, scoreLim)
    
    # Check if any matches were found
    if matches[0][0]:  # If there's at least one valid match
        hyb_col = True
        
        # Process each match
        for match in matches:
            _, query_start, query_end, subject_len, score = match
            match_check = mmp9_dna[query_start:query_end]
            
            # Calculate start and end positions for extracted slice with basePAD buffer
            if query_start <= basePAD:
                q_st = 0
                q_ed = query_end + basePAD
            elif query_end >= len(mmp9_dna) - basePAD:
                q_st = query_start - basePAD
                q_ed = query_end + abs(len(mmp9_dna) - query_end)
            else:
                q_st = query_start - basePAD
                q_ed = query_end + basePAD
            
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
        # If no matches met the score limit, select a random segment from mmp9_dna
        random_start = random.randint(0, len(mmp9_dna) - len(seq) - basePAD*2)
        random_end = random_start + len(seq) + basePAD*2
        random_segment = mmp9_dna[random_start:random_end]
        
        # Append data for the unmatched sequence, including a random segment of mmp9_dna
        dataMATCH.append({
            'hyb_col': False,
            'mmp9_seqs': random_segment,  # Randomly selected segment
            'search_seqs': seq,
            'query_start': random_start,
            'query_end': random_end,
            'subject_len': len(seq),
            'score': 0,  # Default score for random segment
            'Name': dna_row['Name'],
            'Sequence': dna_row['Sequence'],
            'Complement': dna_row['Complement'],
            'Length': dna_row['Length'],
            'MW': dna_row['MW'],
            'Tm(°C)': dna_row['Tm(°C)']
        })

# Convert collected data to DataFrame
dfMATCH = pd.DataFrame(dataMATCH)



# Save DataFrame to Excel file
dfMATCH.to_excel(results_dir+'rn_seqDATA.xlsx', index=True)
# Now we can design and generate RNA sequences containing the complementary sequence to the MMP9, plus two complementary sequences with 3 GC sites to promote Hairpin formation.
# After the design, we can perform the modelling to analyse MFEs and pair probabilities of every generated sequence

# Define the nucleotides
nucleotides = ['A', 'C', 'G', 'T']
reverse_complementary_hairpins =[]
seq_complement=''

# Define panda where sequences will be stored
hairpin_pd = pd.DataFrame()

# Generate all possible combinations of nucleotides of length 7 as this is the hairpin length with over 1 possibility with 2,4 'GC' content

sequences = list(itertools.product(nucleotides, repeat=hpLEN))
# Filter for sequences that contain 3 (or 2) 'GC' nucleotides
filtered_sequences = [seq for seq in sequences if ''.join(seq).count('GC') == 2 or ''.join(seq).count('GC') == 4] #WHY?


# Convert list of tuples to list of strings
list_of_strings = [''.join(t) for t in filtered_sequences]

if doTEST:
    list_of_strings = [list_of_strings[0]]   
    

# Get a list of the reverse complemente haipin sequences
for seq in list_of_strings:
    if WrongHP:
        seq_complement = wrong_complement(seq,'DNA') 
    else:
        seq_complement = _complement(seq,'DNA') #for MIKE: printing in RNA should be DNA possibly
    reverse_complementary_hairpins.append(seq_complement)
    seq_complement=''
    
    
    
# What we want to do now is test whether the hairpin or the hairpin's complement hybridises with MMP9 - this is bad and will help filter out further
 # List to store results for each hairpin-complement combination
# List to store results for each hairpin-complement combination
# Initialize the results list
results = []

# Total number of outer loops (hairpin-complement pairs)


# Correct total_pairs calculation
total_pairs = len(list_of_strings)

print(f"Total pairs: {total_pairs}")



# Loop over each hairpin-complement pair
for hairpin, complement in tqdm(zip(list_of_strings, reverse_complementary_hairpins), total=total_pairs, desc="Hairpin-Complement Pairs Processed"):
    # --- Find best matching segment in MMP9 for hairpin ---
    best_hairpin_segment, bestHPscore = find_best_match(hairpin[::-1], mmp9_dna) # because we're matching we use the complement for the match
    
    # Define strands for hairpin and best matching MMP9 segment
    hp_st = Strand(hairpin, name='hp_st')
    mmp9_st_hairpin = Strand(best_hairpin_segment, name='mmp9_st_hairpin')
    
    # Define complex and set for hairpin
    hyb_complex_hairpin = Complex([hp_st, mmp9_st_hairpin])  
    hyb_set_hairpin = ComplexSet(strands={hp_st: concentration, mmp9_st_hairpin: concentration},
                                 complexes=SetSpec(max_size=0, include=[hyb_complex_hairpin]))
    
    # Perform analysis on the best match for hairpin
    hyb_result_hairpin = complex_analysis(hyb_set_hairpin, compute=['pfunc', 'pairs', 'mfe'],
                                          options={'energy_gap': gap}, model=mmp9_model)
    structure_hairpin = hyb_result_hairpin[hyb_complex_hairpin].mfe[0].structure
    max_prob_hairpin = structure_probability([hp_st, mmp9_st_hairpin], structure_hairpin, mmp9_model)

    # --- Find best matching segment in MMP9 for complement ---
    best_complement_segment, bestCscore = find_best_match(complement[::-1], mmp9_dna)
    
    # Define strands for complement and best matching MMP9 segment
    comp_st = Strand(complement, name='comp_st')
    mmp9_st_complement = Strand(best_complement_segment, name='mmp9_st_complement')
    
    # Define complex and set for complement
    hyb_complex_complement = Complex([comp_st, mmp9_st_complement])
    hyb_set_complement = ComplexSet(strands={comp_st: concentration, mmp9_st_complement: concentration},
                                    complexes=SetSpec(max_size=0, include=[hyb_complex_complement]))
    
    # Perform analysis on the best match for complement
    hyb_result_complement = complex_analysis(hyb_set_complement, compute=['pfunc', 'pairs', 'mfe'],
                                             options={'energy_gap': gap}, model=mmp9_model)
    structure_complement = hyb_result_complement[hyb_complex_complement].mfe[0].structure
    max_prob_complement = structure_probability([comp_st, mmp9_st_complement], structure_complement, mmp9_model)

    # Store hairpin, complement, and both Max Probabilities for plotting and filtering
    results.append((hairpin, complement, max_prob_hairpin, max_prob_complement))

# Convert results to a DataFrame for easier manipulation
results_df = pd.DataFrame(results, columns=['Hairpin', 'Complement', 'Max Prob Hairpin', 'Max Prob Complement'])

# # Set threshold based on the lowest percentage of hairpin and complement probabilities
# threshold_hairpin = results_df['Max Prob Hairpin'].quantile(lowest_probability)
# threshold_complement = results_df['Max Prob Complement'].quantile(lowest_probability)

# Filter results where both hairpin and complement probabilities are below their respective thresholds

if doTEST == 1:
    filtered_results = results_df
else:
    filtered_results = results_df[(results_df['Max Prob Hairpin'] <= lowest_probability) &
                              (results_df['Max Prob Complement'] <= lowest_probability)]

# Create the plot
plt.figure(figsize=(12, 8))
plt.scatter(results_df['Hairpin'], results_df['Max Prob Hairpin'], color='blue', alpha=0.6, label='Hairpin')
plt.scatter(results_df['Hairpin'], results_df['Max Prob Complement'], color='orange', alpha=0.6, label='Complement')
plt.xlabel("Hairpin Sequence")
plt.ylabel("Max Probability")
plt.title("Max Probability for Each Hairpin and Complement Combination")
plt.xticks(rotation=90)

# Highlight the filtered combinations
filtered_hairpins = filtered_results['Hairpin']
filtered_probs_hairpin = filtered_results['Max Prob Hairpin']
filtered_probs_complement = filtered_results['Max Prob Complement']
plt.scatter(filtered_hairpins, filtered_probs_hairpin, color='darkblue', label='Filtered Hairpin', edgecolor='k')
plt.scatter(filtered_hairpins, filtered_probs_complement, color='darkorange', label='Filtered Complement', edgecolor='k')
# Add the threshold line
plt.axhline(y=lowest_probability, color='red', linestyle='--', label=f'Threshold ({lowest_probability})')
plt.legend()
plt.tight_layout()

# Save the plot as both PNG and SVG
plt.savefig(results_dir + "max_probability_hairpin_complement.png", bbox_inches='tight', pad_inches=0)
plt.savefig(results_dir + "max_probability_hairpin_complement.svg", bbox_inches='tight', pad_inches=0)

# Display the plot
plt.show()

# now we need to only select the hairpins that have made it through threshold


    
# We construct the hairpin sequences by adding the hairpin complementary bp to each end of the siRNA sequence
hairpins = []
seq_name = []
hairpin_length = []
hairpin_mw = []
hairpin_mt = []
mmp9_snippets = []

# Initialize a set to track unique sequences
unique_sequences = set()

for sequence in dfMATCH['Complement']:
    # Check if this sequence has already been processed
    if sequence not in unique_sequences:
        unique_sequences.add(sequence)  # Add the unique sequence to the set
        
        # Extract all rows matching the current `sequence` and get `MMP9_Seq` values as a list
        matching_rows = dfMATCH[dfMATCH['Complement'] == sequence]
        mmp9_snippets_for_sequence = matching_rows['mmp9_seqs'].tolist()  # List of all MMP9_Seq values for the sequence
        
        # Fetch the `name` from the first matching row
        name = matching_rows.iloc[0]['Name']
        
        # Process each `mmp9_snippet` along with the `seq` and `rev` pair
        for mmp9_snippet in mmp9_snippets_for_sequence:
            for seq, rev in zip(filtered_results['Hairpin'], filtered_results['Complement']):
                hairpin = seq + _complement(sequence,'DNA') + rev  # Construct the hairpin sequence (the original is the complement of the complement :))
                hairpins.append(hairpin)
                seq_name.append(name)
                hairpin_length.append(len(hairpin))
                hairpin_mw.append(molecular_weight(hairpin, 'DNA'))
                hairpin_mt.append(mt.Tm_Wallace(hairpin))
                mmp9_snippets.append(mmp9_snippet)  # Append each mmp9_snippet per hairpin
            
            
hairpin_pd['Name'] = seq_name
hairpin_pd['Hairpin'] = hairpins
hairpin_pd['Length'] = hairpin_length
hairpin_pd['MW'] = hairpin_mw
hairpin_pd['Tm(°C)'] = hairpin_mt
hairpin_pd['MMP9_Seq'] = mmp9_snippets # This column contains the fragement of MMP9 sequence the hairpin should hybridise with
hairpin_pd['MMP9_Seq_length'] = hairpin_pd['MMP9_Seq'].apply(len)
# print(hairpin_pd)


# At this stage, we have constructed a database with all the possible harpin sequences which can hybridise with MMP9
# Now, we need to predict their looping free energies, as well as, the free energies of the hybridisation process between hairpin and MMP9 to see if the process would be spontaneuos and promoted



#a = mmp9_model.alphabet()

# Calculate hairpin sequences stability analysis and complex(hairpin+MMP9) sequences stability analysis
# Define lists to store dG and MFE values and bp probabilites
hairpin_mfe = []
hairpin_dg = []
hairpin_p = []
hairpin_ep = []

hyb_mfe = []
hyb_dg = []
hyb_p = []
hyb_ep = []
hTCK=0
# hairpin_pd = hairpin_pd.dropna(subset=['MMP9_Seq'])





# Define the iterable based on the flag
if doFullGene == 1:
    iterable = hairpin_pd['Hairpin']
    desc = "Processing Hairpins and full MMP9"
else:
    
    iterable = zip(hairpin_pd['Hairpin'], hairpin_pd['MMP9_Seq'])
    desc = "Processing Hairpin-Gene Pairs"

# Iterate over the dynamically defined iterable with the shared code
hTCK = 0
for item in tqdm(iterable, total=len(hairpin_pd), desc=desc):
    if doFullGene == 1:
        hairpin = item
        gene = mmp9_dna  # Use None or a default value for gene
    else:
        hairpin, gene = item       

    hTCK += 1
    print(f"{hTCK} {hairpin}")
    
    # Define the strands
    hp_st = Strand(hairpin, name='hp_st')
    mmp9_st = Strand(gene, name='mmp9_st') if gene else None
    
    # Define the complex of interest
    hp_complex = Complex([hp_st])  # Secondary structure of hairpin folding
    hyb_complex = Complex([hp_st, mmp9_st]) if gene else None
    
    # Define the complex set
    hp_set = ComplexSet(strands={hp_st: 1e-6}, complexes=SetSpec(max_size=0, include=[hp_complex]))
    hyb_set = (ComplexSet(strands={hp_st: 1e-6, mmp9_st: 1e-6}, 
                         complexes=SetSpec(max_size=0, include=[hyb_complex])) if gene else None)
    
    # Analyze the complex
    hp_result = complex_analysis(hp_set, compute=['pfunc', 'pairs', 'mfe'],
                                 options={'energy_gap': gap}, model=mmp9_model)
    
    if gene:
        hyb_result = complex_analysis(hyb_set, compute=['pfunc', 'pairs', 'mfe'],
                                      options={'energy_gap': gap}, model=mmp9_model)

    # Store results for hairpin
    hairpin_result = hp_result[hp_complex]
    hairpin_dg.append(hairpin_result.free_energy)
    hairpin_mfe.append(hairpin_result.mfe[0][1])
    hairpin_p.append(hairpin_result.pairs.to_array())
    structure = hairpin_result.mfe[0].structure
    hairpin_ep.append(structure_probability([hp_st], structure, mmp9_model))
    
    # Store results for hybridization (if applicable)
    if gene:
        hyb_result = hyb_result[hyb_complex]
        hyb_dg.append(hyb_result.free_energy)
        hyb_mfe.append(hyb_result.mfe[0][1])
        
        temp = hyb_result.pairs.to_array()
        if doTEST:
            plt.figure(figsize=(10, 8))
            sns.heatmap(temp, cmap="viridis", cbar=True, square=True)
            plt.xlabel("Base Index", fontsize=12)
            plt.ylabel("Base Index", fontsize=12)
            plt.title(f"Hairpin: {hairpin} MMP9: {gene}", fontsize=14, fontweight='bold')
            plt.savefig(f"{results_dir}{hairpin}_probability_matrix.png", bbox_inches="tight", pad_inches=0.1)
            plt.savefig(f"{results_dir}{hairpin}_probability_matrix.svg", bbox_inches="tight", pad_inches=0.1)
            plt.show()
            plt.close()
        
        # Set lower triangle to NaN and process
        lower_triangle_indices = np.tril_indices_from(temp, k=-1)
        temp[lower_triangle_indices] = np.nan
        np.fill_diagonal(temp, 0)
        hyb_p.append(temp)
        structure = hyb_result.mfe[0].structure
        hyb_ep.append(structure_probability([hp_st, mmp9_st], structure, mmp9_model))


# # Iterate over the hairpin and MMP9 sequence with a progress bar
# for hairpin, gene in tqdm(zip(hairpin_pd['Hairpin'], hairpin_pd['MMP9_Seq']),
#                           total=len(hairpin_pd),
#                           desc="Processing Hairpin-Gene Pairs"):
# for hairpin in hairpin_pd['Hairpin']:
#     hTCK=hTCK+1
#     # print(str(hTCK) + ' ' + hairpin + ' ' + str(gene))
#     print(str(hTCK) + ' ' + hairpin)
#     # Define the strands
    
#     hp_st = Strand(hairpin, name='hp_st')
#     mmp9_st = Strand(gene, name='mmp9_st')
#     # Define the complex of interest
#     hp_complex = Complex([hp_st]) # Secondary structure of hairpin folding
   
#     hyb_complex = Complex([hp_st, mmp9_st]) # Hybridisation between hairpin and gene of interest
#     #not_incl_complex = Complex([mmp9_st,hp_st]) # Hybridisation between hairpin and gene of interest
#     # Define the complex set to contain only one complex for both cases
#     hp_set = ComplexSet(strands={hp_st: 1e-6}, complexes=SetSpec(max_size=0, include=[hp_complex]))
#     hyb_set = ComplexSet(strands={hp_st: 1e-6, mmp9_st: 1e-6}, complexes=SetSpec(max_size=0, include=[hyb_complex]))
    
 
#     # Analyze the complex 
#     # Calculate pfunc, pairs, mfe, subopt
#     # Analyze the complex 
#    # Calculate pfunc, pairs, mfe, subopt
#     hp_result = complex_analysis(hp_set, compute=['pfunc', 'pairs', 'mfe'],
#                              options={'energy_gap': gap}, model=mmp9_model)
#     hyb_result = complex_analysis(hyb_set, compute=['pfunc', 'pairs', 'mfe'], 
#                              options={'energy_gap': gap}, model=mmp9_model)
#     # Store predicted dG and MFE values and bp probabilites
#     hairpin_result = hp_result[hp_complex]
#     hyb_result = hyb_result[hyb_complex]
#     hairpin_dg.append(hairpin_result.free_energy)
#     hairpin_mfe.append(hairpin_result.mfe[0][1])
#     temp=hairpin_result.pairs.to_array()
   
        
#     # temp[lower_triangle_indices] = np.nan
#     # np.fill_diagonal(temp,0)
#     hairpin_p.append(temp)
#     # print(hairpin,temp.shape)
#     # lets get the equilibrium probability of the secondary structure
#     # Get the secondary structure
#     structure = hairpin_result.mfe[0].structure
#     # Calculate the probability
#     eqP = structure_probability([hp_st], structure, mmp9_model)
#     hairpin_ep.append(eqP)
    
      
       
#     hyb_dg.append(hyb_result.free_energy)
#     hyb_mfe.append(hyb_result.mfe[0][1])
#     temp=hyb_result.pairs.to_array()
    
#     if doTEST:
#         # Create the figure and axis
#         plt.figure(figsize=(10, 8))
        
#         # Plot the heatmap
#         # g=sns.heatmap(
#         #     temp, 
#         #     cmap="viridis", 
#         #     cbar=True, 
#         #     square=True,
#         #     xticklabels=list(hairpin) + list(hairpin)[::-1],  # Use slicing for reverse
#         #     yticklabels=list(gene) + list(gene)[::-1],        # Use slicing for reverse
#         # )
#         g=sns.heatmap(
#             temp, 
#             cmap="viridis", 
#             cbar=True, 
#             square=True,           
#         )
        
#         g.set_yticklabels(g.get_xticklabels(), rotation = 90, fontsize = 8)
#         g.set_yticklabels(g.get_yticklabels(), rotation = 0, fontsize = 8)

#         # Add labels and title
#         plt.xlabel("Base Index", fontsize=12)
#         plt.ylabel("Base Index", fontsize=12)
#         plt.title("Hairpin: "+hairpin +' MMP9: '+gene, fontsize=14, fontweight='bold')
        
#         # Save the figure as PNG and SVG
#         plt.savefig(results_dir + hairpin + "_probability_matrix.png", bbox_inches="tight", pad_inches=0.1)
#         plt.savefig(results_dir + hairpin + "_probability_matrix.svg", bbox_inches="tight", pad_inches=0.1)
#         plt.show()
#         plt.close()
        
#     #Set the lower triangle to NaN
#     lower_triangle_indices = np.tril_indices_from(temp, k=-1)
#     temp[lower_triangle_indices] = np.nan
#     np.fill_diagonal(temp,0)
#     hyb_p.append(temp)    
#     # lets get the equilibrium probability of the secondary structure
#     # Get the secondary structure
#     structure = hyb_result.mfe[0].structure
#     # Calculate the probability
#     eqP = structure_probability([hp_st,mmp9_st], structure, mmp9_model)
#     hyb_ep.append(eqP)
     
       
        

#MFE is how much energy needed to break apart
hairpin_pd['Hairpin_dG (kcal/mol)'] = hairpin_dg
hairpin_pd['Hairpin_MFE (kcal/mol)'] = hairpin_mfe
hairpin_pd['Hairpin_Pb'] = hairpin_p
hairpin_pd['Hairpin_equiPb'] = hairpin_ep

hairpin_pd['Hyb_dG (kcal/mol)'] = hyb_dg
hairpin_pd['Hyb_MFE (kcal/mol)'] = hyb_mfe
hairpin_pd['Hyb_Pb'] = hyb_p
hairpin_pd['Hyb_equiPb'] = hyb_ep

hairpin_pd.to_excel(results_dir+'hairpin_pdDATA_preCLEAN.xlsx', index=True)
# Now, we need to clean the rows in which either the MFE of hybridisation or hairpin are positive
# We also need to clean the rows in which the MFE for hybridisation is higher than the hairpin one - as once inside the cell they will not open to bind to the mRNA of MMP9
hairpin_pd = hairpin_pd[(hairpin_pd['Hairpin_MFE (kcal/mol)'] < 0) & (hairpin_pd['Hyb_MFE (kcal/mol)'] < 0) & (hairpin_pd['Hairpin_MFE (kcal/mol)'] > hairpin_pd['Hyb_MFE (kcal/mol)'])]
# print(hairpin_pd)

# Save DataFrame to Excel file
hairpin_pd.to_excel(results_dir+'hairpin_pdDATA.xlsx', index=True)
# Data visualisation and statistics to choose between prefered sequences
# Plot histogram with the distribution of energies

# fig, ax = plt.subplots()
# fig.set_size_inches(13, 6)
# ax.set(xlabel='MFE (kcal/mol)')
# sns.set_style("darkgrid")
# sns.set(color_codes=True)
# sns.set(palette="muted")
# sns.histplot([hairpin_pd['Hairpin_MFE (kcal/mol)'], hairpin_pd['Hyb_MFE (kcal/mol)']], 
#                 stat="density", bins=250, common_norm=False, kde=True, fill=True, ax=ax)
# plt.show()
# plt.savefig(results_dir + "hairpin_energies.svg")
# plt.savefig(results_dir + "hairpin_energies.png")
# Data visualisation and statistics to choose between prefered sequences
# Plot histogram with the distribution of energies and the difference between them

diff_mfe = np.subtract(hairpin_pd['Hyb_MFE (kcal/mol)'], hairpin_pd['Hairpin_MFE (kcal/mol)']) # Calculate the difference between complexes 
hairpin_pd['MFE_Diff (kcal/mol)'] = diff_mfe
if doTEST:
    hairpin_pd['MFE_Diff_Derivative'] = np.nan # Append the list of differences into the dataframe

else:    
    hairpin_pd['MFE_Diff_Derivative'] = np.gradient(diff_mfe)+np.min(hairpin_pd['Hairpin_MFE (kcal/mol)']) # Append the list of differences into the dataframe

# fig, ax = plt.subplots()
# fig.set_size_inches(13, 6)
# ax.set(xlabel='MFE (kcal/mol)')
# sns.set_style("darkgrid")
# sns.set(color_codes=True)
# sns.set(palette="muted")
# sns.histplot([hairpin_pd['Hairpin_MFE (kcal/mol)'], hairpin_pd['Hyb_MFE (kcal/mol)'], hairpin_pd['MFE_Diff_Derivative']], 
#              stat="density", bins=250, common_norm=False, kde=True, fill=True, ax=ax)

# plt.show()
# plt.savefig(results_dir + "hairpin_energiesSubtraction.svg")
# plt.savefig(results_dir + "hairpin_energiesSubtraction.png")
# Let's remove the rows where the values for both Hairpin and Hybridisation complex match, as well as those lying 1-std around that point
# Take the difference between distributions and build a a new distribution with those values, then calculate the statistics and dicriminate accordingly

diff_mean = np.mean(diff_mfe)
diff_std = np.std(diff_mfe)
#sig_diff = abs(diff_mean) + 2.3*diff_std
#MIKES addition: Nuria says that any difference >5 kcal/mol we keep
sig_diff = 5
# Remove rows from the database where the difference between the hairpin and hybridised form not on the range mean+std
hairpin_pd = hairpin_pd[((hairpin_pd['Hairpin_MFE (kcal/mol)'] - hairpin_pd['Hyb_MFE (kcal/mol)']) > sig_diff)]
# print(hairpin_pd)



# # Plot a histogram to see the difference in energy from the cleaned dataset and check that the overlapping section has been removed

# fig, ax = plt.subplots()
# fig.set_size_inches(13, 6)
# ax.set(xlabel='MFE (kcal/mol)')
# sns.set_style("darkgrid")
# sns.set(color_codes=True)
# sns.set(palette="muted")
# sns.histplot([hairpin_pd['Hairpin_MFE (kcal/mol)'], hairpin_pd['Hyb_MFE (kcal/mol)']], 
#              stat="density", shrink=2, bins=250, fill=True, kde=True, ax=ax)
# plt.show()
# plt.savefig(results_dir + "hairpin_energiesOverlap.svg")
# plt.savefig(results_dir + "hairpin_energiesOverlap.png")
# Data visualisation and statistics to choose between prefered sequences
# Plot swarm plot with the distribution of energies and the difference between them
# To achieve this, we need to arrange the data into a preliminary dataset so seaborn can read it

hairpin_pd_plot = pd.DataFrame() # Define panda and input the complementary data
variable = ['Hairpin MFE (kcal/mol)']*len(hairpin_pd['Hairpin_MFE (kcal/mol)']) + ['Hyb MFE (kcal/mol)']*len(hairpin_pd['Hairpin_MFE (kcal/mol)']) + ['MFE Diff. (kcal/mol)']*len(hairpin_pd['Hairpin_MFE (kcal/mol)'])

values = (hairpin_pd['Hairpin_MFE (kcal/mol)'].to_list() + hairpin_pd['Hyb_MFE (kcal/mol)'].to_list() + hairpin_pd['MFE_Diff (kcal/mol)'].to_list())
hairpin_pd_plot[' '] = variable
hairpin_pd_plot['MFE (kcal/mol)'] = values

# fig, ax = plt.subplots() # plot the data to see the differences in MFE between hybridisation and hairpin structures
# fig.set_size_inches(20, 10)
# sns.swarmplot(x =' ' , y ='MFE (kcal/mol)', hue=' ', data = hairpin_pd_plot, size=2, palette='deep')
# ax.set_title('MFE Distributions for the Hairpin, the Hybridisation Complex and their Difference ', weight='bold', fontsize=12)
# plt.show()
# plt.savefig(results_dir + "hairpin_energiesSubtractionSwarm.svg")
# plt.savefig(results_dir + "hairpin_energiesSubtractionSwarm.png")
# # Plot a swarm to see the difference in energy from the cleaned dataset and check that the overlapping section has been removed
# # To achieve this, we need to arrange the data into a preliminary dataset so seaborn can read it

hairpin_pd_plot = pd.DataFrame() # Panda definition and data input
variable = ['Hairpin MFE (kcal/mol)']*len(hairpin_pd['Hairpin_MFE (kcal/mol)']) + ['Hyb MFE (kcal/mol)']*len(hairpin_pd['Hairpin_MFE (kcal/mol)'])

values = (hairpin_pd['Hairpin_MFE (kcal/mol)'].to_list() + hairpin_pd['Hyb_MFE (kcal/mol)'].to_list())
hairpin_pd_plot[' '] = variable
hairpin_pd_plot['MFE (kcal/mol)'] = values

# fig, ax = plt.subplots() # plot the data
# fig.set_size_inches(13, 6)
# sns.swarmplot(x =' ' , y ='MFE (kcal/mol)', hue=' ', data = hairpin_pd_plot, size=6, palette='deep', linewidth = 0.8)
# ax.set_title('MFE Distributions for the Hairpin and the Hybridisation Complex after energy overlap removal', weight='bold', fontsize=10)
# plt.show()
# plt.savefig(results_dir + "hairpin_energiesCleanSwarm.svg")
# plt.savefig(results_dir + "hairpin_energiesCleanSwarm.png")
# #At this stage, we can plot again both MFE values to see their difference

# # Create a figure and axis with a specific size
# fig, ax = plt.subplots(figsize=(12.8, 6))
# width1 = 0.3
ar = np.arange(len(hairpin_pd['Hairpin'])) # Range values for the set of Hairpin Sequences
# # We plot the data from the dataset for both the Hairpin and the Hybrid complex predicted energies
# ax.bar(ar, abs(hairpin_pd['Hairpin_MFE (kcal/mol)']), width=0.3, align='edge', label='Hairpin_MFE (kcal/mol)', color=sns.xkcd_rgb['windows blue'], alpha=0.5)
# ax.bar(ar + width1, abs(hairpin_pd['Hyb_MFE (kcal/mol)']), width=0.3, align='edge', label='Hyb_MFE (kcal/mol)', color=sns.xkcd_rgb['orangish'], alpha=0.6)
# ax.set_ylabel('ABS MFE (kcal/mol)', fontsize=14)
# ax.set_title('Predicted ABS values of MFE for both Hairpin and Complex Structures', fontsize=10, weight='bold')
# ax.set_xticks(ar, hairpin_pd['Hairpin'], rotation=70, fontweight='bold', fontsize='5', horizontalalignment='right', minor=False)
# ax.legend()
# plt.show()
# plt.savefig(results_dir + "predicted_energiesClean.svg")
# plt.savefig(results_dir + "predicted_energiesClean.png")

# Plot the probability matrices for visualisation. Hairpin sequences
hairpin_pd.to_excel(results_dir+'hairpin_pdDATAcleaned.xlsx', index=True)
if doPLOT:
    num_rows = len(hairpin_pd)
    num_subplots = int(np.ceil(np.sqrt(num_rows)))
    fig, axs = plt.subplots(nrows=num_subplots, ncols=num_subplots, figsize=(25, 55))
    fig.subplots_adjust(wspace=0.5)
    for ax, pb, hairpin in zip(axs.flat, hairpin_pd['Hairpin_Pb'], hairpin_pd['Hairpin']):
        im = ax.imshow(pb, interpolation='nearest', cmap='viridis', aspect='auto', extent=[0,len(hairpin),len(hairpin),0])
        #lets try find best hairpins
        #these need strong connection at ends and lowest probability at other places
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cax.xaxis.set_ticks_position('default')
        cbar = fig.colorbar(im, cax=cax, orientation='vertical', ticks=[0.2, 0.4, 0.6, 0.8])
        cbar.ax.tick_params(labelsize=8)
        
        ax.set_title('Hairpin Sequence: ' + hairpin, fontsize=7, weight='bold')
        ax.set_xlabel('Base index', fontsize=9)
        ax.set_ylabel('Base index', fontsize=9)
        ax.grid(False)
        ax.xaxis.set_tick_params(labelsize=9)
        ax.yaxis.set_tick_params(labelsize=9)
        ax.grid(False)
        
        # Set x-tick labels as hairpin letters
        ax.set_xticks(range(len(hairpin)))
        ax.set_xticklabels(hairpin, fontsize=9)
        
        # Set y-tick labels as hairpin letters
        ax.set_yticks(range(len(hairpin)))
        ax.set_yticklabels(hairpin, fontsize=9)
     
    
    plt.tight_layout()
    plt.show()
    plt.savefig(results_dir + "hairpinProbability.svg")
    plt.savefig(results_dir + "hairpinProbability.png")


if doPLOT:
    # Plot the probability matrices for visualisation. Hybrid complexes
    num_rows = len(hairpin_pd)
    num_subplots = int(np.ceil(np.sqrt(num_rows)))
    fig, axs = plt.subplots(nrows=num_subplots, ncols=num_subplots, figsize=(25, 55))
    fig, axs = plt.subplots(nrows=11, ncols=5, figsize=(25, 55))
    fig.subplots_adjust(wspace=0.5)
    tck=0
    for ax, pb, hairpin, mmp9 in zip(axs.flat, hairpin_pd['Hyb_Pb'], hairpin_pd['Hairpin'], hairpin_pd['MMP9_Seq']):
        tck=tck+1
        im = ax.imshow(pb, interpolation='nearest', cmap='viridis', aspect='auto', extent=[0,len(hairpin),len(mmp9),0])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cax.xaxis.set_ticks_position('default')
        cbar = fig.colorbar(im, cax=cax, orientation='vertical', ticks=[0.2, 0.4, 0.6, 0.8])
        cbar.ax.tick_params(labelsize=8)
        
        ax.set_title('Hairpin Sequence: ' + hairpin, fontsize=7, weight='bold')
        ax.set_xlabel('Base index', fontsize=9)
        ax.set_ylabel('Base index', fontsize=9)
        ax.xaxis.set_tick_params(labelsize=9)
        ax.yaxis.set_tick_params(labelsize=9)
        ax.grid(False)
        
        # Set x-tick labels as hairpin letters
        ax.set_xticks(range(len(hairpin)))
        ax.set_xticklabels(hairpin, fontsize=9)
        
        # Set y-tick labels as hairpin letters
        ax.set_yticks(range(len(mmp9)))
        ax.set_yticklabels(mmp9, fontsize=9)
        
          
    print(tck)
    plt.show()
    plt.tight_layout()
    plt.show()
    plt.savefig(results_dir + "hybridisationProbability.svg")
    plt.savefig(results_dir + "hybridisationProbability.png")
# To clean the data further, we can add the diagonal values from the probability base pair matrix to assess which sequences will have a higher predicted 
# base pairing and hence will hybridise and group into a hairpin at higher rates while mitigating cross-talk


hp_pb_meanDIAG = []
max_sum_hp_Pb = []
non_nan_count_hp_Pb = []
mean_hp_Pb = []
max_sum_hyb_Pb = []
non_nan_count_hyb_Pb = []
mean_hyb_Pb = []
mean_hyb_Pb_weighted = []
hyb_Pb_meanDIAG = []
D_Pb_Hairpin_D_Pb_Hyb_mean_DIAG = []

# Calculate diagonal sum of probability matrix values
for pb_hair, pb_hyb in zip(hairpin_pd['Hairpin_Pb'], hairpin_pd['Hyb_Pb']):
    pb_ha, pb_hy = np.flipud(np.asarray(pb_hair)), np.flipud(np.asarray(pb_hyb))
    
    # first want to get the hairpins that definitely occur
    # this means you need to look at first 7 elements of the diag and get mean and sort descending
    # then need to get 
    temp=(pb_ha)
    innerV,outerV=seqIsolate(temp,hpLEN)
    # this is from diagonal_data[f"Diagonal {i}"] = (diagonal_sum, non_nan_count)
    # which means that diagonal_data[0] is the sum of a diagonal
    # this means that the diagonals with largest number of matches will be rewarded
    # and diagonal_data[1] is how many not nans just so we have the correct length to then normalise by
    diagonal_data, max_sum_index, max_sum, max_mean = sum_diagonals(innerV) 
    max_sum_hp_Pb.append(diagonal_data[f'Diagonal {max_sum_index}'][0])
    non_nan_count_hp_Pb.append(diagonal_data[f'Diagonal {max_sum_index}'][1])
    mean_hp_Pb = diagonal_data[f'Diagonal {max_sum_index}'][0]/diagonal_data[f'Diagonal {max_sum_index}'][1]
    temp=np.nanmean(np.diag(outerV))      
    hp_pb_meanDIAG.append(temp)
    
    temp=(pb_hy)
    innerV,outerV=seqIsolate(temp,hpLEN)
    diagonal_data, max_sum_index, max_sum, max_mean = sum_diagonals(innerV) 
    max_sum_hyb_Pb.append(diagonal_data[f'Diagonal {max_sum_index}'][0])
    non_nan_count_hyb_Pb.append(diagonal_data[f'Diagonal {max_sum_index}'][1])
    mean_hyb_Pb_weighted.append(diagonal_data[f'Diagonal {max_sum_index}'][0]/diagonal_data[f'Diagonal {max_sum_index}'][1])
    mean_hyb_Pb.append(max_mean)
    hyb_Pb_meanDIAG.append(np.nanmean(np.diag(innerV)))
   
    
    
    # now let's work out best way to identify hairpin seq combo for mmp9
    
       
# Add the calculated values in our cleaned database
hairpin_pd.insert(14, 'D_Pb_Hairpin_mean_DIAG', hp_pb_meanDIAG, True)
hairpin_pd.insert(15, 'D_Pb_Hairpin_max_sum', max_sum_hp_Pb, True)
hairpin_pd.insert(16, 'D_Pb_Hairpin_non_nan_count', non_nan_count_hp_Pb, True)
hairpin_pd.insert(17, 'D_Pb_Hyb_max_sum', max_sum_hyb_Pb, True)
hairpin_pd.insert(18, 'D_Pb_Hyb_non_nan_count', non_nan_count_hyb_Pb, True)
hairpin_pd.insert(19, 'D_Pb_Hairpin_mean', mean_hp_Pb, True)
hairpin_pd.insert(20, 'D_Pb_Hyb_mean', mean_hyb_Pb, True)
hairpin_pd.insert(21, 'D_Pb_Hyb_Weighted_Mean', mean_hyb_Pb_weighted, True)
hairpin_pd.insert(22, 'D_Pb_Hyb_mean_DIAG', hyb_Pb_meanDIAG, True)
hairpin_pd.insert(23, 'D_Pb_Hairpin_D_Pb_Hyb_mean_DIAG', np.array(hyb_Pb_meanDIAG)/np.array(hp_pb_meanDIAG), True)


# Let's sort each sequence by the best hairpins

seqIDX = hairpin_pd['Name'].unique()
plt.close('all')
#weighted mean by number of elements in diagonal
for seqSTR in seqIDX:
    temp = hairpin_pd[hairpin_pd['Name'].str.contains(seqSTR)].reset_index(drop=True)
    sequence = globals().get(seqSTR).sequence  # or use locals().get(seqSTR) if needed
    temp.sort_values(by='D_Pb_Hairpin_D_Pb_Hyb_mean_DIAG', ascending=False).reset_index(drop=True)
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(temp['D_Pb_Hairpin_mean_DIAG'], temp['D_Pb_Hyb_Weighted_Mean'], 
                          marker='o', s=50, c=temp.index, cmap='viridis')
    plt.colorbar(scatter, label='Index')
    # Set labels and title
    plt.xlabel('D_Pb_Hairpin_mean_DIAG')
    plt.ylabel('D_Pb_Hyb_Weighted_mean')
    plt.title('Unfiltered results: '+ seqSTR + '  ' + sequence)
    
    # Add a grid for better readability
    plt.grid(True, linestyle='--', alpha=0.7)
    # Set the y-axis limits
    plt.ylim(0.4, 1)
    # plt.axvline(x=0.9, color='r', linestyle='--', alpha=0.5)
    # Show the plot
    plt.tight_layout()
    
    plt.savefig(results_dir + seqSTR +"_meanHPdiagVShybWEIGHTEDmeanDIAG.svg", bbox_inches='tight', pad_inches=0)
    plt.savefig(results_dir + seqSTR +"_meanHPdiagVShybWEIGHTEDmeanDIAG.png",  bbox_inches='tight', pad_inches=0)
    plt.show()
    
#non-weighted mean
for seqSTR in seqIDX:
    temp = hairpin_pd[hairpin_pd['Name'].str.contains(seqSTR)].reset_index(drop=True)
    temp.sort_values(by='D_Pb_Hairpin_D_Pb_Hyb_mean_DIAG', ascending=False).reset_index(drop=True)
    sequence = globals().get(seqSTR).sequence  # or use locals().get(seqSTR) if needed
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(temp['D_Pb_Hairpin_mean_DIAG'], temp['D_Pb_Hyb_mean'], 
                          marker='o', s=50, c=temp.index, cmap='viridis')
    plt.colorbar(scatter, label='Index')
    # Set labels and title
    plt.xlabel('D_Pb_Hairpin_mean_DIAG')
    plt.ylabel('D_Pb_Hyb_mean')
    plt.title('Unfiltered results: '+ seqSTR + '  ' + sequence)
   
    
    # Add a grid for better readability
    plt.grid(True, linestyle='--', alpha=0.7)
    # Set the y-axis limits
    plt.ylim(0.4, 1)
    # plt.axvline(x=0.9, color='r', linestyle='--', alpha=0.5)
    # Show the plot
    plt.tight_layout()
    plt.show()
    plt.savefig(results_dir + seqSTR +"_meanHPdiagVShybMEANdiag.svg", bbox_inches='tight', pad_inches=0)
    plt.savefig(results_dir + seqSTR +"_meanHPdiagVShybMEANdiag.png",  bbox_inches='tight', pad_inches=0)

plt.close('all')
#weighted mean by number of elements in largest diagonal
for seqSTR in seqIDX:
    temp = hairpin_pd[hairpin_pd['Name'].str.contains(seqSTR)].reset_index(drop=True)
    temp.sort_values(by='D_Pb_Hairpin_D_Pb_Hyb_mean_DIAG', ascending=False).reset_index(drop=True)
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(temp['D_Pb_Hairpin_mean_DIAG'], temp['D_Pb_Hyb_mean']/np.nanmax(non_nan_count_hyb_Pb), 
                          marker='o', s=50, c=temp.index, cmap='viridis')
    plt.colorbar(scatter, label='Index')
    # Set labels and title
    plt.xlabel('D_Pb_Hairpin_mean_DIAG')
    plt.ylabel('D_Pb_Hyb_mean')
    plt.title('Unfiltered results: '+ seqSTR + '  ' + sequence)
    
    # Add a grid for better readability
    plt.grid(True, linestyle='--', alpha=0.7)
    # Set the y-axis limits
    plt.ylim(0, 1)
    # plt.axvline(x=0.9, color='r', linestyle='--', alpha=0.5)
    # Show the plot
    plt.tight_layout()
    
    plt.savefig(results_dir + seqSTR +"_meanHPdiagVShybMEANdiagWEIGHTEDbyMAXdiagNUMelements.svg", bbox_inches='tight', pad_inches=0)
    plt.savefig(results_dir + seqSTR +"_meanHPdiagVShybMEANdiagEIGHTEDbyMAXdiagNUMelements.png",  bbox_inches='tight', pad_inches=0)
    plt.show()


# # Sort by 'D_Pb_Hairpin_mean_DIAG' in descending order
hairpin_pd = hairpin_pd.sort_values(by='D_Pb_Hairpin_D_Pb_Hyb_mean_DIAG', ascending=False).reset_index(drop=True)

if np.isnan(numTOP):
    hairpin_pd_top = hairpin_pd
else:
   # Find the top 100 unique hairpin prefixes based on the first 7 characters
   top_hairpin_prefixes = hairpin_pd['Hairpin'].str[:7].drop_duplicates().head(numTOP)

   # Filter the DataFrame to include all rows that match the top hairpin prefixes
   hairpin_pd_top = hairpin_pd[hairpin_pd['Hairpin'].str[:7].isin(top_hairpin_prefixes)].reset_index(drop=True)

unique_hairpin = hairpin_pd_top['Hairpin'].str[:7].unique()


# # Calculate a combined score for each hairpin
# hairpin_pd_top['Combined_Score'] = hairpin_pd_top['D_Pb_Hyb_Weighted_Mean']- hairpin_pd_top['D_Pb_Hairpin_mean_DIAG']

# # Rank the hairpins based on the combined score in descending order
hairpin_pd_top['Rank'] = hairpin_pd_top['D_Pb_Hairpin_D_Pb_Hyb_mean_DIAG'].rank(ascending=False, method='dense').astype(int)

# Get unique hairpins and their ranking
ranked_hairpins = hairpin_pd_top.drop_duplicates(subset='Hairpin', keep='first') #.sort_values(by='Combined_Score', ascending=False)
unique_hairpin = ranked_hairpins['Hairpin'].str[:7].unique()

# Create a color map based on the number of unique hairpins
cmap = cm.get_cmap('viridis', len(unique_hairpin))
color_mapping = {hairpin: cmap(i / len(unique_hairpin)) for i, hairpin in enumerate(unique_hairpin)}

# Calculate the number of rows and columns for the subplot matrix
num_plots = len(unique_hairpin)
num_rows = int(np.ceil(np.sqrt(num_plots)))
num_cols = int(np.ceil(num_plots / num_rows))

# Create a figure with subplots
fig, axs = plt.subplots(num_rows, num_cols, figsize=(20, 20))
fig.subplots_adjust(hspace=0.5, wspace=0.5)  # Adjust space between plots

# Flatten the axs array for easier iteration
axs = axs.flatten()

# Calculate global min and max for the x-axis (D_Pb_Hairpin_mean_DIAG) and y-axis (D_Pb_Hyb_Weighted_Mean)
x_min = hairpin_pd_top['D_Pb_Hairpin_mean_DIAG'].min()
x_max = hairpin_pd_top['D_Pb_Hairpin_mean_DIAG'].max()
y_min = hairpin_pd_top['D_Pb_Hyb_Weighted_Mean'].min()
y_max = hairpin_pd_top['D_Pb_Hyb_Weighted_Mean'].max()

# Loop through each unique hairpin and create a plot
for i, key in enumerate(unique_hairpin):
    # Filter rows that match the current hairpin
    rows = hairpin_pd_top[hairpin_pd_top['Hairpin'].str.startswith(key.split()[0])]
    
    # Get the rank and color for the current hairpin
    rank = ranked_hairpins[ranked_hairpins['Hairpin'].str[:7] == key].iloc[0]['Rank']
    color = color_mapping[key]

    # Plot the data on the subplot with the same color for all points
    axs[i].scatter(
        rows['D_Pb_Hairpin_mean_DIAG'], 
        rows['D_Pb_Hyb_Weighted_Mean'], 
        color=color, 
        alpha=0.7
    )
    
    # Set the title with rank and hairpin name
    axs[i].set_title(f'Rank {rank}: {key}', fontsize=8)
    axs[i].set_xlabel('D_Pb_Hairpin_mean_DIAG', fontsize=6)
    axs[i].set_ylabel('D_Pb_Hyb_Weighted_Mean', fontsize=6)
    axs[i].tick_params(axis='both', which='major', labelsize=6)
    axs[i].set_xlim(x_min, x_max)  # Set x-axis limits to global min and max
    axs[i].set_ylim(y_min, y_max)  # Set y-axis limits to global min and max

# Hide any unused subplots
for j in range(i + 1, len(axs)):
    fig.delaxes(axs[j])

# Display the plot matrix and save it
plt.tight_layout()
plt.show()

# Save the figure in .svg and .png formats
fig.savefig(results_dir + "top_unique_hairpin_ordinal_colored_plots.svg", bbox_inches='tight', pad_inches=0)
fig.savefig(results_dir + "top_unique_hairpin_ordinal_colored_plots.png", bbox_inches='tight', pad_inches=0)



# Sort hairpin_pd_top by the combined score and select the top 50 entries
hairpin_pd_top_sorted = hairpin_pd_top #.sort_values(by='Combined_Score', ascending=False).reset_index(drop=True)
top_50_hairpins = hairpin_pd_top_sorted.head(50)

# Create a new column for combined x-axis labels with both hairpin and MMP strings
top_50_hairpins['Combined_Label'] = top_50_hairpins['Hairpin'] + '\n' + top_50_hairpins['MMP9_Seq']

# Create a figure and axis with a specific size
fig, ax = plt.subplots(figsize=(18, 16))
width1 = 0.3
ar = np.arange(len(top_50_hairpins['Hairpin']))  # Range values for the top 50 Hairpin Sequences

# Plot the data for both the Hairpin and the Hybrid complex predicted energies
ax.bar(ar, top_50_hairpins['D_Pb_Hyb_Weighted_Mean'], width=0.3, align='edge', label='D_Pb_Hyb_Weighted_Mean', color=sns.xkcd_rgb['windows blue'], alpha=0.5)
ax.bar(ar + width1, abs(top_50_hairpins['D_Pb_Hairpin_mean_DIAG']), width=0.3, align='edge', label='D_Pb_Hairpin_mean_DIAG', color=sns.xkcd_rgb['orangish'], alpha=0.6)

# Set the y-axis label and title
ax.set_ylabel('Mean probability of diag weighted by remaining probabilities', fontsize=8)
ax.set_title('Weighted probabilities of Hairpin and Complex Structures (Top 50)', fontsize=12, weight='bold')

# Set the x-ticks with combined labels
ax.set_xticks(ar + width1 / 2)
ax.set_xticklabels(top_50_hairpins['Combined_Label'], rotation=70, fontweight='bold', fontsize='8', horizontalalignment='right')

# Add legend and adjust layout
ax.legend()
plt.tight_layout()

# Show the plot and save it
plt.show()
plt.savefig(results_dir + "weighted_probabilities_SuperClean_top50.svg", bbox_inches='tight', pad_inches=0)
plt.savefig(results_dir + "weighted_probabilities_SuperClean_top50.png", bbox_inches='tight', pad_inches=0)


# Group by each sequence identifier and get the top 50 by combined score
grouped_top_50 = hairpin_pd_top_sorted.groupby('Name').apply(lambda x: x.nlargest(50, 'D_Pb_Hairpin_D_Pb_Hyb_mean_DIAG')).reset_index(drop=True)

# Iterate over each unique sequence and create a separate plot for each
for seq_identifier in grouped_top_50['Name'].unique():
    # Filter the top 50 hairpins for the current sequence
    top_50_hairpins = grouped_top_50[grouped_top_50['Name'] == seq_identifier]

    # Calculate the number of rows and columns for the subplot grid
    num_plots = len(top_50_hairpins)
    num_rows = int(np.ceil(np.sqrt(num_plots)))
    num_cols = int(np.ceil(num_plots / num_rows))

    # Create a figure for the subplots
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(20, 20))
    fig.subplots_adjust(hspace=0.5, wspace=0.5)  # Adjust spacing between plots
    axs = axs.flatten()  # Flatten axs for easier iteration

    # Loop through each top hairpin to create the subplot
    for i, (_, row) in enumerate(top_50_hairpins.iterrows()):
        # Retrieve the Hairpin_Pb matrix, hairpin sequence, MMP9 sequence, and combined score
        pb_matrix = row['Hairpin_Pb']
        hairpin_seq = row['Hairpin']
        mmp9_seq = row['MMP9_Seq']

        # Plot the matrix as an imagesc-like heatmap
        # im = axs[i].imshow(pb_matrix, interpolation='nearest', cmap='viridis', aspect='auto', extent=[0, len(mmp9_seq), len(hairpin_seq), 0])
        im = axs[i].imshow(pb_matrix, interpolation='nearest', cmap='viridis', aspect='auto')
        
        # Add a color bar for each subplot
        divider = make_axes_locatable(axs[i])
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')

        # Set x and y tick labels
        axs[i].set_xticks(np.arange(len(mmp9_seq)))
        axs[i].set_xticklabels(mmp9_seq, fontsize=6, rotation=90)
        axs[i].set_yticks(np.arange(len(hairpin_seq)))
        axs[i].set_yticklabels(hairpin_seq, fontsize=6)
        
        # Set title for each subplot with rank and hairpin name
        axs[i].set_title(f'{seq_identifier}: {hairpin_seq[:7]} vs MMP9', fontsize=8)

    # Hide any remaining empty subplots if there are fewer than num_rows * num_cols
    for j in range(i + 1, len(axs)):
        fig.delaxes(axs[j])

    plt.tight_layout()
    plt.show()

    # Save each figure with a unique filename per sequence
    fig.savefig(results_dir + f"top50_hairpins_{seq_identifier}_imagesc_plots.svg", bbox_inches='tight', pad_inches=0)
    fig.savefig(results_dir + f"top50_hairpins_{seq_identifier}_imagesc_plots.png", bbox_inches='tight', pad_inches=0)

# Optionally save the combined top 50 data for all sequences
grouped_top_50.to_excel(results_dir + 'hairpin_pdDATA_FinalResults_Top50PerSequence.xlsx', index=True)




