#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 13:18:32 2024

@author: michaelbh23
"""

# hairpin_generator.py

import itertools
import pandas as pd
from sequence_utils import complement

def generate_hairpins(dna_seq, hpLEN):
    """
    Generate hairpin sequences based on the input DNA sequences.
    
    Parameters:
    dna_seq (DataFrame): The DataFrame containing DNA sequences.
    hpLEN (int): Length of the hairpin sequence.

    Returns:
    pd.DataFrame: DataFrame with hairpin sequences and associated data.
    """
    nucleotides = ['A', 'C', 'G', 'T']
    reverse_complementary_hairpins = []
    
    # Generate all possible combinations of nucleotides of specified length
    sequences = list(itertools.product(nucleotides, repeat=hpLEN))
    filtered_sequences = [seq for seq in sequences if ''.join(seq).count('GC') in [2, 4]]
    list_of_strings = [''.join(t) for t in filtered_sequences]
    
    # Get reverse complement for each generated hairpin sequence
    for seq in list_of_strings:
        seq_complement = complement(seq, 'DNA')
        reverse_complementary_hairpins.append(seq_complement)
    
    # Construct hairpins and store in DataFrame
    hairpins = []
    seq_name = []
    hairpin_length = []
    
    for name, sequence in zip(dna_seq['Name'], dna_seq['Sequence']):
        for seq, rev in zip(list_of_strings, reverse_complementary_hairpins):
            hairpin = seq + sequence + rev  # Construct the hairpin sequence
            hairpins.append(hairpin)
            seq_name.append(name)
            hairpin_length.append(len(hairpin))

    # Create a DataFrame with the generated hairpin data
    hairpin_pd = pd.DataFrame({
        'Name': seq_name,
        'Hairpin': hairpins,
        'Length': hairpin_length
    })
    
    return hairpin_pd