#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 13:02:42 2024

@author: michaelbh23
"""

# sequence_utils.py

def complement(sequence, nucleotide_type):
    '''Return reverse complement of DNA or RNA strand'''
    if nucleotide_type == "DNA":
        seq = sequence.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c").upper()[::-1]
    elif nucleotide_type == "RNA":
        seq = sequence.replace("A", "u").replace("C", "g").replace("U", "a").replace("G", "c").upper()[::-1]
    else:
        raise ValueError("Unrecognised nucleotide type.")
    return seq

def dna_to_mrna(inputDNA):
    '''Convert DNA sequence to mRNA'''
    return inputDNA.replace("T", "U")