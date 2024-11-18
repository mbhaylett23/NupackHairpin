#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 13:23:15 2024

@author: michaelbh23
"""

# thermo_analysis.py

import numpy as np
from nupack import Model, Strand, Complex, ComplexSet, SetSpec, structure_probability
import pandas as pd

def analyse_thermodynamics(hairpin_pd, mmp9_dna, temperature=37, sodium=0.5, magnesium=0.2, gap=1.1):
    """
    Analyse thermodynamic properties of hairpins and their hybridisation with MMP9.
    
    Parameters:
    hairpin_pd (DataFrame): DataFrame containing hairpin sequences.
    mmp9_dna (str): DNA sequence of MMP9.
    temperature (float): Temperature in Celsius.
    sodium (float): Sodium concentration in M.
    magnesium (float): Magnesium concentration in M.
    gap (float): Energy gap for suboptimal structures in kcal/mol.
    
    Returns:
    pd.DataFrame: Updated DataFrame with thermodynamic properties.
    """
    # Define physical model
    mmp9_model = Model(material='dna04', celsius=temperature, sodium=sodium, magnesium=magnesium)
    
    # Lists to store thermodynamic properties
    hairpin_dg = []
    hairpin_mfe = []
    hairpin_p = []
    hairpin_ep = []

    hyb_dg = []
    hyb_mfe = []
    hyb_p = []
    hyb_ep = []

    # Perform analysis for each hairpin sequence and its hybrid with MMP9
    for hairpin, mmp9_seq in zip(hairpin_pd['Hairpin'], hairpin_pd['MMP9_Seq']):
        # Define strands
        hp_st = Strand(hairpin, name='hp_st')
        mmp9_st = Strand(mmp9_seq, name='mmp9_st')
        
        # Define complexes
        hp_complex = Complex([hp_st])
        hyb_complex = Complex([hp_st, mmp9_st])

        # Define complex sets for hairpin and hybrid structures
        hp_set = ComplexSet(strands={hp_st: 1e-6}, complexes=SetSpec(max_size=0, include=[hp_complex]))
        hyb_set = ComplexSet(strands={hp_st: 1e-6, mmp9_st: 1e-6}, complexes=SetSpec(max_size=0, include=[hyb_complex]))

        # Hairpin complex analysis
        hp_result = hp_set.compute('pfunc', 'pairs', 'mfe', 'subopt', model=mmp9_model, energy_gap=gap)
        hairpin_dg.append(hp_result[hp_complex].free_energy)
        hairpin_mfe.append(hp_result[hp_complex].mfe[0][1])
        hp_pair_matrix = hp_result[hp_complex].pairs.to_array()
        np.fill_diagonal(hp_pair_matrix, 0)  # Zero out diagonal for probability matrix
        hairpin_p.append(hp_pair_matrix)
        eqP = structure_probability([hp_st], hp_result[hp_complex].mfe[0].structure, mmp9_model)
        hairpin_ep.append(eqP)

        # Hybridisation complex analysis
        hyb_result = hyb_set.compute('pfunc', 'pairs', 'mfe', 'subopt', model=mmp9_model, energy_gap=gap)
        hyb_dg.append(hyb_result[hyb_complex].free_energy)
        hyb_mfe.append(hyb_result[hyb_complex].mfe[0][1])
        hyb_pair_matrix = hyb_result[hyb_complex].pairs.to_array()
        np.fill_diagonal(hyb_pair_matrix, 0)
        hyb_p.append(hyb_pair_matrix)
        eqP = structure_probability([hp_st, mmp9_st], hyb_result[hyb_complex].mfe[0].structure, mmp9_model)
        hyb_ep.append(eqP)

    # Append thermodynamic properties to hairpin DataFrame
    hairpin_pd['Hairpin_dG (kcal/mol)'] = hairpin_dg
    hairpin_pd['Hairpin_MFE (kcal/mol)'] = hairpin_mfe
    hairpin_pd['Hairpin_Pb'] = hairpin_p
    hairpin_pd['Hairpin_equiPb'] = hairpin_ep
    hairpin_pd['Hyb_dG (kcal/mol)'] = hyb_dg
    hairpin_pd['Hyb_MFE (kcal/mol)'] = hyb_mfe
    hairpin_pd['Hyb_Pb'] = hyb_p
    hairpin_pd['Hyb_equiPb'] = hyb_ep

    return hairpin_pd