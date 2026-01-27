# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 13:10:51 2025

@author: m.roska
"""
# %% packages
import pandas as pd
import re
import numpy as np
import periodictable
from collections import defaultdict


# %% support functions

# count atom number
def parse_chemical_formula(formula):
    # Regex pattern to capture elements and their counts
    pattern = r'([A-Z][a-z]*)(\d*)'
    # Find all matches
    matches = re.findall(pattern, formula)
    # Dictionary to store the counts of each element
    element_counts = defaultdict(int)
    for (element, count) in matches:
        if count == '':
            count = 1
        else:
            count = int(count)
        element_counts[element] += count
    return element_counts

# apply count number of atoms
def atom_num(formula, atom):
    element_counts = parse_chemical_formula(formula)
    return element_counts[atom]

# helper fucntion for differnt naimgngs
def get_best_count(counts, *possible_keys):
    """Return the first available value for given keys."""
    for key in possible_keys:
        if key in counts:
            return counts[key]
    return 0  # fallback if none found


# %% main functions 

# calcualte volatility based on formula
# Li et al. 2016
# implemented with suport of Alexandra Tsimpidi
def calc_vol(spc_eqn):
    # checks if carbon and hydrogen are in formula
    if all(element in spc_eqn for element in ['C', 'H']):
        # checks if Br Cl I and F are not in formula
        if all(element not in spc_eqn for element in ['Br', 'Cl', 'I', 'F']):
            # calcualte number of atoms
            nC = atom_num(spc_eqn, 'C')
            nH = atom_num(spc_eqn, 'H')
            nO = atom_num(spc_eqn, 'O')
            nN = atom_num(spc_eqn, 'N')
            nS = atom_num(spc_eqn, 'S')
            # Apply equation of Li et al. 2016 (Molecular corridors and parameterizations of volatility in the chemical evolution of organic aerosols, ACP, 2016)
            # log10 C0
            n0C = 23.80
            bC  = 0.4861
            bO  = 0.0
            bCO = 0.0
            bN  = 0.0
            bS  = 0.0
            # choose set of variables based on atom types contained
            if 'O' in spc_eqn :
                n0C = 22.66
                bC  = 0.4481
                bO  = 1.656
                bCO = -0.7790
                bN  = 0.0
                bS  = 0.0
            if 'N' in spc_eqn :
                n0C = 24.59
                bC  = 0.4066
                bO  = 0.0
                bCO = 0.0
                bN  = 0.9619
                bS  = 0.0
            if 'O' in spc_eqn and 'N' in spc_eqn :
                n0C = 24.13
                bC  = 0.3667
                bO  = 0.7732
                bCO = -0.07790
                bN  = 1.114
                bS  = 0.0
            if 'O' in spc_eqn and 'S' in spc_eqn:
                n0C = 24.06
                bC  = 0.3637
                bO  = 1.327
                bCO = -0.3988
                bN  = 0.0
                bS  = 0.7579
            if 'O' in spc_eqn and 'N' in spc_eqn and 'S' in spc_eqn:
                n0C = 28.50
                bC  = 0.3848
                bO  = 1.011
                bCO = 0.2921
                bN  = 1.053
                bS  = 1.316
            # calculate voaltility as log(C0) of formula
            vol = ((n0C-nC)*bC-nO*bO-2*(nC*nO/(nC+nO))*bCO-nN*bN-nS*bS)
        else:
            print('Br, Cl, I or F in forumla')
            vol = np.nan
    else:
        print('no C and/or H in forumla')
        vol = np.nan
    return vol


# calcuates molecular mass
# potentailly add a reagiont ion etc to the formula
def calculate_mass(formula, added = None):
    matches = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    total_mass = 0.0
    for element, count in matches:
        count = int(count) if count else 1
        try:
            atomic_mass = getattr(periodictable, element).mass
            total_mass += atomic_mass * count
        except AttributeError:
            print(f"Warning: Unknown element '{element}' in formula")
    if added == None:
        added_mass = 0.0
    else:
        matches = re.findall(r'([A-Z][a-z]*)(\d*)', added)
        added_mass = 0.0
        for element, count in matches:
            count = int(count) if count else 1
            try:
                atomic_mass = getattr(periodictable, element).mass
                added_mass += atomic_mass * count
            except AttributeError:
                print(f"Warning: Unknown element '{element}' in added")
                
    return total_mass + added_mass


# calcualte soa yield basesd on formula
# very rought estiamtion
def estimate_soayield_from_formula(counts):
    """
    counts: dict or row containing 'C','H','O' counts (integers)
    Returns: approximate SOA yield fraction (0–1)
    """
    C = get_best_count(counts, 'C', 'carbon_count')
    O = get_best_count(counts, 'O', 'oxygen_count')
    # basic class-based default
    if O == 0:
        # assume hydrocarbon: if C>=5 assume aromatic → yield ~0.25, else alkane ~0.00
        return 0.25 if C >= 5 else 0.00
    else:
        # oxygenated: yield increases with O/C ratio
        ocr = O / max(C, 1)
        return min(ocr * 0.5, 0.5)  # cap at ~50%

def add_soa_yield_column(df):
    def row_to_counts(row):
        return {
            'C': get_best_count(row, 'C', 'carbon_count'),
            'H': get_best_count(row, 'H', 'hydrogen_count'),
            'O': get_best_count(row, 'O', 'oxygen_count'),
        }

    df['estimated_SOAyield'] = df.apply(lambda row: estimate_soayield_from_formula(row_to_counts(row)), axis=1)
    return df