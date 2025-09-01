# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 13:10:51 2025

@author: m.roska
"""
# %% packages
import pandas as pd
import re


# %% functions

# helper fucntion for differnt naimgngs
def get_best_count(counts, *possible_keys):
    """Return the first available value for given keys."""
    for key in possible_keys:
        if key in counts:
            return counts[key]
    return 0  # fallback if none found


# calcualte soa yield basesd on formula
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