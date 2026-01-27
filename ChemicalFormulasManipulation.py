# -*- coding: utf-8 -*-
"""
Created on Thu Jul  3 11:20:51 2025

@author: m.roska
"""

# %% packages
import periodictable as pt
import re
from collections import defaultdict

# %% fucntions


# function to drop NH4+ from formula reducing n by 1 h by 4 and strip the +
def strip_nh4_plus(formula: str):
    # Remove + sign
    formula = formula.rstrip('+')

    # Extract elements and their counts
    matches = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    element_counts = {}
    for elem, count in matches:
        count = int(count) if count else 1
        element_counts[elem] = element_counts.get(elem, 0) + count

    # Adjust counts
    element_counts['N'] = max(0, element_counts.get('N', 0) - 1)
    element_counts['H'] = max(0, element_counts.get('H', 0) - 4)

    # Reconstruct formula with alphabetical order
    stripped_formula = ''.join(
        f"{elem}{element_counts[elem] if element_counts[elem] > 1 else ''}"
        for elem in sorted(element_counts) if element_counts[elem] > 0
    )
    return stripped_formula


# fucntion to spereate abse formula and the NH4+ ion
def extract_nh4_plus(formula: str):
    # Parse chemical formula into element counts
    def parse_formula(formula):
        tokens = re.findall(r'[A-Z][a-z]?|\d+|[+.-]', formula)
        elements = defaultdict(int)
        i = 0
        charge = ''
        while i < len(tokens):
            token = tokens[i]
            if re.fullmatch(r'[A-Z][a-z]?', token):
                elem = token
                count = 1
                if i + 1 < len(tokens) and tokens[i + 1].isdigit():
                    count = int(tokens[i + 1])
                    i += 1
                elements[elem] += count
            elif token in '+-.':
                charge += token
            i += 1
        return dict(elements), charge

    elements, charge = parse_formula(formula)

    # Check for NH4+
    if elements.get('N', 0) >= 1 and elements.get('H', 0) >= 4:
        elements['N'] -= 1
        elements['H'] -= 4
        nh4 = 'NH4+'
    else:
        return formula, None

    # Rebuild remaining compound
    remaining = ''.join(
        f"{el}{'' if count == 1 else count}"
        for el, count in elements.items() if count > 0
    )
    return remaining, nh4


# fucntion to seperate formula and the NH4+ ion
# keep them together with a dot
def extract_nh4_plus_cdot(formula: str):
    # Parse chemical formula into element counts
    def parse_formula(formula):
        tokens = re.findall(r'[A-Z][a-z]?|\d+|[+.-]', formula)
        elements = defaultdict(int)
        i = 0
        charge = ''
        while i < len(tokens):
            token = tokens[i]
            if re.fullmatch(r'[A-Z][a-z]?', token):
                elem = token
                count = 1
                if i + 1 < len(tokens) and tokens[i + 1].isdigit():
                    count = int(tokens[i + 1])
                    i += 1
                elements[elem] += count
            elif token in '+-.':
                charge += token
            i += 1
        return dict(elements), charge

    elements, charge = parse_formula(formula)

    # Check for NH4+
    if elements.get('N', 0) >= 1 and elements.get('H', 0) >= 4:
        elements['N'] -= 1
        elements['H'] -= 4
        nh4 = 'NH4+'
    else:
        return formula, None

    # Rebuild remaining compound
    remaining = ''.join(
        f"{el}{'' if count == 1 else count}"
        for el, count in elements.items() if count > 0
    )

    updated_formula = f'{remaining} · {nh4}'
    return updated_formula

# function to count all elements in a formula to a list
# read out output like this count_elements(formula).get('C', 0)
def count_elements(formula):
    """Return dict of element counts from a formula like 'C7H14O2'."""
    matches = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    counts = {}
    for elem, count in matches:
        count = int(count) if count else 1
        counts[elem] = counts.get(elem, 0) + count
    return counts

# calcaultes exact mass for a formula
def mass_from_formula(formula):
    """Compute molecular mass from formula like 'C7H14O'."""
    matches = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    total_mass = 0
    for elem, count in matches:
        count = int(count) if count else 1
        try:
            total_mass += getattr(pt, elem).mass * count
        except AttributeError:
            return None  # invalid element
    return total_mass


# clasify group of compound CH, CHO, CHON, CHON2, CHON3
def classify_group(formula):
    counts = count_elements(formula)
    has_C = 'C' in counts
    has_H = 'H' in counts
    has_O = 'O' in counts
    N = counts.get('N', 0)

    if has_C and has_H:
        if N > 0:
            return f'CHON{N}' if N <= 3 else 'CHON>3'
        elif has_O:
            return 'CHO'
        else:
            return 'CH'
    return 'Other'