import argparse
import json
import sys
import time
from functools import wraps
from multiprocessing import Pool
from pathlib import Path
from typing import Any, Dict, Generator, List, Tuple

from Bio.SeqIO.FastaIO import SimpleFastaParser
from pro2codon import pn2codon

GENETIC_TABLE = """
{
    "1": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG",
            "TGA"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "2": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG",
            "AGA",
            "AGG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGA",
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG"
        ],
        "I": [
            "ATT",
            "ATC"
        ],
        "M": [
            "ATA",
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "3": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGA",
            "TGG"
        ],
        "T": [
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC"
        ],
        "M": [
            "ATA",
            "ATG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "ATT",
            "ATC"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "4": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGA",
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "5": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC",
            "AGA",
            "AGG"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGA",
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG"
        ],
        "I": [
            "ATT",
            "ATC"
        ],
        "M": [
            "ATA",
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "6": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "Q": [
            "TAA",
            "TAG",
            "CAA",
            "CAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "*": [
            "TGA"
        ],
        "W": [
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "TAA",
            "TAG",
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "9": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC",
            "AGA",
            "AGG"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGA",
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC",
            "AAA"
        ],
        "K": [
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "AAA",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "10": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG"
        ],
        "C": [
            "TGT",
            "TGC",
            "TGA"
        ],
        "W": [
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "11": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG",
            "TGA"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "12": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "CTG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG",
            "TGA"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "13": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGA",
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG"
        ],
        "I": [
            "ATT",
            "ATC"
        ],
        "M": [
            "ATA",
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "G": [
            "AGA",
            "AGG",
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "14": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC",
            "AGA",
            "AGG"
        ],
        "Y": [
            "TAT",
            "TAC",
            "TAA"
        ],
        "*": [
            "TAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGA",
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC",
            "AAA"
        ],
        "K": [
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "AAA",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "15": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TGA"
        ],
        "Q": [
            "TAG",
            "CAA",
            "CAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "TAG",
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "16": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "TAG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TGA"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "TAG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "21": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC",
            "AGA",
            "AGG"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGA",
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG"
        ],
        "I": [
            "ATT",
            "ATC"
        ],
        "M": [
            "ATA",
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC",
            "AAA"
        ],
        "K": [
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "AAA",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "22": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "TAG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCG",
            "AGT",
            "AGC"
        ],
        "*": [
            "TCA",
            "TAA",
            "TGA"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "TAG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "23": {
        "F": [
            "TTT",
            "TTC"
        ],
        "*": [
            "TTA",
            "TAA",
            "TAG",
            "TGA"
        ],
        "L": [
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "24": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC",
            "AGA"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGA",
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG",
            "AGG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "25": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "G": [
            "TGA",
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "W": [
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "26": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TAG",
            "TGA"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGG"
        ],
        "A": [
            "CTG",
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "27": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "Q": [
            "TAA",
            "TAG",
            "CAA",
            "CAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGA",
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "TAA",
            "TAG",
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "28": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "Q": [
            "TAA",
            "TAG",
            "CAA",
            "CAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGA",
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "TAA",
            "TAG",
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "29": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC",
            "TAA",
            "TAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "*": [
            "TGA"
        ],
        "W": [
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "30": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "E": [
            "TAA",
            "TAG",
            "GAA",
            "GAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "*": [
            "TGA"
        ],
        "W": [
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "TAA",
            "TAG",
            "GAA",
            "GAG",
            "CAA",
            "CAG"
        ]
    },
    "31": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "E": [
            "TAA",
            "TAG",
            "GAA",
            "GAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGA",
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "TAA",
            "TAG",
            "GAA",
            "GAG",
            "CAA",
            "CAG"
        ]
    },
    "32": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC"
        ],
        "Y": [
            "TAT",
            "TAC"
        ],
        "*": [
            "TAA",
            "TGA"
        ],
        "W": [
            "TAG",
            "TGG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG",
            "AGA",
            "AGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    },
    "33": {
        "F": [
            "TTT",
            "TTC"
        ],
        "L": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG"
        ],
        "S": [
            "TCT",
            "TCC",
            "TCA",
            "TCG",
            "AGT",
            "AGC",
            "AGA"
        ],
        "Y": [
            "TAT",
            "TAC",
            "TAA"
        ],
        "*": [
            "TAG"
        ],
        "C": [
            "TGT",
            "TGC"
        ],
        "W": [
            "TGA",
            "TGG"
        ],
        "P": [
            "CCT",
            "CCC",
            "CCA",
            "CCG"
        ],
        "H": [
            "CAT",
            "CAC"
        ],
        "Q": [
            "CAA",
            "CAG"
        ],
        "R": [
            "CGT",
            "CGC",
            "CGA",
            "CGG"
        ],
        "I": [
            "ATT",
            "ATC",
            "ATA"
        ],
        "M": [
            "ATG"
        ],
        "T": [
            "ACT",
            "ACC",
            "ACA",
            "ACG"
        ],
        "N": [
            "AAT",
            "AAC"
        ],
        "K": [
            "AAA",
            "AAG",
            "AGG"
        ],
        "V": [
            "GTT",
            "GTC",
            "GTA",
            "GTG"
        ],
        "A": [
            "GCT",
            "GCC",
            "GCA",
            "GCG"
        ],
        "D": [
            "GAT",
            "GAC"
        ],
        "E": [
            "GAA",
            "GAG"
        ],
        "G": [
            "GGT",
            "GGC",
            "GGA",
            "GGG"
        ],
        "X": [],
        "B": [
            "AAT",
            "AAC",
            "GAT",
            "GAC"
        ],
        "J": [
            "TTA",
            "TTG",
            "CTT",
            "CTC",
            "CTA",
            "CTG",
            "ATT",
            "ATC",
            "ATA"
        ],
        "Z": [
            "CAA",
            "CAG",
            "GAA",
            "GAG"
        ]
    }
}
"""

DICT_TABLES = json.loads(GENETIC_TABLE)


def return_aligned_paths(
    glob_paths_taxa: List[Path],
    glob_paths_genes: List[Path],
    path_aligned: Path,
    d,
) -> Generator[Path, Any, Any]:
    for path_nt, path_aa in zip(glob_paths_taxa, glob_paths_genes):
        if not path_nt.is_file() or not path_aa.is_file():
            continue

        stem_taxon = Path(path_nt.stem).stem
        stem_gene = Path(path_aa.stem).stem

        if stem_taxon != stem_gene:
            continue

        yield (
            path_aa,
            path_nt,
            path_aligned.joinpath(Path(f"{stem_taxon}.nt.fa")),
            d,
        )


def prepare_taxa_and_genes(input: str, d) -> Tuple[Generator[
    Tuple[Path, Path, Path],
    Any,
    Any
], int]:
    input_path = Path(input)

    joined_mafft = input_path.joinpath(Path("mafft"))
    joined_nt = input_path.joinpath(Path("nt"))
    joined_nt_aligned = input_path.joinpath(Path("nt_aligned"))

    if not joined_nt_aligned.exists():
        joined_nt_aligned.mkdir()

    glob_genes = sorted(list(joined_mafft.glob("*aa.fa")))
    glob_taxa = sorted(list(joined_nt.glob("*nt.fa")))

    out_generator = return_aligned_paths(
        glob_taxa,
        glob_genes,
        joined_nt_aligned,
        d,
    )

    return out_generator, len(glob_genes)

from pprint import pprint

DICT_OPP = {'NT Seqs': "AA Seqs", "AA Seqs": 'NT_Seqs'}
from copy import copy
def read_and_convert_fasta_files(
    aa_file: str,
    nt_file: str,
) -> Dict[str, 
        Tuple[
            List[Tuple[str, str]],
            List[Tuple[str, str]]
        ]
]:
    nt_seqs = SimpleFastaParser(open(nt_file))
    aa_seqs = SimpleFastaParser(open(aa_file))


    aas = []
    nts = {}


    for aa, nt in zip(aa_seqs, nt_seqs):
        aa_header, aa_seq = aa
        nt_header, nt_seq = nt

        aa_header, aa_seq = aa_header.strip(), aa_seq.strip()
        nt_header, nt_seq = nt_header.strip(), nt_seq.strip()

        if aa_header[-2:] == '.a':
            print("- ", i)
            break

        nts[nt_header] = (nt_header, nt_seq)
        aas.append((aa_header, aa_seq))

        

    


    ret = {}
    
    i = -1

    for aa in aas:
        try:
            if aa[0] not in ret:
                i += 1

            ret[aa[0]] = (aa, (i, *nts[aa[0]]))
            

        except:
            print("ERROR CAUGHT: There is a single header in PEP sequence FASTA file that does not exist in NUC sequence FASTA file")

            sys.exit()
    return ret


def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(

        usage="pal2nal.py [OPTIONS] [FILES]...",

        description="Batch-Convert Amino Acids to Codons with error checking using the default NCBI table, or your specified table."
    )

    parser.add_argument('-i', '--input', type=str, default='Parent',
                        help='Parent input path.')
    parser.add_argument('-p', '--processes', type=int, default=4,
                        help='Number of threads used to call processes.')
    parser.add_argument('-t', '--table', type=int, default=1,
                        help='Table ID.')
    return parser


def worker(tup: Tuple[Path, Path, Path, Dict]):
    aa_file, nt_file, out_file, d = tup
    seqs = read_and_convert_fasta_files(
        aa_file,
        nt_file
    )

    stem = out_file.stem
    res = pn2codon(stem, d, seqs)

    out_file.write_text(res)    



def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(
            f'Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds')
        return result
    return timeit_wrapper


@timeit
def run_batch_threaded(num_threads: int, ls: List[
        List[
            List[Tuple[Tuple[Path, Path, Path], Dict]]
        ]]):


    with Pool(num_threads) as pool:
        list(pool.map(worker, ls, chunksize=100))


if __name__ == "__main__":
    arg_parser = init_argparse()
    args = arg_parser.parse_args()

    d = DICT_TABLES[str(args.table)]

    generator, _ = prepare_taxa_and_genes(args.input, d)

    run_batch_threaded(num_threads=args.processes, ls=generator)
