from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def dna_rna(record):
    return str(record.seq).upper().replace('T', 'U')

def rna_to_protein(record):
    codon_table = {
        "UUU": "F", #Phenylalanine 
        "UUC": "F",
        "UUA": "L", #Leucine
        "UUG": "L",
        "UCU": "S", #Serine
        "UCC": "S",
        "UCA": "S",
        "UCG": "S",
        "UAU": "Y", #Tyrosine
        "UAC": "Y",
        "UAA": "#", #Ochre
        "UAG": "#", #Amber
        "UGU": "C", #Cysteine
        "UGC": "C",
        "UGA": "#", #Opal
        "UGG": "W", #Tryptophan
        "CUU": "L", # Leucine
        "CUC": "L",
        "CUA": "L",
        "CUG": "L",
        "CCU": "P", #Proline
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAU": "H", #Histidine
        "CAC": "H",
        "CAA": "Q", #Glutamine
        "CAG": "Q", 
        "CGU": "R", #Argenine
        "CGC": "R", 
        "CGA": "R",
        "CGG": "R",
        "AUU": "I", #Isoleucine
        "AUC": "I",
        "AUA": "I",
        "AUG": "M", #Methionine
        "ACU": "T", #Threonine
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAU": "N", #Asparagine
        "AAC": "N",
        "AAA": "K", #Lysine
        "AAG": "K",
        "AGU": "S", #Serine
        "AGC": "S",
        "AGA": "R", #Argenine
        "AGG": "R",
        "GUU": "V", #Valine
        "GUC": "V",
        "GUA": "V", 
        "GUG": "V",
        "GCU": "A", #Alanine
        "GCG": "A",
        "GCA": "A",
        "GCG": "A",
        "GAU": "D", #Aspartic acid
        "GAC": "D",
        "GAA": "E" #Glutamic acid
        "GAG": "E",
        "GGU": "G", #Glycine
        "GGC": "G",
        "GGA": "G",
        "GGG": "G"
    }

    str_record = str(record.seq)
    result = []
     for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = codon_table.get(codon, "?")  # "?" for unknown codons
        if aa == "#":
            break  # stop at stop codon
        result.append(aa)

    return "".join(result)