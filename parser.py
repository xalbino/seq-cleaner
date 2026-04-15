from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path 

def parse_fasta(filepath):
    return list(SeqIO.parse(filepath, "fasta"))

