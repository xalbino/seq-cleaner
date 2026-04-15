from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def dna_rna(record):
    return str(record.seq).upper().replace('T', 'U')