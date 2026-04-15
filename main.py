from Bio import SeqIO
from parser import parse_fasta
from cleaner import clean_all
from pathlib import Path 
from cleaner import clean_all
from parser import parse_fasta

records = parse_fasta(Path("data/example.fasta"))
cleaned = clean_all(records)
SeqIO.write(cleaned, Path("output/output.fasta"), "fasta")
print(f"{len(cleaned)} Sequenzen bereinigt.")

