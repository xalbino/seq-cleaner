from Bio import SeqIO
from parser import parse_fasta
from cleaner import clean_all
from pathlib import Path 
from analyzer import n_content
from analyzer import gc_content
from analyzer import records_stats

# records = parse_fasta(Path("data/example.fasta"))
# cleaned = clean_all(records)
# SeqIO.write(cleaned, Path("output/output.fasta"), "fasta")
# print(f"{len(cleaned)} Sequenzen bereinigt.")

records = parse_fasta(Path("output/output.fasta"))
# for record in records:
  #  print(gc_content(record))
  #   print(n_content(record))

print(records_stats(records))

