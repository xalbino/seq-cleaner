from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

VALID = set("ACGTN")

def clean_record(record):
    cleaned_seq = "".join(c if c in VALID else "N" for c in str(record.seq).upper())
    return SeqRecord(Seq(cleaned_seq), id=record.id, description=record.description)

def clean_all(records):
    return [clean_record(r) for r in records]
