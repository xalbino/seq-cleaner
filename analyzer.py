from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def gc_content(record):
    """
    Function calculates the length of the shortest sequence in a fasta, the length of the longest sequence in a fasta and the average length of records inside of a fasta.
    
    Parameters:
        records: fasta to analyze

    Requirements: 
        Must already be cleaned.

    Returns:
        list: [min, max, avg]

    Time Complexity: O(n)
    Space Complexity: O(n)
    """
    seq_string = str(record.seq)
    gc_counter = 0
    for c in seq_string.upper():
        if c == 'G' or c == 'C':
            gc_counter += 1
    return gc_counter/len(seq_string)

def records_stats(records):
    """
    Function calculates the length of the shortest sequence in a fasta, the length of the longest sequence in a fasta and the average length of records inside of a fasta.
    
    Parameters:
        records: fasta to analyze

    Requirements: 
        Must already be cleaned.

    Returns:
        list: [min, max, avg]

    Time Complexity: O(n)
    Space Complexity: O(n)
    """
    max_length = 0
    min_length = float('inf')
    total_length = 0
    total_records = 0
    for record in records:
        seq_len = len(str(record.seq))
        if seq_len > max_length:
            max_length = seq_len
        if seq_len < min_length:
            min_length = seq_len
        total_length += seq_len
        total_records += 1
    return [min_length, max_length, int(total_length/total_records)]

def n_content(record):
    """
    Function compares two sequences and shows in how many chars they differ.
    
    Parameters:
        record_a (record): the first record 
        record_b (record): the second record

    Requirements: 
        Sequences must be of equal length.

    Returns:
        int: The amount of nucleotids that differ in both sequences.

    Time Complexity: O(n^2): n = len(record_a) 
    Space Complexity: O(n)
    """
    n_counter = 0
    for c in str(record.seq).upper():
        if c == 'N':
            n_counter += 1
    return n_counter/len(record.seq)

def seq_len(record)
    return len(str(record.seq))

def point_mutations(record_a, record_b)
    """
    Function compares two sequences and shows in how many chars they differ.
    
    Parameters:
        record_a (record): the first record 
        record_b (record): the second record

    Requirements: 
        Sequences must be of equal length.

    Returns:
        int: The amount of nucleotids that differ in both sequences.

    Time Complexity: O(n^2): n = len(record_a) 
    Space Complexity: O(n)
    """
    str_a = str(record_a.seq)
    str_b = str(record_b.seq)

    if len(str_a) != len(str_b):
        raise ValueError('Sequences must be of equal length')

    counter = 0
    for c in str_a:
        for d in str_b:
            if c != d:
                counter += 1
    return counter

def check_motif_in_sequence(sequence: str, motif: str) -> int:
    """
    Function provides the index of the first occurrence of a provided motif within a given sequence.
    If the motif is not found, the function returns -1.
    
    Parameters:
        sequence (str): The sequence in which to search for the motif.
        motif (str): The motif to search for within the sequence.

    Requirements: 
        Only one occurrence of the motif is expected in the sequence.

    Returns:
        int: The index of the first occurrence of the motif in the sequence, or -1 if not found.

    Time Complexity: O(n*m), n = length of sequence, m = length of motif
    Space Complexity: O(1)
    """

    return sequence.upper().find(motif.upper())

def check_motif_positions(seq, motif):
    """
    Function provides a list of all occurrences of a provided motif within a given sequence.
    
    Parameters:
        sequence (str): The sequence in which to search for the motif.
        motif (str): The motif to search for within the sequence.

    Requirements: 
        None.

    Returns:
        list: A list of starting indices of all occurrences of the motif in the sequence.

    Time Complexity: O(n*m), n = length of sequence, m = length of motif

    Space Complexity: O(n)
    """

    pos = []

    seq_str = str(seq.seq).upper()
    mot = str(motif.seq).upper()


    strt = 0
    
    while True:
        idx = seq_str.find(mot, strt)
        if idx == -1:
            break
        pos.append(idx + 1)
        strt = idx + 1
    return pos
            