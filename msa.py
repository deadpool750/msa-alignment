from Bio import SeqIO

from needleman_wunsch import generate_array, analyze_alignment


def alignment_score(aligned1, aligned2, match_score, mismatch_penalty, gap_penalty):
    """
    Calculates the total alignment score based on aligned sequences.
    """

    #comparing the aligned strands letter by letter and adding scores based on it
    total_score = 0
    for x, y in zip(aligned1, aligned2):
        if x == y and x != '-' and y != '-':
            total_score += match_score
        elif x == '-' or y == '-':
            total_score += gap_penalty
        else:
            total_score += mismatch_penalty
    return total_score

#new method with aligning two sequences
def align_two_sequences(seq1, seq2, match, mismatch, gap):

    #generating a needleman wunsch matrix
    matrix = generate_array(seq1, seq2, match, mismatch, gap)

    #analyzing the alignment of the strands
    aligned1, aligned2 = analyze_alignment(matrix, match, mismatch, gap, seq1, seq2)

    #calculating the alignment score
    score = alignment_score(aligned1, aligned2, match, mismatch, gap)
    return score, ''.join(aligned1), ''.join(aligned2)

def get_from_fasta_file(file_name):
    """
    Reads at least two sequences from a fasta file with the name provided and returns them as a tuple.
    (I am using Biopython library here)
    """

    #initializing the list for the strands
    sequences = []

    #extracting all the strands from the file
    for record in SeqIO.parse(file_name, "fasta"):
        sequences.append(str(record.seq))

    #checking if there are at least 2 values in the file
    if len(sequences) <= 2:
        raise ValueError("At least two values in the fasta file expected!")

    return sequences

def create_distance_matrix(sequences, match_score, mismatch_penalty, gap_penalty):
    """
    Builds a distance/similarity matrix for a list of sequences.
    Each element [i][j] contains the alignment score between sequence i and j.
    """

    #initializing the array length according to strand length
    n = len(sequences)
    distance_matrix = [[0 for _ in range(n)] for _ in range(n)]

    #setting all zeros for the diagonal because the previous function set it to +n
    for i in range(len(sequences)):
        distance_matrix[i][i] = 0

    #diagonal must be all 0 because distance
    #only computing the upper 3 values past the 0 diagonal (bottom ones are the same)
    for i in range(n):
        for j in range(i + 1, n):
            score, _, _ = align_two_sequences(sequences[i], sequences[j], match_score, mismatch_penalty, gap_penalty)

            #values are the same because the matrix is symmetrical
            distance_matrix[i][j] = score
            distance_matrix[j][i] = score

    return distance_matrix

def find_center_sequence(distance_matrix):
    """
    Finds the index of the center sequence (the one most similar to all others).
    """
    #initializing the n length of a distance matrix for iterative operations
    n = len(distance_matrix)

    #list to hold the scores
    total_scores = [0] * n

    #finding the center sequence based on the highest row value from the distance matrix
    for i in range(n):
        for j in range(n):
            total_scores[i] += distance_matrix[i][j]

    #picking the row max value
    center_index = total_scores.index(max(total_scores))
    return center_index




