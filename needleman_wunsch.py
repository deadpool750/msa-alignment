import numpy as np
from Bio import SeqIO

def maxvalue_needleman(matrix, x_pos, y_pos, match_score, mismatch_penalty, gap_penalty):
    """
    Calculate and set maximum score for a given cell in the matrix
    based on match, mismatch, and gap values provided.
    """

    #first we are checking the diagonal we have a condition based on the matching or mismatching parts of the sequence
    if matrix[x_pos][0] == matrix[0][y_pos]: #match
        diagonal = matrix[x_pos - 1][y_pos - 1] + match_score
    else:
        diagonal = matrix[x_pos - 1][y_pos - 1] + mismatch_penalty #mismatch

    #we calculate the left value
    left = matrix[x_pos - 1][y_pos] + gap_penalty

    #we calculate the top value
    top = matrix[x_pos][y_pos - 1] + gap_penalty

    #the final value set in the matrix is the maximum value out of the three
    matrix[x_pos][y_pos] = max(left, top, diagonal)

def match_or_mismatch(matrix, x_pos, y_pos, match_score, mismatch_penalty):
    """
    Returns match or mismatch value based on the nucleotides provided in the 0th row and column.
    """

    #if the values are the same return match score else return mismatch penalty
    return match_score if matrix[x_pos][0] == matrix[0][y_pos] else mismatch_penalty

def traceback_needleman(matrix, match_score, mismatch_penalty, gap_penalty):
    """
    Traces back the path through the matrix and returns the path values
    and their coordinates.
    """

    #initializing lists for the coordinates and the values
    path_list = []
    path_coordinates = []

    #getting the length and width of the matrix for iteratign through it
    max_length = matrix.shape[0]
    max_width = matrix.shape[1]

    #subtracting 1 from it since we are operating on matrices
    position_x = max_length - 1
    position_y = max_width - 1

    #while loop that stops after we reach coordinates x < 1 or y < 1
    while position_x > 1 or position_y > 1:

        #adding the results to the list
        path_list.append(matrix[position_x, position_y])
        path_coordinates.append((position_x, position_y))

        #first we check the diagonal
        if matrix[position_x][position_y] == matrix[position_x - 1][position_y - 1] + match_or_mismatch(
                matrix, position_x, position_y, match_score, mismatch_penalty):
            position_x -= 1
            position_y -= 1

        #if not the diagonal then top value
        elif matrix[position_x][position_y] == matrix[position_x - 1][position_y] + gap_penalty:
            position_x -= 1

        #if not diagonal and top then left value
        else:
            position_y -= 1

    #adding the [1,1] and [0,0] coordinates and values manually because they are the same everytime
    path_list.append(matrix[1, 1])
    path_coordinates.append((1, 1))

    path_list.append(matrix[0, 0])
    path_coordinates.append((0, 0))

    #coordinates given in the
    return path_list, path_coordinates

import numpy as np

def generate_array(strand1, strand2, match_value, mismatch_value, gap_value):
    """
    Generates the Needleman-Wunsch alignment matrix (numeric only)
    based on input strands and scoring values.
    Includes debug output for invalid characters.
    """

    #clean up sequences: uppercase and remove non-letter characters
    strand1 = ''.join(filter(str.isalpha, strand1.upper()))
    strand2 = ''.join(filter(str.isalpha, strand2.upper()))

    #allowed nucleotide bases
    bases = {'A', 'C', 'G', 'T', 'U'}

    #check for invalid bases
    invalid1 = set(strand1) - bases
    invalid2 = set(strand2) - bases

    if invalid1 or invalid2:
        if invalid1:
            print(f"[ERROR] Invalid characters in Strand 1: {invalid1}")
        if invalid2:
            print(f"[ERROR] Invalid characters in Strand 2: {invalid2}")
        print("Invalid bases detected!")
        return None

    #initialize the alignment matrix
    rows = len(strand1) + 1
    cols = len(strand2) + 1
    matrix = np.zeros((rows, cols), dtype=int)

    #initialize first row and first column with gap penalties
    for i in range(1, rows):
        matrix[i][0] = matrix[i - 1][0] + gap_value
    for j in range(1, cols):
        matrix[0][j] = matrix[0][j - 1] + gap_value

    #fill the matrix based on dynamic programming rules
    for i in range(1, rows):
        for j in range(1, cols):
            if strand1[i - 1] == strand2[j - 1]:
                diagonal = matrix[i - 1][j - 1] + match_value
            else:
                diagonal = matrix[i - 1][j - 1] + mismatch_value
            up = matrix[i - 1][j] + gap_value
            left = matrix[i][j - 1] + gap_value
            matrix[i][j] = max(diagonal, up, left)

    return matrix

def analyze_alignment(matrix, match_score, mismatch_penalty, gap_penalty, strand1, strand2):
    #ensure sequences are clean (strip non-letters and uppercase)
    strand1 = ''.join(filter(str.isalpha, strand1.upper()))
    strand2 = ''.join(filter(str.isalpha, strand2.upper()))

    #initialize lists to store the aligned sequences
    aligned1 = []
    aligned2 = []

    #start from the bottom-right cell of the matrix
    i = len(strand1)
    j = len(strand2)

    #trace back through the matrix until reaching the top-left corner
    while i > 0 or j > 0:
        #check for a diagonal move (match or mismatch)
        if i > 0 and j > 0:
            c1 = strand1[i - 1]
            c2 = strand2[j - 1]
            match = match_score if c1 == c2 else mismatch_penalty
            diag = matrix[i - 1][j - 1]

            #if the current cell equals the diagonal cell + match/mismatch, move diagonally
            if matrix[i][j] == diag + match:
                aligned1.append(c1)
                aligned2.append(c2)
                i -= 1
                j -= 1
                continue

        #check for a vertical move (gap in strand2)
        if i > 0 and matrix[i][j] == matrix[i - 1][j] + gap_penalty:
            aligned1.append(strand1[i - 1])
            aligned2.append('-')
            i -= 1
            continue

        #check for a horizontal move (gap in strand1)
        if j > 0 and matrix[i][j] == matrix[i][j - 1] + gap_penalty:
            aligned1.append('-')
            aligned2.append(strand2[j - 1])
            j -= 1

    #reverse the aligned sequences to correct the order
    aligned1.reverse()
    aligned2.reverse()

    #return the aligned sequences
    return aligned1, aligned2


