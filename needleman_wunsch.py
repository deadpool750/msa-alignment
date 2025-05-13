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

def generate_array(strand1, strand2, match_value, mismatch_value, gap_value):
    bases = {'A', 'C', 'G', 'T', 'U'}
    if not set(strand1.upper()).issubset(bases) or not set(strand2.upper()).issubset(bases):
        print("Invalid bases detected!")
        return None

    strand1 = strand1.upper()
    strand2 = strand2.upper()

    rows = len(strand1) + 1
    cols = len(strand2) + 1
    matrix = np.zeros((rows, cols), dtype=int)

    for i in range(1, rows):
        matrix[i][0] = matrix[i - 1][0] + gap_value
    for j in range(1, cols):
        matrix[0][j] = matrix[0][j - 1] + gap_value

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
    aligned1 = []
    aligned2 = []
    i = len(strand1)
    j = len(strand2)

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            c1 = strand1[i - 1]
            c2 = strand2[j - 1]
            match = match_score if c1 == c2 else mismatch_penalty
            diag = matrix[i - 1][j - 1]
            if matrix[i][j] == diag + match:
                aligned1.append(c1)
                aligned2.append(c2)
                i -= 1
                j -= 1
                continue

        if i > 0 and matrix[i][j] == matrix[i - 1][j] + gap_penalty:
            aligned1.append(strand1[i - 1])
            aligned2.append('-')
            i -= 1
            continue

        if j > 0 and matrix[i][j] == matrix[i][j - 1] + gap_penalty:
            aligned1.append('-')
            aligned2.append(strand2[j - 1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()

    return aligned1, aligned2


#need adjustments
def save_alignment_to_file(filename, aligned1, aligned2, match, mismatch, gap, alignment_length, match_count, gap_count, identity_percentage):
    """
    saves all relevant alignment results to a text file
    """

    #structure how it should be saved into a .txt file
    with open(filename, "w") as f:
        f.write("Needleman-Wunsch Alignment Result\n")
        f.write("===============================\n\n")
        f.write("Aligned Sequences:\n")
        f.write("Strand 1: " + ''.join(aligned1) + "\n")
        f.write("Strand 2: " + ''.join(aligned2) + "\n\n")
        f.write(f"Scoring:\n  Match = {match}, Mismatch = {mismatch}, Gap = {gap}\n\n")
        f.write(f"Alignment Length: {alignment_length}\n")
        f.write(f"Number of Matches: {match_count}\n")
        f.write(f"Number of Gaps: {gap_count}\n")
        f.write(f"Identity Percentage: {identity_percentage}%\n")
