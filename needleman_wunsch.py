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
    """
    Generates the Needleman-Wunsch alignment matrix based on the input strands and scoring values.
    """

    #creating appropriate length and width for the matrix taking into account the space for gap penalties and 0 values
    row_length = len(strand1) + 2
    column_length = len(strand2) + 2

    #validating the strand nucleotides
    bases = {'A', 'C', 'G', 'T', 'U'}
    if not set(strand1).issubset(bases) or not set(strand2).issubset(bases):
        print("Invalid bases detected!")
        return None

    #creating the array originally all with zeros
    matrix = np.zeros((row_length, column_length), dtype='object')

    #placing the two strands into the matrix
    matrix[0, 2:2 + len(strand2)] = list(strand2)
    matrix[2:2 + len(strand1), 0] = list(strand1)

    #placing the zero manually as it will be always in the same place
    matrix[1, 1] = 0

    #placing the gap values in every cell of the 0th row
    n = gap_value
    for i in range(2, row_length):
        matrix[i, 1] = n
        n += gap_value

    #placing the gap values in every cell of the 0th column
    n = gap_value
    for i in range(2, column_length):
        matrix[1, i] = n
        n += gap_value

    #filling in the values using previously created maxvalue_needleman function
    for i in range(2, row_length):
        for j in range(2, column_length):
            maxvalue_needleman(matrix, i, j, match_value, mismatch_value, gap_value)

    return matrix

def analyze_alignment(matrix, match_score, mismatch_penalty, gap_penalty, strand1, strand2):
    #initializing lists to store the aligned sequences
    aligned1 = []
    aligned2 = []

    #setting starting positions at the bottom-right of the matrix
    i = len(strand1) + 1
    j = len(strand2) + 1

    #tracing back through the matrix until reaching the top-left corner
    while i > 1 or j > 1:
        #checking for a diagonal move (match or mismatch)
        if i > 1 and j > 1:
            match = match_score if strand1[i - 2] == strand2[j - 2] else mismatch_penalty
            if matrix[i][j] == matrix[i - 1][j - 1] + match:
                aligned1.append(strand1[i - 2])
                aligned2.append(strand2[j - 2])
                i -= 1
                j -= 1
                continue

        #checking for a move from the top (gap in strand2)
        if i > 1 and matrix[i][j] == matrix[i - 1][j] + gap_penalty:
            aligned1.append(strand1[i - 2])
            aligned2.append("-")
            i -= 1
            continue

        #checking for a move from the left (gap in strand1)
        if j > 1 and matrix[i][j] == matrix[i][j - 1] + gap_penalty:
            aligned1.append("-")
            aligned2.append(strand2[j - 2])
            j -= 1
            continue

    #reversing the sequences to get the final alignment
    aligned1.reverse()
    aligned2.reverse()

    #returning the aligned sequences
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
