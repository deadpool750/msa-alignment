from Bio import SeqIO
from needleman_wunsch import generate_array, analyze_alignment

def alignment_score(aligned1, aligned2, match_score, mismatch_penalty, gap_penalty):
    """
    Calculate the total alignment score between two aligned sequences.
    Applies the scoring scheme based on matches, mismatches, and gaps.
    """
    #calculates the total alignment score based on a scoring scheme
    total_score = 0
    for x, y in zip(aligned1, aligned2):
        #match case
        if x == y and x != '-' and y != '-':
            total_score += match_score
        #gap case
        elif x == '-' or y == '-':
            total_score += gap_penalty
        #mismatch case
        else:
            total_score += mismatch_penalty
    return total_score

def align_two_sequences(seq1, seq2, match, mismatch, gap):
    """
    Aligns two sequences using the Needleman-Wunsch algorithm and returns
    the alignment score and aligned sequences.
    """
    #aligns two sequences using needleman-wunsch and returns score and alignments
    matrix = generate_array(seq1, seq2, match, mismatch, gap)
    if matrix is None:
        raise ValueError("Invalid input sequences.")
    aligned1, aligned2 = analyze_alignment(matrix, match, mismatch, gap, seq1, seq2)
    score = alignment_score(aligned1, aligned2, match, mismatch, gap)
    return score, ''.join(aligned1), ''.join(aligned2)

def get_from_fasta_file(file_name):
    """
    Reads sequences from a FASTA file and returns them in a list.
    Requires at least two sequences.
    """
    #reads sequences from a fasta file and returns them in a list
    sequences = []
    for record in SeqIO.parse(file_name, "fasta"):
        sequences.append(str(record.seq))
    if len(sequences) < 2:
        raise ValueError("At least two sequences required.")
    return sequences

def create_distance_matrix(sequences, match, mismatch, gap):
    """
    Creates a distance matrix of pairwise alignment scores
    between all sequences using the Needleman-Wunsch algorithm.
    """
    #builds a pairwise distance matrix of alignment scores
    n = len(sequences)
    distance_matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            score, _, _ = align_two_sequences(sequences[i], sequences[j], match, mismatch, gap)
            distance_matrix[i][j] = score
            distance_matrix[j][i] = score
    return distance_matrix

def find_center_sequence(distance_matrix):
    """
    Finds and returns the index of the center sequence,
    i.e., the sequence with the highest total similarity to others.
    """
    #selects the center sequence index with highest total similarity
    scores = [sum(row) for row in distance_matrix]
    return scores.index(max(scores))

def merge_alignments(aligned_seqs, center_index, seq_index, new_center, new_seq):
    """
    Merges a new aligned sequence into the existing set of aligned sequences,
    inserting gaps as necessary to maintain global alignment.
    """
    #merges aligned sequences into center alignment preserving gaps
    merged = []
    pointer_old = 0
    for i in range(len(new_center)):
        #if there's a new gap in the center, insert gaps in all aligned sequences
        if new_center[i] == '-':
            for j in range(len(aligned_seqs)):
                if aligned_seqs[j]:
                    aligned_seqs[j] = aligned_seqs[j][:pointer_old] + '-' + aligned_seqs[j][pointer_old:]
            pointer_old += 1
        else:
            pointer_old += 1
        #ensure the list is long enough
        if len(merged) <= seq_index:
            merged.append('')
        merged_seq = aligned_seqs[seq_index] if seq_index < len(aligned_seqs) else ''
        merged_seq += new_seq[i]
        aligned_seqs[seq_index] = merged_seq
    aligned_seqs[center_index] = new_center
    return aligned_seqs

def align_all_to_center(sequences, center_index, match, mismatch, gap):
    """
    Aligns all sequences to the chosen center sequence using the center-star method.
    Returns the final aligned sequence list.
    """
    #aligns all sequences to the center using center-star method
    n = len(sequences)
    aligned_seqs = [''] * n
    aligned_seqs[center_index] = sequences[center_index]
    current_center = sequences[center_index]
    for i in range(n):
        if i == center_index:
            continue
        _, new_center, new_seq = align_two_sequences(current_center, sequences[i], match, mismatch, gap)
        aligned_seqs = merge_alignments(aligned_seqs, center_index, i, new_center, new_seq)
        current_center = new_center
    return aligned_seqs

def compute_alignment_stats(aligned_sequences):
    """
    Computes alignment statistics across all aligned sequences:
    number of matches, mismatches, gaps, and identity percentage.
    """
    #computes match/mismatch/gap stats and identity percentage
    match_count = 0
    mismatch_count = 0
    gap_count = 0
    total_cols = len(aligned_sequences[0])
    num_seqs = len(aligned_sequences)

    for i in range(total_cols):
        column = [seq[i] for seq in aligned_sequences]
        for j in range(num_seqs):
            for k in range(j + 1, num_seqs):
                a, b = column[j], column[k]
                if a == '-' or b == '-':
                    gap_count += 1
                elif a == b:
                    match_count += 1
                else:
                    mismatch_count += 1

    total_comparisons = (num_seqs * (num_seqs - 1) // 2) * total_cols
    identity = (match_count / total_comparisons) * 100 if total_comparisons > 0 else 0

    return {
        "identity_percent": round(identity, 2),
        "matches": match_count,
        "mismatches": mismatch_count,
        "gaps": gap_count
    }

def save_msa_to_fasta(filename, aligned_sequences):
    """
    Saves aligned sequences to a file in FASTA format.
    """
    #saves aligned sequences to a fasta-format file
    with open(filename, 'w') as f:
        for i, seq in enumerate(aligned_sequences):
            f.write(f">Sequence_{i+1}\n{seq}\n")

def save_stats_to_file(filename, stats, match, mismatch, gap):
    """
    Saves alignment statistics and scoring parameters to a text file.
    """
    #saves alignment statistics to a separate file
    with open(filename, 'w') as f:
        f.write("MSA Alignment Statistics\n")
        f.write("========================\n")
        f.write(f"Scoring:\n  Match = {match}, Mismatch = {mismatch}, Gap = {gap}\n\n")
        f.write(f"Matches: {stats['matches']}\n")
        f.write(f"Mismatches: {stats['mismatches']}\n")
        f.write(f"Gaps: {stats['gaps']}\n")
        f.write(f"Identity %: {stats['identity_percent']}%\n")

def save_full_report(filename, aligned_sequences, stats, match, mismatch, gap):
    """
    Saves both the aligned sequences and the alignment statistics to a single report file.
    """
    #saves both alignment and stats to a combined text report
    with open(filename, 'w') as f:
        f.write("Multiple Sequence Alignment Report\n")
        f.write("=================================\n\n")

        f.write("Aligned Sequences:\n")
        for i, seq in enumerate(aligned_sequences):
            f.write(f">Sequence_{i + 1}\n{seq}\n")

        f.write("\nScoring:\n")
        f.write(f"  Match = {match}\n")
        f.write(f"  Mismatch = {mismatch}\n")
        f.write(f"  Gap = {gap}\n\n")

        f.write("Alignment Statistics:\n")
        f.write(f"  Matches: {stats['matches']}\n")
        f.write(f"  Mismatches: {stats['mismatches']}\n")
        f.write(f"  Gaps: {stats['gaps']}\n")
        f.write(f"  Identity %: {stats['identity_percent']}%\n")

def center_star_msa(sequences, match, mismatch, gap):
    """
    Performs the center-star multiple sequence alignment process.
    Returns the final aligned sequences.
    """
    #wraps the full MSA process using center-star strategy
    distance_matrix = create_distance_matrix(sequences, match, mismatch, gap)
    center_index = find_center_sequence(distance_matrix)
    return align_all_to_center(sequences, center_index, match, mismatch, gap)

def main():
    """
    Main function that interacts with the user, performs MSA,
    and saves the result to a report file.
    """
    #get scoring scheme from user
    match = int(input("Match score: "))
    mismatch = int(input("Mismatch penalty: "))
    gap = int(input("Gap penalty: "))

    #choose data source
    source = input("Enter 'fasta' to load from file or 'manual' to enter sequences: ").strip().lower()
    if source == "fasta":
        fasta_file = input("Enter FASTA file name: ")
        sequences = get_from_fasta_file(fasta_file)
    elif source == "manual":
        sequences = []
        while True:
            seq = input("Enter sequence (or press Enter to finish): ").strip().upper()
            if not seq:
                break
            sequences.append(seq)
    else:
        print("Invalid input.")
        return

    #perform center-star MSA
    msa = center_star_msa(sequences, match, mismatch, gap)

    #print aligned sequences to console
    print("\nFinal Multiple Sequence Alignment:")
    for i, seq in enumerate(msa):
        print(f"Sequence {i+1}: {seq}")

    #compute stats and save full report
    stats = compute_alignment_stats(msa)
    save_full_report("msa_report.txt", msa, stats, match, mismatch, gap)

    #confirm saving
    print("\nMSA + statistics saved to 'msa_report.txt'")

#if the script is run directly, call main
if __name__ == "__main__":
    main()
