from Bio import SeqIO
from needleman_wunsch import generate_array, analyze_alignment

def alignment_score(aligned1, aligned2, match_score, mismatch_penalty, gap_penalty):
    total_score = 0
    for x, y in zip(aligned1, aligned2):
        if x == y and x != '-' and y != '-':
            total_score += match_score
        elif x == '-' or y == '-':
            total_score += gap_penalty
        else:
            total_score += mismatch_penalty
    return total_score

def align_two_sequences(seq1, seq2, match, mismatch, gap):
    matrix = generate_array(seq1, seq2, match, mismatch, gap)
    if matrix is None:
        raise ValueError("Invalid input sequences.")
    aligned1, aligned2 = analyze_alignment(matrix, match, mismatch, gap, seq1, seq2)
    score = alignment_score(aligned1, aligned2, match, mismatch, gap)
    return score, ''.join(aligned1), ''.join(aligned2)

def get_from_fasta_file(file_name):
    sequences = []
    for record in SeqIO.parse(file_name, "fasta"):
        sequences.append(str(record.seq))
    if len(sequences) < 2:
        raise ValueError("At least two sequences required.")
    return sequences

def create_distance_matrix(sequences, match, mismatch, gap):
    n = len(sequences)
    distance_matrix = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            score, _, _ = align_two_sequences(sequences[i], sequences[j], match, mismatch, gap)
            distance_matrix[i][j] = score
            distance_matrix[j][i] = score
    return distance_matrix

def find_center_sequence(distance_matrix):
    scores = [sum(row) for row in distance_matrix]
    return scores.index(max(scores))

def merge_alignments(aligned_seqs, center_index, seq_index, new_center, new_seq):
    merged = []
    pointer_old = 0
    for i in range(len(new_center)):
        if new_center[i] == '-':
            for j in range(len(aligned_seqs)):
                if aligned_seqs[j]:
                    aligned_seqs[j] = aligned_seqs[j][:pointer_old] + '-' + aligned_seqs[j][pointer_old:]
            pointer_old += 1
        else:
            pointer_old += 1
        if len(merged) <= seq_index:
            merged.append('')
        merged_seq = aligned_seqs[seq_index] if seq_index < len(aligned_seqs) else ''
        merged_seq += new_seq[i]
        aligned_seqs[seq_index] = merged_seq
    aligned_seqs[center_index] = new_center
    return aligned_seqs

def align_all_to_center(sequences, center_index, match, mismatch, gap):
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
    with open(filename, 'w') as f:
        for i, seq in enumerate(aligned_sequences):
            f.write(f">Sequence_{i+1}\n{seq}\n")

def save_stats_to_file(filename, stats, match, mismatch, gap):
    with open(filename, 'w') as f:
        f.write("MSA Alignment Statistics\n")
        f.write("========================\n")
        f.write(f"Scoring:\n  Match = {match}, Mismatch = {mismatch}, Gap = {gap}\n\n")
        f.write(f"Matches: {stats['matches']}\n")
        f.write(f"Mismatches: {stats['mismatches']}\n")
        f.write(f"Gaps: {stats['gaps']}\n")
        f.write(f"Identity %: {stats['identity_percent']}%\n")


def save_full_report(filename, aligned_sequences, stats, match, mismatch, gap):
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
    distance_matrix = create_distance_matrix(sequences, match, mismatch, gap)
    center_index = find_center_sequence(distance_matrix)
    return align_all_to_center(sequences, center_index, match, mismatch, gap)

def main():
    match = int(input("Match score: "))
    mismatch = int(input("Mismatch penalty: "))
    gap = int(input("Gap penalty: "))

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

    msa = center_star_msa(sequences, match, mismatch, gap)
    print("\nFinal Multiple Sequence Alignment:")
    for i, seq in enumerate(msa):
        print(f"Sequence {i+1}: {seq}")

    stats = compute_alignment_stats(msa)
    save_full_report("msa_report.txt", msa, stats, match, mismatch, gap)

    print("\nMSA + statistics saved to 'msa_report.txt'")


if __name__ == "__main__":
    main()
