import numpy as np


def smith_waterman_match_matrix(seq1, seq2, match=1, mismatch=-1):
    # Matrix for checking if the letters are the same
    matrix_is_match = np.zeros((len(seq1), len(seq2)))
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if seq1[i] == seq2[j]:
                matrix_is_match[i][j] = match
            else:
                matrix_is_match[i][j] = mismatch
    return matrix_is_match


def smith_waterman_initialization(seq1, seq2):
    # Matrix for the Smith-Waterman algorithm
    matrix = np.zeros((len(seq1) + 1, len(seq2) + 1))
    matrices = []  # List of matrices for the animation

    # Initialization - all zeros
    matrices.append(matrix.copy())
    return matrix, matrices


def smith_waterman_filling(seq1, seq2, matrix, matrix_is_match, gap=-2):
    matrices = []  # List of matrices for the animation

    # Filling of the matrix
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            matrix[i][j] = max(0,
                               matrix[i - 1][j - 1] + matrix_is_match[i - 1][j - 1],
                               matrix[i - 1][j] + gap,
                               matrix[i][j - 1] + gap)
            matrices.append(matrix.copy())

    return matrix, matrices


def smith_waterman_traceback(seq1, seq2, matrix, matrix_is_match, gap=-2):
    matrices = []  # List of matrices for the animation

    # Find the highest scoring cell
    max_score = matrix.max()
    max_position = np.where(matrix == max_score)
    ti, tj = max_position[0][0], max_position[1][0]

    # Traceback
    aligned_1 = ""
    aligned_2 = ""

    while matrix[ti][tj] > 0:
        if ti > 0 and tj > 0 and matrix[ti][tj] == matrix[ti - 1][tj - 1] + matrix_is_match[ti - 1][tj - 1]:
            aligned_1 = seq1[ti - 1] + aligned_1
            aligned_2 = seq2[tj - 1] + aligned_2
            ti -= 1
            tj -= 1
        elif ti > 0 and matrix[ti][tj] == matrix[ti - 1][tj] + gap:
            aligned_1 = seq1[ti - 1] + aligned_1
            aligned_2 = "-" + aligned_2
            ti -= 1
        elif tj > 0 and matrix[ti][tj] == matrix[ti][tj - 1] + gap:
            aligned_1 = "-" + aligned_1
            aligned_2 = seq2[tj - 1] + aligned_2
            tj -= 1

        matrices.append(matrix.copy())

    return aligned_1, aligned_2, matrices


def smith_waterman_algorithm(seq1, seq2, gap=-2):
    matrix_is_match = smith_waterman_match_matrix(seq1, seq2)
    all_matrices = []  # List of all matrices for the animation

    matrix, init_matrices = smith_waterman_initialization(seq1, seq2)
    all_matrices.extend(init_matrices)

    matrix, filling_matrices = smith_waterman_filling(seq1, seq2, matrix, matrix_is_match, gap)
    all_matrices.extend(filling_matrices)

    aligned_1, aligned_2, traceback_matrices = smith_waterman_traceback(seq1, seq2, matrix, matrix_is_match, gap)
    all_matrices.extend(traceback_matrices)

    return aligned_1, aligned_2, all_matrices


# Example usage
seq1 = 'ACTG'
seq2 = 'ACTTG'
aligned_1, aligned_2, all_matrices = smith_waterman_algorithm(seq1, seq2)
print("Aligned Sequences:")
print(aligned_1)
print(aligned_2)
# Optionally, you can print the matrices or analyze them as needed.
for matrix in all_matrices:
    print(matrix)
    print()
