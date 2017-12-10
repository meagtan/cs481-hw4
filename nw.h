/*
 * CS 481
 * Ata Deniz Aydin
 * 21502637
 *
 * Header file for Needleman-Wunsch algorithm.
 */

#include <stdio.h> // for FILE

// bit representation of nucleotides
// A - 0x41, C - 0x43, G - 0x47, T - 0x54 (0x20 added for lowercase)
// / 2 0x20,     0x21,     0x23,     0x2A
// % 4   00,       01,       11,       10
#define BITS(c) (((c) >> 1) & 3)

// match, mismatch and gap score
#define MATCH 4
#define MISMATCH -3
#define GAP -2

#define SUB(c1, c2) (((c1) == (c2)) ? MATCH : MISMATCH)

// Smith-Waterman algorithm using linear gap penalties and affine gap penalties
int naivegap(char *seq1, char *seq2, int n, int m);

// write traceback to file starting from (besti, bestj) using scoring and traceback matrices
// void traceback(FILE *out, char *seq1, char *seq2, int besti, int bestj, int **scores, int **bestis, int **bestjs);
