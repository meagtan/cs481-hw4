/*
 * CS 481
 * Ata Deniz Aydin
 * 21502637
 *
 * Implementation of Needleman-Wunsch algorithm as given in nw.h.
 */

#include <stdlib.h>

#include "nw.h"

// no need for traceback
int nw(char *seq1, char *seq2, int n, int m)
{
	// scoring matrix
	// the parent of (i,j) is (bestis[i][j], bestjs[i][j])
	int **scores, i, j;

	// allocate matrix
	scores = calloc(n+1, sizeof(int *));
	for (i = 0; i < n+1; ++i)
		scores[i] = calloc(m+1, sizeof(int));

	// fill scoring matrix row by row
	for (i = 0; i < n; ++i) {
		for (j = 0; j < m; ++j) {
			// find best parent for scores[i+1][j+1]
			scores[i+1][j+1] = scores[i][j] + SUB(seq1[i], seq2[j]);
			if (scores[i+1][j+1] < scores[i+1][j] + GAP)
				scores[i+1][j+1] = scores[i+1][j] + GAP;
			if (scores[i+1][j+1] < scores[i][j+1] + GAP)
				scores[i+1][j+1] = scores[i][j+1] + GAP;
		}
	}

	int res = scores[n+1][m+1];

	// free matrices
	for (i = 0; i < n+1; ++i)
		free(scores[i]);
	free(scores);

	return res;
}
