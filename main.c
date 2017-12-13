#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "nw.h"
#include "upgma.h"

#define INPUT "sequences.fasta"
#define OUTPUT "upgma_tree.txt"

// default max number of sequences and size of string buffers
#define MAXLEN 16
#define SEQLEN 256

int main()
{
	double **dists;
	char **names, **seqs;
	int *lens;
	int n = 0, maxlen = MAXLEN;
	FILE *f;

	int i, j;

	// allocate arrays
	names = malloc(MAXLEN * sizeof(char *));
	seqs  = malloc(MAXLEN * sizeof(char *));
	lens  = calloc(MAXLEN, sizeof(int));

	// read file
	f = fopen(INPUT, "r");
	if (!f) {
		fprintf(stderr, "Error: file %s does not exist.\n", INPUT);
		return 1;
	}

	// read each sequence
	fscanf(f, ">"); // skip initial >
	while (!feof(f)) {
		// read name
		names[n] = malloc(SEQLEN); // sizeof(char) == 1
		fscanf(f, "%s\n", names[n]);

		// read sequence from file, skipping whitespace
		seqs[n] = malloc(SEQLEN);
		int siz = SEQLEN;
		do {
			seqs[n][lens[n]] = fgetc(f);		  // write each character read
			if (isalpha(seqs[n][lens[n]])) ++lens[n]; // override whitespace read
			if (lens[n] >= siz) {
				siz <<= 1;
				seqs[n] = realloc(seqs[n], siz);
			}
		} while (!feof(f) && seqs[n][lens[n]] != '>');
		seqs[n][lens[n]] = '\0';

		// extend arrays if necessary
		if (++n >= maxlen) {
			maxlen <<= 1;
			names = realloc(names, maxlen * sizeof(char *));
			seqs  = realloc(seqs,  maxlen * sizeof(char *));
			lens  = realloc(lens,  maxlen * sizeof(int));
		}
	}

	fclose(f);

	// calculate alignment distances
	dists = malloc(n * sizeof(double *));
	for (i = 0; i < n; ++i) {
		dists[i] = malloc(n * sizeof(double));
		dists[i][i] = 0;
		for (j = 0; j < i; ++j)
			dists[i][j] = dists[j][i] = nw(seqs[i], seqs[j], lens[i], lens[j]);
	}

	// output UPGMA

	f = fopen(OUTPUT, "w");
	if (!f) {
		fprintf(stderr, "Error: could not open file %s.\n", OUTPUT);
		return 1;
	}

	upgma(f, dists, names, n);

	fclose(f);

	// deallocate arrays
	for (i = 0; i < n; ++i) {
		free(dists[i]);
		free(names[i]);
		free(seqs[i]);
	}
	free(dists);
	free(names);
	free(seqs);
	free(lens);
}
