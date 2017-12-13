/*
 * CS 481
 * Ata Deniz Aydin
 * 21502637
 *
 * Construct phylogenetic tree from distance matrix using UPGMA.
 */

#include <stdio.h> // write tree to file

// dist nxn distance matrix, names array of n strings
// writes tree to file
void upgma(FILE *f, double **dist, char **names, int n);
