#include <stdlib.h>

#include "upgma.h"

// output upgma recursively
void upgmaout(FILE *f, int i, int m, int *inodes, int *jnodes, int *iprevs, int *jprevs, double *heights, char **names);

// modifies dist
void upgma(FILE *f, double **dists, char **names, int n)
{
	int i, j, k, l;
	int mini, minj;

	// each internal node obtained from leaves i, j, i < j, may be represented by i
	// n-1 internal nodes in total
	// when an internal node (i,j) is formed, (i,j) is added to the following lists
	// TODO also keep links to last cluster including i, j, and a dynamic index of when each i was last clustered
	int *sizes = malloc(n * sizeof(int));	// size of each node, initially all 1
	int *cluster = malloc(n * sizeof(int)); // cluster[i] : the last time i was clustered
	int *inodes = malloc((n-1) * sizeof(int)); // (inodes[m], jnodes[m]) clustered at m
	int *jnodes = malloc((n-1) * sizeof(int));
	int *iprevs = malloc((n-1) * sizeof(int)); // iprevs[m] : last time inodes[m] was clustered
	int *jprevs = malloc((n-1) * sizeof(int)); // jprevs[m] : last time jnodes[m] was clustered
	double *heights = malloc((n-1) * sizeof(double)); // height of each internal node
	int m = 0; // size of above arrays, number of clusters formed

	for (i = 0; i < n; ++i) {
		sizes[i] = 1;
		cluster[i] = -1;
	}

	// maintain minimum
	mini = 0; minj = 1;
	for (j = 0; j < n; ++j) {
		for (i = 0; i < j; ++i) {
			if (dists[i][j] < dists[mini][minj] || (dists[i][j] == dists[mini][minj] && i <= mini && minj <= j)) {
				mini = i;
				minj = j;
			}
		}
	}

	// printf("minimum %lf at (%d,%d)\n", dists[mini][minj], mini, minj);

	// combine until one cluster is left
	for (m = 0; m < n-1; ++m) {
		// find (i,j) pair with shortest distance, s.t. i<j, skipping 0 distances
		i = mini;
		j = minj;

		// combine (i,j)
		heights[m] = dists[i][j] / 2;
		inodes[m]  = i;
		jnodes[m]  = j;
		iprevs[m]  = cluster[i];
		jprevs[m]  = cluster[j];
		cluster[i] = cluster[j] = m;

		// printf("(%d,%d) in new cluster with size %d and height %lf\n", i, j, sizes[i] + sizes[j], heights[m]);

		// average the distance vectors of i and j, update the rest
		// here, should (mini, minj) be the smallest possible or the largest possible? would it introduce any inconsistency with cluster, inodes, jnodes?
		// TODO is it necessary to go through entire matrix, or just to average ith and jth rows? need consistent matrix in order to calculate minimum
		for (k = 0; k < n; ++k) {
			for (l = k+1; l < n; ++l) {
				// average the distances from and to i, j
				// TODO also take into account the elements previously clustered with i or j // fixed by replacing k == i by dists[k][i] == 0
				// since i will be visited before j, after updating i's value just set j's value to i's
				if ((dists[k][i] == 0 && dists[l][j] == 0) || (dists[k][j] == 0 && dists[l][i] == 0) || dists[k][l] == 0) // hack
					dists[k][l] = 0; // guarantees that (i,j) will be skipped next time
				else if (k == i)
					dists[k][l] = (sizes[i] * dists[i][l] + sizes[j] * dists[j][l]) / (sizes[i] + sizes[j]);
				else if (dists[k][j] == 0)
					dists[k][l] = dists[i][l];
				else if (l == i || (i < k && l == j)) // first calculation
					dists[k][l] = (sizes[i] * dists[k][i] + sizes[j] * dists[k][j]) / (sizes[i] + sizes[j]);
				else if (dists[l][j] == 0)
					dists[k][l] = dists[k][i];
			}
		}

		mini = 0; minj = 1; // may have zero distance
		for (k = 0; k < n; ++k) {
			for (l = k+1; l < n; ++l) {
				dists[l][k] = dists[k][l];
				// maintain minimum, skipping members of same cluster
				if (dists[k][l] != 0 && (dists[mini][minj] == 0 || dists[k][l] < dists[mini][minj ]|| (dists[k][l] == dists[mini][minj] && k <= mini && minj <= l))) {
					mini = k;
					minj = l;
				}
			}
		}

		/*
		for (k = 0; k < n; ++k) {
			for (l = 0; l < n; ++l)
				printf("%f\t", dists[k][l]);
			printf("\n");
		}
		printf("minimum %lf at (%d,%d)\n", dists[mini][minj], mini, minj);
		*/

		// update cluster sizes
		sizes[i] += sizes[j];
		sizes[j]  = sizes[i];
	}

	/*
	printf("inodes\tjnodes\tiprevs\tjprevs\theights\n");
	for (m = 0; m < n-1; ++m)
		printf("%d\t%d\t%d\t%d\t%g\n", inodes[m], jnodes[m], iprevs[m], jprevs[m], heights[m]);
	*/

	// output clusters from heights, inodes and jnodes recursively
	upgmaout(f, 0, n-2, inodes, jnodes, iprevs, jprevs, heights, names);
	fprintf(f, ":0.0\n");

	free(sizes);
	free(cluster);
	free(inodes);
	free(jnodes);
	free(iprevs);
	free(jprevs);
	free(heights);
}

void upgmaout(FILE *f, int i, int m, int *inodes, int *jnodes, int *iprevs, int *jprevs, double *heights, char **names)
{
	// printf("printing %d at cluster %d\n", i, m);
	if (m == -1) {
		// output node name
		fprintf(f, "%s", names[i]);
	} else {
		int i = inodes[m], iprev = iprevs[m],
		    j = jnodes[m], jprev = jprevs[m];

		// subtree with smaller height comes first
		if (iprev != -1 && (jprev == -1 || heights[jprev] < heights[iprev])) {
			i = jnodes[m]; iprev = jprevs[m];
			j = inodes[m]; jprev = iprevs[m];
		}

		// "remove" cluster from tree, recursively print each subtree
		fprintf(f, "[");
		upgmaout(f, i, iprev, inodes, jnodes, iprevs, jprevs, heights, names);
		fprintf(f, ":%g-", heights[m] - (iprev != -1) * heights[iprev]);
		upgmaout(f, j, jprev, inodes, jnodes, iprevs, jprevs, heights, names);
		fprintf(f, ":%g]", heights[m] - (jprev != -1) * heights[jprev]);
	}
}
