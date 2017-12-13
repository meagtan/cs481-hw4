#include "upgma.h"

// output upgma recursively
void upgmaout(FILE *f, int i, int m, int *inodes, int *jnodes, double *heights, char **names);

// modifies dist
void upgma(FILE *f, double **dist, char **names, int n)
{
	int i, j, k, l;
	double min = 0;
	int mini, minj;

	// each internal node obtained from leaves i, j, i < j, may be represented by i
	// n-1 internal nodes in total
	// when an internal node (i,j) is formed, (i,j) is added to the following lists
	// TODO also keep links to last cluster including i, j, and a dynamic index of when each i was last clustered
	int *cluster = malloc((n-1) * sizeof(int)); // cluster[i] : the last time i was clustered
	int *inodes = malloc((n-1) * sizeof(int));
	int *jnodes = malloc((n-1) * sizeof(int));
	double *heights = malloc((n-1) * sizeof(double)); // height of each internal node
	int m = 0; // size of above arrays, number of clusters formed

	// size of each node, initially all 1
	int *sizes = malloc(n * sizeof(int));
	for (i = 0; i < n; ++i) {
		sizes[i] = 1;
		cluster[i] = -1;
	}

	// maintain minimum
	for (i = 0, max = 0; i < n; ++i) {
		for (j = i+1; j < n; ++j) {
			if (dists[i][j] < min) {
				min = dists[i][j];
				mini = i;
				minj = j;
			}
		}
	}

	// combine until one cluster is left
	for (m = 0; m < n-1; ++m) {
		// find (i,j) pair with shortest distance, s.t. i<j, skipping 0 distances
		i = mini;
		j = minj;

		// combine (i,j)
		heights[m] = dist[i][j] / 2;
		/*
		inodes[m]  = i;
		jnodes[m]  = j;
		*/
		inodes[m]  = cluster[i];
		jnodes[m]  = cluster[j];
		cluster[i] = cluster[j] = m;

		// average the distance vectors of i and j, update the rest
		// here, should (mini, minj) be the smallest possible or the largest possible? would it introduce any inconsistency with cluster, inodes, jnodes?
		min = 0;
		for (k = 0; k < n; ++k) {
			for (l = k+1; l < n; ++l) {
				// average the distances from and to i, j
				// since i will be visited before j, after updating i's value just set j's value to i's
				if (k == i && l == j || k == j && l == i)
					dists[k][l] = 0; // guarantees that (i,j) will be skipped last time
				else if (k == i)
					dists[k][l] = (sizes[i] * dists[i][l] + sizes[j] * dists[j][l]) / (sizes[i] + sizes[j]);
				else if (k == j)
					dists[k][l] = dists[i][l];
				else if (l == i)
					dists[k][l] = (sizes[i] * dists[k][i] + sizes[j] * dists[k][j]) / (sizes[i] + sizes[j]);
				else if (l == j)
					dists[k][l] = dists[k][i];
				// maintain minimum, skipping members of same cluster
				if (dists[k][l] != 0 && dists[k][l] < min) {
					min = dists[k][l];
					mini = k;
					minj = l;
				}
			}
		}

		// update cluster sizes
		sizes[i] += sizes[j];
		sizes[j]  = sizes[i];
	}

	// output clusters from heights, inodes and jnodes recursively
	upgmaout(f, 0, m, inodes, jnodes, heights);
}

void upgmaout(FILE *f, int i, int m, int *inodes, int *jnodes, double *heights, char **names)
{
	if (m == 0) {
		// output node name
		fprintf(f, "%s", names[i]);
	} else {
		// "remove" cluster from tree
		fprintf(f, "[");
		upgmaout(f, i, inodes[i], inodes, jnodes, heights);
		fprintf(f, ":%lf-", heights[m]);
		upgmaout(f, j, jnodes[j], inodes, jnodes, heights);
		upgmaout(f, ":%lf]", heights[m]);
	}
}
