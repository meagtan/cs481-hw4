#include <stdlib.h>

#include "upgma.h"

// output upgma recursively
void upgmaout(FILE *f, int i, int m, int *inodes, int *jnodes, int *iprevs, int *jprevs, double *heights, char **names);

// modifies dist
void upgma(FILE *f, double **dists, char **names, int n)
{
	int i, j, k, l;
	int mini, minj;

	// found a way to implement UPGMA without actual trees
	// (but still effectively treelike structures are stored in arrays)
	// n-1 clusters are formed in total
	// for each cluster m, store the smallest and largest index that is part of the cluster
	//  as well as the previous clusters containing each index
	// for each index i, store and update the largest cluster containing it
	// after clustering, each element in the cluster is represented by its minimum element
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

		// average the distance vectors of i and j, update the rest
		// is it necessary to go through entire matrix, or just to average ith and jth rows? need consistent matrix in order to calculate minimum
		for (k = 0; k < n; ++k) {
			for (l = k+1; l < n; ++l) {
				// average the distances from and to i, j
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
				// the last disjunct guarantees that mini is the minimum possible cluster member and minj is the maximum possible
				if (dists[k][l] != 0 && (dists[mini][minj] == 0 || dists[k][l] < dists[mini][minj] || (dists[k][l] == dists[mini][minj] && k <= mini && minj <= l))) {
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
	if (m == -1) { // no cluster, output node name
		fprintf(f, "%s", names[i]);
	} else { // if cluster m contains (i,j), print each cluster
		int i = inodes[m], iprev = iprevs[m],
		    j = jnodes[m], jprev = jprevs[m];

		// subtree with smaller height comes first (has been clustered earlier)
		if (iprev != -1 && (jprev == -1 || heights[jprev] < heights[iprev])) {
			i = jnodes[m]; iprev = jprevs[m];
			j = inodes[m]; jprev = iprevs[m];
		}

		// distance of i and j to parent node
		double iheight = heights[m] - (iprev != -1) * heights[iprev],
		       jheight = heights[m] - (jprev != -1) * heights[jprev];

		// "remove" cluster from tree, recursively print each subtree
		fprintf(f, "[");
		upgmaout(f, i, iprev, inodes, jnodes, iprevs, jprevs, heights, names);
		fprintf(f, ":%g%s-", iheight, iheight == (int) iheight ? ".0" : ""); // pad integer height with decimal
		upgmaout(f, j, jprev, inodes, jnodes, iprevs, jprevs, heights, names);
		fprintf(f, ":%g%s]", jheight, jheight == (int) jheight ? ".0" : "");
	}
}
