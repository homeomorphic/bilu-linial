/* vim: set noet ts=8 sw=8 sts=8: */

#define MAXN WORDSIZE
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include <nauty24r2/nauty.h>
#include <nauty24r2/gtools.h>

#include <lapacke.h>

#define PANIC(msg) \
	do { \
		fputs(msg "\n", stderr); \
		abort(); \
	} while (0)


static void usage(const char *cmd)
{
	fprintf(stderr,
		"%s: verify the Bilu-Linial conjecture for graphs.\n"
		"\n"
		"%s <infile>\n"
		"\n"
		"  infile - a file containing graphs in graph6 format\n",
		cmd, cmd);
}

static int regularity_degree(graph *g, int n)
{
	if (n == 0) {
		return 0;
	}
	int deg0 = POPCOUNT(*GRAPHROW(g, 0, MAXM));
	for (int i = 1; i < n; i++) {
		set *gv = GRAPHROW(g, i, MAXM);
		int deg = POPCOUNT(*gv);
		if (deg != deg0) {
			return -1; /* not regular */
		}
	}
	return deg0;
}

static void dump_graph(graph *g, int n)
{
	for (int i = 0; i < n; i++) {
		set *gv = GRAPHROW(g, i, MAXM);
		for (int j = 0; j < n; j++) {
			fputs(ISELEMENT(gv, j) ? "+" : ".", stderr);
		}
		fputs("\n", stderr);
	}

}

static bool has_small_eigenvalues(graph *restrict g, graph *restrict subg, int n, int d)
{
	double *mat = malloc(n * n * sizeof *mat);
	double *eig = malloc(n * sizeof *eig);
	if (!mat || !eig) {
		PANIC("malloc failure");
	}

	for (int i = 0; i < n; i++) {
		set *gv = GRAPHROW(g, i, MAXM);
		set *subgv = GRAPHROW(subg, i, MAXM);
		for (int j = 0; j < n; j++) {
			if (!ISELEMENT(gv, j)) {
				mat[i + j * n] = 0.0;
			} else if (ISELEMENT(subgv, j)) {
				mat[i + j * n] = 1.0;
			} else {
				mat[i + j * n] = -1.0;
			}
		}
		eig[i] = 0.0/0.0;
	}

	lapack_int info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 
		'N', /* eigenvalues only */
		'U', /* upper triangular is stored */
		n, /* order of the matrix */
		mat, /* symmetric matrix */
		n, /* leading dimension */
		eig /* eigenvalues */);

	if (info > 0) {
		PANIC("LAPACKE_dsyev failed");
	}


	/* Is abs(lambda) <= 2 sqrt(d-1) for each
	   eigenvalue lambda? */
	double max_eig = 2 * sqrt(d-1);
	bool result = true;
	for (int i = 0; i < n; i++) {
		if (fabs(eig[i]) > max_eig) {
			result = false;
			break;
		}
	}

	free(eig);
	free(mat);
	return result;
}

static bool go(graph *restrict g, graph *restrict subg, int n, int d, int mine)
{
	/* Find next edge. */
	for (int e = mine; e < n * n; e++) {
		int i = e / n;
		int j = e % n;
		if (j <= i) {
			continue;
		}
		set *gv = GRAPHROW(g, i, MAXM);
		if (!ISELEMENT(gv, j)) {
			continue;
		}

		set *subgv;
		/* add (i, j) and recurse */
		subgv = GRAPHROW(subg, i, MAXM);
		ADDELEMENT(subgv, j);
		subgv = GRAPHROW(subg, j, MAXM);
		ADDELEMENT(subgv, i);
		if (go(g, subg, n, d, e + 1)) {
			return true;
		}

		/* remove (i, j) and recurse */
		subgv = GRAPHROW(subg, i, MAXM);
		DELELEMENT(subgv, j);
		subgv = GRAPHROW(subg, j, MAXM);
		DELELEMENT(subgv, i);
		if (go(g, subg, n, d, e + 1)) {
			return true;
		}

		return false;
	}

	/* No edges found. */
	return has_small_eigenvalues(g, subg, n, d);
}

static void test_bilu_linial(graph *g, int n)
{
	int d = regularity_degree(g, n);	
	if (d < 0) {
		PANIC("input contained a graph which is not regular");
	}

	if (d <= 1) {
		/* the Bilu-Linial conjecture is vacuously
		   true for graphs of degree <= 1. */
		return;
	}

	graph subg[MAXN];
	for (int i = 0; i < n; i++) {
		set *gv = GRAPHROW(subg, i, MAXM);
		EMPTYSET(gv, MAXM);
	}

	bool satisfied = go(g, subg, n, d, 0);

	if (!satisfied) {
		fputs(" ** COUNTEREXAMPLE **\n", stderr);
		dump_graph(g, n);
		/* We can call this a success, I guess. Of course
		   rounding errors could still have occurred. */
		exit(EXIT_SUCCESS);
	}
}


int main(int argc, char **argv)
{
	if (MAXM != 1) {
		/* This might not be a disaster, but this would shatter my
		   world view. Let's check it. */
		PANIC("MAXM != 1 ?!");
	}
	nauty_check(WORDSIZE, MAXM, MAXN, NAUTYVERSIONID);

	if (argc != 2) {
		usage(argv[0]);
		return EXIT_FAILURE;
	}

	FILE *f = fopen(argv[1], "r");
	if (!f) {
		PANIC("could not open file");
	}

	int i = 0;
	for (;;) {
		graph *g;
		int m, n;
		g = readg(f, NULL, 0, &m, &n);
		if (!g) {
			break;
		}
		if (n > MAXN) {
			PANIC("encountered an input graph which is too large");
		}

		test_bilu_linial(g, n);

		free(g);
		i++;
		if (i % 100 == 0) {
			fprintf(stderr, "Processed %d graphs now.\n", i);
		}
	}

	fclose(f);

	fprintf(stderr, "%d graphs processed. No counterexamples found.\n", i);

	return EXIT_SUCCESS;
}

