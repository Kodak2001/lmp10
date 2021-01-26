#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#define delta 1.0e-6

int silnia(int n)
{
	if (n == 0)
		return 1;
	if (n == 1)
		return 1;
	return n*silnia(n-1);
}
double laguerr(int n, int alfa, double x)
{
	if (n==0)
		return 1;
	if (n==1)
		return 1 + (double)alfa - x;
	return (((2.0 * n - 1 + alfa - x)  * laguerr(n-1, alfa, x) - (n - 1 + alfa) * laguerr(n - 2, alfa ,x)) / n);
}

double pochodna(int k, int n, int alfa, double x)
{
	int a = k % 2 == 0 ? 1 : -1;
	if (k <= n)
		return a * laguerr(n - k, alfa + k, x);
	return 0;
}

void make_spl(points_t * pts, spline_t * spl, int baza)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	
  
	eqs = make_matrix(nb, nb + 1);

#ifdef DEBUG
#define TESTBASE 500
	{
		FILE           *tst = fopen("debug_base_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for (i = 0; i < TESTBASE; i++) {
			fprintf(tst, "%g", a + i * dx);
			for (j = 0; j < baza; j++) {
				fprintf(tst, " %g", fi  (j, a + i * dx));
				fprintf(tst, " %g", dfi (j, a + i * dx));
				fprintf(tst, " %g", d2fi(j, a + i * dx));
				fprintf(tst, " %g", d3fi(j, a + i * dx));
			}
			fprintf(tst, "\n");
		}
		fclose(tst);
	}
#endif

	for (j = 0; j < baza; j++) {
		for (i = 0; i < baza; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi(i, x[k]) * fi( j, x[k]));
		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, baza, y[k] * fi(j, x[k]));
	}

#ifdef DEBUG
	write_matrix(eqs, stdout);
#endif

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}
#ifdef DEBUG 
	write_matrix(eqs, stdout);
#endif

	if (alloc_spl(spl, baza) == 0) {
		double xx;
		double ck;
		for (i = 0; i < spl->n; i++) {
			xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < baza; k++) {
				ck = get_entry_matrix(eqs, k, baza);
				spl->f[i]  += ck * laguerr (k, 0, xx);
				spl->f1[i] += ck * pochodna (1, k, 0, xx); 
				spl->f2[i] += ck * pochodna (2, k, 0, xx); 
				spl->f3[i] += ck * pochodna (3, k, 0, xx); 
			}
		}
	}

#ifdef DEBUG 
	{
		FILE           *tst = fopen("debug_spline_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		double yi= 0;
		double dyi= 0;
		double d2yi= 0;
		double d3yi= 0;
		double xi;
		for (i = 0; i < TESTBASE; i++) {

			xi= a + i * dx;
			for( k= 0; k < baza; k++ ) {
							yi += get_entry_matrix(eqs, k, baza) *   laguerr (k, 0, xi);
							dyi += get_entry_matrix(eqs, k, baza) *  pochodna (1, k, 0, xi); 
							d2yi += get_entry_matrix(eqs, k, baza) * pochodna (2, k, 0, xi); 
							d3yi += get_entry_matrix(eqs, k, baza) * pochodna (3, k, 0, xi); 
			}
			fprintf(tst, "%g %g %g %g %g\n", xi, yi, dyi, d2yi, d3yi );
		}
		fclose(tst);
	}
#endif

}
