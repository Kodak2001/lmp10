#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#define _USE_MATH_DEFINES
#include <math.h>

/* UWAGA: liczbę używanych f. bazowych można ustawić przez wartość
          zmiennej środowiskowej APPROX_BASE_SIZE
*/

/*
 * Funkcje bazowe:
 * n - liczba funkcji 
 * a,b - granice przedzialu aproksymacji 
 * i- numer funkcji 
 * x - wspolrzedna dla ktorej obliczana jest wartosc funkcji
 */
double dx = 0.01;
double laguerr(int n, double x)
{
	if (n==0)
		return 1;
	if (n==1)
		return 1 - x;
	return ((2.0 * n - 1 -x)  * laguerr(n-1, x) - (n-1)*(n-1) * laguerr(n - 2,x));
}
double calka_trapez(double a, double b, double ya, double yb)
{
	return (ya + yb) * (b - a) / 2;
}
double pochodna(double val1, double val2, double dx)
{
	return (val2 - val1) / dx;
}
double fwagowa(double x)
{
	return pow(M_E, -x);	
}
double fi(double a, double b, int n, int i, double x)
{
	
	double h = (b - a) / (n - 1);
	double hx[5];
	int hi[5] = {i-2, i-1, i, i+1, i+2};
	int j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[1])
		return calka_trapez(hx[0], hx[1], laguerr(i, hx[0]), laguerr(i, hx[1]));
	else if (x > hx[1] && x <= hx[2])
		return calka_trapez(hx[1], hx[2], laguerr(i, hx[1]), laguerr(i, hx[2]));
	else if (x > hx[2] && x <= hx[3])
		return calka_trapez(hx[2], hx[3], laguerr(i, hx[2]), laguerr(i, hx[3]));   
	else			/* if (x > hx[3]) && (x <= hx[4]) */
		return calka_trapez(hx[3], hx[4], laguerr(i, hx[3]), laguerr(i, hx[4]));

}

/* Pierwsza pochodna fi */
double dfi(double a, double b, int n, int i, double x)
{
	double h = (b - a) / (n - 1);
	int hi [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx [5];
	int j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[4])
		return pochodna(laguerr(i, x), laguerr(i, x + dx), dx);
}

/* Druga pochodna fi */
double d2fi(double a, double b, int n, int i, double x)
{
	double h = (b - a) / (n - 1);
	int hi [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx [5];
	int j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[4])
		return pochodna(dfi(a, b, n, i, x), dfi(a, b, n, i, x + dx), dx);
}

/* Trzecia pochodna fi */
double d3fi(double a, double b, int n, int i, double x)
{
	double h = (b - a) / (n - 1);
	int hi [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx [5];
	int j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[4])
		return pochodna(d2fi(a, b, n, i, x), d2fi(a, b, n, i, x + dx), dx);

}

void
make_spl(points_t * pts, spline_t * spl)
{

	matrix_t       *eqs= NULL;
	double         *x = pts->x;
	double         *y = pts->y;
	double		a = x[0];
	double		b = x[pts->n - 1];
	int		i, j, k;
	int		nb = pts->n - 3 > 10 ? 10 : pts->n - 3;
  char *nbEnv= getenv( "APPROX_BASE_SIZE" );

	if( nbEnv != NULL && atoi( nbEnv ) > 0 )
		nb = atoi( nbEnv );

	eqs = make_matrix(nb, nb + 1);

#ifdef DEBUG
#define TESTBASE 500
	{
		FILE           *tst = fopen("debug_base_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for( j= 0; j < nb; j++ )
			xfi( a, b, nb, j, tst );
		for (i = 0; i < TESTBASE; i++) {
			fprintf(tst, "%g", a + i * dx);
			for (j = 0; j < nb; j++) {
				fprintf(tst, " %g", fi  (a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", dfi (a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d2fi(a, b, nb, j, a + i * dx));
				fprintf(tst, " %g", d3fi(a, b, nb, j, a + i * dx));
			}
			fprintf(tst, "\n");
		}
		fclose(tst);
	}
#endif

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi(a,b, nb, i, x[k]) * fi(a, b, nb, j, x[k]));

		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(a, b, nb, j, x[k]));
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

	if (alloc_spl(spl, nb) == 0) {
		for (i = 0; i < spl->n; i++) {
			double xx = spl->x[i] = a + i*(b-a)/(spl->n-1);
			xx+= 10.0*DBL_EPSILON;  // zabezpieczenie przed ulokowaniem punktu w poprzednim przedziale
			spl->f[i] = 0;
			spl->f1[i] = 0;
			spl->f2[i] = 0;
			spl->f3[i] = 0;
			for (k = 0; k < nb; k++) {
				double ck = get_entry_matrix(eqs, k, nb);
				spl->f[i]  += ck * fi  (a, b, nb, k, xx);
				spl->f1[i] += ck * dfi (a, b, nb, k, xx);
				spl->f2[i] += ck * d2fi(a, b, nb, k, xx);
				spl->f3[i] += ck * d3fi(a, b, nb, k, xx);
			}
		}
	}

#ifdef DEBUG
	{
		FILE           *tst = fopen("debug_spline_plot.txt", "w");
		double		dx = (b - a) / (TESTBASE - 1);
		for (i = 0; i < TESTBASE; i++) {
			double yi= 0;
			double dyi= 0;
			double d2yi= 0;
			double d3yi= 0;
			double xi= a + i * dx;
			for( k= 0; k < nb; k++ ) {
							yi += get_entry_matrix(eqs, k, nb) * fi(a, b, nb, k, xi);
							dyi += get_entry_matrix(eqs, k, nb) * dfi(a, b, nb, k, xi);
							d2yi += get_entry_matrix(eqs, k, nb) * d2fi(a, b, nb, k, xi);
							d3yi += get_entry_matrix(eqs, k, nb) * d3fi(a, b, nb, k, xi);
			}
			fprintf(tst, "%g %g %g %g %g\n", xi, yi, dyi, d2yi, d3yi );
		}
		fclose(tst);
	}
#endif

}
