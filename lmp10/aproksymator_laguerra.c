#include "makespl.h"
#include "piv_ge_solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#define DEBUG
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
double delta = 1.0e-6;
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
double calka_trapez(double a, double b, double ya, double yb)
{
	return (ya + yb) * (b - a) / 2;
}
double pochodna(int k, int n, int alfa, double x)
{
	int a = k % 2 == 0 ? 1 : -1;
	if (k <= n)
		return a * laguerr(n - k, alfa + k, x);
	return 0;
}
double fi(int i, double x)
{	
	//wydaje sie zbyt proste by byc realne ale chuj z tym
	return laguerr(i, 0, x);
}

/* Pierwsza pochodna fi */
double dfi(int i, double x)
{
	/*i double h = (b - a) / (n - 1);
	int hi [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx [5];
	int j; */
	return pochodna(1, i, 0, x);
}
//chwilowo wyjebalem reszte pochodnych zeby sie kompilowalo

double d2fi(int i, double x)
{ 
	/*double h = (b - a) / (n - 1);
	int hi [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx [5];
	int j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[4]) 
		return pochodna(dfi(a, b, n, i, x), dfi(a, b, n, i, x + dx), dx);
		*/
	return pochodna(2, i, 0, x);
		
}


double d3fi(int i, double x)
{ 
/*	double h = (b - a) / (n - 1);
	int hi [5] = {i - 2, i - 1, i, i + 1, i + 2};
	double hx [5];
	int j;

	for (j = 0; j < 5; j++)
		hx[j] = a + h * hi[j];

	if ((x < hx[0]) || (x > hx[4]))
		return 0;
	else if (x >= hx[0] && x <= hx[4])
		return pochodna(d2fi(i, x), d2fi(a, b, n, i, x + dx), dx);
*/
	return pochodna(3, i, 0, x);
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
		for (i = 0; i < TESTBASE; i++) {
			fprintf(tst, "%g", a + i * dx);
			for (j = 0; j < nb; j++) {
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

	for (j = 0; j < nb; j++) {
		for (i = 0; i < nb; i++)
			for (k = 0; k < pts->n; k++)
				add_to_entry_matrix(eqs, j, i, fi(i, x[k]) * fi( j, x[k]));
					//poprawilem w fi na 2 arguemnty zamiast 5
		for (k = 0; k < pts->n; k++)
			add_to_entry_matrix(eqs, j, nb, y[k] * fi(j, x[k]));
	}

#ifdef DEBUG
	printf("1\n");
	write_matrix(eqs, stdout);
#endif

	if (piv_ge_solver(eqs)) {
		spl->n = 0;
		return;
	}
#ifdef DEBUG 
	printf("2\n");
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
				spl->f[i]  += ck * fi  (k, xx);
				spl->f1[i] += ck * dfi (k, xx);
				spl->f2[i] += ck * d2fi(k, xx);
				spl->f3[i] += ck * d3fi(k, xx); 
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
							yi += get_entry_matrix(eqs, k, nb) * fi(k, xi);
							dyi += get_entry_matrix(eqs, k, nb) * dfi(k, xi);
							d2yi += get_entry_matrix(eqs, k, nb) * d2fi(k, xi);
							d3yi += get_entry_matrix(eqs, k, nb) * d3fi(k, xi);
			}
			fprintf(tst, "%g %g %g %g %g\n", xi, yi, dyi, d2yi, d3yi );
		}
		fclose(tst);
	}
#endif

}
