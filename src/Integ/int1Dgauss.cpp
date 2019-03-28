#include "stdafx.h"
#include "int1Dgauss.h"

//!
//! Konstruktor. deg - stopien kwadratury Int1DGauss
//!
Int1DGauss::Int1DGauss(size_t deg) : m_deg(deg)
{
	m_w.resize(m_deg);
	m_t.resize(m_deg);

	util::gauleg(-1, 1, m_t, m_w, m_deg);
}

//!
//! Destruktro
//!
Int1DGauss::~Int1DGauss(void)
{
}

//!
//! Zwraca wartosc calki $\int_a^b f(t) dt$
//! Calka obliczana jest adaptacyjnie z wykorzystaniem kwadratur Gaussa.
//! Precyzja obliczen wynosi $absTol$
//!
double Int1DGauss::Adapt(const util::Fun1D& f, double a, double b, double errAbs) const
{
static const double ZERO = 1e-16;
double vl, vr, v;

	assert(a < b);

	vl = Calc(f, a, 0.5 * (a + b));
	vr = Calc(f, 0.5 * (a + b), b);
	if(fabs(vl + vr) < ZERO)
		return 0;

	v = Calc(f, a, b);
	if(fabs(vl + vr - v) < errAbs)
		return vl + vr;

	vl = Adapt(f, a, 0.5 * (a + b), errAbs);
	vr = Adapt(f, 0.5 * (a + b), b, errAbs);

	return vl + vr;
}

//!
//! Oblicza wartosc funkcji "f" na przdziale [a, b].
//!
double Int1DGauss::Calc(const util::Fun1D& f, double a, double b) const
{
double p, q, sum;

	// Przeskalowanie do przedzialu [-1, 1].
	q = 0.5 * (a + b);
	p = 0.5 * (b - a);
	
	sum = 0;
	for(size_t i = 0; i < m_deg; i++)
		sum += m_w[i] * f.Get(p * m_t[i] + q);

	return p * sum;
}

/*
//! 
//! Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
//! arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Int1DGauss-
//! Legendre n-point quadrature formula.
//! eps - dokladnosc obliczania wspolrzednych kwadratur Gaussa
//! 
void Int1DGauss::gauleg(double x1, double x2, Vec& x, Vec& w, size_t n)
{
const double eps = 100.0 * DBL_MIN;
size_t m, j, i;
double z1, z, xm, xl, pp, p3, p2, p1; // High precision is a good idea for this routine.
const double pi = 4 * atan(1.);

	assert(x.size() == n);
	assert(w.size() == n);

	m = (n + 1) / 2; // The roots are symmetric in the interval, so
	xm = 0.5 * (x2 + x1); //we only have to find half of them.
	xl = 0.5 * (x2 - x1);
	for(i = 1; i <= m; i++) // Loop over the desired roots.
	{ 
		z = cos(pi * (i - 0.25) / (n + 0.5));
		// Starting with the above approximation to the ith root, we enter the main loop of
		// refinement by Newton's method.
		do
		{
			p1 = 1.0;
			p2 = 0.0;
			for(j = 1; j <= n; j++) //Loop up the recurrence relation to get the
			{
				p3 = p2; // Legendre polynomial evaluated at z.
				p2 = p1;
				p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
			}
			// p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
			// by a standard relation involving also p2, the polynomial of one lower order.
			pp = n * (z * p1 - p2) / (z * z - 1.0);
			z1 = z;
			z = z1 - p1 / pp; // Newton's method.
		}
		while(fabs(z - z1) > eps);

		x[i - 1] = xm - xl * z; // Scale the root to the desired interval,
		x[n - i] = xm + xl * z; // and put in its symmetric counterpart.
		w[i - 1] = 2.0 * xl / ((1.0 - z * z) * pp * pp); //Compute the weight
		w[n - i] = w[i - 1]; //and its symmetric counterpart.
	}
}
*/

