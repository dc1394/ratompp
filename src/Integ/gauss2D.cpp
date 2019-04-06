#include "stdafx.h"
#include "gauss2D.h"


//!
//! Konstruktor
//!
Gauss2D::Gauss2D(double, size_t order)
{
std::vector<double> w1D(order), x1D(order);
size_t i, j, k;

	m_recNo = 0;
	m_funCall = 0;
	GauLeg(-1., 1., x1D, w1D);

	m_w.resize(order * order);
	m_x.resize(order * order);
	m_y.resize(order * order);

	for(k = 0, i = 0; i < order; i++)
	{
		for(j = 0; j < order; j++)
		{
			m_w[k] = w1D[i] * w1D[j];
			m_x[k] = x1D[i];
			m_y[k] = x1D[j];
			k++;
		}
	}
}



//!
//! Zwraca wartosc calki z funkcji "f" po prostokacie. 
//! Algorytm adaptacyjny dzieli prostkoat na CZTERY mniejsze prostokaty.
//!
double Gauss2D::CalcAdapt(const Fun2D& f, double x0, double x1, double y0, double y1)
{
	m_funCall = 0;
	m_recNo = 0;

	return Run(f, x0, x1, y0, y1);
}


double Gauss2D::Run(const Fun2D& f, double x0, double x1, double y0, double y1)
{
double s1, s2, s3, s4, s, sum, xMid, yMid;

	assert(x1 > x0);
	assert(y1 > y0);

	m_recNo += 1;

	s = Calc(f, x0, x1, y0, y1);

	yMid = 0.5 * (y0 + y1);
	xMid = 0.5 * (x0 + x1);

	// Dzile prostokat na cztery mniejsze (podobne) prostokaty
	s1 = Calc(f, x0, xMid, y0, yMid);
	s2 = Calc(f, xMid, x1, y0, yMid);

	s3 = Calc(f, x0, xMid, yMid, y1);
	s4 = Calc(f, xMid, x1, yMid, y1);

	// Suma kwadratur z 4 prostokatow
	sum = s1 + s2 + s3 + s4;

//	// Dokladnosc wzgledna
//	if(fabs(s - sum) < m_epsRel * fabs(s))
//		return sum;

	// Dokladnosc absolutna
	if(fabs(s - sum) < m_epsAbs)
		return sum;

	m_recNo += (4 - 1); // Cztery nowe prostokaty na miejsce jednego starego

	// Wywolanie rekurencyjne
	s1 = Run(f, x0, xMid, y0, yMid);
	s2 = Run(f, xMid, x1, y0, yMid);

	s3 = Run(f, x0, xMid, yMid, y1);
	s4 = Run(f, xMid, x1, yMid, y1);

	return s1 + s2 + s3 + s4;
}


//!
//! Zwraca wartosc calki z funkcji "f" po prostokacie. Wspolrzedne prostokata 
//! podane sa na rysunku ponizej. 
//! Calka liczona jest na prostokacie uzywajac kwadratur Gaussa.
//!
//!     +--------+ <-- y1
//!     |        |
//!     |        |
//!     |        |
//!     |        |
//!     +--------+ <-- y0
//!     ^        ^
//!     x0       x1
//!
double Gauss2D::Calc(const Fun2D& f, double x0, double x1, double y0, double y1)
{
double ax, bx, ay, by;
double t, s, v = 0.;

	// Wspolczynniki przekstalcenia x(t) = ax * t + bx, -1 <= t <= 1
	ax = 0.5 * (x1 - x0);
	bx = 0.5 * (x1 + x0);

	// Wspolczynniki przekstalcenia y(s) = ay * s + by, -1 <= s <= 1
	ay = 0.5 * (y1 - y0);
	by = 0.5 * (y1 + y0);

	for(size_t i = 0; i < m_w.size(); i++)
	{
		// Przejscie do nowych wspolrzednych
		t = ax * m_x[i] + bx;
		s = ay * m_y[i] + by;
		v += m_w[i] * f.Get(t, s);
	}

	m_funCall += m_w.size();

	// Mnoze przez jakobian przejscia
	return ax * ay * v;
}

//!
//! Zwraca wartosc calki z funkcji "f" na kwadracie jednostkowym.
//!
double Gauss2D::Gauss(const Fun2D& f) const
{
double v = 0.;

	for(size_t i = 0; i < m_w.size(); i++)
	{
		v += m_w[i] * f.Get(m_x[i], m_y[i]);
	}
	return v;
}


//!
//! Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
//! arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
//! Legendre n-point quadrature formula.
//!
//! Przepisane z Numerical Recipes in C
//! 
void Gauss2D::GauLeg(double x1, double x2, std::vector<double>& x, std::vector<double>& w)
{
// const double pi = 3.1415926535897932384626433832795;
const double pi = 4 * atan(1.);
const double EPS = 3.0e-11;
int m, j, i, n;
double z1, z, xm, xl, pp, p3, p2, p1; // High precision is a good idea for this routine.

	assert(x.size() == w.size());
	n = static_cast<int>(x.size());

	m = (n + 1) / 2; // The roots are symmetric in the interval, so we only have to find half of them.
	xm = 0.5 * (x2 + x1);
	xl = 0.5 * (x2 - x1);
	for(i = 0; i < m; i++) // Loop over the desired roots.
	{ 
		z = cos(pi * (i + 1 - 0.25) / (n + 0.5));
		// Starting with the above approximation to the ith root, we enter the main loop of
		// refinement by Newton’s method.
		do
		{
			p1 = 1.0;
			p2 = 0.0;
			// Loop up the recurrence relation to get the Legendre polynomial evaluated at z.
			for(j = 1; j <= n; j++)
			{ 				
				p3 = p2;
				p2 = p1;
				p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
			}
			// p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
			// by a standard relation involving also p2, the polynomial of one lower order.
			pp = n * (z * p1 - p2) / (z * z - 1.0);
			z1 = z;
			z = z1 - p1 / pp; // Newton’s method.
		}
		while(fabs(z - z1) > EPS);

		x[i] = xm - xl * z;			// Scale the root to the desired interval,
		x[n - 1 - i] = xm + xl * z; // and put in its symmetric counterpart.
		w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp); // Compute the weight
		w[n - 1 - i] = w[i];						// and its symmetric counterpart.
	}
}




double Gauss2D::ErrAbs() const
{
	// assert(0);
	return 0.;
}


double Gauss2D::ErrRel() const
{
	// assert(0);
	return 0.;
}

size_t Gauss2D::RecNo() const
{
	return m_recNo;
}

size_t Gauss2D::FunCall() const
{
	return m_funCall;
}

