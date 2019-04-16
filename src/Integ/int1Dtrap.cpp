#include "stdafx.h"
#include "int1Dtrap.h"


//!
//! Konstruktor.
//! \param "ext" maksymalna liczba eksttrapolacji.
//! \param "eps" dokladnosc obliczen.
//!
Int1DTrap::Int1DTrap(double eps, size_t ext) : m_ext(ext), m_eps(eps)
{
size_t i;

	assert(m_ext > 0);

	// Przedzielenie pamieci
	m_denom = new double[m_ext];
	assert(m_denom != NULL);

	m_p = new double*[m_ext];
	assert(m_p != NULL);
	for(i = 0; i < m_ext; i++)
	{
		m_p[i] = new double[m_ext];
		assert(m_p[i] != NULL);
	}


	// Przypisanie wartosci mianownikowi
	for(i = 0; i < m_ext; i++)
		m_denom[i] = 1. / ((2 << i) - 1.);
}

//!
//! Destruktor.
//!
Int1DTrap::~Int1DTrap()
{
	// Zwolnienie pamieci
	delete [] m_denom;

	for(unsigned long int i = 0; i < m_ext; i++)
		delete [] m_p[i];
	delete [] m_p;
}


//!
//! Oblicza wartosc funkcji "f" na przdziale [a, b].
//! Zwraca wartosc calki obliczona metoda trapezow z wykorzystaniem ektrapolacji Romberga.
//! Dane do ekrapolacji gromadzone sa w macierzy.
double Int1DTrap::Calc(util::Fun1D const& f, double a, double b) const
{
unsigned long int i, k, divNo;
// double v1, v2;

	assert(a <= b);

	if(a == b)
		return 0.;

	divNo = 1; 
	m_p[0][0] = Val(f, a, b, divNo);

	for(i = 1; i < m_ext; i++)
	{
		divNo = divNo << 1;

		// p[i][0] = Val(f, a, b, divNo);
		m_p[i][0] = Fast(f, a, b, divNo, m_p[i - 1][0]);

		// Extrapoluje wartosc calki algorytmem Neville'a korzytsjac z danych zgromadzonych
		// w macierzy "m_p". Extrapolacja zrobiona jest dla wiersza "i".
		for(k = 1; k <= i; k++)
		{
			m_p[i][k] = m_p[i][k - 1] + (m_p[i][k - 1] - m_p[i - 1][k - 1]) * m_denom[k];
		}

		if(i > 5) // Piersze oszacowania moga byc bledne
		{
			// Aby napisy byly czytelniejsze
			const double v1 = m_p[i][i];
			const double v2 = m_p[i - 1][i - 1];

			// if((::fabs(v1 - v2) < m_eps * ::fabs(v1)) || (v1 == 0 && v2 == 0.))
			if((::fabs(v1 - v2) < m_eps) || (v1 == 0 && v2 == 0.))
				return v1;
		}

	}
	// printf("Za mala liczba extrapolacji. Zwieksz maxExt = %ld\n", m_ext);
	return m_p[m_ext - 1][m_ext - 1];

}


//!
//! Oblicza wartosc calki z funkcji "f" metoda trapezow z zadana liczba podzialow odcinka [a, b]. 
//!
double Int1DTrap::Val(util::Fun1D const& f, double a, double b, long int divNo) const
{
long int i;
double res, h, x;

	assert(b > a);

	h = (b - a) / divNo;
	res = (f.Get(a) + f.Get(b)) * 0.5;

	x = a;
	for(i = 0; i < divNo - 1; i++)
	{
		x += h;
		res += f.Get(x);
	}
	return h * res;
}

//!
//! Oblicza wartosc calki z funkcji "f" metoda trapezow z zadana liczba podzialow odcinka [a, b].
//! Aby przyspieszyc obliczenia wykorzystywana jest wartosc calki obliczona na mniejszej liczbie
//! podzialow odcinka [a, b]. Wartosc ta przekazywana jest na zmiennej "resOld".
//!
double Int1DTrap::Fast(util::Fun1D const & f, double a, double b, long int divNo, double resOld) const
{
long int i;
double res, h, x;

	assert(b > a);

	h = (b - a) / divNo;
	x = a + h;
	h *= 2.;
	res = 0.;
	divNo = divNo >> 1;
	for(i = 0; i < divNo; i++)
	{
		res += f.Get(x);
		x += h;
	}

	return 0.5 * (resOld + h * res);
}


