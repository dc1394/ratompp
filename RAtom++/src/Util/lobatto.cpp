#include "stdafx.h"
#include "lobatto.h"
#include "helpfun.h"

//
// Constructor
//
Lobatto::Lobatto(void)
{
	SetGauss();
	CalcMemiMesi();
}

//
// Destructor
//
Lobatto::~Lobatto(void)
{
}

//
// Calculates MESI and MEMI matrices
//
void Lobatto::CalcMemiMesi()
{
size_t i, j;

	// Calculate the master element stiffness integrals MESI
	for(i = 0; i < MAXP + 1; i++)
	{
		for(j = i; j < MAXP + 1; j++)
		{
			m_mesi[i][j] = CalcMesi(i, j);
			m_mesi[j][i] = m_mesi[i][j];

			m_memi[i][j] = CalcMemi(i, j);
			m_memi[j][i] = m_memi[i][j];

                        // In order to check analytical and numerical evaluation of integrals
			//if(fabs(CalcMemi2(i, j) - CalcMemi(i, j)) > 0.0001)
			//	printf("i=%ld  j=%ld\n", i, j);
		}
	}
}

//
// Returns integral $\int_{-1}^1 \psi_i'(x) \psi_j'(x) dx$.
// Integrals evaluated analiticaly.
//
double Lobatto::CalcMesi(size_t i, size_t j) const
{
	// Upper (triangular) symmetric part only
	assert(j >= i);

	if((i == 0 && j == 0) || (i == 1 && j == 1))
		return 0.5;

	if(i == 0 && j == 1)
		return -0.5;

	if(i == j)
		return 1.;

	return 0;
}


//
// Returns integral $\int_{-1}^1 \psi_i(x) \psi_j(x) dx$.
// Integral is evaluated based on the Gaussian quadratures
//
double Lobatto::CalcMemi2(size_t i, size_t j) const
{
double w, x, v = 0;

	assert(m_xGauss.size() == m_wGauss.size());
	assert(m_xGauss.size() > 0);

	for(size_t n = 0; n < m_xGauss.size(); n++)
	{
		x = m_xGauss[n];
		w = m_wGauss[n];
		v += w * Basis(i, x) * Basis(j, x);
	}
	return v;
}


//
// Returns integral $\int_{-1}^1 \psi_i(x) \psi_j(x) dx$.
// Integrals evaluated analiticaly.
//
double Lobatto::CalcMemi(size_t i, size_t j) const
{
const double c0[4]  = {2./3  , 1/3.       , -1/sqrt(6.)         , 1./(3*sqrt(10.))};
const double c1[3]  = {2./3  , -1/sqrt(6.), -1/(3*sqrt(10.))    };
const double c2[3]  = {2./5  , 0          , -1/(5*sqrt(21.))    };
const double c3[3]  = {2./21 , 0          , -1/(21*sqrt(5.))    };
const double c4[3]  = {2./45 , 0          , -1/(9*sqrt(77.))    };
const double c5[3]  = {2./77 , 0          , -1/(33*sqrt(13.))   };
const double c6[3]  = {2./117, 0          , -1/(13*sqrt(165.))  };
const double c7[3]  = {2./165, 0          , -1/(15*sqrt(221.))  };
const double c8[3]  = {2./221, 0          , -1/(17*sqrt(285.))  };
const double c9  = 2./285;
const double c10 = 2./357;
double v;

	// Upper (triangular) symmetric part only
	assert(j >= i);

	const size_t k = j - i;

	v = 0;
	switch(i)
	{
		case 0: if(k < 4) v = c0[k]; break;
		case 1: if(k < 3) v = c1[k]; break;
		case 2: if(k < 3) v = c2[k]; break;
		case 3: if(k < 3) v = c3[k]; break;
		case 4: if(k < 3) v = c4[k]; break;
		case 5: if(k < 3) v = c5[k]; break;
		case 6: if(k < 3) v = c6[k]; break;
		case 7: if(k < 3) v = c7[k]; break;
		case 8: if(k < 3) v = c8[k]; break;
		case 9: if(j == 9) v = c9; break;
		case 10: if(j == 10) v = c10; break;
		default:
			assert(0);
		break;
	}
	
	return v;
}

//
// Returns pre-computed MEMI
//
double Lobatto::Memi(size_t i, size_t j) const
{
	assert(i < MAXP + 1);
	assert(j < MAXP + 1);

	return m_memi[i][j];
}

//
// Returns pre-computed MESI
//
double Lobatto::Mesi(size_t i, size_t j) const
{
	assert(i < MAXP + 1);
	assert(j < MAXP + 1);

	return m_mesi[i][j];
}


//
// Returns the value of basis function $\psi_i(s)$
// First eleven (for i = 0, 1, ..., 10) Lobatto hierarchic shape functions.
//
double Lobatto::Basis(size_t i, double s) const
{
static const double c[11] = {0.5, 0.5, sqrt(3./2.)/2., sqrt(5./2)/2., sqrt(7./2.)/8., sqrt(9./2.)/8.,
	sqrt(11./2.)/16., sqrt(13./2.)/16., sqrt(15./2.)/128., sqrt(17./2.)/128., sqrt(19./2.)/256.};
const double s2 = s * s;
double v = 0;

	assert(s >= -1 && s <= 1);

	switch(i)
	{
		case 0: v = 1 - s; break;

		case 1:	v = 1 + s; break;

		case 2: v = s2 - 1; break;

		case 3: v = (s2 - 1) * s; break;

		case 4: v = (s2 - 1) * (5 * s2 - 1); break;

		case 5: v = (s2 - 1) * (7 * s2 - 3) * s; break;

		case 6: v = (s2 - 1) * (21 * s2 * s2 - 14 * s2 + 1); break;

		case 7: v = (s2 - 1) * (33 * s2 * s2 - 30 * s2 + 5) * s; break;

		case 8: v = (s2 - 1) * (429 * s2 * s2 * s2 - 495 * s2 * s2 + 135 * s2 - 5); break;

		case 9: v = (s2 - 1) * (715 * s2 * s2 * s2 - 1001 * s2 * s2 + 385 * s2 - 35) * s; break;

		case 10: 
		{
			const double s4 = s2 * s2, s6 = s4 * s2, s8 = s6 * s2;
			v = (s2 - 1) * (2431 * s8 - 4004 * s6 + 2002 * s4 - 308 * s2 + 7); 
		}
		break;

		default:
			assert(0);
		break;
	}

	return c[i] * v;
}

//
// Returns the value of first derivative of basis function $\psi_i''(s)$
// First eleven (for i = 0, 1, ..., 10) second derivative of Lobatto hierarchic shape functions.
//
double Lobatto::BasisD1(size_t i, double s) const
{
static const double c[11] = {-0.5, 0.5, sqrt(3./2.), sqrt(5./2)/2., sqrt(7./2.)/2., sqrt(9./2.)/8.,
	sqrt(11./2.)/8., sqrt(13./2.)/16., sqrt(15./2.)/16., sqrt(17./2.)/128., sqrt(19./2.)/256.};
const double s2 = s * s;
double v = 0;

	assert(s >= -1 && s <= 1);

	switch(i)
	{
		case 0: v = 1; break;

		case 1:	v = 1; break;

		case 2: v = s; break;

		case 3: v = 3 * s2 - 1; break;

		case 4: v = (5 * s2 - 3) * s; break;

		case 5: v = 35 * s2 * s2 - 30 * s2 + 3; break;

		case 6: v = (63 * s2 * s2 - 70 * s2 + 15) * s; break;

		case 7: v = 231 * s2 * s2 * s2 - 315 * s2 * s2 + 105 * s2 - 5; break;

		case 8: v = (429 * s2 * s2 * s2 - 693 * s2 * s2 + 315 * s2 - 35) * s; break;

		case 9: v = 6435 * s2 * s2 * s2 * s2 - 12012 * s2 * s2 * s2 + 6930 * s2 * s2 - 1260 * s2 + 35; break;

		case 10: 
		{
			const double s4 = s2 * s2, s6 = s4 * s2, s8 = s6 * s2;
			v = 12155 * s8 - 25740 * s6 + 18018 * s4 - 4620 * s2 + 315; 
		}
		break;

		default:
			assert(0);
		break;
	}

	return c[i] * v;
}


//
// Returns the value of second derivative of basis function $\psi_i''(s)$
// First eleven (for i = 0, 1, ..., 10) second derivative of Lobatto hierarchic shape functions.
//
double Lobatto::BasisD2(size_t i, double s) const
{
static const double c[11] = {0, 0, sqrt(3./2.), 3*sqrt(5./2), (3./2.)*sqrt(7./2.), (5./2.)*sqrt(9./2.),
	(15./8.)*sqrt(11./2.), (21./8.)*sqrt(13./2.), (7./16.)*sqrt(15./2.), (9./16.)*sqrt(17./2.), (45./128.)*sqrt(19./2.)};
const double s2 = s * s;
double v = 0;

	assert(s >= -1 && s <= 1);

	switch(i)
	{
		case 0: v = 0; break;

		case 1:	v = 0; break;

		case 2: v = 1; break;

		case 3: v = s; break;

		case 4: v = 5 * s2 - 1; break;

		case 5: v = (7 * s2 - 3) * s; break;

		case 6: v = 21 * s2 * s2 - 14 * s2 + 1; break;

		case 7: v = (33 * s2 * s2 - 30 * s2 + 5) * s; break;

		case 8: v = 429 * s2 * s2 * s2 - 495 * s2 * s2 + 135 * s2 - 5; break;

		case 9: v = (715 * s2 * s2 * s2 - 1001 * s2 * s2 + 385 * s2 - 35) * s; break;

		case 10: 
		{
			const double s4 = s2 * s2, s6 = s4 * s2, s8 = s6 * s2;
			v = 2431 * s8 - 4004 * s6 + 2002 * s4 - 308 * s2 + 7; 
		}
		break;

		default:
			assert(0);
		break;
	}

	return c[i] * v;
}



//
// Calculates the gaussian node coorinates and weights
//
void Lobatto::SetGauss()
{
const size_t n = 3*MAXP; // Is it OK?

	m_xGauss.resize(n);
	m_wGauss.resize(n);
	::gauleg(-1, 1, m_xGauss, m_wGauss, n);
}

