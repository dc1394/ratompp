#include "stdafx.h"
#include "corrlyp.h"



const double CorrLyp::m_a = 0.04918;
const double CorrLyp::m_b = 0.132;
const double CorrLyp::m_c = 0.2533;
const double CorrLyp::m_d = 0.349;
const double CorrLyp::m_cf = (3. / 10) * pow(3 * util::HelpFun::M_PI * util::HelpFun::M_PI, 2. / 3.);

CorrLyp::CorrLyp(void)
{
}

CorrLyp::~CorrLyp(void)
{
}

//!
//! Zwraca gestosc energii wymiany przypdajaca na jeden elektron
//! rhoa - gestosc elektronowa alpha
//! rhob - gestosc elektronowa beta
//! gaa  - gamma_{\alpha \alpha}
//! gab  - gamma_{\alpha \beta};
//! gbb  - gamma_{\beta \beta}
//!
double CorrLyp::E(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
const double rho = rhoa + rhob;
const double g = gaa + 2 * gab + gbb;
const double rhos = pow(rho, -1./3.); // $\rho^*$
const double w = (rhoa * gaa / rho + rhob * gbb / rho) / 9.;
const double mul = rhoa * rhob; 

const double v1 = -4 * m_a / (1 + m_d * rhos) * mul / rho;
const double v2 = -m_a * m_b * OmegaT(rhos);
const double v3 = pow(2., 11./3.) * m_cf * mul * (pow(rhoa, 8./3.) + pow(rhob, 8./3.));
const double v4 = DeltaT(rhos) * mul * (1./18. * (gaa + gbb) - 7./18. * g - w);
const double v5 = mul * (47./18. * g - 5./2. * (gaa + gbb) + 11. * w);
const double v6 = 2./3. * rho * rho * (gaa + gbb - g) - rhoa * rhoa * gbb - rhob * rhob * gaa;
 

	return v1 + v2 * (v3 + v4 + v5 + v6);
}

//!
//! Funkcja pomocnicza
//! f = -(1/9) * a * b * \omega(rho)
//!
double CorrLyp::A(double f, double delta, double rhoa, double rhob) const
{
const double rho = rhoa + rhob;
const double v1 = 1. - 3. * delta - (delta - 11.) * rhoa / rho;

	return f * (rhoa * rhob * v1 - 9. * rhob * rhob);
}

//!
//! Funkcja pomocnicza
//! f = -(1/9) * a * b * \omega(rho)
//!
double CorrLyp::B(double f, double delta, double rhoa, double rhob) const
{
const double rho = rhoa + rhob;
const double v1 = 47. - 7. * delta;

	return f * (rhoa * rhob * v1 - 12. * rho * rho);

}

//!
//! Funkcja pomocnicza
//! f = -(1/9) * a * b * \omega(rho)
//!
double CorrLyp::C(double f, double delta, double rhoa, double rhob) const
{
	// Wystarczy zmienic argumenty wywolania w funcki A()
	return A(f, delta, rhob, rhoa);
}

/*
//!
//! Zwraca wartosc funkcji $\omega(rho)$
//!
double CorrLyp::Omega(double rho) const
{
const double rho13 = pow(rho, -1. / 3.);

	return exp(-m_c * rho13) * pow(rho, -11. / 3.) / (1. + m_d * rho13);
}
*/

//!
//! Zwraca wartosc funkcji $\omega(rho)$ przyjmujac, ze x = rho^{-1/3}
//!
double CorrLyp::OmegaT(double x) const
{
// Szymkie obliczanie pow(x, 11);
double y = x * x; // x^2

	y = y * y; // x^4;
	y = y * y; // x^8
	y = y * x * x * x; // x^11

	// return exp(-m_c * x) * pow(x, 11.) / (1. + m_d * x);
	return exp(-m_c * x) * y / (1. + m_d * x);
}

//!
//! Zwraca wartosc funkcji $\Tilde{\omega}'(x)$
//!
double CorrLyp::OmegaDerT(double x) const
{
const double v1 = -(1./3.) * x * x * x * x * OmegaT(x);
const double v2 = 11. / x - m_c - m_d / (1. + m_d * x);

	return v1 * v2;
}


/*
//!
//! Zwraca wartosc funkcji $\delta(rho)$
//!
double CorrLyp::Delta(double rho) const
{
const double rho13 = pow(rho, -1. / 3.);

	return m_c * rho13 + m_d * rho13 / (1 + m_d * rho13);
}
*/

//!
//! Zwraca wartosc funkcji $\delta(rho)$ przyjmujac, ze x = rho^{-1/3}
//!
double CorrLyp::DeltaT(double x) const
{
	return m_c * x + m_d * x / (1. + m_d * x);
}

//!
//! Zwraca wartosc funkcji $\Tilde{\delta}'(x)$
//!
double CorrLyp::DeltaDerT(double x) const
{
const double v0 = 1. + m_d * x;
const double v1 = m_d * m_d * x * x * x * x * x;
const double v2 = x * x * x * DeltaT(x);

	return (v1 / (v0 * v0) - v2) / 3.;
}


//!
//! Zwraca pochodna funkcjonalu po $\gamma_{\alpha \alpha}$.
//! rhoa - gestosc elektronowa alpha
//! rhob - gestosc elektronowa beta
//! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
//! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
//!
double CorrLyp::Vgaa(double rhoa, double rhob, double /* gaa */, double /* gab */, double /* gbb */) const
{
const double rho = rhoa + rhob;
const double x = pow(rho, -1. / 3.); // \rho^*
const double delta = DeltaT(x);
const double omega = OmegaT(x);
const double f = - (1./9.) * m_a * m_b * omega;

	return A(f, delta, rhoa, rhob);
}

//!
//! Zwraca pochodna funkcjonalu po $\gamma_{\beta \beta}$.
//! rhoa - gestosc elektronowa alpha
//! rhob - gestosc elektronowa beta
//! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
//! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
//!
double CorrLyp::Vgbb(double rhoa, double rhob, double /* gaa */, double /* gab */, double /* gbb */) const
{
const double rho = rhoa + rhob;
const double x = pow(rho, -1. / 3.); // \rho^*
const double delta = DeltaT(x);
const double omega = OmegaT(x);
const double f = - (1./9.) * m_a * m_b * omega;

	return C(f, delta, rhoa, rhob);
}

//!
//! Zwraca pochodna funkcjonalu po $\gamma_{\beta \beta}$.
//! rhoa - gestosc elektronowa alpha
//! rhob - gestosc elektronowa beta
//! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
//! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
//!
double CorrLyp::Vgab(double rhoa, double rhob, double /* gaa */, double /* gab */, double /* gbb */) const
{
const double rho = rhoa + rhob;
const double x = pow(rho, -1. / 3.); // \rho^*
const double delta = DeltaT(x);
const double omega = OmegaT(x);
const double f = - (1./9.) * m_a * m_b * omega;

	return B(f, delta, rhoa, rhob);
}



//!
//! Zwraca wartosc pochodnej czastkowej po funkcjonale gestosci korelacjosci.
//! Pochodna czastkwa jest po $\rho_{\alpha}$
//!
//! \frac{\partial \text{LYP}_c}{\partial \rho_{\rho_{\alpha}}}
//!
double CorrLyp::Vrhoa(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
const double rho = rhoa + rhob;
const double x = pow(rho, -1./3.); // \rho^*

const double omega = OmegaT(x);
const double omegaDer = OmegaDerT(x);
const double delta = DeltaT(x);
const double deltaDer = DeltaDerT(x);

const double y = 1. / x; // pow(rho, 1./3.)
const double rho53 = y * y * y * y * y; // pow(rho, 5./3.)
const double ra = -4. * m_a  * rhob;
const double rb = 3. * rhob * y + m_d * (rhoa + 3. * rhob);
const double rc = 3. * rho53 * (m_d + y) * (m_d + y);
const double r = ra * rb / rc;


const double rhoa83 = pow(rhoa, 8./3.);
const double rhob83 = pow(rhob, 8./3.);

const double wa = omegaDer * rhoa * rhob * (rhoa83 + rhob83);
const double wb = omega * rhob * ((11./3.) * rhoa83 + rhob83);
const double w = -pow(2., 11./3.) * m_cf * m_a * m_b * (wa + wb);

const double f = - (1./9.) * m_a * m_b * omega;
const double om = omegaDer / omega;

	return r + w + 
		gaa * Ader(f, om, delta, deltaDer, rhoa, rhob) + 
		gab * Bder(f, om, delta, deltaDer, rhoa, rhob) + 
		gbb * Cder(f, om, delta, deltaDer, rhoa, rhob);
}

//!
//! Zwraca wartosc pochodnej czastkowej po funkcjonale gestosci korelacjosci.
//! Pochodna czastkwa jest po $\rho_{\rho_{\beta}}$
//!
//! \frac{\partial \text{LYP}_c}{\partial \rho_{\beta}}
//!
double CorrLyp::Vrhob(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
	// Wystarczy zamienic \alpha na \beta
	return Vrhoa(rhob, rhoa, gbb, gab, gaa);
}


//!
//! Pochodna czastkowa $\frac{\partial A(\rho_{\alpha}, \rho_{\beta})}{\partial \rho_{\alpha}}$
//!
double CorrLyp::Ader(double f, double om, double delta, double deltaDer, double rhoa, double rhob) const
{
const double rho = rhoa + rhob;
const double v1 = 1. - 3. * delta;
const double v2 = -(delta - 11.) * rhoa / rho * (1. + rhob / rho);
const double v3 = -(3. + rhoa / rho) * deltaDer * rhoa;

	return om * A(f, delta, rhoa, rhob) + f * rhob * (v1 + v2 + v3);
}

//!
//! Pochodna czastkowa $\frac{\partial B(\rho_{\alpha}, \rho_{\beta})}{\partial \rho_{\alpha}}$
//!
double CorrLyp::Bder(double f, double om, double delta, double deltaDer, double rhoa, double rhob) const
{
const double rho = rhoa + rhob;
const double v1 = 47. - 7. * delta - 7. * rhoa * deltaDer; 

	return om * B(f, delta, rhoa, rhob) + f * (rhob * v1 - 24. * rho);
}

//!
//! Pochodna czastkowa $\frac{\partial C(\rho_{\alpha}, \rho_{\beta})}{\partial \rho_{\alpha}}$
//!
double CorrLyp::Cder(double f, double om, double delta, double deltaDer, double rhoa, double rhob) const
{
const double rho = rhoa + rhob;
const double v1 = 1. - 3. * delta;
const double v2 = -(delta - 11.) * rhob / rho * (1. - rhoa / rho);
const double v3 = -(3. + rhob / rho) * deltaDer * rhoa;
const double v4 = rhob * (v1 + v2 + v3) - 18. * rhoa;

	return om * C(f, delta, rhoa, rhob) + f * v4;

}

