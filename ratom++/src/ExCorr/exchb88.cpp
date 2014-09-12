#include "stdafx.h"
#include "exchb88.h"

// March 7th, 2014	Added by dc1394
#include <utility>

const double ExchB88::m_rhoMin = 1e-20;

//!
//! Konstruktor
//!
// Modified by dc1394 
//ExchB88::ExchB88(void)
ExchB88::ExchB88(std::function<double(double)> rhoTilde,
				 std::function<double(double)> rhoTildeDeriv,
				 std::function<double(double)> rhoTildeLapl)
	:	Xc(std::move(rhoTilde), std::move(rhoTildeDeriv), std::move(rhoTildeLapl))
{
	// March 7th, 2014	Added by dc1394
	xc_func_init(pxcfunc_.get(), XC_GGA_X_B88, XC_UNPOLARIZED);
}

//!
//! Destruktor
//!
ExchB88::~ExchB88(void)
{
}

//!
//! Zwraca gestosc energii wymiany przypdajaca na jeden elektron
//! rhoa - gestosc elektronowa alpha
//! rhob - gestosc elektronowa beta
//! gaa  - gamma_{\alpha \alpha}
//! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb  - gamma_{\beta \beta}
//!
double ExchB88::E(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
	// return (Ehelp(rhoa, gaa) + Ehelp(rhob, gbb)) / (rhoa + rhob);
	return Ehelp(rhoa, gaa) + Ehelp(rhob, gbb);
}

//!
//! Zwraca "kawalek" gestosc energii wymiany dla gestsci "rho" i gradiencie "g"
//!
double ExchB88::Ehelp(double rho, double g) const
{
	if(rho < m_rhoMin)
		return 0.;

const double rho43 = pow(rho, 4. / 3.);
const double x = sqrt(g) / rho43;

	return rho43 * G(x);
}

//!
//! Zwraca wartosc funkcji pomocniczej 
//!
double ExchB88::G(double x) const
{
const double q = 3. / 2. * pow(3. / (4 * M_PI), 1. / 3.);
const double b = 0.0042; // Stala z pracy B88

	return -q - b * x * x / (1 + 6 * b * x * asinh(x));
}

//!
//! Zwraca wartosc pochodnej funkcji pomocniczej $G'(x)$
//!
double ExchB88::Gprim(double x) const
{
const double b = 0.0042; // Stala z pracy B88
const double q = asinh(x); // Aby raz liczyc $asinh(x)$
double v1, v2;

	v1 = x / sqrt(x * x + 1) - q;
	v2 = 1 + 6 * b * x * q;

	return (6 * b * b * x * x * v1 - 2 * b * x) / (v2 * v2);
}


//!
//! Zwraca pochodna funkcjonalu po $\rho_{\alpha}$.
//! rhoa - gestosc elektronowa alpha
//! rhob - gestosc elektronowa beta; NIE UZYWANE.
//! gaa  - gamma_{\alpha \alpha}
//! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
//!
double ExchB88::Vrhoa(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
	return Vhelp(rhoa, gaa);
}

//!
//! Zwraca pochodna funkcjonalu po $\rho_{\beta}$.
//! rhoa - gestosc elektronowa alpha; NIE UZYWANE.
//! rhob - gestosc elektronowa beta
//! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
//! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb  - gamma_{\beta \beta}
//!
double ExchB88::Vrhob(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
	return Vhelp(rhob, gbb);
}

//!
//! Funkcja pomocnicza do obliczania pochodnej po 
//!
double ExchB88::Vhelp(double rho, double g) const
{
	if(rho < m_rhoMin)
		return 0.;

const double rho13 = pow(rho, 1. / 3.);
const double x = sqrt(g) / (rho * rho13);

	return (4. / 3.) * rho13 * (G(x) - x * Gprim(x));
}

//!
//! Zwraca pochodna funkcjonalu po $\gamma_{\alpha \alpha}$.
//! rhoa - gestosc elektronowa alpha
//! rhob - gestosc elektronowa beta; NIE UZYWANE.
//! gaa  - gamma_{\alpha \alpha}
//! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
//!
double ExchB88::Vgaa(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
	if(rhoa < m_rhoMin || gaa < m_rhoMin)
		return 0.;

const double q = sqrt(gaa);
const double xa = q / pow(rhoa, 4. / 3.);

	return Gprim(xa) / (2 * q);
}

//!
//! Zwraca pochodna funkcjonalu po $\gamma_{\beta \beta}$.
//! rhoa - gestosc elektronowa alpha; NIE UZYWANE.
//! rhob - gestosc elektronowa beta
//! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
//! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb  - gamma_{\beta \beta}
//!
double ExchB88::Vgbb(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
	if(rhob < m_rhoMin || gbb < m_rhoMin)
		return 0.;

const double q = sqrt(gbb);
const double xb = q / pow(rhob, 4. / 3.);

	return Gprim(xb) / (2 * q);
}

//!
//! Zwraca pochodna funkcjonalu po $\gamma_{\alpha \beta}$. To jest zawsze rowne zero.
//! rhoa - gestosc elektronowa alpha; NIE UZYWANE.
//! rhob - gestosc elektronowa beta; NIE UZYWANE.
//! gaa - gamma_{\alpha \alpha}; NIE UZYWANE.
//! gab - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb - gamma_{\beta \beta}; NIE UZYWANE.
//!
double ExchB88::Vgab(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
	return 0;
}

// March 8th, 2014	Added by dc1394
double ExchB88::V(double r) const
{
	return my_xc_gga_vxc(r);
}

// March 8th, 2014	Added by dc1394
double ExchB88::E(double r) const
{
	return my_xc_gga_exc(r);
}
