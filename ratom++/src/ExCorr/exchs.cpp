#include "stdafx.h"
#include "exchs.h"



const double ExchS::m_alpha = 2./3;
const double ExchS::m_c = pow(3. / (4. * util::HelpFun::M_PI), 1./3.);

//
// Construktor
//
ExchS::ExchS(void)
{
}

//
// Destructor
//
ExchS::~ExchS(void)
{
}

//!
//! Zwraca gestosc energii wymiany przypdajaca na jeden elektron
//! rhoa - gestosc elektronowa alpha
//! rhob - gestosc elektronowa beta
//! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
//! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
//!
double ExchS::E(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
const double f = -(9./4.) * m_alpha * m_c;

	// return f * (pow(rhoa, 4./3.) + pow(rhob, 4./3.)) / (rhoa + rhob);
	return f * (pow(rhoa, 4./3.) + pow(rhob, 4./3.));
}



//!
//! Zwraca pochodna funkcjonalu po $\rho_{\alpha}$.
//! rhoa - gestosc elektronowa alpha
//! rhob - gestosc elektronowa beta; NIE UZYWANE.
//! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
//! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
//!
double ExchS::Vrhoa(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
const double f = -3. * m_alpha * m_c;

	return f * pow(rhoa, 1./3.);
}

//!
//! Zwraca pochodna funkcjonalu po $\rho_{\beta}$.
//! rhoa - gestosc elektronowa alpha; NIE UZYWANE.
//! rhob - gestosc elektronowa beta
//! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
//! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
//!
double ExchS::Vrhob(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
const double f = -3. * m_alpha * m_c;

	return f * pow(rhob, 1./3.);
}


//!
//! Zwraca pochodna funkcjonalu po $\gamma_{\alpha \alpha}$.
//! rhoa - gestosc elektronowa alpha; NIE UZYWANE.
//! rhob - gestosc elektronowa beta; NIE UZYWANE.
//! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
//! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
//!
double ExchS::Vgaa(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
	return 0.;
}

//!
//! Zwraca pochodna funkcjonalu po $\gamma_{\beta \beta}$.
//! rhoa - gestosc elektronowa alpha; NIE UZYWANE.
//! rhob - gestosc elektronowa beta; NIE UZYWANE.
//! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
//! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
//!
double ExchS::Vgbb(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
	return 0.;
}

//!
//! Zwraca pochodna funkcjonalu po $\gamma_{\alpha \beta}$. To jest zawsze rowne zero.
//! rhoa - gestosc elektronowa alpha; NIE UZYWANE.
//! rhob - gestosc elektronowa beta; NIE UZYWANE.
//! gaa - gamma_{\alpha \alpha}; NIE UZYWANE.
//! gab - gamma_{\alpha \beta}; NIE UZYWANE.
//! gbb - gamma_{\beta \beta}; NIE UZYWANE.
//!
double ExchS::Vgab(double rhoa, double rhob, double gaa, double gab, double gbb) const
{
	return 0;
}



