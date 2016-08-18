//#include "stdafx.h"
//#include "exchpbe.h"
//
//// March 8th, 2014	Added by dc1394
//#include <utility>
//
//
//const double ExchPbe::m_kpTable[2] = {0.804, 1.245};
//// const double ExchPbe::m_beta = 0.066725;
//// const double ExchPbe::m_mu = m_beta * util::HelpFun::M_PI * util::HelpFun::M_PI / 3.;
//const double ExchPbe::m_mu = 0.2195149727645171;
//const ExchS ExchPbe::m_slater;
//const double ExchPbe::m_coef = 1. / (2 * pow(3. * util::HelpFun::M_PI * util::HelpFun::M_PI, 1./3.));
//
//
////!
////! Konstruktor
////! rev - jezeli "true", to "revPBE"
////!
////ExchPbe::ExchPbe(bool rev)
////{
////	m_kp = m_kpTable[rev ? 1 : 0];
////}
//// March 7th, 2014	Modified by dc1394
//ExchPbe::ExchPbe(bool rev, std::function<double(double)> rhoTilde,
//				 std::function<double(double)> rhoTildeDeriv,
//				 std::function<double(double)> rhoTildeLapl)
//	:	Xc(std::move(rhoTilde), std::move(rhoTildeDeriv), std::move(rhoTildeLapl))
//{
//	if (rev)
//		xc_func_init(pxcfunc_.get(), XC_GGA_X_RPBE, XC_UNPOLARIZED);
//	else
//		xc_func_init(pxcfunc_.get(), XC_GGA_X_PBE, XC_UNPOLARIZED);
//}
//
////!
////! Destruktor
////!
//ExchPbe::~ExchPbe(void)
//{
//}
//
////!
////! Zwraca gestosc energii wymiany przypdajaca na jeden elektron
////! rhoa - gestosc elektronowa alpha
////! rhob - gestosc elektronowa beta
////! gaa  - gamma_{\alpha \alpha}
////! gab  - gamma_{\alpha \beta}, NIE UZYWANE.
////! gbb  - gamma_{\beta \beta}
////!
//double ExchPbe::E(double rhoa, double rhob, double gaa, double gab, double gbb) const
//{
//	return 0.5 * (Ehelp(rhoa, gaa) + Ehelp(rhob, gbb));
//}
//
//double ExchPbe::Ehelp(double rho, double g) const
//{
//	return m_slater.E(rho, rho, 0., 0., 0.) * Fx(S(2. * rho, 4. * g));
//}
//
//double ExchPbe::Fx(double s) const
//{
//	return 1. + m_kp - m_kp / (1. + m_mu * s * s / m_kp);
//}
//
////!
////! \frac{ \partial F_x(s) } { \partial s}
////!
//double ExchPbe::FxDer(double s) const
//{
//const double v1 = 2. * m_kp * m_kp * m_mu * s;
//const double v2 = m_kp + m_mu * s * s;
//
//	return v1 / (v2 * v2);
//}
//
////!
////! grho = $|\nabla \rho|^2$
////!
//double ExchPbe::S(double rho, double grho) const
//{
//	return sqrt(grho) * m_coef * pow(rho, -4./3.);
//}
//
////!
////! \frac{ \partial s }{ \partial \rho }
////!
//double ExchPbe::SRho(double rho, double grho) const
//{
//	return -(4./3.) * sqrt(grho) * m_coef * pow(rho, -7./3.);
//}
//
//
////!
////! \frac{ \partial s }{ \partial |\rho| }
////!
//double ExchPbe::SGradRho(double rho, double grho) const
//{
//	return m_coef * pow(rho, -4./3.);
//}
//
//double ExchPbe::Kf(double rho) const
//{
//	return pow(3 * util::HelpFun::M_PI * util::HelpFun::M_PI * rho, 1./3.);
//}
//
// 
////!
////! Zwraca pochodna funkcjonalu po $\rho_{\alpha}$.
////! rhoa - gestosc elektronowa alpha
////! rhob - gestosc elektronowa beta; NIE UZYWANE.
////! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
////! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
////! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
////!
//double ExchPbe::Vrhoa(double rhoa, double rhob, double gaa, double gab, double gbb) const
//{
//	return  pow(rhoa, 1./3.);
//}
//
////!
////! Zwraca pochodna funkcjonalu po $\rho_{\beta}$.
////! rhoa - gestosc elektronowa alpha; NIE UZYWANE.
////! rhob - gestosc elektronowa beta
////! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
////! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
////! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
////!
//double ExchPbe::Vrhob(double rhoa, double rhob, double gaa, double gab, double gbb) const
//{
//
//	return pow(rhob, 1./3.);
//}
//
//
////!
////! Zwraca pochodna funkcjonalu po $\gamma_{\alpha \alpha}$.
////! rhoa - gestosc elektronowa alpha; NIE UZYWANE.
////! rhob - gestosc elektronowa beta; NIE UZYWANE.
////! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
////! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
////! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
////!
//double ExchPbe::Vgaa(double rhoa, double rhob, double gaa, double gab, double gbb) const
//{
//const double v1 = m_slater.E(rhoa, rhoa, 0., 0., 0.);
//const double s = S(2. * rhoa, 4. * gaa);
//const double v2 = FxDer(s);
//const double v3 = SGradRho(2. * rhoa, 4. * gaa);
//
//	return 0.5 * v1 * v2 * v3;
//}
//
////!
////! Zwraca pochodna funkcjonalu po $\gamma_{\beta \beta}$.
////! rhoa - gestosc elektronowa alpha; NIE UZYWANE.
////! rhob - gestosc elektronowa beta; NIE UZYWANE.
////! gaa  - gamma_{\alpha \alpha}; NIE UZYWANE.
////! gab  - gamma_{\alpha \beta}; NIE UZYWANE.
////! gbb  - gamma_{\beta \beta}; NIE UZYWANE.
////!
//double ExchPbe::Vgbb(double rhoa, double rhob, double gaa, double gab, double gbb) const
//{
//	return 0.;
//}
//
////!
////! Zwraca pochodna funkcjonalu po $\gamma_{\alpha \beta}$. To jest zawsze rowne zero.
////! rhoa - gestosc elektronowa alpha; NIE UZYWANE.
////! rhob - gestosc elektronowa beta; NIE UZYWANE.
////! gaa - gamma_{\alpha \alpha}; NIE UZYWANE.
////! gab - gamma_{\alpha \beta}; NIE UZYWANE.
////! gbb - gamma_{\beta \beta}; NIE UZYWANE.
////!
//double ExchPbe::Vgab(double rhoa, double rhob, double gaa, double gab, double gbb) const
//{
//	return 0;
//}
//
//
//// March 8th, 2014	Added by dc1394
//double ExchPbe::V(double r) const
//{
//	return my_xc_gga_vxc(r);
//}
//
//// March 8th, 2014	Added by dc1394
//double ExchPbe::E(double r) const
//{
//	return my_xc_gga_exc(r);
//}
