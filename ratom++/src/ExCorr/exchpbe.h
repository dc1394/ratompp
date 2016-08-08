//#ifndef __RATOM_EXCHPBE_H__
//#define __RATOM_EXCHPBE_H__
//
//
///** \brief PBE approximation.
//*
//* \author Zbigniew Romanowski [ROMZ@wp.pl]
//*
//*/
//
//// March 7th, 2014	Modified by dc1394
//
//#include "exchs.h"
//
//
//// March 8th, 2014	Added by dc1394
//#include "xc.h"
//
//class ExchPbe : public Xc
//{
//public:
//	// March 8th, 2014	Modified by dc1394
//	//ExchPbe(bool rev);
//	ExchPbe(bool rev, std::function<double(double)> rhoTilde,
//			std::function<double(double)> derivRhoTilde,
//			std::function<double(double)> laplRhoTilde);
//	virtual ~ExchPbe(void);
//
//	virtual double E(double rhoa, double rhob, double gaa, double gab, double gbb) const;
//
//	virtual double Vrhoa(double rhoa, double rhob, double gaa, double gab, double gbb) const;
//	virtual double Vrhob(double rhoa, double rhob, double gaa, double gab, double gbb) const;
//
//	virtual double Vgaa(double rhoa, double rhob, double gaa, double gab, double gbb) const;
//	virtual double Vgbb(double rhoa, double rhob, double gaa, double gab, double gbb) const;
//	virtual double Vgab(double rhoa, double rhob, double gaa, double gab, double gbb) const;
//
//	// March 8th, 2014	Added by dc1394
//	virtual double V(double r) const;
//	virtual double E(double r) const;
//
//	virtual const char* Name() const
//	{
//		return "(rev) PBE Exchange";
//	}
//
//private:
//	double Ehelp(double rho, double g) const;
//
//	double Kf(double rho) const;
//	// double Uni(double rho) const;
//
//	double S(double rho, double grho) const;
//	double SRho(double rho, double grho) const;
//	double SGradRho(double rho, double grho) const;
//
//	double Fx(double s) const;
//	double FxDer(double s) const;
//
//private:
//	//! Dwie wartosci $\kappa$
//	static const double m_kpTable[2];
//
//	//! $\mu$
//	static const double m_mu;
//
//	//! $\beta$
//	// static const double m_beta;
//
//	static const ExchS m_slater;
//
//	//! $\frac{1}{ 2 (3 \pi)^{1/3} }$
//	static const double m_coef;
//
//	//! $\kappa$
//	double m_kp;
//};
//
//#endif
//
