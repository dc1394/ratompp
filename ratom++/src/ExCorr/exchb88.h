//#ifndef __RATOM_EXCHB88_H__
//#define __RATOM_EXCHB88_H__
//
//
///** \brief Becke 1988 approximation.
//*
//* \author Zbigniew Romanowski [ROMZ@wp.pl]
//*
//*/
//
//// March 7th, 2014 Modified by dc1394
//
//#include "xc.h"
//
//class ExchB88 : public Xc
//{
//public:
//	// March 7th, 2014 Modified by dc1394
//	//ExchB88(void);
//	ExchB88(std::function<double(double)> rhoTilde,
//		    std::function<double(double)> rhoTildeDeriv,
//			std::function<double(double)> rhoTildeLapl);
//
//	virtual ~ExchB88(void);
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
//		return pxcfunc_->info->name;
//	}
//
//private:
//	double Ehelp(double rho, double g) const;
//	double Vhelp(double rho, double g) const;
//
//	double G(double x) const;
//	double Gprim(double x) const;
//
//	static double asinh(double x)
//	{
//		return log(x + sqrt(x * x + 1.));
//	}
//
//
//protected:
//	static const double m_rhoMin;
//};
//
//
//#endif
//
