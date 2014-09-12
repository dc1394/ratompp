#ifndef __RATOM_CORRLYP_H__
#define __RATOM_CORRLYP_H__


/** \brief LYP approximation.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/



#include "xc.h"


class CorrLyp //: public Xc
{
public:
	CorrLyp(void);
	virtual ~CorrLyp(void);

	virtual double E(double rhoa, double rhob, double gaa, double gab, double gbb) const;

	virtual double Vrhoa(double rhoa, double rhob, double gaa, double gab, double gbb) const;
	virtual double Vrhob(double rhoa, double rhob, double gaa, double gab, double gbb) const;

	virtual double Vgaa(double rhoa, double rhob, double gaa, double gab, double gbb) const;
	virtual double Vgbb(double rhoa, double rhob, double gaa, double gab, double gbb) const;
	virtual double Vgab(double rhoa, double rhob, double gaa, double gab, double gbb) const;


private:
	double A(double f, double delta, double rhoa, double rhob) const;
	double B(double f, double delta, double rhoa, double rhob) const;
	double C(double f, double delta, double rhoa, double rhob) const;

	// double Omega(double rho) const;
	double OmegaT(double rho) const;
	double OmegaDerT(double rho) const;

	// double Delta(double rho) const;
	double DeltaT(double rho) const;
	double DeltaDerT(double rho) const;

	double Ader(double omega, double omegaDer, double delta, double deltaDer, double rhoa, double rhob) const;
	double Bder(double omega, double omegaDer, double delta, double deltaDer, double rhoa, double rhob) const;
	double Cder(double omega, double omegaDer, double delta, double deltaDer, double rhoa, double rhob) const;

private:
	// stale
	static const double m_a;
	static const double m_b;
	static const double m_c;
	static const double m_d;
	static const double m_cf;

};


#endif

