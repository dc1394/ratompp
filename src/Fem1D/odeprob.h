#ifndef __RATOM_ODEPROB_H__
#define __RATOM_ODEPROB_H__


/** \brief Ordinary Differential equation, ODE, with form
*     -\gamma u''(x) + g(x) u(x) = f(x)
* Functions f(x), g(x) and constant \gamma are given.
* Function u(x) is searched.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/


#include "prob.h"
#include "eltinfo.h"

class OdeProb : public Prob
{
public:
	OdeProb();
	OdeProb(Bndr left, Bndr right, double gamma, const util::Fun1D* g, const util::Fun1D* f);
	~OdeProb(void);

	void Define(Bndr left, Bndr right, double gamma, const util::Fun1D* g, const util::Fun1D* f);

	void Solve();
	void SolveAdapt(double absMaxCoef);

	double GetSol(double x) const;
	void WriteSol(const char* path, size_t pointNo) const;

private:
	void Assemble();
	void Malloc(void);
	void MaxMinCoef(EltInfo& eltInfo) const;

private:
	//! Stifness matrix
	ClpMtxBand* m_s;

        //! Right hand side vektor for equation Sy = b
	Vec* m_b;

        //! Coefficient vector y.
	Vec* m_y;
};

#endif

