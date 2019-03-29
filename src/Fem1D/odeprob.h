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
	OdeProb() = default;
	OdeProb(Bndr left, Bndr right, double gamma, std::shared_ptr<util::Fun1D> const & g, std::shared_ptr<util::Fun1D> const & f);
	~OdeProb() = default;

	void Define(Bndr left, Bndr right, double gamma, std::shared_ptr<util::Fun1D> const & g, std::shared_ptr<util::Fun1D> const & f);

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
	std::shared_ptr<ClpMtxBand> m_s;

    //! Right hand side vektor for equation Sy = b
	std::shared_ptr<Vec> m_b;

    //! Coefficient vector y.
	std::shared_ptr<Vec> m_y;
};

#endif

