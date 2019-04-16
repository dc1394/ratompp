#ifndef __RATOM_EIGPROB_H__
#define __RATOM_EIGPROB_H__


/** \brief Differential eigenvalue problem, with form
*     -\gamma u''(x) + g(x) u(x) = \lambda f(x)
* Function g(x) and constant \gamma are given.
* Eigenfunctions u(x) and eigenvalues \lambda are searched.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/


#include "prob.h"
#include "eltinfo.h"

class EigProb : public Prob
{
public:
	EigProb() = default;
	~EigProb() override = default;

	void Define(double gamma, std::shared_ptr<util::Fun1D const> const & g);

	void Solve(size_t eigNo, double abstol);
	void SolveAdapt(size_t eigNo, double abstol, double absMaxCoef);

	double GetEigVal(size_t eig) const;
	double GetEigFun(size_t eig, double x) const;
	double GetEigFunD1(size_t eig, double x) const;
	double GetEigFunD2(size_t eig, double x) const;

	void WriteEigFun(const char* path, size_t eig, size_t pointNo) const;
	void WriteEigFunEx(const char* path, size_t ll,  size_t eig, size_t pointNo) const;

	void WriteEigCoef(const char* path, size_t eigNo) const;


private:
	void Malloc(void);
	void Assemble();
	void MaxMinCoef(std::vector<EltInfo>& eltInfo) const;


private:
	// Stifness matrix
	std::shared_ptr<ClpMtxBand> m_s;

	// Vector with calculated eigenvalues
	std::shared_ptr<Vec> m_w;

    // Matrix with calculated eignvectors
	std::shared_ptr<ClpMtx> m_z;

	// Overlap matrix
	std::shared_ptr<ClpMtxBand> m_o;
};

#endif

