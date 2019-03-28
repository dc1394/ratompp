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
	EigProb();
	EigProb(double gamma, const util::Fun1D* g);
	~EigProb(void);

	void Define(double gamma, const util::Fun1D* g);

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
	ClpMtxBand* m_s;

	// Vector with calculated eigenvalues
	Vec* m_w;

        // Matrix with calculated eignvectors
	ClpMtx* m_z;

	// Overlap matrix
	ClpMtxBand* m_o;
};

#endif

