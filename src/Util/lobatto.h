#ifndef __RATOM_LOBATTO_H__
#define __RATOM_LOBATTO_H__


/** \brief Lobatto functions and integral evaluation
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/


//! Maximal element degree
#define MAXP 10

#include "helpfun.h"

class Lobatto
{
public:
	Lobatto(void);
	~Lobatto(void);

	double Basis(size_t i, double s) const;
	double BasisD1(size_t i, double s) const;
	double BasisD2(size_t i, double s) const;

protected:
	double Memi(size_t i, size_t j) const;
	double Mesi(size_t i, size_t j) const;

private:
	void CalcMemiMesi();
	double CalcMesi(size_t i, size_t j) const;
	double CalcMemi(size_t i, size_t j) const;
	double CalcMemi2(size_t i, size_t j) const;	

	
	void SetGauss();
	

protected:
	//! Coordinates of Gauss quadrature
	Vec m_xGauss;

	//! Weights of Gauss quadrature
	Vec m_wGauss;

private:
	//! Master elements stiffness integrals, MESI
	double m_mesi[MAXP + 1][MAXP + 1];

	//! Master elements mass integrals, MEMI
	double m_memi[MAXP + 1][MAXP + 1];



};

#endif

