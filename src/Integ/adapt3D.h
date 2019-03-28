#ifndef __RATOM_ADAPT3D_H__
#define __RATOM_ADAPT3D_H__


/**
*	\brief Oblicza calke w trzech wymiarach po prostopadloscianie. 
*		Prostopadloscian ma sciany rownolegle do plaszczyzn ukladu wspolrzednych. 
*
*	\author Zbigniew Romanowski, ROMZ
*
*	\version 07-Lip-2007 [romz]
*	\version 12-Sty-2009 [romz]
*
*/


#include "eltset.h"
#include "brick.h"
#include "pool.h"


class Adapt3D
{
public:
	Adapt3D(double epsAbs, size_t initSize);
	~Adapt3D(void);

	size_t BrickNo() const;
	size_t FunCall() const;

	double GetEpsAbs(void) const { return m_epsAbs; }
	void SetEpsAbs(double epsAbs) { m_epsAbs = epsAbs; }

	double Calc(const Fun3D& f, double x0, double x1, double y0, double y1, double z0, double z1);

private:
//	double Rule(const Fun3D& f, double x0, double x1, double y0, double y1, double z0, double z1, 
//		double* errAbs, double* errRel, int* divAx);

//	void RuleEx(const Fun3D& fun, Brick* rec);

	void Rule(const Fun3D& fun, Brick* rec);

	void WeightAbsci();

	Brick* NewBrick();

	void Clear();

protected:
	// Liczba obliczen funkcji
	size_t m_funCall;

	// Dokladnosc z jaka ma byc obliczona calka
	double m_epsAbs;

	// Wspolrzedne punktow kwadratury. W oryginalnej pracy nazywaja sie $\lambda$
	double m_q2, m_q3, m_q4, m_q5;

	// Wagi kwadratur 7-ego rzedu
	double m_w7[5];

	// Wagi kwadratur 5-ego rzedu
	double m_w5[4];

	// Zarzadzanie pamiecia
	Pool<Brick> m_pool;
	
	// Zbior prostokatow dzielacych prostokat oryginalny
	EltSet<Brick> m_set;

};


#endif

