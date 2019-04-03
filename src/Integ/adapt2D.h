#ifndef __RATOM_ADAPT2D_H__
#define __RATOM_ADAPT2D_H__


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
#include "rec.h"
#include "pool.h"

enum Adapt2DRule
{
	Adapt2DRule_57,
	Adapt2DRule_79,
	Adapt2DRule_911,
	Adapt2DRule_913
};



class Adapt2D
{
public:
	Adapt2D(double epsAbs, size_t initSize, Adapt2DRule rule = Adapt2DRule_79);
	~Adapt2D(void);

	size_t RecNo() const;
	size_t FunCall() const;

	double GetEpsAbs(void) const { return m_epsAbs; }
	void SetEpsAbs(double epsAbs) { m_epsAbs = epsAbs; }

	double Calc(const Fun2D& f, double x0, double x1, double y0, double y1);

private:
	void Rule(const Fun2D& fun, Rec* rec);
	void Rule57(const Fun2D& fun, Rec* rec);
	void Rule79(const Fun2D& fun, Rec* rec);
	void Rule911(const Fun2D& fun, Rec* rec);
	void Rule913(const Fun2D& fun, Rec* rec);

	double Rule13(const Fun2D& fun, Rec* rec) const;
	double Rule11(const Fun2D& fun, Rec* rec) const;
	double Rule9 (const Fun2D& fun, Rec* rec) const;

	void WeightAbsci();

	Rec* NewRec();

	void Clear();

protected:
	// Liczba obliczen funkcji
	size_t m_funCall;

	// Dokladnosc z jaka ma byc obliczona calka
	double m_epsAbs;

	// Zbior prostokatow dzielacych prostokat oryginalny
	EltSet<Rec> m_set;
	// Zarzadzanie pamiecia
	Pool<Rec> m_pool;

	// Sposob obliczania kwadratur
	const Adapt2DRule m_rule;
};


#endif

