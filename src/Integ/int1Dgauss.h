#ifndef __RATOM_INT1DGAUSS_H__
#define __RATOM_INT1DGAUSS_H__


/**
* \brief Reprezentuje algorytm calkowania numerycznego funkcji jednej zmiennej za pomoca kwadratur Gaussa.
*
* \note Calkowanie wykonywane jest na przedziale [a, b]. 
*	Kwadratury Gaussa sa przeznaczone do calkowania na przedziale [-1, 1]. Dlatego
*	nalezy wykonac zmiane zmiennych, aby calkowac na dowolnym przedziale [a, b].
*	te podstawienie jest wykonywane wewnatrz skladowej Calc().
*
*	Stopien kwadratury Gaussa okresla liczbe punktow na przedziale [-1, 1], dla ktorych
*	obliczana jest wartosc funkcji. Mowiac nieprecyzyjnie, liczba wezlow determinuje
*	dokladnosc obliczanej calki.
*
*	*** Komentarz numeryczny ***
*	Kwadratury Gaussa zdefiniowane sa na przedziale [-1, 1], wiec nalezy zrobic 
*	podstawienie zmieniajace granice calkowania z [a, b] na [-1, 1]. Mozna to zrobic
*	za pomoca przeksztalcenia liniowego:
*
*      x(t) = P * t + Q    ==> dx = P * dt
*
*	Przeksztalcenie to musi spelniac warunki:
*		x(-1) = a   oraz   x(1) = b
*	co jest rownowazne ukladowi rownan na P i Q
*		a = -P + Q
*		b =  P + Q
*	Z powyzszego otrzymujemy Q = (a + b) / 2,  P = (b - a) / 2
*	Ostatecznie otrzymujemy:
*	\f[
*		I = \int_a^b f(x) dx = P * \int_{-1}^1 f(P * t + Q) dt
*	\f]
*	Ten wzor jest wykorzystany do obliczania calki
*
*	*** Uwaga ***
*	Wagi kwadratur sa zapisane w tablicy $w$.
*	Wspolrzedne kwadratur zapisane sa w tablicy $t$.
*
* \author Zbigniew Romanowski [ROMZ]
*
* \version 10-Mar-2005 (romz) - utworzenie
*/

#include "int1D.h"

class Int1DGauss final : public Int1D
{
public:
	Int1DGauss(size_t deg);
	~Int1DGauss() override = default;

	double Calc(util::Fun1D const& f, double a, double b) const override;
	double Adapt(util::Fun1D const& f, double a, double b, double errAbs) const;

	// static void gauleg(double x1, double x2, Vec& x, Vec& w, size_t n);
	
private:
	double Run(util::Fun1D const& f, double p, double q) const;

private:
	// Stopien kwadratury
	size_t m_deg;

	// Wagi kwadratury
	Vec m_w;

	// Wspolrzedne wezlow
	Vec m_t;
};

#endif

