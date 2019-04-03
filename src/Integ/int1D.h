#ifndef __RATOM_INT1D_H__
#define __RATOM_INT1D_H__


/**
*	\brief Oblicza calke w jednym wymiarze po odcinku. Klasa abstrakcyjna.
*
*	\author Zbigniew Romanowski, ROMZ
*
*	\version 12-Cze-2007 [romz]
*
*/

class Int1D
{
public:
	Int1D() = default;
	virtual ~Int1D() = default;

	// Zwraca wartosc calki z funkcji $f$ na prezdziale $[a, b]$.
	virtual double Calc(const util::Fun1D& f, double a, double b) const = 0;
};

#endif

