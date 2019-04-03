#ifndef __RATOM_INT1DTRAP_H__
#define __RATOM_INT1DTRAP_H__


/**
* \brief Reprezentuje metode trapezow (z ekrapolacja Romberga) calkowania numerycznego 
*	funkcji jednej zmiennej na przedziale [a,b]. 
*
* \note Czasmi ten algorytm nazywa sie metoda Romberga.
*
* \author Zbigniew Romanowski [ROMZ]
*
* \version 09-Mar-2005 (romz) - przepisanie (plus mala modyfikacja) ze starej wersji
*
* \todo Dodac referencje.
*/

#include "../Util/fun1D.h"
#include "int1D.h"


class Int1DTrap : public Int1D
{
public:
	// ext = 20, to z powodzeniem w wiekszosci przypadkow wystrczy, gdyz jest 
	// to 2^{20} podzialow odcinka [a,b]
	Int1DTrap(double eps, size_t ext = 20);
	virtual ~Int1DTrap();

	virtual double Calc(const util::Fun1D &f, double a, double b) const;

private:
	double Val(const util::Fun1D& f, double a, double b, long int divNo) const;
	double Fast(const util::Fun1D& f, double a, double b, long int divNo, double resOld) const;

private:
	// Maksymalna liczba ektrapolacji
	size_t m_ext;

	// Dokladnosc obliczenia calki 
	double m_eps;

	// Macierz kwadratowa (ext x ext) przechowujaca dane wykorzystywane do ekstrapolacji
	double** m_p;

	// Wartosci mianownika. Przechowywane w tablicy (o rozmiarze ext), aby bylo szybciej.
	double* m_denom;
};

#endif

