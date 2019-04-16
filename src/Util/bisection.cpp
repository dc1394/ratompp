#include "stdafx.h"
#include "bisection.h"

//!
//! Zwraca miejsce zerowe funkcji "f" na przedziale [a, b].
//! Obliczenia wykonane z dokladnoscia "eps".
//!
double Bisection::Zero(util::Fun1D const* f, double a, double b, double eps)
{
double va, vb, vc, c;

	assert(a < b);
//	if(a > b)
//		throw std::exc("Error! Bisection::Zero. Niepoprawny przedzial [a,b].");

	va = f->Get(a);
	vb = f->Get(b);

	assert(va * vb <= 0);
//	if(va * vb > 0)
//		throw Except("Error! Bisection::Zero. Nie jest spelniony warunek f(a) * f(b) < 0.");

	if(va == 0)
		return a;

	if(vb == 0)
		return b;


	while(b - a > eps)
	{
		c = 0.5 * (a + b);
		vc = f->Get(c);

		if(vc == 0)
			break;

		if(vc * vb < 0)
		{
			a = c;
			va = vc;
		}
		else 
		{
			b = c;
			vb = vc;
		}
	}
	return 0.5 * (a + b);
}

//!
//! Oblicza kilka miejsc zerowych. Liczba obliczanych miejsc zerowych okreslona jest przez
//! dlugosci wektora "zero".
//! xi - ograniczenie dolne miejsca zerowego
//! delta - ograniczenie dolne odleglosci miedzy miejscami zerowymi
//! xMax - maksymalny argument. Poszukiwania sa przerywane, jezeli argument jest wiekszy niz "xMax".
//!
void Bisection::ZeroSeq(util::Fun1D const* f, Vec& zero, double xMin, double, double delta, double eps)
{
double a, va, vb;
size_t i = 0;

	delta = fabs(xMin) / 50;

	a = xMin;
	va = f->Get(a);
	while(i < zero.size())
	{
		vb = f->Get(a + delta);
		// Sprawdzam czy jest miejsce zerowe na [a, b]
		if(va * vb <= 0)
		{
			zero[i] = Zero(f, a, a + delta, eps);
			a = zero[i] + 2 * eps; // Klejne miejsce zerowe moze byc bardzo blisko!
			// printf("%lf  ", a);
			va = f->Get(a);

			delta = fabs(zero[i]) / 50;
			i++;
			continue;
		}
		a += delta;
		va = vb;

//		if(a > xMax)
//			throw Except("Error! Bisection::ZeroSeq. Argument przekroczyl limit.");
	}
}

//!
//! Ogranicza (okracza) miejsce zerowe funkcji "f". Miejsce zerowe znajduje sie na przedziale [a, b]
//! xi - oszcowanie miejsca zerowego
//! a - poszukany lewy kraniec przedzialu
//! b - poszukany prawy kraniec przedzialu
//! del - krok poszukiwan
//! delMax - najwiekszy krok poszukiwan
//! alpha - wspolczynnik zwiekszania kroku poszukiwan
//!
void Bisection::Braket(util::Fun1D const* f, double xi, double del, double delMax, double alpha, double* a, double* b)
{
	assert(del > 0);
	assert(alpha > 1);

	while(true)
	{
		*a = xi - del;
		*b = xi + del;
		if(f->Get(*a) * f->Get(*b) <= 0)
			return;

		del *= alpha;

		assert(del < delMax);
//		if(del > delMax)
//			throw Except("Error! Bisection::Braket. Delta przekroczyla limit.");
	}
}


