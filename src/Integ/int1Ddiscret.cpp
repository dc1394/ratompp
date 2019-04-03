#include "stdafx.h"
#include "int1Ddiscret.h"



Int1DDiscret::Int1DDiscret(void)
{
}

Int1DDiscret::~Int1DDiscret(void)
{
}

/*
//!
//! Zwraca wartosc calki liczona w sferycznym ukladzie wspolrzednych z funkcji "f".
//! Wartosci funkcji sa zapisane w wektorze "f".
//!
double Int1DDiscret::Spher(const Vec& f, const Mesh* mesh)
{
Vec tmp(f.size());

	assert(mesh);
	assert(mesh->NodeNo() == f.size());

	for(size_t i = 0; i < f.size(); i++)
		tmp[i] = mesh->R2(i) * f[i];

	return 4. * util::HelpFun::M_PI * Order4a(tmp, mesh->Dr());
}

//!
//! Zwraca wartosc calki liczona w sferycznym ukladzie wspolrzednych z iloczynu "f * g".
//! Wartosci funkcji sa zapisane w wektorze "f".
//!
double Int1DDiscret::Spher(const Vec& f, const Vec& g, const Mesh* mesh)
{
Vec tmp(f.size());

	assert(mesh);
	assert(mesh->NodeNo() == f.size());
	assert(mesh->NodeNo() == g.size());

	for(size_t i = 0; i < f.size(); i++)
		tmp[i] = mesh->R2(i) * f[i] * g[i];

	return 4. * util::HelpFun::M_PI * Order4a(tmp, mesh->Dr());

}
*/


//!
//! Zwraca wartosc calki z funkcji o  wartosciach zapisanych w wektorze "f".
//! Forth order integration routine
//! Nie wiem skad ten wzor pochodzi
//! "h" - odleglosc miedzy wezlami
//!
double Int1DDiscret::Order4a(const Vec& f, double h)
{
const double c[4] = {17./48., 59./48., 43./48., 49./48.};
const size_t n = f.size();

	assert(n >= 8);

double sum = c[0] * f[0] + c[1] * f[1] + c[2] * f[2] + c[3] * f[3];

	for(size_t i = 4; i < n - 4; i++)
		sum += f[i];

	// double qq = f[n - 1];
	sum += c[0] * f[n - 1] + c[1] * f[n - 2] + c[2] * f[n - 3] + c[3] * f[n - 4];
	return sum * h;

}

//!
//! Zwraca wartosc calki na odcinku.
//! Wartosci funkcji sa zapisane w wektorze "f".
//! Wzor (4.1.12) z Numerical Recipes
//! "h" - odleglosc miedzy wezlami siatki
//!
double Int1DDiscret::Order3(const Vec& f, double h)
{
const double c[4] = {5./12., 13./12., 13./12., 5./12.};
const size_t n = f.size();

	assert(n > 4);

double sum = c[0] * f[0] + c[1] * f[1] + c[2] * f[n - 2] + c[3] * f[n - 1];

	for(size_t i = 2; i < n - 2; i++)
		sum += f[i];

	return sum * h;
}

//!
//! Zwraca wartosc calki z funkcji o  wartosciach zapisanych w wektorze "f".
//! Wzor (4.1.14) z Numerical Recipes
//! "h" - odleglosc miedzy wezlami
//!
double Int1DDiscret::Order4(const Vec& f, double h)
{
const double c[3] = {3./8., 7./6., 23./24.};
const size_t n = f.size();
double sum = (c[0] * f[0] + c[1] * f[1] + c[2] * f[2]);

	sum += (c[0] * f[n - 1] + c[1] * f[n - 2] + c[2] * f[n - 3]);
	for(size_t i = 3; i < n - 3; i++)
		sum += f[i];

	return sum * h;
}

//!
//! Zwraca wartosc calki z funkcji o  wartosciach zapisanych w wektorze "f".
//! Jankowscy "Przeglad algorytmow i metod numerycznych"
//! "h" - odleglosc miedzy wezlami
//!
double Int1DDiscret::NewtonCotes3(const Vec& f, double h)
{
const double c[4] = {1./8., 3./8., 3./8., 1./8.};
const size_t n = f.size();
double sum = 0;

	assert(n >= 4);
	assert((n - 1) % 3 == 0);

	// sum = (c[3] * f[n - 1] + c[2] * f[n - 2] + c[1] * f[n - 3] + c[0] * f[n - 4]);
	for(size_t i = 0; i < n - 1; i += 3)
		sum += (c[0] * f[i] + c[1] * f[i + 1] + c[2] * f[i + 2] + c[3] * f[i + 3]);

	return 3 * h * sum;
}

//!
//! Zwraca wartosc calki z funkcji o  wartosciach zapisanych w wektorze "f".
//! Jankowscy "Przeglad algorytmow i metod numerycznych"
//! "h" - odleglosc miedzy wezlami
//!
double Int1DDiscret::NewtonCotes4(const Vec& f, double h)
{
const double c[5] = {7./90., 32./90., 12./90., 32./90., 7./90.};
const size_t n = f.size();
double sum = 0;

	assert(n >= 5);
	assert((n - 1) % 4 == 0);

	for(size_t i = 0; i < n - 1; i += 4)
		sum += (c[0] * f[i] + c[1] * f[i + 1] + c[2] * f[i + 2] + c[3] * f[i + 3] + c[4] * f[i + 4]);

	return 4 * h * sum;
}

//!
//! Zwraca wartosc calki z funkcji o  wartosciach zapisanych w wektorze "f".
//! Jankowscy "Przeglad algorytmow i metod numerycznych"
//! "h" - odleglosc miedzy wezlami
//!
double Int1DDiscret::NewtonCotes5(const Vec& f, double h)
{
const double c[6] = {19./288., 75./288., 50./288., 50./288., 75./288., 19./288.};
const size_t n = f.size();
double sum = 0;
size_t i, j;

	assert(n >= 6);
	assert((n - 1) % 5 == 0);

	for(i = 0; i < n - 1; i += 5)
		for(j = 0; j < 6; j++)
			sum += c[j] * f[i + j];

	return 5 * h * sum;
}


//!
//! Zwraca wartosc calki z funkcji o  wartosciach zapisanych w wektorze "f".
//! Jankowscy "Przeglad algorytmow i metod numerycznych"
//! "h" - odleglosc miedzy wezlami
//!
double Int1DDiscret::NewtonCotes6(const Vec& f, double h)
{
const double c[7] = {41./840., 216./840., 27./840., 272./840., 27./840., 216./840., 41./840.};
const size_t n = f.size();
double sum = 0;
size_t i, j;

	assert(n >= 7);
	assert((n - 1) % 6 == 0);

	for(i = 0; i < n - 1; i += 6)
		for(j = 0; j < 7; j++)
			sum += c[j] * f[i + j];

	return 6 * h * sum;
}

//!
//! Zwraca wartosc calki z funkcji o  wartosciach zapisanych w wektorze "f".
//! Liczba elementow w wektorze "f" dowolna.
//!
double Int1DDiscret::NewtonCotes6Ex(const Vec& f, double h)
{
const double c[7] = {41./840., 216./840., 27./840., 272./840., 27./840., 216./840., 41./840.};
const size_t n = f.size();
double sumA, sumB;
size_t i, j, q;


	assert(n >= 7);

	q = (n - 1) % 6;
	assert((n - 1 - q) % 6 == 0);

	sumA = 0;
	for(i = 0; i < n - 1 - q; i += 6)
		for(j = 0; j < 7; j++)
			sumA += c[j] * f[i + j];

	// Dla reszy wezlow metoda trapezow
	sumB = 0;
	for(i = n - q - 1; i < n - 1; i++)
		sumB += (f[i] + f[i + 1]);

	return (6 * h * sumA + 0.5 * h * sumB);
}


