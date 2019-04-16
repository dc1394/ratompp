#ifndef __RATOM_BISECTION_H__
#define __RATOM_BISECTION_H__

#pragma once

/** \brief Bisection algorithm with some addtional functionality.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
* \note It searches zero of continouous function.
*/


#include "fun1D.h"
#include "vec.h"

class Bisection final
{
public:
	Bisection() = default;
	~Bisection() = default;

	static double Zero(util::Fun1D const* f, double a, double b, double eps);

	static void ZeroSeq(util::Fun1D const* f, Vec& zero, double xMin, double, double delta, double eps);

	static void Braket(util::Fun1D const* f, double xi, double del, double delMax, double alpha, double* a, double* b);
};

#endif

