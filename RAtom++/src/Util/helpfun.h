#ifndef __RATOM_HELPFUN_H__
#define __RATOM_HELPFUN_H__



/** \brief Helper functions
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/


#include "point3D.h"
#include "vec.h"

// Sometimes M_PI is defined as a macro
#ifdef M_PI
# undef M_PI
#endif
const double M_PI = 4 * atan(1.);  // \pi
const double M_4PI = 4 * M_PI;		// 4 * \pi



template <typename T>
T Max(const T& a, const T& b)
{
	return (a > b) ? a : b;
}

template <typename T>
T Min(const T& a, const T& b)
{
	return (a < b) ? a : b;
}

void Sph2Car(double r, double theta, double phi, double* x, double* y, double* z);
void Car2Sph(double x, double y, double z, double* r, double*theta, double*phi);

void RotMatT(double* rotMtxT, const Point3D* a, const Point3D* b);
void RotMatT(double* rotMtxT, const double* u, double alpha);
void RotMatTsimple(double* rotMtxT, const double* u, double alpha);

double Pow(double x, int n);
double Factor(int n);
double Acosh(double x);

void Spline(const std::vector<double>& x, const std::vector<double>& y, double yp1, double ypn, std::vector<double>& y2);
double SplInt(const std::vector<double>& xa, const std::vector<double>& ya, const std::vector<double>& y2a, double x);
double SplIntD2(const std::vector<double>& xa, const std::vector<double>& y2a, double x);

void gauleg(double x1, double x2, Vec& x, Vec& w, size_t n);


#endif

