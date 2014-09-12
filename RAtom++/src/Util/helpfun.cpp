#include "stdafx.h"
#include "helpfun.h"
#include <float.h>


//
// Transformation of spherical coordinates (r, theta, phi)
// into cartesian coordinates (x, y, z)
//
void Sph2Car(double r, double theta, double phi, double* x, double* y, double* z)
{
const double st = sin(theta);

	assert(x && y && z);

	*x = r * st * cos(phi);
	*y = r * st * sin(phi);
	*z = r * cos(theta);
}


//
// Transformation of cartesian coordinates (x, y, z)
// into spherical coordinates (r, theta, phi)
//
void Car2Sph(double x, double y, double z, double* r, double*theta, double*phi)
{
	*r = sqrt(x * x + y * y + z * z);
	*theta = atan2(sqrt(x * x + y * y), z);
	*phi = atan2(y, x);
}

//
// Creation  of transposed rotation matrix "rotMtxT".
// Rotation angle and rotation axis is determined by two points "a" and "b".
// After rotation vector from "a" to "b" is position along versor (0,0,1)
//
void RotMatT(double* rotMtxT, const Point3D* a, const Point3D* b)
{
double u[3], qx, qy, qz, q;

	qx = b->m_x - a->m_x;
	qy = b->m_y - a->m_y;
	qz = b->m_z - a->m_z;

	q = sqrt(qx * qx + qy * qy);

	if(q < FLT_MIN) // Rotation matrix is identity matrix
	{
		rotMtxT[0] = 1; rotMtxT[1] = 0; rotMtxT[2] = 0;
		rotMtxT[3] = 0; rotMtxT[4] = 1; rotMtxT[5] = 0;
		rotMtxT[6] = 0; rotMtxT[7] = 0; rotMtxT[8] = 1;
		return;
	}

	u[0] =  qy / q;
	u[1] = -qx / q;
	u[2] = 0;

	
	q = sqrt(qx * qx + qy * qy + qz * qz);

	const double alpha = acos(qz / q);

	// RotMatT(rotMtxT, u, alpha);
	RotMatTsimple(rotMtxT, u, alpha);
}


//
// Creation  of transposed rotation matrix "rotMtxT".
// "alpha" is a rotation angle and "u[3]" is a rotation axis.
// One-dimensional vector "rotMtx" represents matrix in following way
//           |r0 r3 r6|
//  rotMtx = |r1 r4 r7|
//           |r2 r5 r8|
// 
void RotMatT(double* rotMtxT, const double* u, double alpha)
{
const double ca = cos(alpha), sa = sin(alpha);
const double ux = u[0], uy = u[1], uz = u[2]; // For short-hand notation

	rotMtxT[0] = ux * ux + ca * (1 - ux * ux);
	rotMtxT[3] = ux * uy * (1 - ca) - uz * sa;
	rotMtxT[6] = uz * ux * (1 - ca) + uy * sa;

	rotMtxT[1] = ux * uy * (1 - ca) + uz * sa;
	rotMtxT[4] = uy * uy + ca * (1 - uy * uy);
	rotMtxT[7] = uy * uz * (1 - ca) - ux * sa;

	rotMtxT[2] = uz * ux * (1 - ca) - uy * sa;
	rotMtxT[5] = uy * uz * (1 - ca) + ux * sa;
	rotMtxT[8] = uz * uz + ca * (1 - uz * uz);
}


//
// Creation  of transposed rotation matrix "rotMtxT".
// "alpha" is a rotation angle.
// "u[2]" is a rotation axis it is assumed that "u[3] = 0"
// One-dimensional vector "rotMtx" represents matrix in following way
//           |r0 r3 r6|
//  rotMtx = |r1 r4 r7|
//           |r2 r5 r8|
// It is simplified version of function RotMatT
//
void RotMatTsimple(double* rotMtxT, const double* u, double alpha)
{
const double ca = cos(alpha), sa = sin(alpha);
const double ux = u[0], uy = u[1]; // Short-hand notation

	
	rotMtxT[0] = ux * ux + ca * (1 - ux * ux);
	rotMtxT[1] = ux * uy * (1 - ca);
	rotMtxT[2] = -uy * sa;
	
	rotMtxT[3] = ux * uy * (1 - ca);
	rotMtxT[4] = uy * uy + ca * (1 - uy * uy);
	rotMtxT[5] = ux * sa;

	rotMtxT[6] = uy * sa;
	rotMtxT[7] = -ux * sa;
	rotMtxT[8] = ca;
}


//
// Returns $x^n$ for integer "n" and real "x".
//
double Pow(double x, int n)
{
double v = 1.;

	switch(labs(n))
	{
		case 0: v = 1.;				break;
		case 1: v = x;				break;
		case 2: v = x * x;			break;
		case 3: v = x * x * x;			break;
		case 4: v = x * x * x * x;		break;
		case 5: v = x * x * x * x * x;		break;
		case 6: v = x * x * x * x * x * x;	break;

		default:
			for(int i = 1; i <= labs(n); i++)
				v *= x;
		break;
	}
	if(n < 0)
		return 1. / v;
	return v;
}


//
// Returns "n!".
//
double Factor(int n)
{
double res = 1.;

	assert(n >= 0);
	// 0! = 1
	// 1! = 1
	// 2! = 2
	for(int i = 2; i <= n; i++)
		res *= i;

	return res;
}

//
// Returns value of arccosh (inverse hyperbolic cosine)
//
double Acosh(double x)
{
	return log(x + sqrt(x * x - 1));
}


// 
// Given arrays $x[1..n]$ and $y[1..n]$ containing a tabulated function, i.e., $y_i = f(x_i)$, with
// $x_1 < x_2 < .. . < x_N$, and given values $yp_1$ and $yp_n$ for the first derivative of the interpolating
// function at points 1 and n, respectively, this routine returns an array $y2[1..n]$ that contains
// the second derivatives of the interpolating function at the tabulated points $x_i$. If $yp_1$ and/or
// $yp_n$ are equal to $1  10^30$ or larger, the routine is signaled to set the corresponding boundary
// condition for a natural spline, with zero second derivative on that boundary.
//
// Procedure from Numerical Recipes in C with std::vector<double>
//
void Spline(const std::vector<double>& x, const std::vector<double>& y, double yp1, double ypn, std::vector<double>& y2)
{
std::vector<double> u(x.size());
size_t i, k;
double p, qn, sig, un;
const size_t n = x.size();

	assert(x.size() == y.size());
	assert(x.size() == y2.size());


	// The lower boundary condition is set to be "natural"...
	if(yp1 > 0.99e30)
	{
		y2[0] = 0.;
		u[0] = 0.;
	}
	else // ... or else to have a specific first derivative
	{
		y2[0] = -0.5;
		u[0] = (3. / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
	}

	// This is the decomposition loop of the tridiagonal algorithm.
	// y2 and u are used for temporary storage of the decomposed factors.
	for(i = 1; i < n - 1; i++)
	{
		sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
		p = sig * y2[i - 1] + 2;
		y2[i] = (sig - 1) / p;
		u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
		u[i] = (6 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
	}

	// The upper boundary condition is set to be "natural"...
	if(ypn > 0.99e30)
	{
		qn = 0.0;
		un = 0.0;
	}
	else // ... or else to have a specific first derivative
	{
		qn = 0.5;
		un = (3 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
	}
	
	y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);

	// This is the backsubstitution loop of the tridiagonal algorithm. 
	for(k = n - 1; k >= 1; k--) 
		y2[k - 1] = y2[k - 1] * y2[k] + u[k - 1];
}


//
// Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xa's in order),
// and given the array y2a[1..n], which is the output from spline above, and given a value of
// x, this routine returns a cubic-spline interpolated value y.
//
// Procedure from Numerical Recipes in C with std::vector<double>
//
double SplInt(const std::vector<double>& xa, const std::vector<double>& ya, const std::vector<double>& y2a, double x)
{
size_t klo, khi, k;
double h, b, a;
const size_t n = xa.size();

	assert(xa.size() == ya.size());
	assert(xa.size() == y2a.size());

	// We will find the right place in the table by means of
	// bisection. This is optimal if sequential calls to this
	// routine are at random values of x. If sequential calls
	// are in order, and closely spaced, one would do better
	// to store previous values of klo and khi and test if
	// they remain appropriate on the next call.

	klo = 0;
	khi = n - 1;
	while(khi - klo > 1)
	{
		k = (khi + klo) >> 1;
		if(xa[k] > x)
			khi = k;
		else
			klo = k;
	}
	// klo and khi now bracket the input value of x.
	// March 19th, 2014 Modified by dc1394
	//assert(xa[klo] <= x && x <= xa[khi]);

	h = xa[khi] - xa[klo];
	// if (h == 0.0) nrerror("Bad xa input to routine splint"); The xas must be distinct.
	assert(h != 0.);

	a = (xa[khi] - x) / h;
	// b = (x - xa[klo]) / h; // Aby zaoszczedzic jedno dzielenie
	b = 1 - a;

	// Cubic spline polynomial is now evaluated.
	return a * ya[klo] + b * ya[khi] + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
}

//
// Returns value of second derivative of spline for "x"
// Procedure from Numerical Recipes in C with std::vector<double>
//
double SplIntD2(const std::vector<double>& xa, const std::vector<double>& y2a, double x)
{
size_t klo, khi, k;
double h, b, a;
const size_t n = xa.size();

	assert(xa.size() == y2a.size());

	// We will find the right place in the table by means of
	// bisection. This is optimal if sequential calls to this
	// routine are at random values of x. If sequential calls
	// are in order, and closely spaced, one would do better
	// to store previous values of klo and khi and test if
	// they remain appropriate on the next call.

	klo = 0;
	khi = n - 1;
	while(khi - klo > 1)
	{
		k = (khi + klo) >> 1;
		if(xa[k] > x)
			khi = k;
		else
			klo = k;
	}
	// klo and khi now bracket the input value of x.
	assert(xa[klo] < x && x <= xa[khi]);

	h = xa[khi] - xa[klo];
	// if (h == 0.0) nrerror("Bad xa input to routine splint"); The xa's must be distinct.
	assert(h != 0.);

	a = (xa[khi] - x) / h;
	// b = (x - xa[klo]) / h; // In order to avoid one division
	b = 1 - a;

	// Interpolowana jest druga pochodna
	return a * y2a[klo] + b * y2a[khi];
}



// 
// Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
// arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Int1DGauss-
// Legendre n-point quadrature formula.
// "eps" - accuracy of Gauss coordinates calculations
// 
void gauleg(double x1, double x2, Vec& x, Vec& w, size_t n)
{
// const double eps = 100.0 * DBL_MIN;
const double eps = 1E-14;
size_t m, j, i;
double z1, z, xm, xl, pp, p3, p2, p1; // High precision is a good idea for this routine.
const double pi = 4 * atan(1.);

	assert(x.size() == n);
	assert(w.size() == n);

	m = (n + 1) / 2; // The roots are symmetric in the interval, so
	xm = 0.5 * (x2 + x1); //we only have to find half of them.
	xl = 0.5 * (x2 - x1);
	for(i = 1; i <= m; i++) // Loop over the desired roots.
	{ 
		z = cos(pi * (i - 0.25) / (n + 0.5));
		// Starting with the above approximation to the ith root, we enter the main loop of
		// refinement by Newton's method.
		do
		{
			p1 = 1.0;
			p2 = 0.0;
			for(j = 1; j <= n; j++) //Loop up the recurrence relation to get the
			{
				p3 = p2; // Legendre polynomial evaluated at z.
				p2 = p1;
				p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
			}
			// p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
			// by a standard relation involving also p2, the polynomial of one lower order.
			pp = n * (z * p1 - p2) / (z * z - 1.0);
			z1 = z;
			z = z1 - p1 / pp; // Newton's method.
		}
		while(fabs(z - z1) > eps);

		x[i - 1] = xm - xl * z; // Scale the root to the desired interval,
		x[n - i] = xm + xl * z; // and put in its symmetric counterpart.
		w[i - 1] = 2.0 * xl / ((1.0 - z * z) * pp * pp); //Compute the weight
		w[n - i] = w[i - 1]; //and its symmetric counterpart.
	}
}


