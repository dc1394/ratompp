#include "stdafx.h"
#include "corrvwn.h"

// March 18th, 2014	Added by dc1394
#include <utility>

// March 18th, 2014	Modified by dc1394
//CorrVwn::CorrVwn(void) : 
CorrVwn::CorrVwn(std::function<double(double)> rhoTilde)
	:	Xc(std::move(rhoTilde), nullptr, nullptr)
{
	
}

CorrVwn::~CorrVwn(void)
{
}

//!
//! Potencjal
//!
double CorrVwn::V(double rho, double gRho) const
{
double ec, vc;

	Help(rho, &ec, &vc);
	return vc;
}


// March 18th, 2014 Added by dc1394
//!
//! Potential
//!
double CorrVwn::V(double r) const
{

	return my_xc_lda_vxc(r);
}


//!
//! Gestosc energii
//!
double CorrVwn::E(double rho, double gRho) const
{
double ec, vc;

	Help(rho, &ec, &vc);
	return ec;
}


// March 18th, 2014 Added by dc1394
//!
//! Gestosc energii
//!
double CorrVwn::E(double r) const
{
	return my_xc_lda_exc(r);
}


//!
//! Roznica gestosci elektronowej i potencjalu. Wersja zoptymalizowana :-)
//!
double CorrVwn::EdiffV(double rho, double gRho) const
{
double ec, vc;

	Help(rho, &ec, &vc);
	return ec - vc;
}


//!
//! Funkcja pomocnicza
//!
void CorrVwn::Help(double rho, double* ec, double* vc) const
{
	if(rho == 0)
	{
		*ec = 0;
		*vc = 0;
		return;
	}

const double a = 0.0621814, b = 3.72744, c = 12.9352, x0 = -0.10498;
const double rs = Rs(rho);

const double x = sqrt(rs);
const double q = sqrt(4. * c - b * b);
const double f1 = 2. * b / q;
const double f2 = -b * x0 / (x0 * x0 + b * x0 + c);
const double f3 = 2. * (b + 2. * x0) * f2 / q;

const double px = rs + b * x + c;
const double dx = rs / px;
const double yx = atan(q / (2. * x + b));
const double wx = (x - x0) * (x - x0) / px;

	*ec = 0.5 * a * (log(dx) + f1 * yx + f2 * log(wx) + f3 * yx);
	
const double v1 = c * (x - x0) - b * x * x0;
const double v2 = (x - x0) * px;

	*vc = (*ec) - a * v1 / (6. * v2);
}


