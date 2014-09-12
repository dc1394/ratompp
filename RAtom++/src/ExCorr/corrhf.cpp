#include "stdafx.h"
#include "corrHF.h"

//
// Constructor
//
CorrHf::CorrHf()
	: Xc(nullptr, nullptr, nullptr)
{
}


//
// Destructor
//
CorrHf::~CorrHf(void)
{
}


//
// Potencial
//
double CorrHf::V(double r) const
{
	return 0.0;
}


//
// Energy density
//
double CorrHf::E(double r) const
{
	return 0.0;
}


