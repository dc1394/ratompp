#include "stdafx.h"
#include "vec.h"
#include <stdio.h>



//
// Constructor
//
Vec::Vec(const Vec& v) : std::vector<double, util::mkl_allocator<double> >(v.size())
{
	// resize(v.size());
	for(size_t i = 0; i < v.size(); i++)
		(*this)[i] = v[i];
}

//
// Constructor
//
Vec::Vec(size_t n) : std::vector<double, util::mkl_allocator<double> >(n)
{
}

//
// Set all elements to zero
//
void Vec::Zero()
{
	for(size_t i = 0; i < size(); i++)
		(*this)[i] = 0.;
}

//
// Sets value of "the first" element of vector based on "linear extrapolation"
//
void Vec::Extrap()
{
	assert(size() > 2);
	(*this)[0] = 2 * (*this)[1] - (*this)[2];
}

//
// Saves vector into file
//
void Vec::Write(const char* path) const
{
FILE* out = fopen(path, "wt");

	assert(out);
	for(size_t i = 0; i < size(); i++)
		fprintf(out, "%4lu %lf\n", static_cast<unsigned long>(i), Get(i));

	fclose(out);
}
