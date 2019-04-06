#include "stdafx.h"
#include "vec.h"
#include <cstdio>   // for std::fopen
#include <memory>   // for std::unique_ptr


//
// Constructor
//
Vec::Vec(const Vec& v) : std::vector<double>(v.size())
{
	// resize(v.size());
	for(size_t i = 0; i < v.size(); i++)
		(*this)[i] = v[i];
}

//
// Constructor
//
Vec::Vec(size_t n) : std::vector<double>(n)
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
    auto out = std::unique_ptr<FILE, decltype(&std::fclose)>(std::fopen(path, "wt"), std::fclose);

	assert(out);
    for (size_t i = 0; i < size(); i++)
    {
        fprintf(out.get(), "%4lu %lf\n", static_cast<unsigned long>(i), Get(i));
    }
}
