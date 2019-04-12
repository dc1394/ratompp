#include "stdafx.h"
#include "vec.h"
#include <cstdint>  // for std::size_t
#include <memory>   // for std::unique_ptr


//
// Constructor
//
#ifdef USE_MKL
Vec::Vec(const Vec& v) : std::vector<double, util::mkl_allocator<double> >(v.size())
#else
Vec::Vec(const Vec& v) : std::vector<double>(v.size())
#endif
{
	// resize(v.size());
	for(size_t i = 0; i < v.size(); i++)
		(*this)[i] = v[i];
}

//
// Constructor
//

#ifdef USE_MKL
Vec::Vec(std::size_t n) : std::vector<double, util::mkl_allocator<double> >(n)
#else
Vec::Vec(std::size_t n) : std::vector<double>(n)
#endif
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
        std::fprintf(out.get(), "%4lu %lf\n", static_cast<unsigned long>(i), Get(i));
    }
}

