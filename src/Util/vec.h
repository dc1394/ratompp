#ifndef __RATOM_VEC_H__
#define __RATOM_VEC_H__


/** \brief Vector of real numbers
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/

#ifdef USE_MKL
    #include "mkl_allocator.h"
#endif
#include <vector>           // for std::vector

#ifdef USE_MKL
class Vec final : public std::vector<double, util::mkl_allocator<double> >
#else
class Vec final : public std::vector<double>
#endif
{
public:
	Vec() = default;
	Vec(const Vec& v);
	Vec(size_t n);
    ~Vec() = default;

	double  Get(size_t i) const { return (*this)[i]; }
	double& Set(size_t i) { return (*this)[i]; }

	void Write(const char* path) const;

	void Zero();
	void Extrap();
};


#endif

