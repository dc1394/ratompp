#ifndef __RATOM_VEC_H__
#define __RATOM_VEC_H__


/** \brief Vector of real numbers
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/



class Vec : public std::vector<double>
{
public:
	Vec();
	Vec(const Vec& v);
	Vec(size_t n);
        ~Vec();

	double  Get(size_t i) const { return (*this)[i]; }
	double& Set(size_t i) { return (*this)[i]; }

	void Write(const char* path) const;

	void Zero();
	void Extrap();
};


#endif

