//#ifndef __RATOM_EXCHHF_H__
//#define __RATOM_EXCHHF_H__
//
//
///** \brief HF Exchange
//*
//* \author @dc1394
//* March 28th, 2014
//*/
//
//
//#include "xc.h"
//
//
//
//class ExchHf : public Xc
//{
//public:
//	// March 7th, 2014	Modified by dc1394
//	ExchHf(std::function<double(double)> Vh, double Z);
//	virtual ~ExchHf(void);
//
//	// March 7th, 2014	Modified by dc1394
//	virtual double V(double r) const;
//	virtual double E(double r) const;
//
//	virtual const char* Name() const
//	{
//		return "hartree-fock";
//	}
//
//
//private:
//	// Debug only
//#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__)
//	static constexpr double ZERO = 1.0E-12;
//#else
//	static const double ZERO;
//#endif
//	const std::function<double(double)> Vh_;
//	const double Z_;
//};
//
//#endif	// __RATOM_EXCHHF_H__
