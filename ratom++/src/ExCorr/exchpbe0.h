//#ifndef __RATOM_EXCHPBE0_H__
//#define __RATOM_EXCHPBE0_H__
//
//
///** \brief PBE0 Exchange
//*
//* \author @dc1394
//* April 1st, 2014
//*/
//
//
//#include "xc.h"
//
//
//
//class ExchPbe0 : public excorr::Xc
//{
//public:
//	ExchPbe0(std::function<double(double)> rhoTilde,
//			 std::function<double(double)> rhoTildeDeriv,
//			 std::function<double(double)> rhoTildeLapl,
//			 std::function<double(double)> Vh, double Z);
//	virtual ~ExchPbe0(void);
//
//	// March 7th, 2014	Modified by dc1394
//	virtual double V(double r) const;
//	virtual double E(double r) const;
//
//	virtual const char* Name() const
//	{
//		return "PBE0 Exchange";
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
//
//	const std::function<double(double)> Vh_;
//	const double Z_;
//};
//
//#endif	// __RATOM_EXCHPBE0_H__
