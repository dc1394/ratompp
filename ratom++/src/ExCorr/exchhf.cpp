//#include "stdafx.h"
//#include "exchHF.h"
//#include <utility>
//#include <cmath>
//
//#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
//	const double ExchHf::ZERO = ZERO;
//#endif
//
////
//// Constructor
////
//ExchHf::ExchHf(std::function<double(double)> Vh, double Z)
//	: Xc(nullptr, nullptr, nullptr),
//	  Vh_(std::move(Vh)),
//	  Z_(Z)
//{
//}
//
//
////
//// Destructor
////
//ExchHf::~ExchHf(void)
//{
//}
//
//
////
//// Potencial
////
//double ExchHf::V(double r) const
//{
//	if (std::fabs(Z_ - 1.0) < ZERO)
//		return -Vh_(r);
//	else if (std::fabs(Z_ - 2.0) < ZERO)
//		return -0.5 * Vh_(r);
//	else
//		throw std::runtime_error("–¢ŽÀ‘•");
//}
//
//
////
//// Energy density
////
//double ExchHf::E(double r) const
//{
//	if (std::fabs(Z_ - 1.0) < ZERO)
//		return -0.5 * Vh_(r);
//	else if (std::fabs(Z_ - 2.0) < ZERO)
//		return - 0.25 * Vh_(r);
//	else
//		throw std::runtime_error("–¢ŽÀ‘•");
//}
//
