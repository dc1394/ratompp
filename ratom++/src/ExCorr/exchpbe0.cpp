//#include "stdafx.h"
//#include "exchPBE0.h"
//
//#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
//	const double ExchPbe0::ZERO = 1.0E-12;
//#endif
//
////
//// Constructor
////
//// March 8th, 2014	Modified by dc1394
////ExchPbe0::ExchPbe0(void) : m_c(-pow(1.5 / M_PI, 2. / 3.))
//ExchPbe0::ExchPbe0(std::function<double(double)> rhoTilde,
//				   std::function<double(double)> rhoTildeDeriv,
//				   std::function<double(double)> rhoTildeLapl,
//				   std::function<double(double)> Vh, double Z)
//	: Xc(std::move(rhoTilde), std::move(rhoTildeDeriv), std::move(rhoTildeLapl)),
//	  Vh_(std::move(Vh)),
//	  Z_(Z)
//{
//	xc_func_init(pxcfunc_.get(), XC_GGA_X_PBE, XC_UNPOLARIZED);
//}
//
//
//// March 8th, 2014	Modified by dc1394
////
//// Destructor
////
//ExchPbe0::~ExchPbe0(void)
//{
//}
//
//
////
//// Potencial
////
//double ExchPbe0::V(double r) const
//{
//	const double exchvpbe = my_xc_gga_vxc(r);
//
//	if (std::fabs(Z_ - 1.0) < ZERO)
//		return 0.75 * exchvpbe - 0.25 * Vh_(r);
//	else if (std::fabs(Z_ - 2.0) < ZERO)
//		return 0.75 * exchvpbe - 0.125 * Vh_(r);
//	else
//		throw std::runtime_error("–¢ŽÀ‘•");
//}
//
//
////
//// Energy density
////
//double ExchPbe0::E(double r) const
//{
//	const double exchepbe = my_xc_gga_exc(r);
//
//	if (std::fabs(Z_ - 1.0) < ZERO)
//		return 0.75 * exchepbe - 0.125 * Vh_(r);
//	else if (std::fabs(Z_ - 2.0) < ZERO)
//		return 0.75 * exchepbe - 0.0625 * Vh_(r);
//	else
//		throw std::runtime_error("–¢ŽÀ‘•");
//}
//
