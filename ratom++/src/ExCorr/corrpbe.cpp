//#include "stdafx.h"
//#include "corrpbe.h"
//
//#include <utility>
//
//CorrPbe::CorrPbe(std::function<double(double)> rhoTilde, std::function<double(double)> rhoTildeDeriv, std::function<double(double)> rhoTildeLapl)
//	:	Xc(std::move(rhoTilde), std::move(rhoTildeDeriv), std::move(rhoTildeLapl))
//{
//	xc_func_init(pxcfunc_.get(), XC_GGA_C_PBE, XC_UNPOLARIZED);
//}
//
//CorrPbe::~CorrPbe(void)
//{
//}
//
//double CorrPbe::V(double r) const
//{
//	return my_xc_gga_vxc(r);
//}
//
//double CorrPbe::E(double r) const
//{
//	return my_xc_gga_exc(r);
//}
