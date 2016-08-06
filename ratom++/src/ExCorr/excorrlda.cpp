/*! \file excorrlda.cpp
    \brief Represents any LDA Exchange-Correlation potential.

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#include "excorrlda.h"
#include <array>        // for std::array

namespace excorr {
    ExCorrLDA::ExCorrLDA(std::function<std::pair<double, double> (double)> && rhoTilde, std::uint32_t xc_type)
        : rhoTilde_(std::move(rhoTilde)),
          pxcfunc_(new xc_func_type, xcfunc_deleter)
    {
        xc_func_init(pxcfunc_.get(), xc_type, XC_POLARIZED);
    }
    
    std::string ExCorrLDA::Name() const
    {
        return std::string(pxcfunc_->info->name);
    }

    double ExCorrLDA::xc_exc(double r) const
    {
        const std::array<double, 2> rho = { rhoTilde_(r).first, rhoTilde_(r).second };
        std::array<double, 1> zk;
        xc_lda_exc(pxcfunc_.get(), 1, rho.data(), zk.data());

        return zk[0];
    }

    std::pair<double, double> ExCorrLDA::xc_vxc_impl(double r) const
    {
        std::array<double, 2> const rho = { rhoTilde_(r).first, rhoTilde_(r).second };
        std::array<double, 2> zk;
        xc_lda_vxc(pxcfunc_.get(), 1, rho.data(), zk.data());

        return std::make_pair(zk[0], zk[1]);
    }
}
