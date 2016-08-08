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

    std::pair<double, double> ExCorrLDA::xc_exc_impl(double r) const
    {
        std::array<double, 2> const rho = { rhoTilde_(r).first, rhoTilde_(r).second };
        std::array<double, 2> zk;
        xc_lda_exc(pxcfunc_.get(), 1, rho.data(), zk.data());

        return std::make_pair(zk[0], zk[1]);
    }

    std::pair<double, double> ExCorrLDA::xc_vxc_impl(double r) const
    {
        std::array<double, 2> const rho = { rhoTilde_(r).first, rhoTilde_(r).second };
        std::array<double, 2> zk;
        xc_lda_vxc(pxcfunc_.get(), 1, rho.data(), zk.data());

        return std::make_pair(zk[0], zk[1]);
    }

    // #region templateメンバ関数の実体化

    template <>
    double ExCorrLDA::xc_exc<util::Spin::Alpha>(double r) const
    {
        return xc_exc_impl(r).first;
    }

    template <>
    double ExCorrLDA::xc_exc<util::Spin::Beta>(double r) const
    {
        return xc_exc_impl(r).second;
    }

    template <>
    double ExCorrLDA::xc_vxc<util::Spin::Alpha>(double r) const
    {
        return xc_vxc_impl(r).first;
    }

    template <>
    double ExCorrLDA::xc_vxc<util::Spin::Beta>(double r) const
    {
        return xc_vxc_impl(r).second;
    }

    // #endregion templateメンバ関数の実体化
}
