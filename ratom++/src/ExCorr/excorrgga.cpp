/*! \file excorrgga.cpp
    \brief Represents GGA Exchange-Correlation potential.

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#include "excorrgga.h"
#include <array>        // for std::array

namespace excorr {
    // #region コンストラクタ

    ExCorrGGA::ExCorrGGA(std::function<std::pair<double, double>(double)> && rhoTilde, std::function<std::pair<double, double>(double)> && rhoTildeDeriv, std::function<std::pair<double, double>(double)> && rhoTildeLapl, std::uint32_t xc_type)
        :   pxcfunc_(new xc_func_type, xcfunc_deleter),
            rhoTilde_(std::move(rhoTilde)),
            rhoTildeDeriv_(std::move(rhoTildeDeriv)),
            rhoTildeLapl_(std::move(rhoTildeLapl))
    {
        xc_func_init(pxcfunc_.get(), xc_type, XC_POLARIZED);
    }
    
    // #endregion コンストラクタ

    // #region publicメンバ関数

    double ExCorrGGA::xc_exc(double r) const
    {
        std::array<double, 2> const rho = { rhoTilde_(r).first, rhoTilde_(r).second };

        auto const rhotdaa = rhoTildeDeriv_(r).first * rhoTildeDeriv_(r).first;
        auto const rhotdab = rhoTildeDeriv_(r).first * rhoTildeDeriv_(r).second;
        auto const rhotdbb = rhoTildeDeriv_(r).second * rhoTildeDeriv_(r).second;

        std::array<double, 3> const sigma = { rhotdaa, rhotdab, rhotdbb };
        std::array<double, 1> zk;
        
        xc_gga_exc(pxcfunc_.get(), 1, rho.data(), sigma.data(), zk.data());

        return zk[0];
    }

    // #endregion publicメンバ関数

    // #region privateメンバ関数
    
    std::pair<double, double> ExCorrGGA::xc_vxc_impl(double r) const
    {
        std::array<double, 2> const rho = { rhoTilde_(r).first, rhoTilde_(r).second };
        std::array<double, 2> zk;
        xc_lda_vxc(pxcfunc_.get(), 1, rho.data(), zk.data());

        return std::make_pair(zk[0], zk[1]);
    }

    // #endregion privateメンバ関数

    // #region templateメンバ関数の実体化

    template <>
    double ExCorrGGA::xc_vxc<util::Spin::Alpha>(double r) const
    {
        return xc_vxc_impl(r).first;
    }

    template <>
    double ExCorrGGA::xc_vxc<util::Spin::Beta>(double r) const
    {
        return xc_vxc_impl(r).second;
    }

    // #endregion templateメンバ関数の実体化
}
