/*! \file excorrlda.cpp
    \brief Represents LDA Exchange-Correlation potential.

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#include "excorrlda.h"
#include "xcfunc_deleter.h"
#include <array>            // for std::array

namespace excorr {
    // #region コンストラクタ

    ExCorrLDA::ExCorrLDA(std::function<std::pair<double, double> (double)> && rhoTilde, std::uint32_t xc_type)
        :   ExCorr(),
            pxcfunc_(new xc_func_type, xcfunc_deleter),
            rhoTilde_(std::move(rhoTilde))
    {
        xc_func_init(pxcfunc_.get(), xc_type, XC_POLARIZED);
    }
    
    // #endregion コンストラクタ

    // #region publicメンバ関数

    double ExCorrLDA::xc_exc(double r) const
    {
        std::array<double, 2> const rho = { rhoTilde_(r).first, rhoTilde_(r).second };
        std::array<double, 1> zk;
        xc_lda_exc(pxcfunc_.get(), 1, rho.data(), zk.data());

        return zk[0];
    }

    // #endregion publicメンバ関数

    // #region privateメンバ関数
    
    ExCorr::dpair ExCorrLDA::xc_vxc_impl(double r) const
    {
        std::array<double, 2> const rho = { rhoTilde_(r).first, rhoTilde_(r).second };
        std::array<double, 2> zk;
        xc_lda_vxc(pxcfunc_.get(), 1, rho.data(), zk.data());

        return std::make_pair(zk[0], zk[1]);
    }

    // #endregion privateメンバ関数
}
