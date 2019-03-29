/*! \file exchpbe0.h
    \brief Represents PBE0 Exchange potential.

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/


#include "exchpbe0.h"
#include <cmath>        // for std::fabs

namespace excorr {
    // #region コンストラクタ

    ExchPbe0::ExchPbe0(std::function<dpair(double)> && rhoTilde,
                       std::function<dpair(double)> && rhoTildeDeriv,
                       std::function<dpair(double)> && rhoTildeLapl,
                       std::function<double(double)> const & Vh,
                       std::uint32_t xc_type,
                       double Z)
        :   ExchHf(Vh, Z),
            ExCorrGGA(std::move(rhoTilde), std::move(rhoTildeDeriv), std::move(rhoTildeLapl), xc_type)
    {
    }

    // #endregion コンストラクタ

    // #region publicメンバ関数

    double ExchPbe0::xc_exc(double r) const
    {
        if (std::fabs(Z_ - 1.0) < ZERO) {
            return 0.75 * ExCorrGGA::xc_exc(r) - 0.125 * Vh_(r);
        }
        else if (std::fabs(Z_ - 2.0) < ZERO) {
            return 0.75 * ExCorrGGA::xc_exc(r) - 0.0625 * Vh_(r);
        }
        else {
            return NAN;
        }
    }

    // #endregion publicメンバ関数

    // #region privateメンバ関数

    ExCorr::dpair ExchPbe0::xc_vxc_impl(double r) const
    {
        if (std::fabs(Z_ - 1.0) < ZERO) {
            return std::make_pair(0.75 * ExCorrGGA::xc_vxc_impl(r).first - 0.25 * Vh_(r), 0.75 * ExCorrGGA::xc_vxc_impl(r).second);
        }
        else if (std::fabs(Z_ - 2.0) < ZERO) {
            return std::make_pair(0.75 * ExCorrGGA::xc_vxc_impl(r).first - 0.125 * Vh_(r), 0.75 * ExCorrGGA::xc_vxc_impl(r).second - 0.125 * Vh_(r));
        }
        else {
            return std::make_pair(NAN, NAN);
        }
    }

    // #endregion privateメンバ関数
}
