/*! \file excorrgga.cpp
    \brief Represents GGA Exchange-Correlation potential.

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#include "excorrgga.h"
#include "xcfunc_deleter.h"
#include <array>            // for std::array

namespace excorr {
    // #region コンストラクタ

    ExCorrGGA::ExCorrGGA(std::function<dpair(double)> && rhoTilde, std::function<dpair(double)> && rhoTildeDeriv, std::function<dpair(double)> && rhoTildeLapl, std::uint32_t xc_type)
        :   ExCorr(),
            pxcfunc_(new xc_func_type, xcfunc_deleter),
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

        auto const rhodaa = rhoTildeDeriv_(r).first * rhoTildeDeriv_(r).first;
        auto const rhodab = rhoTildeDeriv_(r).first * rhoTildeDeriv_(r).second;
        auto const rhodbb = rhoTildeDeriv_(r).second * rhoTildeDeriv_(r).second;

        std::array<double, 3> const sigma = { rhodaa, rhodab, rhodbb };
        std::array<double, 1> zk;
        
        xc_gga_exc(pxcfunc_.get(), 1, rho.data(), sigma.data(), zk.data());

        return zk[0];
    }

    // #endregion publicメンバ関数

    // #region privateメンバ関数
    
    ExCorr::dpair ExCorrGGA::xc_vxc_impl(double r) const
    {
        std::array<double, 2> const rho = { rhoTilde_(r).first, rhoTilde_(r).second };
        std::array<double, 2> const rhod = { rhoTildeDeriv_(r).first, rhoTildeDeriv_(r).second };

        auto const rhodaa = rhod[0] * rhod[0];
        auto const rhodab = rhod[0] * rhod[1];
        auto const rhodbb = rhod[1] * rhod[1];

        std::array<double, 3> const sigma = { rhodaa, rhodab, rhodbb };

        std::array<double, 2> vrho;
        std::array<double, 3> deds;

        xc_gga_vxc(pxcfunc_.get(), 1, rho.data(), sigma.data(), vrho.data(), deds.data());

        std::array<double, 3> dedgrad = { 2.0 * deds[0] * rhod[0] + deds[1] * rhod[1], 2.0 * deds[2] * rhod[1] + deds[1] * rhod[0] };
        std::array<double, 3> d2edn2;
        std::array<double, 6> d2ednds, d2eds2;
        xc_gga_fxc(pxcfunc_.get(), 1, rho.data(), sigma.data(), d2edn2.data(), d2ednds.data(), d2eds2.data());

        std::array<double, 4> d2edrhodgrad;
        std::array<double, 3> d2edgrad2;

        d2edrhodgrad[0] = 2.0 * rhod[0] * d2ednds[0] + rhod[1] * d2ednds[1];
        d2edrhodgrad[1] = 2.0 * rhod[0] * d2ednds[3] + rhod[1] * d2ednds[4];
        d2edrhodgrad[2] = 2.0 * rhod[1] * d2ednds[2] + rhod[0] * d2ednds[1];
        d2edrhodgrad[3] = 2.0 * rhod[1] * d2ednds[5] + rhod[0] * d2ednds[4];

        d2edgrad2[0] = 2.0 * deds[0] + 4.0 * sigma[0] * d2eds2[0] + 4.0 * sigma[1] * d2eds2[1] + sigma[2] * d2eds2[3];
        d2edgrad2[1] = deds[1] + 4.0 * sigma[1] * d2eds2[2] + 2.0 * sigma[0] * d2eds2[1] + 2.0 * sigma[2] * d2eds2[4] + sigma[1] * d2eds2[3];
        d2edgrad2[2] = 2.0 * deds[2] + 4.0 * sigma[2] * d2eds2[5] + 4.0 * sigma[1] * d2eds2[4] + sigma[0] * d2eds2[3];

        vrho[0] += - 2.0 / r * dedgrad[0]
                   - d2edgrad2[0] * (rhoTildeLapl_(r).first - 2.0 / r * rhod[0])
                   - d2edgrad2[1] * (rhoTildeLapl_(r).second - 2.0 / r * rhod[1])
                   - d2edrhodgrad[0] * rhod[0]
                   - d2edrhodgrad[1] * rhod[1];

        vrho[1] += - 2.0 / r * dedgrad[1]
                   - d2edgrad2[2] * (rhoTildeLapl_(r).second - 2.0 / r * rhod[1])
                   - d2edgrad2[1] * (rhoTildeLapl_(r).first - 2.0 / r * rhod[0])
                   - d2edrhodgrad[2] * rhod[0]
                   - d2edrhodgrad[3] * rhod[1];

        return std::make_pair(vrho[0], vrho[1]);
    }
}
