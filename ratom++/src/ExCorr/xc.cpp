#include "stdafx.h"
#include "xc.h"

// use Debug only
//#include "XC_PBE.h"

namespace excorr {
    /*! rでの交換相関ポテンシャルを返す
    @param[in]      r    原点からの距離（極座標）
    @return         rでの交換相関ポテンシャル
    @exception      none
    */
    template <util::Spin Spin>
    double Xc<Spin>::V(double r) const
    {
        return V_(r);
    }

    /*! rでの交換相関エネルギー密度を返す
    @param[in]      r    原点からの距離（極座標）
    @return         rでの交換相関エネルギー密度
    @exception      none
    */
    template <util::Spin Spin>
    double Xc<Spin>::E(double r) const
    {
        return E_(r);
    }

    //
    // Returns difference betwee energy density and potential
    // This function always works. Somtimes it is possible to give
    // more efficient version of this function
    //
    //double Xc::EdiffV(double rho, double gRho) const
    //{
    //    return E(rho, gRho) - V(rho, gRho);
    //}

    /*! Returns difference betwee energy density and potential for radius "r"
    @param[in]      r    原点からの距離（極座標）
    @return         difference betwee energy density and potential for radius "r"
    @exception      none
    */
    template <util::Spin Spin>
    double Xc<Spin>::EdiffV(double r) const
    {
        return E_(r) - V_(r);
    }

    //********************************************************
    // GGA exchange-correlation energy density.
    // This routine dose not take into account spin-poralized
    // states.
    //
    // March 8th, 2014	Released by dc1394
    //
    //*********************************************************/
    //double Xc::my_xc_gga_exc(double r) const
    //{
    //    const std::array<double, 1> rho = { rhoTilde_(r) };
    //    const double rhotid = derivRhoTilde_(r);
    //    const std::array<double, 1> sigma = { rhotid * rhotid };
    //    std::array<double, 1> zk;

    //    xc_gga_exc(pxcfunc_.get(), 1, rho.data(), sigma.data(), zk.data());

    //    return zk[0];
    //}

    ////********************************************************
    //// GGA exchange-correlation potential.
    //// This routine dose not take into account spin-poralized
    //// states.
    ////
    //// March 8th, 2014	Released by dc1394
    ////
    ////*********************************************************/
    //double Xc::my_xc_gga_vxc(double r) const
    //{
    //    const std::array<double, 1> rho = { rhoTilde_(r) };
    //    const double rhotid = derivRhoTilde_(r);
    //    const std::array<double, 1> sigma = { rhotid * rhotid };
    //    std::array<double, 1> vrho, deds;

    //    xc_gga_vxc(pxcfunc_.get(), 1, rho.data(), sigma.data(), vrho.data(), deds.data());

    //    std::array<double, 1> d2edn2, d2ednds, d2eds2;
    //    xc_gga_fxc(pxcfunc_.get(), 1, rho.data(), sigma.data(), d2edn2.data(), d2ednds.data(), d2eds2.data());

    //    const double dedgrad = 2.0 * deds[0] * rhotid;
    //    const double d2edrhodgrad = 2.0 * rhotid * d2ednds[0];
    //    const double d2edgrad2 = 2.0 * deds[0] + 4.0 * sigma[0] * d2eds2[0];

    //    const double pot = -2.0 / r * dedgrad -
    //        d2edgrad2 * (rhoTildeLapl_(r) - 2.0 / r * rhotid) -
    //        d2edrhodgrad * rhotid;

    //    return vrho[0] + pot;
    //}

    /*! 交換相関汎関数の名前を返す
    @return         交換相関汎関数の名前
    @exception      none
    */
    template <util::Spin Spin>
    char const * Xc<Spin>::Name() const
    {
        return Name_();
    }

    template class Xc<util::Spin::Alpha>;
    template class Xc<util::Spin::Beta>;
}
