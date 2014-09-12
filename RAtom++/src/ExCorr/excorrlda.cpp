#include "excorrlda.h"

namespace excorr {
    ExCorrLDA::ExCorrLDA(std::function<std::pair<double, double> (double)> rhoTilde,
                         std::uint32_t xc_type)
        : rhoTilde_(std::move(rhoTilde)),
          deleter_([](xc_func_type * xcfunc)
                   { xc_func_end(xcfunc);
                     delete xcfunc; }),
          pxcfunc_(new xc_func_type, deleter_)
    {
        xc_func_init(pxcfunc_.get(), xc_type, XC_POLARIZED);
    }

    ExCorrLDA::~ExCorrLDA()
    {
    }

    /*! rでの交換相関エネルギー密度（LDA）を返す
    @param[in]      r    原点からの距離（極座標）
    @return         rでの交換相関エネルギー密度（LDA）
    @exception      none
    */
    double ExCorrLDA::xc_exc(double r) const
    {
        const std::array<double, 2> rho = { rhoTilde_(r).first, rhoTilde_(r).second };
        std::array<double, 1> zk;
        xc_lda_exc(pxcfunc_.get(), 1, rho.data(), zk.data());

        return zk[0];
    }

    /*! rでの交換相関ポテンシャル（LDA）を返す
    @param[in]      r    原点からの距離（極座標）
    @return         rでの交換相関ポテンシャル（LDA）
    @exception      none
    */
    std::array<double, 2> ExCorrLDA::xc_vxc_impl(double r) const
    {
        const std::array<double, 2> rho = { rhoTilde_(r).first, rhoTilde_(r).second };
        std::array<double, 2> zk;
        xc_lda_vxc(pxcfunc_.get(), 1, rho.data(), zk.data());

        return std::move(zk);
    }

    /*! rでのαスピンに対する交換相関ポテンシャル（LDA）を返す
    @param[in]      r    原点からの距離（極座標）
    @return         rでのαスピンに対する交換相関ポテンシャル（LDA）
    @exception      none
    */
    template <util::IsSpin Spin, alpha_enabler<Spin> T>
    double ExCorrLDA::xc_vxc(double r) const
    {
        return xc_vxc_impl(r)[0];
    }

    /*! rでのβスピンに対する交換相関ポテンシャル（LDA）を返す
    @param[in]      r    原点からの距離（極座標）
    @return         rでのβスピンに対する交換相関ポテンシャル（LDA）
    @exception      none
    */
    template <util::IsSpin Spin, beta_enabler<Spin> T>
    double ExCorrLDA::xc_vxc(double r) const
    {
        return xc_vxc_impl(r)[1];
    }

    /*! 交換相関汎関数の名前を返す
    @return         交換相関汎関数の名前
    @exception      none
    */
    char const * ExCorrLDA::Name() const
    {
        return pxcfunc_->info->name;
    }

    template double ExCorrLDA::xc_vxc<util::IsSpin::Alpha>(double r) const;
    template double ExCorrLDA::xc_vxc<util::IsSpin::Beta>(double r) const;
}
