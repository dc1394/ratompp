#include "excorr.h"

namespace excorr {
    // #region templateメンバ関数の実体化

    template <>
    double ExCorr::xc_vxc<util::Spin::Alpha>(double r) const
    {
        return xc_vxc_impl(r).first;
    }

    template <>
    double ExCorr::xc_vxc<util::Spin::Beta>(double r) const
    {
        return xc_vxc_impl(r).second;
    }

    // #endregion templateメンバ関数の実体化
}

