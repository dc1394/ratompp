/*! \file exchhf.cpp
    \brief Represents Hartree-Fock Exchange potential.

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#include "exchhf.h"
#include <cmath>

namespace excorr {
    // #region コンストラクタ

    ExchHf::ExchHf(std::function<double(double)> const & Vh, double Z)
        :   Vh_(Vh),
            Z_(Z)
    {
    }
    
    // #endregion コンストラクタ

    // #region publicメンバ関数
    
    double ExchHf::xc_exc(double r) const
    {
        if (std::fabs(Z_ - 1.0) < ZERO) {
            return - 0.5 * Vh_(r);
        }
        else if (std::fabs(Z_ - 2.0) < ZERO) {
            return -0.25 * Vh_(r);
        }
        else {
            return NAN;
        }
    }

    // #endregion publicメンバ関数

    // #region privateメンバ関数

    ExCorr::dpair ExchHf::xc_vxc_impl(double r) const
    {
        if (std::fabs(Z_ - 1.0) < ZERO) {
            return std::make_pair(- Vh_(r), 0.0);
        }
        
        if (std::fabs(Z_ - 2.0) < ZERO) {
            return std::make_pair(-0.5 * Vh_(r), -0.5 * Vh_(r));
        }
        
        return std::make_pair(NAN, NAN);
    }

    // #endregion privateメンバ関数
}