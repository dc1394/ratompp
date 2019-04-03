﻿/*! \file corrhf.h
    \brief Represents Hartree-Fock Correlation potential.

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#ifndef __RATOM_CORRHF_H__
#define __RATOM_CORRHF_H__

#include "excorr.h"
#include "../Util/spin.h"

namespace excorr {
    //! A class.
    /*!
        Represents Hartree-Fock Correlation potential
    */
    class CorrHf final : public ExCorr {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            デフォルトコンストラクタ
        */
        CorrHf() = default;

        //! A copy constructor.
        /*!
        デフォルトコピーコンストラクタ
        */
        CorrHf(CorrHf const &) = default;

        //! A destructor.
        /*!
        デフォルトデストラクタ
        */
        ~CorrHf() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function (const). 
        /*!
        交換相関汎関数の名前を返す
        \return 交換相関汎関数の名前
        */
        std::string name() const override
        {
            return "dummy";
        }

        //! A public member function (const).
        /*!
            rでの交換相関エネルギー密度を返す関数
            \param r 原点からの距離（極座標）
            \return rでの交換エネルギー密度
        */
        double xc_exc(double r) const override
        {
            return 0.0;
        }

        // #endregion publicメンバ関数

        // #region privateメンバ関数

        //!  private member function (const).
        /*! rでの交換相関ポテンシャル（Hartree-Fock）を返す
        \param r 原点からの距離（極座標）
        \return rでの交換相関ポテンシャル（Hartree-Fock）のαスピンとβスピンのstd::pair
        */
        dpair xc_vxc_impl(double r) const override
        {
            return std::make_pair(0.0, 0.0);
        }

        // #endregion privateメンバ関数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        CorrHf& operator=(CorrHf const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif	// __RATOM_EXCHHF_H__