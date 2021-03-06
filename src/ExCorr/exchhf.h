﻿/*! \file exchhf.h
    \brief Represents Hartree-Fock Exchange potential.

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#ifndef __RATOM_EXCHHF_H__
#define __RATOM_EXCHHF_H__

#include "excorr.h"
#include "../Util/spin.h"
#include <functional>       // for std::function
#include <string>           // for std::string
#include <utility>          // for std::pair

namespace excorr {
    //! A class.
    /*!
        Represents Hartree-Fock Exchange potential
    */
    class ExchHf : public virtual ExCorr {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param rhoTilde 電子密度
            \param rhoTildeDeriv 電子密度の勾配
            \param rhoTildeLapl 電子密度のラプラシアン
            \param xc_type 汎関数の種類
        */
        ExchHf(std::function<double(double)> const & Vh, double Z);

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~ExchHf() override = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function (const). 
        /*!
            交換相関汎関数の名前を返す
            \return 交換相関汎関数の名前
        */
        virtual std::string name() const override
        {
            return "hartree-fock Exchange";
        }

        //! A public member function (const).
        /*!
            rでの交換相関エネルギー密度（Hartree-Fock）を返す関数
            \param r 原点からの距離（極座標）
            \return rでの交換エネルギー密度（Hartree-Fock）
        */
        double xc_exc(double r) const override;

        // #endregion publicメンバ関数

        // #region privateメンバ関数

    private:
        //! A private member function (const).
        /*! rでの交換相関ポテンシャル（Hartree-Fock）を返す
            \param r 原点からの距離（極座標）
            \return rでの交換相関ポテンシャル（Hartree-Fock）のαスピンとβスピンのstd::pair
        */
        dpair xc_vxc_impl(double r) const override;

        // #endregion privateメンバ関数

    protected:
        // #region protectedメンバ変数

        //! A protected member variable (constant expression).
        /*!
            ゼロ判定用の定数
        */
        static auto constexpr ZERO = 1.0E-12;

        //! A protected member variable (constant).
        /*!
            Hartreeポテンシャル
        */
        std::function<double(double)> const Vh_;

        //! A protected member variable (constant).
        /*!
            原子番号
        */
        double const Z_;

        // #endregion protectedメンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

    public:
        //! A default constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        ExchHf() = delete;

        //! A copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト（未使用）
        */
        ExchHf(ExchHf const & dummy) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        ExchHf & operator=(ExchHf const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif	// __RATOM_EXCHHF_H__
