/*! \file exchhf.h
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
    class ExchHf final : public ExCorr {
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

        //! A copy constructor.
        /*!
            デフォルトコピーコンストラクタ
        */
        ExchHf(ExchHf const &) = default;

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~ExchHf() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function (const). 
        /*!
            交換相関汎関数の名前を返す
            \return 交換相関汎関数の名前
        */
        std::string name() const override
        {
            return "hartree-fock";
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

        //!  private member function (const).
        /*! rでの交換相関ポテンシャル（Hartree-Fock）を返す
            \param r 原点からの距離（極座標）
            \return rでの交換相関ポテンシャル（Hartree-Fock）のαスピンとβスピンのstd::pair
        */
        dpair xc_vxc_impl(double r) const override;

        // #endregion privateメンバ関数

        // #region メンバ変数

    private:
        //! A private member variable (constant expression).
        /*!
            ゼロ判定用の定数
        */
        static constexpr double ZERO = 1.0E-12;

        //! A private member variable (constant).
        /*!
            Hartreeポテンシャル
        */
        std::function<double(double)> const Vh_;

        //! A private member variable (constant).
        /*!
            原子番号
        */
        double const Z_;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        ExchHf() = delete;

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
