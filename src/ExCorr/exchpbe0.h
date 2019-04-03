/*! \file exchpbe0.h
    \brief Represents PBE0 Exchange potential.

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/
#ifndef __RATOM_EXCHPBE0_H__
#define __RATOM_EXCHPBE0_H__

#include "exchhf.h"
#include "excorrgga.h"

namespace excorr {
    class ExchPbe0 : public ExchHf, public ExCorrGGA {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param rhoTilde 電子密度
            \param rhoTildeDeriv 電子密度の勾配
            \param rhoTildeLapl 電子密度のラプラシアン
            \param Vh Hartreeポテンシャル
            \param xc_type 汎関数の種類
            \param Z 原子番号
        */
        ExchPbe0(std::function<dpair(double)> && rhoTilde,
                 std::function<dpair(double)> && rhoTildeDeriv,
                 std::function<dpair(double)> && rhoTildeLapl,
                 std::function<double(double)> const & Vh,
                 std::uint32_t xc_type,
                 double Z);
        
        //! A copy constructor.
        /*!
            デフォルトコピーコンストラクタ
        */
        ExchPbe0(ExchPbe0 const &) = default;

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~ExchPbe0() = default;
        
        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function (const). 
        /*!
            交換相関汎関数の名前を返す
            \return 交換相関汎関数の名前
        */
        std::string name() const override
        {
            return "PBE0 Exchange";
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

    };

}

#endif	// __RATOM_EXCHPBE0_H__
