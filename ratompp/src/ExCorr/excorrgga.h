/*! \file excorrgga.h
    \brief Represents GGA Exchange Correlation potential.

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#ifndef __RATOM_EXCORRGGA_H__
#define __RATOM_EXCORRGGA_H__

#include "../Util/spin.h"
#include "xcfunc_deleter.h"
#include <cstdint>          // for std::uint32_t
#include <functional>       // for std::function
#include <memory>           // for std::unique_ptr
#include <string>           // for std::string
#include <utility>          // for std::pair

namespace excorr {
    //! A class.
    /*!
        Represents Slater Exchange potential
    */
    class ExCorrGGA final {
        // #region 型エイリアス

        // A typedef.
        using dpair = std::pair<double, double>;

        // #endregion 型エイリアス

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
        ExCorrGGA(std::function<dpair(double)> && rhoTilde, std::function<dpair(double)> && rhoTildeDeriv, std::function<dpair(double)> && rhoTildeLapl, std::uint32_t xc_type);

        //! A copy constructor.
        /*!
            デフォルトコピーコンストラクタ
        */
        ExCorrGGA(ExCorrGGA const &) = default;

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~ExCorrGGA() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function (const). 
        /*!
            交換相関汎関数の名前を返す
            \return 交換相関汎関数の名前
        */
        std::string ExCorrGGA::Name() const
        {
            return std::string(pxcfunc_->info->name);
        }

        //! A public member function (const).
        /*!
            rでの交換相関エネルギー密度（Slater）を返す関数
            \param r 原点からの距離（極座標）
            \return rでの交換エネルギー密度（Slater）
        */
        double xc_exc(double r) const;

        //! A public member function (const).
        /*!
            rでのα or βスピンに対する交換相関ポテンシャル（Slater）を返す関数
            \param r 原点からの距離（極座標）
            \return rでのα or βスピンに対する交換ポテンシャル（Slater）
        */
        template <util::Spin S> double xc_vxc(double r) const;

        // #endregion publicメンバ関数

        // #region privateメンバ関数

        //!  private member function (const).
        /*! rでの交換相関ポテンシャル（Slater）を返す
            \param r 原点からの距離（極座標）
            \return rでの交換相関ポテンシャル（Slater）のαスピンとβスピンのstd::pair
        */
        dpair ExCorrGGA::xc_vxc_impl(double r) const;

        // #endregion privateメンバ関数

        // #region メンバ変数

    private:
        
        //! A private member variable (constant).
        /*!
            ratomへのスマートポインタ
        */
        std::shared_ptr<xc_func_type> const pxcfunc_;

        //! A private member variable (constant).
        /*!
            電子密度
        */
        std::function<dpair(double)> const rhoTilde_;

        //! A private member variable (constant).
        /*!
            電子密度の勾配
        */
        std::function<dpair(double)> const rhoTildeDeriv_;

        //! A private member variable (constant).
        /*!
            電子密度のラプラシアン
        */
        std::function<dpair(double)> const rhoTildeLapl_;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        ExCorrGGA() = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        ExCorrGGA & operator=(ExCorrGGA const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };        
}

#endif  // __RATOM_EXCORRGGA_H__
