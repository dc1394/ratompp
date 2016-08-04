/*! \file excorrlda.h
    \brief Represents any LDA Exchange-Correlation potential.

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#ifndef __RATOM_EXCORR_LDA_H__
#define __RATOM_EXCORR_LDA_H__

#include "../Util/spin.h"
#include <cstdint>                      // for std::uint32_t
#include <functional>                   // for std::function
#include <memory>                       // for std::unique_ptr
#include <string>                       // for std::string
#include <utility>                      // for std::pair
#include <xc.h>                         // for xc_func_end
#include <boost/checked_delete.hpp>     // for boost::checked_delete

namespace excorr {
    //! A lambda expression.
    /*!
        xc_func_typeへのポインタを解放するラムダ式
        \param xcfunc xc_func_type へのポインタ
    */
    static auto const xcfunc_deleter = [](xc_func_type * xcfunc) {
        xc_func_end(xcfunc);
        boost::checked_delete(xcfunc);
    };

    //! A class.
    /*!
        Represents any LDA Exchange-Correlation potential
    */
    class ExCorrLDA final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param rhoTilde 電子密度
        */
        ExCorrLDA(std::function<std::pair<double, double>(double)> && rhoTilde, std::uint32_t xc_type);

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~ExCorrLDA() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function (const). 
        /*!
            交換相関汎関数の名前を返す
            \return 交換相関汎関数の名前
        */
        std::string Name() const;

        //! A public member function (const).
        /*!
            rでの交換相関エネルギー密度（LDA）を返す関数
            \param r 原点からの距離（極座標）
            \return rでの交換相関エネルギー密度（LDA）
        */
        double xc_exc(double r) const;
        
        //! A public member function (const).
        /*!
            rでのα or βスピンに対する交換相関ポテンシャル（LDA）を返す関数
            \param r 原点からの距離（極座標）
            \return rでのα or βスピンに対する交換相関ポテンシャル（LDA）
        */
        template <util::Spin S> double xc_vxc(double r) const;

        // #endregion publicメンバ関数

        // #region privateメンバ関数

        //!  private member function (const).
        /*! rでの交換相関ポテンシャル（LDA）を返す
            \param r 原点からの距離（極座標）
            \return rでの交換相関ポテンシャル（LDA）のαスピンとβスピンのstd::pair
        */
        std::pair<double, double> ExCorrLDA::xc_vxc_impl(double r) const;

        // #endregion privateメンバ関数

        // #region メンバ変数

    private:

        //! A private member variable (constant).
        /*!
            電子密度
        */
        std::function<std::pair<double, double> (double)> const rhoTilde_;
        
        //! A private member variable (constant).
        /*!
            ratomへのスマートポインタ
        */
        std::unique_ptr<xc_func_type, decltype(xcfunc_deleter)> const pxcfunc_;

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        ExCorrLDA() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
        */
        ExCorrLDA(ExCorrLDA const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        ExCorrLDA & operator=(ExCorrLDA const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };

    // #region templateメンバ関数の実体化

    template <>
    double ExCorrLDA::xc_vxc<util::Spin::Alpha>(double r) const
    {
        return xc_vxc_impl(r).first;
    }

    template <>
    double ExCorrLDA::xc_vxc<util::Spin::Beta>(double r) const
    {
        return xc_vxc_impl(r).second;
    }
    
    // #endregion templateメンバ関数の実体化
}

#endif  // __RATOM_EXCORR_LDA_H__
