/*! \file excorrlda.h
    \brief Represents LDA Exchange-Correlation potential.

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#ifndef __RATOM_EXCORRLDA_H__
#define __RATOM_EXCORRLDA_H__

#include "excorr.h"
#include "../Util/spin.h"
#include "xcfunc_deleter.h"
#include <cstdint>          // for std::uint32_t
#include <functional>       // for std::function
#include <memory>           // for std::shared_ptr
#include <string>           // for std::string

namespace excorr {
    //! A class.
    /*!
        Represents LDA Exchange-Correlation potential.
    */
    class ExCorrLDA final : public ExCorr {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param rhoTilde 電子密度
            \param xc_type 汎関数の種類
        */
        ExCorrLDA(std::function<std::pair<double, double>(double)> && rhoTilde, std::uint32_t xc_type);

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~ExCorrLDA() override = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function (const). 
        /*!
            交換相関汎関数の名前を返す
            \return 交換相関汎関数の名前
        */
        std::string name() const override
        {
            return std::string(pxcfunc_->info->name);
        }

        //! A public member function (const).
        /*!
            rでの交換相関エネルギー密度（LDA）を返す関数
            \param r 原点からの距離（極座標）
            \return rでの交換相関エネルギー密度（LDA）
        */
        double xc_exc(double r) const override;

        // #endregion publicメンバ関数

        // #region privateメンバ関数

        //!  private member function (const).
        /*! 
            rでの交換相関ポテンシャル（LDA）を返す
            \param r 原点からの距離（極座標）
            \return rでの交換相関ポテンシャル（LDA）のαスピンとβスピンのstd::pair
        */
        dpair xc_vxc_impl(double r) const override;

        // #endregion privateメンバ関数

        // #region メンバ変数

    private:
        //! A private member variable (constant).
        /*!
            libxcへのスマートポインタ
        */
        std::unique_ptr<xc_func_type, decltype(xcfunc_deleter)> const pxcfunc_;

        //! A private member variable (constant).
        /*!
            電子密度を与える関数オブジェクト
        */
        std::function<std::pair<double, double> (double)> const rhoTilde_;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

    public:
        //! A default constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        ExCorrLDA() = delete;

        //! A copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト（未使用）
        */
        ExCorrLDA(ExCorrLDA const & dummy) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        ExCorrLDA & operator=(ExCorrLDA const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };        
}

#endif  // __RATOM_EXCORRLDA_H__
