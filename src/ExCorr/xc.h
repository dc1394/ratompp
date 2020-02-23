/*! \file xc.h
    \brief Represents any exchange-coreelation potential.

    Copyright ©  2016 Zbigniew Romanowski [ROMZ@wp.pl] and @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#ifndef __RATOM_XC_H__
#define __RATOM_XC_H__

#include "excorr.h"
#include "../Util/spin.h"
#include <memory>       // for std::shared_ptr    
#include <string>       // for std::string

namespace excorr {
    //! A template class.
    /*!
        Represents any exchange-coreelation potential
    */
    template <util::Spin S>
    class Xc final {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param obj Excorrオブジェクトへのスマートポインタ
        */
        Xc(std::shared_ptr<ExCorr> const & obj)
            :   obj_(obj)
        {}

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~Xc() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function (const).
        /*!
            rでの交換相関エネルギーの値を返す
            \param r 原点からの距離（極座標）
            \return rでの交換相関エネルギー
        */
        double E(double r) const
        {
            return obj_->xc_exc(r);
        }

        //! A public member function (const).
        /*!
            Returns difference betwee energy density and potential for radius "r"
            rでの交換相関エネルギーと交換相関ポテンシャルの差の値を返す
            \param r 原点からの距離（極座標）
            \return rでの交換相関エネルギーと交換相関ポテンシャルの差の値
        */
        double EdiffV(double r) const
        {
            return obj_->xc_exc(r) - obj_->template xc_vxc<S>(r);
        }

        //! A public member function (const).
        /*!
            交換相関汎関数の名前を返す
            \return 交換相関汎関数の名前
        */
        std::string Name() const
        {
            return obj_->name();
        }

        //! A public member function (const).
        /*!
            rでの交換相関ポテンシャルの値を返す
            \param r 原点からの距離（極座標）
            \return rでの交換相関ポテンシャル
        */
        double V(double r) const
        {
            return obj_->template xc_vxc<S>(r);
        }

        // #endregion publicメンバ関数

        // #region メンバ変数

    private:
        //! A private member variable (const).
        /*!
            ExCorrオブジェクトへのスマートポインタ
        */
        std::shared_ptr<ExCorr> const obj_;
        
        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

    public:
        //! A default constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        Xc() = delete;

        //! A copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト（未使用）
        */
        Xc(Xc const & dummy) = delete;

        //! A public member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param dummy コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Xc & operator=(Xc const & dummy) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif

