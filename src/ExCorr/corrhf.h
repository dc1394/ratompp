/*! \file corrhf.h
    \brief Represents Hartree-Fock Correlation potential.

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#ifndef __RATOM_CORRHF_H__
#define __RATOM_CORRHF_H__

#include "excorr.h"

namespace excorr {
    //! A class.
    /*!
        Represents Hartree-Fock Correlation potential
    */
    class CorrHf final : public virtual ExCorr {
        // #region コンストラクタ・デストラクタ

    public:
        //! A default constructor.
        /*!
            デフォルトコンストラクタ
        */
        CorrHf() = default;

        //! A destructor.
        /*!
        デフォルトデストラクタ
        */
        ~CorrHf() override = default;

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
            交換相関エネルギー密度を返す関数
            \return rでの交換エネルギー密度
        */
        double xc_exc(double) const override
        {
            return 0.0;
        }

        // #endregion publicメンバ関数

    private:
        // #region privateメンバ関数

        //!  private member function (const).
        /*! rでの交換相関ポテンシャル（Hartree-Fock）を返す
            \return rでの交換相関ポテンシャル（Hartree-Fock）のαスピンとβスピンのstd::pair
        */
        dpair xc_vxc_impl(double) const override
        {
            return std::make_pair(0.0, 0.0);
        }

        // #endregion privateメンバ関数

    public:
        // #region 禁止されたコンストラクタ・メンバ関数
        
        //! A copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト（未使用）
        */
        CorrHf(CorrHf const & dummy) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param dummy コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        CorrHf& operator=(CorrHf const & dummy) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif	// __RATOM_EXCHHF_H__
