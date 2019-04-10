/*! \file excorr.h
    \brief 交換相関ポテンシャルクラスの基底クラスの宣言

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#ifndef __RATOM_EXCORR_H__
#define __RATOM_EXCORR_H__

#include "../Util/spin.h"
#include <string>   // for std::string
#include <utility>  // for std::pair

namespace excorr {
    //! A class.
    /*!
        交換相関ポテンシャルクラスの基底クラス
    */
    class ExCorr {
        // #region 型エイリアス

    protected:
        // A typedef.
        using dpair = std::pair<double, double>;

        // #endregion 型エイリアス

        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            デフォルトコンストラクタ
        */
        ExCorr() = default;

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        virtual ~ExCorr() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public pure virtual member function (const). 
        /*!
            交換相関汎関数の名前を返す
            \return 交換相関汎関数の名前
        */
        virtual std::string name() const = 0;

        //! A public pure virtual member function (const).
        /*!
            rでの交換相関エネルギー密度を返す関数
            \param r 原点からの距離（極座標）
            \return rでの交換エネルギー密度
        */
        virtual double xc_exc(double r) const = 0;
        
        //! A public member function (const).
        /*!
            rでのα or βスピンに対する交換相関ポテンシャルを返す関数
            \param r 原点からの距離（極座標）
            \return rでのα or βスピンに対する交換ポテンシャル
        */
        template <util::Spin S> double xc_vxc(double r) const;

        // #endregion publicメンバ関数

        // #region privateメンバ関数

    private:
        //! A private pure virtual member function (const).
        /*! rでの交換相関ポテンシャルを返す
            \param r 原点からの距離（極座標）
            \return rでの交換相関ポテンシャルのαスピンとβスピンのstd::pair
        */
        virtual dpair xc_vxc_impl(double r) const = 0;

        // #endregion privateメンバ関数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        ExCorr & operator=(ExCorr const &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif  // __RATOM_EXCORR_H__
