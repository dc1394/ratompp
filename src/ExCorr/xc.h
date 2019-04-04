/*! \file xc.h
    \brief Represents any exchange-coreelation potential.

    Copyright ©  2016 Zbigniew Romanowski [ROMZ@wp.pl] and @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#ifndef __RATOM_XC_H__
#define __RATOM_XC_H__

#include "../Util/spin.h"
#include <functional>       // for std::function
#include <string>           // for std::string

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
            \param obj オブジェクト
        */
        template <typename T>
        Xc(const T & obj)
            :   E_([&obj](double r) { return obj.xc_exc(r); }),
                Name_([&obj]() { return obj.name(); }),
                V_([obj](double r) { return obj.template xc_vxc<S>(r); })
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
            return E_(r);
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
            return E_(r) - V_(r);
        }

        //! A public member function (const).
        /*!
            交換相関汎関数の名前を返す
            \return 交換相関汎関数の名前
        */
        std::string Name() const
        {
            return Name_();
        }

        //! A public member function (const).
        /*!
            rでの交換相関ポテンシャルの値を返す
            \param r 原点からの距離（極座標）
            \return rでの交換相関ポテンシャル
        */
        double V(double r) const
        {
            return V_(r);
        }

        // #endregion publicメンバ関数

        // #region メンバ変数

    private:
        //! A private member variable (const).
        /*!
            交換相関エネルギーの値を返す関数オブジェクト
        */
        std::function<double(double)> const E_;

        //! A private member variable (const).
        /*!
            交換相関汎関数の名前を返す関数オブジェクト
        */
        std::function<std::string()> const Name_;
        
        //! A private member variable (const).
        /*!
            交換相関ポテンシャルの値を返す関数オブジェクト
        */
        std::function<double(double)> const V_;
        
        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

        //! A private constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        Xc() = delete;

        //! A private copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param コピー元のオブジェクト（未使用）
        */
        Xc(Xc const &) = delete;

        //! A private member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        Xc & operator=(const Xc &) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif

