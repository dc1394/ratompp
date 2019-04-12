/*! \file kohnsham.h
    \brief Linear Kohn-Sham equation
    Copyright ©  2016 Zbigniew Romanowski [ROMZ@wp.pl] and @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#ifndef __RATOM_KOHNSHAM_H__
#define __RATOM_KOHNSHAM_H__

#include "paramdb.h"
#include "stateset.h"
#include "../Util/spin.h"

namespace ks {
    //! A template class.
    /*!
        Linear Kohn-Sham equation
    */
    template <util::Spin S>
    class KohnSham final : public util::Fun1D
    {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param db パラメータへのスマートポインタ
            \param stateSet StateSetオブジェクトへのスマートポインタ
        */
        KohnSham(std::shared_ptr<const ParamDb> && db, std::shared_ptr<StateSet> const & stateSet);

        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param db パラメータへのスマートポインタ
            \param stateSet StateSetオブジェクトへのスマートポインタ
        */
        KohnSham(std::shared_ptr<const ParamDb> && db, std::shared_ptr<StateSet> && stateSet);

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~KohnSham() override = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function.
        /*!
            Solver configurations.
            Required parameters are properly set
            \param pot ポテンシャルへのスマートポインタ
        */
        void Config(std::shared_ptr<util::Fun1D> const & pot);
        
        //! A public member function (const).
        /*!
            Returns value of electron density for radius "r"
            \param r 動径方向の値r
        */
        double Get(double r) const override;

        //! A public member function (const).
        /*!
            Solves linear eqigenvalue problem
        */
        void Solve();

        //! A public member function (const).
        /*!
            Writes eigenfunctions into file.
            Lobato coeffictienst are written as well.
        */
        void WriteEigen() const;

        // #ebdregion publicメンバ関数

        // #region インナークラス

    private:
        //! A inner class.
        /*!
            Interaction potential used in radial Kohn-Sham equation
        */
        class PotRad final : public util::Fun1D
        {
            std::shared_ptr<util::Fun1D> m_pot;

        public:
            PotRad() : m_l(0) {}
            PotRad(std::shared_ptr<util::Fun1D> const & pot) : m_pot(pot), m_l(0) {}
            ~PotRad() override = default;

            double Get(double r) const override
            {
                assert(r > 0);
                return m_pot->Get(r) + m_l * (m_l + 1) / (2 * r * r);
            }

            std::shared_ptr<util::Fun1D> & getPot()
            {
                return m_pot;
            }

        public:
            // Angular quantum number
            size_t m_l;
        };

        // #endregion インナークラス

        // #region メンバ変数

        //! A private member variable.
        /*!
            database of parameters
        */
        std::shared_ptr<const ParamDb> const m_db;

        //! A private member variable.
        /*!
            One solver for each angular quantum number L
        */
        std::vector<EigProb> m_eigProb;

        //! A private member variable.
        /*!
            Effective interaction potential used in radial Kohn-Sham equation
        */
        std::shared_ptr<PotRad> m_radPot;

        //! A private member variable.
        /*!
            Number of calculated eigenvalues for each angular quantum number L
        */
        std::vector<size_t> m_eigNo;

        //! A private member variable.
        /*!
            Set of eigenstates
        */
        std::shared_ptr<StateSet> m_stateSet;

        // #endregion メンバ変数

        // #region 禁止されたコンストラクタ・メンバ関数

    public:
        //! A default constructor (deleted).
        /*!
            デフォルトコンストラクタ（禁止）
        */
        KohnSham() = delete;

        //! A default copy constructor (deleted).
        /*!
            コピーコンストラクタ（禁止）
            \param dummy コピー元のオブジェクト（未使用）
        */
        KohnSham(KohnSham const & dummy) = delete;

        //! A public member function (deleted).
        /*!
            operator=()の宣言（禁止）
            \param dummy コピー元のオブジェクト（未使用）
            \return コピー元のオブジェクト
        */
        KohnSham & operator=(KohnSham const & dummy) = delete;

        // #endregion 禁止されたコンストラクタ・メンバ関数
    };
}

#endif
