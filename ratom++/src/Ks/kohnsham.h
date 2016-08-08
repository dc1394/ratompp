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
#include <boost/mpl/int.hpp>    // for boost::mpl::int_

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
            \param m_stateSet 
        */
        KohnSham(std::shared_ptr<ParamDb> const & db, std::shared_ptr<StateSet> const & stateSet) : m_db(db), m_stateSet(stateSet) {}

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~KohnSham() = default;

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
        double Get(double r) const;

        //! A public member function (const).
        /*!
            Returns value of electron density for radius "r" for alpha spin
            \param r 動径方向の値r
            \param オーバーロード用のダミー引数
        */
        double Get(double r, boost::mpl::int_<static_cast<std::int32_t>(util::Spin::Alpha)>) const;
        
        //! A public member function (const).
        /*!
            Returns value of electron density for radius "r" for beta spin
            \param r 動径方向の値r
            \param オーバーロード用のダミー引数
        */
        double Get(double r, boost::mpl::int_<static_cast<std::int32_t>(util::Spin::Beta)>) const;

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

        // #region メンバ変数

    private:
        // Interaction potential used in radial Kohn-Sham equation
        class PotRad : public util::Fun1D
        {
            std::shared_ptr<util::Fun1D> m_pot;

        public:
            PotRad(void) : m_l(0) {}
            PotRad(std::shared_ptr<util::Fun1D> const & pot) : m_pot(pot), m_l(0) {}
            virtual ~PotRad() = default;

            virtual double Get(double r) const
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

        //! A private member variable.
        /*!
            One solver for each angular quantum number L
        */
        std::vector<EigProb> m_eigProb;

        //! A private member variable.
        /*!
            Effective interaction potential used in radial Kohn-Sham equation
        */
        PotRad m_radPot;

        //! A private member variable.
        /*!
            Number of calculated eigenvalues for each angular quantum number L
        */
        std::vector<size_t> m_eigNo;

        //! A private member variable.
        /*!
            database of parameters
        */
        std::shared_ptr<ParamDb> const m_db;

        //! A private member variable.
        /*!
            Set of eigenstates
        */
        std::shared_ptr<StateSet> m_stateSet;
    };
}

#endif
