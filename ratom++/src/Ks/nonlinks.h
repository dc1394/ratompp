/*! \file nonlinks.h
    \brief Nonlinear Kohn-Sham equation

    Copyright ©  2016 Zbigniew Romanowski [ROMZ@wp.pl] and @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#ifndef __RATOM_NONLINKS_H__
#define __RATOM_NONLINKS_H__

#include "energy.h"
#include "kohnsham.h"
#include "paramdb.h"
#include "pot.h"
#include "rho.h"
#include "../Util/spin.h"
#include <chrono>           // for std::chrono::high_resolution_clock::time_point

namespace ks {
    //! A class.
    /*!
        Nonlinear Kohn-Sham equation
    */
    class NonLinKs final
    {
        // #region コンストラクタ・デストラクタ

    public:
        //! A constructor.
        /*!
            唯一のコンストラクタ
            \param path パス
        */
        NonLinKs(char const * path);
        
        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~NonLinKs() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region publicメンバ関数

        //! A public member function.
        /*!
            Read in input parameters. Stars the calculation. Writes results to files.
            \return 関数が成功したかどうかの値 
        */
        int Run();

        // #endregion publicメンバ関数

        // #region privateメンバ関数

    private:
        
        //! A private member function (const).
        /*!
            Returns "true", if required accuarcy reached.
            \return 収束したかどうか
        */
        bool IsFinished() const;

        //! A private member function.
        /*!
            Iterative solution of nonlinear Kohn-Sham equation, SCF loop
        */
        void Scf();

        //! A private member function (const).
        /*!
            Writes energy and eigenvalues into files.
            "sec" - time of calcuylations in seconds.
            \param sec タイムポイント
        */
        void WriteInfo(std::chrono::duration<double> const & sec) const;
        
        //! A private member function (const).
        /*!
            Write results into files
            \param sec タイムポイント
        */
        void WriteRes(std::chrono::duration<double> const & sec) const;

        // #endregion privateメンバ関数
        
        // #region メンバ変数

        //! A private member variable.
        /*!
            Interaction potential（αスピンとβスピンのstd::pair）
        */
        std::pair<std::shared_ptr<Pot<util::Spin::Alpha>>, std::shared_ptr<Pot<util::Spin::Beta>>> m_pot;

        //! A private member variable.
        /*!
            Solver for LINER Kohna-Shama equation（αスピンとβスピンのstd::pair）
        */
        std::pair<std::shared_ptr<KohnSham<util::Spin::Alpha>>, std::shared_ptr<KohnSham<util::Spin::Beta>>> m_ks;

        //! A private member variable.
        /*!
            Set of all electronic states (eigenfunctions)
        */
        std::shared_ptr<StateSet> m_ss_alpha;

        //! A private member variable.
        /*!
            Set of all electronic states (eigenfunctions)
        */
        std::shared_ptr<StateSet> m_ss_beta;
        
        //! A private member variable.
        /*!
            Database of parameters
        */
        std::shared_ptr<ParamDb> m_db;

        //! A private member variable.
        /*!
            Calculates required energy of atom
        */
        std::unique_ptr<Energy> m_energy;
        
        //! A private member variable.
        /*!
            Electron density（αスピンとβスピンのstd::pair）
        */
        std::pair<std::shared_ptr<Rho>, std::shared_ptr<Rho>> m_rho;
        
        // #endregion メンバ変数 

    private:
        //
        // Linear mixing of states neaded to aobtain SCF convergence
        //
        class RhoMix : public util::Fun1D
        {
        public:
            RhoMix(std::shared_ptr<ParamDb> const & db)
            {
                m_scfMix = atof(db->Get("Scf_Mix"));
            }
            virtual ~RhoMix() = default;
            // May23rd, 2014 Modified by dc1394
            //void SetRho(util::Fun1D* rhoCur, util::Fun1D* rhoOld)
            void SetRho(std::shared_ptr<util::Fun1D> const & rhoCur, std::shared_ptr<util::Fun1D> const & rhoOld)
            {
                m_rhoCur = rhoCur;
                m_rhoOld = rhoOld;
            }
            virtual double Get(double r) const
            {
                return m_scfMix * m_rhoCur->Get(r) + (1 - m_scfMix) * m_rhoOld->Get(r);
            }
        private:
            // May 23rd, 2014 Modified by dc1394
            //util::Fun1D* m_rhoCur;
            std::shared_ptr<util::Fun1D> m_rhoCur;
            //util::Fun1D* m_rhoOld;
            std::shared_ptr<util::Fun1D> m_rhoOld;
            double m_scfMix;
        };
    };
}

#endif

