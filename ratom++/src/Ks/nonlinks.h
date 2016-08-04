#ifndef __RATOM_NONLINKS_H__
#define __RATOM_NONLINKS_H__


/** \brief Nonlinear Kohn-Sham equation
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/
// March 16th, 2014 Modified by dc1394

#include "paramdb.h"
#include "rho.h"
#include "pot.h"
#include "kohnsham.h"
#include "energy.h"
#include "../Util/spin.h"

namespace ks {
    // May 26th, 2014 Added by dc1394
    class NonLinKs final
    {
    public:
        NonLinKs(const char* path);
        ~NonLinKs(void);

        int Run(void);

    private:
        void Scf(void);
        bool IsFinished(void) const;

        void WriteRes(time_t sec) const;
        void WriteInfo(time_t sec) const;

    private:
        // Interaction potential
        //Pot *m_pot;
        std::pair<std::shared_ptr<Pot<util::Spin::Alpha>>, std::shared_ptr<Pot<util::Spin::Beta>>> m_pot;

        // Solver for LINER Kohna-Shama equation
        // May 25th, 2014 Modified by dc1394
        //KohnSham *m_ks;
        std::pair <std::shared_ptr<KohnSham<util::Spin::Alpha>>, std::shared_ptr<KohnSham<util::Spin::Beta>>> m_ks;

        // Set of all electronic states (eigenfunctions)
        StateSet *m_ss;

        // Database of parameters
        ParamDb m_db;

        // Calculates required energy of atom
        Energy *m_energy;

        // Electron density
        // May 23rd, 2014 Modified by dc1394
        //Rho *m_rho;
        std::pair<std::shared_ptr<Rho>, std::shared_ptr<Rho>> m_rho;



    private:
        //
        // Linear mixing of states neaded to aobtain SCF convergence
        //
        class RhoMix : public util::Fun1D
        {
        public:
            RhoMix(const ParamDb* db) : m_rhoCur(NULL), m_rhoOld(NULL)
            {
                m_scfMix = atof(db->Get("Scf_Mix"));
            }
            virtual ~RhoMix() { }
            // May23rd, 2014 Modified by dc1394
            //void SetRho(util::Fun1D* rhoCur, util::Fun1D* rhoOld)
            void SetRho(const std::shared_ptr<util::Fun1D> & rhoCur, const std::shared_ptr<util::Fun1D> & rhoOld)
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

