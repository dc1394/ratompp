#ifndef __RATOM_POT_H__
#define __RATOM_POT_H__


/** \brief Interaction potential (without part ddependent on L) in Kohn-Sham equation
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/

#include "paramdb.h"
#include "../ExCorr/xc.h"
#include <memory>
#include <utility>

extern void* enabler;
namespace ks {
    template <util::Spin S>
    class Pot final : public util::Fun1D
    {
        friend class Energy;
        friend class FunEner;

    private:
        // May 23rd, 2014 Added by dc1394
        class RhoHelp : public util::Fun1D
        {
        public:
            RhoHelp() /*: m_rho(NULL)*/ {}
            virtual ~RhoHelp() {}
            virtual double Get(double r) const
            {
                //assert(m_rho); 
                // return 4 * M_PI * r * m_rho->Get(r); 
                //return m_rho->Get(r) / r;
                // alpha spin + beta spin
                return (m_rho.first->Get(r) + m_rho.second->Get(r)) / r;
            }
        public:
            //util::Fun1D* m_rho;
            std::pair<std::shared_ptr<util::Fun1D>, std::shared_ptr<util::Fun1D>> m_rho;
        };


    public:
        Pot(const ParamDb* db);
        virtual ~Pot(void);

        // May 23rd, 2014 Modified by dc1394
        //void SetRho(util::Fun1D* rho);
        void SetRho(std::pair<std::shared_ptr<util::Fun1D>, std::shared_ptr<util::Fun1D>> const & rho);
        void SolvePoisson(void);

        virtual double Get(double r) const;
        template <util::Spin S>
        double GetRho(double r) const;
        template <util::Spin S>
        double GetRho(double r) const;

        // March 31st, 2014 Added by @dc1394
        void Write() /*const*/;

    private:
        double Vn(double r) const;
        //double Vx(double rho, double rhoDer) const;
        // March 7th, 2014	Added by dc1394
        double Vx(double r) const;
        //double Vc(double rho, double rhoDer) const;
        // March 8th, 2014	Added by dc1394
        double Vc(double r) const;
        double Vh(double r) const;

        //double Ex(double rho, double rhoDer) const;
        // March 7th, 2014	Added by dc1394
        double Ex(double r) const;
        //double Ec(double rho, double rhoDer) const;
        // March 7th, 2014	Added by dc1394
        double Ec(double r) const;
        double XcEdiffV(double r) const;

        void SetXc(const char* exch, const char* corr);

        // April 3rd, 2014 Added by dc1394
        std::pair<double, double> GetRhoTilde(double r);
        std::pair<double, double> GetRhoTildeDeriv(double r);
        std::pair<double, double> GetRhoTildeLapl(double r);

    private:
        // Electron density
        // May 23rd, 2014 Modified by dc1394
        //util::Fun1D* m_rho;
        std::pair<std::shared_ptr<util::Fun1D>, std::shared_ptr<util::Fun1D>> m_rho;

        // March 7th, 2014	Modified by dc1394 
        // Exchenge potential
        //Xc* m_exch;
        std::unique_ptr<excorr::Xc<S>> m_exch;
        // Correlation potential
        //Xc* m_corr;
        std::unique_ptr<excorr::Xc<S>> m_corr;

        // Number of protons in atom
        double m_z;

        // Maximal radius of atom
        double m_rc;

        // Hartree potential
        OdeProb m_hart;

        // Helper funtion required to solve Poisson equation
        RhoHelp m_rhoHelp;

        // Database of parameters
        const ParamDb* m_db;

#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__)
        static constexpr double MAX = 10.0;
        static constexpr double DR = 1.0E-3;
#else
        static const double MAX;
        static const double DR;
#endif
    };
}
#endif

