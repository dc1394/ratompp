#ifndef __RATOM_POT_H__
#define __RATOM_POT_H__


/** \brief Interaction potential (without part ddependent on L) in Kohn-Sham equation
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/

#include "paramdb.h"
#include "../Util/property.h"
#include "../ExCorr/xc.h"
#include <memory>
#include <utility>
#include <boost/mpl/int.hpp>    // for boost::mpl::int_

namespace ks {
    template <util::Spin S>
    class Pot final : public util::Fun1D
    {
        friend class Energy;

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
                // return 4 * util::HelpFun::M_PI * r * m_rho->Get(r); 
                //return m_rho->Get(r) / r;
                // alpha spin + beta spin
                return (m_rho.first->Get(r) + m_rho.second->Get(r)) / r;
            }
        public:
            //util::Fun1D* m_rho;
            std::pair<std::shared_ptr<util::Fun1D>, std::shared_ptr<util::Fun1D> > m_rho;
        };


    public:
        Pot(std::shared_ptr<ParamDb> const & db);
        virtual ~Pot() = default;

        // May 23rd, 2014 Modified by dc1394
        //void SetRho(util::Fun1D* rho);
        void SetRho(std::pair<std::shared_ptr<util::Fun1D>, std::shared_ptr<util::Fun1D>> const & rho);
        void SetRho(std::pair<std::shared_ptr<util::Fun1D>, std::shared_ptr<util::Fun1D>> const & rho, boost::mpl::int_<static_cast<std::int32_t>(util::Spin::Alpha)>);
        void SetRho(std::pair<std::shared_ptr<util::Fun1D>, std::shared_ptr<util::Fun1D>> const & rho, boost::mpl::int_<static_cast<std::int32_t>(util::Spin::Beta)>);
        void SolvePoisson(void);

        virtual double Get(double r) const;
        double GetRho(double r) const;
        double GetRho(double r, boost::mpl::int_<static_cast<std::int32_t>(util::Spin::Alpha)>) const;
        double GetRho(double r, boost::mpl::int_<static_cast<std::int32_t>(util::Spin::Beta)>) const;

        // March 31st, 2014 Added by @dc1394
        void Write() /*const*/;

        double Vn(double r) const;
        double Vx(double r) const;
        double Vc(double r) const;
        double Vh(double r) const;

        double Ex(double r) const;
        double Ec(double r) const;
        double XcEdiffV(double r) const;

        void SetXc(std::string const & exch, std::string const & corr);

        // April 3rd, 2014 Added by dc1394
        std::pair<double, double> GetRhoTilde(double r);
        std::pair<double, double> GetRhoTildeDeriv(double r);
        std::pair<double, double> GetRhoTildeLapl(double r);

        //! A property.
        /*!
        */
        util::Property< std::shared_ptr<OdeProb> &> Hart;

    private:
        // Electron density
        std::pair<std::shared_ptr<util::Fun1D>, std::shared_ptr<util::Fun1D>> m_rho;

    public:
        // March 7th, 2014	Modified by dc1394 
        // Exchenge potential
        //Xc* m_exch;
        std::unique_ptr<excorr::Xc<S>> m_exch;
        // Correlation potential
        //Xc* m_corr;
        std::unique_ptr<excorr::Xc<S>> m_corr;

    private:
        // Number of protons in atom
        double m_z;

        // Maximal radius of atom
        double m_rc;

        // Hartree potential
        std::shared_ptr<OdeProb> m_hart;

        // Helper funtion required to solve Poisson equation
        RhoHelp m_rhoHelp;

        // Database of parameters
        std::shared_ptr<ParamDb> const m_db;

        static auto constexpr MAX = 10.0;
        static auto constexpr DR = 1.0E-3;
    };
}
#endif

