// modified by dc1394 - March 7th, 2014

#include "../ExCorr/exchhf.h"
#include "../ExCorr/corrhf.h"
#include "../ExCorr/exchpbe0.h"
#include "../ExCorr/excorrgga.h"
#include "../ExCorr/excorrlda.h"
#include "stdafx.h"
#include "pot.h"
#include "rho.h"
#include <cstdint>
#include <boost/cast.hpp>

namespace ks {
    // May 24th, 2014 Modified by dc1394
    //
    // Constructor
    //
    template <util::Spin S>
    Pot<S>::Pot(std::shared_ptr<ParamDb> const & db)
        :   Hart([this] { return std::cref(m_hart); }, [this](std::shared_ptr<OdeProb> const & hart) { m_hart = hart; return std::cref(m_hart); }),
            m_db(db),
            m_hart(std::make_shared<OdeProb>())
    {
        // March 28th, 2014 Modified by dc1394
        m_z = db->GetDouble("Atom_Proton");
        SetXc(db->Get("XC_Exch"), db->Get("XC_Corr"));

        //m_z = db->GetDouble("Atom_Proton");
        m_rc = m_db->GetDouble("Atom_Rc");
    }
    
    template <util::Spin S>
    double Pot<S>::GetRho(double r) const
    {
        return GetRho(r, boost::mpl::int_<static_cast<std::int32_t>(S)>());
    }

    // May 24th, 2014 Modified by dc1394
    //
    // Returns value of electron density for radius "r"
    //
    template <util::Spin S>
    double Pot<S>::GetRho(double r, boost::mpl::int_<static_cast<std::int32_t>(util::Spin::Alpha)>) const
    {
    	//return m_rho->Get(r);
        return m_rho.first->Get(r);
    }


    //
    // Returns value of electron density for radius "r"
    //
    template <util::Spin S>
    double Pot<S>::GetRho(double r, boost::mpl::int_<static_cast<std::int32_t>(util::Spin::Beta)>) const
    {
        //return m_rho->Get(r);
        return m_rho.second->Get(r);
    }


    /*! Returns effective potential for radial Kohn-Sham equation for radius "r"
    @param[in]      r    原点からの距離（極座標）
    @return         effective potential for radial Kohn-Sham equation for radius "r"
    @exception      none
    */
    template <util::Spin S>
    double Pot<S>::Get(double r) const
    {
        double rho, vx, vc, vn, vh;

        assert(r > 0);

        //	// For testing purpose - hydrogen atom
        //	return -1 / r; 

        //assert(m_rho);
        //rho = m_rho->Get(r) / (util::HelpFun::M_4PI * r * r);

        //assert(rho >= 0);

        // March 7th, 2014	Modified by dc1394
        vn = Vn(r);
        vx = Vx(r);
        vc = Vc(r);
        vh = Vh(r);
        // vl = m_l * (m_l + 1) / (2 * r * r);  // Dependency on angular quantum number

        // printf("%lf   %lf  %lf   %lf   %lf\n", rho, vx, vc, vn, ve);
        return vx + vc + vn + vh;
    }

    /*! Defines electron density for solving Poisson equation
    @param[in]      rho    αスピン密度とβスピン密度のpair
    @exception      none
    */
    //void Pot::SetRho(util::Fun1D* rho)
    template <util::Spin S>
    void Pot<S>::SetRho(std::pair<std::shared_ptr<util::Fun1D>, std::shared_ptr<util::Fun1D>> const & rho)
    {
        //const double rc = m_db->GetDouble("Atom_Rc");
        //const size_t psnNode = m_db->GetSize_t("Solver_PsnNode");
        //const size_t psnDeg = m_db->GetSize_t("Solver_PsnDeg");

        //const double gamma = 1;

        //// const Bndr left(BndrType_Dir, 0), right(BndrType_Dir, m_z);
        //// Zero Dirichlet boundary conditions
        //const Bndr left(BndrType_Dir, 0), right(BndrType_Dir, 0);

        ////assert(rho);
        //m_rho = rho;
        //m_rhoHelp.m_rho = rho;

        //m_hart->Define(left, right, gamma, NULL, &m_rhoHelp);
        //m_hart->GenMeshLin(0, rc, psnNode, psnDeg);
        SetRho(rho, boost::mpl::int_<static_cast<std::int32_t>(S)>());
    }

    template <util::Spin S>
    void Pot<S>::SetRho(std::pair<std::shared_ptr<util::Fun1D>, std::shared_ptr<util::Fun1D>> const & rho, boost::mpl::int_<static_cast<std::int32_t>(util::Spin::Alpha)>)
    {
        const double rc = m_db->GetDouble("Atom_Rc");
        const size_t psnNode = m_db->GetSize_t("Solver_PsnNode");
        const size_t psnDeg = m_db->GetSize_t("Solver_PsnDeg");

        const double gamma = 1;

        // const Bndr left(BndrType_Dir, 0), right(BndrType_Dir, m_z);
        // Zero Dirichlet boundary conditions
        const Bndr left(BndrType_Dir, 0), right(BndrType_Dir, 0);

        //assert(rho);
        m_rho = rho;
        m_rhoHelp.m_rho = rho;

        m_hart->Define(left, right, gamma, NULL, &m_rhoHelp);
        m_hart->GenMeshLin(0, rc, psnNode, psnDeg);
    }

    template <util::Spin S>
    void Pot<S>::SetRho(std::pair<std::shared_ptr<util::Fun1D>, std::shared_ptr<util::Fun1D>> const & rho, boost::mpl::int_<static_cast<std::int32_t>(util::Spin::Beta)>)
    {
        m_rho = rho;
        m_rhoHelp.m_rho = rho;
    }


    /*! Solves Poission equation
    @exception      none
    */
    template <util::Spin S>
    void Pot<S>::SolvePoisson(void)
    {
        const double absMaxCoef = m_db->GetDouble("Solver_PsnAbsMaxCoef");
        const bool adapt = m_db->GetBool("Solver_PsnAdapt");

        if (adapt)
            m_hart->SolveAdapt(absMaxCoef);
        else
            m_hart->Solve();

        // m_hart->WriteSol("C:\\romz\\Syf\\OdeProb\\Hartree.dat", 2000);
    }


    //
    // 
    //

    /*! Electrostatic potential of atomic core
    @param[in]      r    原点からの距離（極座標）
    @return         Electrostatic potential of atomic core
    @exception      none
    */
    template <util::Spin S>
    double Pot<S>::Vn(double r) const
    {
        return -m_z / r;
    }


    // May 25th, 2014 comment out by dc1394
    //
    // Exchange potential
    //
    //double Pot::Vx(double rho, double rhoDer) const
    //{
    //	assert(m_exch);
    //	return m_exch->V(rho, rhoDer);
    //}
    /*! Exchange potential for radius "r"
    @param[in]      r    原点からの距離（極座標）
    @return         Exchange potential for radius "r"
    @exception      none
    */
    template <util::Spin S>
    double Pot<S>::Vx(double r) const
    {
        return m_exch->V(r);
    }

    // May 25th, 2014 comment out by dc1394
    //
    // Correlation potential
    //
    //double Pot::Vc(double rho, double rhoDer) const
    //{
    //	assert(m_corr);
    //	return m_corr->V(rho, rhoDer);
    //}

    //
    // Correlation potential for radius "r"
    // March 7th, 2014	Added by dc1394
    //
    template <util::Spin S>
    double Pot<S>::Vc(double r) const
    {
        return m_corr->V(r);
    }

    //
    // Hartree potential
    //
    template <util::Spin S>
    double Pot<S>::Vh(double r) const
    {
        assert(r > 0);
        const double val = m_hart->GetSol(r);

        // Apply non-zero Dirichlet boundary conditions
        const double ua = 0, ub = m_z, a = 0, b = m_rc;
        const double alpha = (ub - ua) / (b - a), beta = ua - a * alpha;

        return (val + alpha * r + beta) / r;
    }

    // May 25th, 2014 comment out by dc1394
    //
    // Density of exchnage energy
    //
    //double Pot::Ex(double rho, double rhoDer) const
    //{
    //	assert(m_exch);
    //	return m_exch->E(rho, rhoDer);
    //}

    //
    // Density of exchnage energy for radius "r"
    // March 7th, 2014	Added by dc1394
    //
    template <util::Spin S>
    double Pot<S>::Ex(double r) const
    {
        return m_exch->E(r);
    }

    // May 25th, 2014 comment out by dc1394
    //
    // Density of correlation energy
    //
    //double Pot::Ec(double rho, double rhoDer) const
    //{
    //	assert(m_corr);
    //	return m_corr->E(rho, rhoDer);
    //}

    //
    // Density of correlation energy for radius "r"
    // March 8th, 2014	Added by dc1394
    //
    template <util::Spin S>
    double Pot<S>::Ec(double r) const
    {
        return m_corr->E(r);
    }

    //
    // Returns difference between density of energy and potential
    //
    template <util::Spin S>
    double Pot<S>::XcEdiffV(double r) const
    {
        return m_exch->EdiffV(r) + m_corr->EdiffV(r);
    }


    //
    // Defines exchange and correlation approximation
    // exch - name of exchange potential
    // corr - name of correlation potential
    //
    template <util::Spin S>
    void Pot<S>::SetXc(std::string const & exch, std::string const & corr)
    {
        // March 7th, 2014	Modified by dc1394
        if (exch == "slater") {
            m_exch.reset(new excorr::Xc<S>(excorr::ExCorrLDA([this](double r) { return GetRhoTilde(r); }, XC_LDA_X)));
        }
        else if (exch == "b88") {
            m_exch.reset(
                new excorr::Xc<S>(
                    excorr::ExCorrGGA(
                        [this](double r) { return GetRhoTilde(r); },
                        [this](double r) { return GetRhoTildeDeriv(r); },
                        [this](double r) { return GetRhoTildeLapl(r); },
                        XC_GGA_X_B88)));
        }
        else if (exch == "pbe") {
            m_exch.reset(
                new excorr::Xc<S>(
                    excorr::ExCorrGGA(
                        [this](double r) { return GetRhoTilde(r); },
                        [this](double r) { return GetRhoTildeDeriv(r); },
                        [this](double r) { return GetRhoTildeLapl(r); },
                        XC_GGA_X_PBE)));
        }
        else if (exch == "pbe0") {
            m_exch.reset(
                new excorr::Xc<S>(
                    excorr::ExchPbe0(
                        [this](double r) { return GetRhoTilde(r); },
                        [this](double r) { return GetRhoTildeDeriv(r); },
                        [this](double r) { return GetRhoTildeLapl(r); },
                        [this](double r) { return Vh(r); },
                        XC_GGA_X_PBE,
                        m_z)));
        }
        else if (exch == "hf") {
            m_exch.reset(new excorr::Xc<S>(excorr::ExchHf([this](double r) { return Vh(r); }, m_z)));
        }
        //else {
        //    throw std::invalid_argument("Unknown exchange type");
        //}

        if (corr == "vwn") {
            m_corr.reset(new excorr::Xc<S>(excorr::ExCorrLDA([this](double r) { return GetRhoTilde(r); }, XC_LDA_C_VWN)));
        }
        else if (corr == "pbe") {
            m_corr.reset(
                new excorr::Xc<S>(
                    excorr::ExCorrGGA(
                        [this](double r) { return GetRhoTilde(r); },
                        [this](double r) { return GetRhoTildeDeriv(r); },
                        [this](double r) { return GetRhoTildeLapl(r); },
                        XC_GGA_C_PBE)));
        }
        //else if (strcmp(corr, "pbe") == 0) {
        //    m_corr.reset(new CorrPbe(std::bind(&Pot::GetRhoTilde, std::ref(*this),
        //        std::placeholders::_1),
        //        std::bind(&Pot::GetRhoTildeDeriv, std::ref(*this),
        //        std::placeholders::_1),
        //        std::bind(&Pot::GetRhoTildeLapl, std::ref(*this),
        //        std::placeholders::_1)));
        //}
        else if (corr == "hf") {
            m_corr.reset(new excorr::Xc<S>(excorr::CorrHf()));
        }
        //else
        //    throw std::invalid_argument("Unknown correlation type");
    }


    //
    // Writes potential int file
    //
    // March 31st, 2014	Added by dc1394
    template <util::Spin S>
    void Pot<S>::Write() /*const*/
    {
        auto out = m_db->OpenFile("pot", "wt");

        //m_exch.reset(new excorr::Xc<S>(excorr::ExCorrLDA([this](double r) { return GetRhoTilde(r); }, XC_LDA_X)));
        //m_corr.reset(new excorr::Xc<S>(excorr::ExCorrLDA([this](double r) { return GetRhoTilde(r); }, XC_LDA_C_VWN)));
        //m_exch.reset(
        //    new excorr::Xc<S>(
        //        excorr::ExCorrGGA(
        //            [this](double r) { return GetRhoTilde(r); },
        //            [this](double r) { return GetRhoTildeDeriv(r); },
        //            [this](double r) { return GetRhoTildeLapl(r); },
        //            XC_GGA_X_PBE)));

        //m_corr.reset(
        //    new excorr::Xc<S>(
        //        excorr::ExCorrGGA(
        //            [this](double r) { return GetRhoTilde(r); },
        //            [this](double r) { return GetRhoTildeDeriv(r); },
        //            [this](double r) { return GetRhoTildeLapl(r); },
        //            XC_GGA_C_PBE)));
        
        /*m_exch.reset(
            new excorr::Xc<S>(
                excorr::ExchPbe0(
                    [this](double r) { return GetRhoTilde(r); },
                    [this](double r) { return GetRhoTildeDeriv(r); },
                    [this](double r) { return GetRhoTildeLapl(r); },
                    [this](double r) { return Vh(r); },
                    XC_GGA_X_PBE,
                    m_z)));

        m_corr.reset(
            new excorr::Xc<S>(
                excorr::ExCorrGGA(
                    [this](double r) { return GetRhoTilde(r); },
                    [this](double r) { return GetRhoTildeDeriv(r); },
                    [this](double r) { return GetRhoTildeLapl(r); },
                    XC_GGA_C_PBE)));*/

        fprintf(out.get(), "%16s \t %16s \t %16s \t %16s \n", "R", "Hartree potential", "Exchange potential", "Correlation potential");

        auto const m = boost::numeric_cast<std::int32_t>(MAX / DR);
        for (auto k = 1; k <= m; ++k)
        {
            auto const r = static_cast<double>(k)* DR;
            fprintf(out.get(), "%16.6E \t %16.6E \t %16.6E \t %16.6E\n", r, Vh(r), Vx(r), Vc(r));
        }
    }

    // April 3rd, 2014 Added by dc1394
    // Returns value of electron density (divided by 4 * pi * r^2) for radius "r"
    //
    template <util::Spin S>
    std::pair<double, double> Pot<S>::GetRhoTilde(double r)
    {
        return std::make_pair(m_rho.first->Get(r) / (util::HelpFun::M_4PI * r * r),
            m_rho.second->Get(r) / (util::HelpFun::M_4PI * r * r));
    }


    // April 3rd, 2014 Added by dc1394
    // Returns value of electron density (divided by 4 * pi * r^2) derivative for radius "r"
    //
    template <util::Spin S>
    std::pair<double, double> Pot<S>::GetRhoTildeDeriv(double r)
    {
        auto const alpha = (m_rho.first->GetDeriv(r) - 2.0 * m_rho.first->Get(r) / r) / (util::HelpFun::M_4PI * r * r);
        auto const beta = (m_rho.second->GetDeriv(r) - 2.0 * m_rho.second->Get(r) / r) / (util::HelpFun::M_4PI * r * r);

        return std::make_pair(alpha, beta);
    };


    //
    // Returns value of electron density (divided by 4 * pi * r^2) laplacian for radius "r"
    // April 3rd, 2014 added by dc1394
    //
    template <util::Spin S>
    std::pair<double, double> Pot<S>::GetRhoTildeLapl(double r)
    {
        auto const first = (m_rho.first->Get2ndDeriv(r) - 2.0 / r * m_rho.first->GetDeriv(r) +
            2.0 / (r * r) * m_rho.first->Get(r)) / (util::HelpFun::M_4PI * r * r);
        auto const second = (m_rho.second->Get2ndDeriv(r) - 2.0 / r * m_rho.second->GetDeriv(r) +
            2.0 / (r * r) * m_rho.second->Get(r)) / (util::HelpFun::M_4PI * r * r);

        return std::make_pair(first, second);
    };

    template class Pot<util::Spin::Alpha>;
    template class Pot<util::Spin::Beta>;
}
