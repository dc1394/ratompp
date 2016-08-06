#include "stdafx.h"
#include "kohnsham.h"

namespace ks {
    template <util::Spin S>
    void KohnSham<S>::Config(std::shared_ptr<util::Fun1D> const & pot)
    {
        const auto rc = m_db->GetDouble("Atom_Rc");
        const auto eigNode = m_db->GetSize_t("Solver_EigNode");
        const auto eigDeg = m_db->GetSize_t("Solver_EigDeg");

        const auto gamma = 0.5; // Paramter for radial Kohn-Sham equation
        const auto Lmax = m_stateSet->GetLmax();

        // May 25th, 2014 Modified by dc1394
        //m_radPot.m_pot = pot;
        m_radPot.getPot() = pot;

        m_eigProb.resize(Lmax);
        for (auto l = 0U; l < Lmax; l++)
        {
            m_eigProb[l].Define(gamma, &m_radPot);
            m_eigProb[l].GenMeshLin(0, rc, eigNode, eigDeg);
        }

        m_eigNo.resize(Lmax);
        for (auto l = 0U; l < Lmax; l++)
            m_eigNo[l] = m_stateSet->GetNmax(l);
    }

    template <util::Spin S>
    void KohnSham<S>::Solve(void)
    {
        const auto abstol = m_db->GetDouble("Solver_EigAbsTol");
        const auto absMaxCoef = m_db->GetDouble("Solver_EigAbsMaxCoef");
        const auto adapt = m_db->GetBool("Solver_EigAdapt");

        double eigVal;

        auto const size = m_eigProb.size();
        for (auto l = 0U; l < size; l++)
        {
            m_radPot.m_l = l;

            if (adapt)
                m_eigProb[l].SolveAdapt(m_eigNo[l], abstol, absMaxCoef);
            else
                m_eigProb[l].Solve(m_eigNo[l], abstol);

            // Sets eigenvalues of states
            for (auto n = 0U; n < m_eigNo[l]; n++)
            {
                eigVal = m_eigProb[l].GetEigVal(n);
                m_stateSet->SetEigVal(l, n, eigVal);
            }
        }
    }

    template <util::Spin S>
    double KohnSham<S>::Get(double r) const
    {
        // スピンに応じた関数を呼び出す
        return Get(r, boost::mpl::int_<static_cast<std::int32_t>(S)>());

        //const double rc = m_db->GetDouble("Atom_Rc");
        //const size_t Lmax = m_stateSet->GetLmax();
        //size_t l, n;
        //double rnl, occ, rhoL, rho = 0;

        //if (r >= rc)
        //    return 0;

        //// For all quantum angular menetum numbers
        //for (l = 0; l < Lmax; l++)
        //{
        //    rhoL = 0;

        //    // For all states for fixed "L"
        //    for (n = 0; n < m_eigNo[l]; n++)
        //    {
        //        occ = m_stateSet->Occ(l, n);
        //        //if (occ > 0)
        //        if (occ > 0.0 && occ <= 2.0 * static_cast<double>(l)+1.0) {
        //            
        //            if (AlphaOrBetaFlag) {
        //                rnl = m_eigProb[l].GetEigFun(n, r); // R_{n, \ell}(r)
        //                rhoL += occ * rnl * rnl;
        //            }
        //            // Beta spin
        //        {
        //        }
        //        }
        //        else if (occ > 2.0 * static_cast<double>(l)+1.0) {
        //            // Alpha spin
        //            if (AlphaOrBetaFlag) {
        //                rnl = m_eigProb[l].GetEigFun(n, r); // R_{n, \ell}(r)
        //                rhoL += (2.0 * static_cast<double>(l)+1.0) * rnl * rnl;
        //            }
        //            // Beta spin
        //            else {
        //                rnl = m_eigProb[l].GetEigFun(n, r); // R_{n, \ell}(r)
        //                rhoL += (occ - (2.0 * static_cast<double>(l)+1.0)) * rnl * rnl;
        //            }
        //        }
        //    }

        //    // Sum up all constituents
        //    rho += rhoL;
        //}

        //return rho;
    }

    template <util::Spin S>
    double KohnSham<S>::Get(double r, boost::mpl::int_<static_cast<std::int32_t>(util::Spin::Alpha)>) const
    {
        const double rc = m_db->GetDouble("Atom_Rc");
        const size_t Lmax = m_stateSet->GetLmax();
        size_t l, n;
        double rnl, occ, rhoL, rho = 0;

        if (r >= rc)
            return 0;

        // For all quantum angular menetum numbers
        for (l = 0; l < Lmax; l++)
        {
            rhoL = 0;

            // For all states for fixed "L"
            for (n = 0; n < m_eigNo[l]; n++)
            {
                occ = m_stateSet->Occ(l, n);
                //if (occ > 0)
                // occにはαスピンの電子のみ
                if (occ > 0.0 && occ <= static_cast<double>(Lmax)) {
                    rnl = m_eigProb[l].GetEigFun(n, r); // R_{n, \ell}(r)
                    // Lmaxまではαスピンの電子数はoccと等しい
                    rhoL += occ * rnl * rnl;
                }
                // occにβスピンの電子が含まれる
                else if (occ > static_cast<double>(Lmax)) {
                    rnl = m_eigProb[l].GetEigFun(n, r); // R_{n, \ell}(r)
                    // αスピンの電子数はLmaxと等しい
                    rhoL += static_cast<double>(Lmax) * rnl * rnl;
                }
            }

            // Sum up all constituents
            rho += rhoL;
        }

        return rho;
    }

    template <util::Spin S>
    double KohnSham<S>::Get(double r, boost::mpl::int_<static_cast<std::int32_t>(util::Spin::Beta)>) const
    {
        const double rc = m_db->GetDouble("Atom_Rc");
        const size_t Lmax = m_stateSet->GetLmax();
        size_t l, n;
        double rnl, occ, rhoL, rho = 0;

        if (r >= rc)
            return 0;

        // For all quantum angular menetum numbers
        for (l = 0; l < Lmax; l++)
        {
            rhoL = 0;

            // For all states for fixed "L"
            for (n = 0; n < m_eigNo[l]; n++)
            {
                occ = m_stateSet->Occ(l, n);
                //if (occ > 0)
                // occがLmaxより大きければ、occにβスピンの電子が含まれる
                if (occ > static_cast<double>(Lmax)) {
                    rnl = m_eigProb[l].GetEigFun(n, r); // R_{n, \ell}(r)
                    // βスピンの電子数は、occからLmax分差し引いたもの
                    rhoL += (occ - static_cast<double>(Lmax)) * rnl * rnl;
                }
            }

            // Sum up all constituents
            rho += rhoL;
        }

        return rho;
    }

    template <util::Spin S>
    void KohnSham<S>::WriteEigen(void) const
    {
        const size_t eigNode = m_db->GetSize_t("Out_EigNode");
        const size_t Lmax = m_stateSet->GetLmax();
        size_t l, n;
        std::string path;
        char buf[10];


        // For each angular quantum number
        for (l = 0; l < Lmax; l++)
        {
            m_db->GetPath(path);
            sprintf(buf, ".L=%lu", static_cast<unsigned long>(l));
            path += buf;
            m_eigProb[l].WriteEigCoef(path.c_str(), m_eigNo[l]);

            // For each state for fixed "L"
            for (n = 0; n < m_eigNo[l]; n++)
            {
                m_db->GetPath(path);
                path += m_stateSet->Name(l, n);
                path += ".nto";

                m_eigProb[l].WriteEigFun(path.c_str(), n, eigNode);

            }
        }
    }

    template class KohnSham<util::Spin::Alpha>;
    template class KohnSham<util::Spin::Beta>;
}
