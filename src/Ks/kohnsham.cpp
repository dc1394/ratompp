#include "stdafx.h"
#include "kohnsham.h"
#include <boost/format.hpp>     // for boost::format.hpp

namespace ks {
    template <util::Spin S>
    KohnSham<S>::KohnSham(std::shared_ptr<const ParamDb> && db, std::shared_ptr<StateSet> const & stateSet)
        :   m_db(std::move(db)),
            m_radPot(std::make_shared<PotRad>()),
            m_stateSet(stateSet)
    {}

    template <util::Spin S>
    KohnSham<S>::KohnSham(std::shared_ptr<const ParamDb> && db, std::shared_ptr<StateSet> && stateSet)
        : m_db(std::move(db)),
        m_radPot(std::make_shared<PotRad>()),
        m_stateSet(std::move(stateSet))
    {}

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
        m_radPot->getPot() = pot;

        m_eigProb.resize(Lmax);
        for (auto l = 0U; l < Lmax; l++)
        {
            m_eigProb[l].Define(gamma, m_radPot);
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
            m_radPot->m_l = l;

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

    template <>
    double KohnSham<util::Spin::Alpha>::Get(double r) const
    {
        auto const rc = m_db->GetDouble("Atom_Rc");
        auto const Lmax = m_stateSet->GetLmax();

        if (r >= rc) {
            return 0.0;
        }

        auto rho = 0.0;
        // For all quantum angular menetum numbers
        for (auto l = 0U; l < Lmax; l++)
        {
            auto rhoL = 0.0;

            // For all states for fixed "L"
            for (auto n = 0U; n < m_eigNo[l]; n++)
            {
                auto const occ = m_stateSet->Occ(l, n);
                //if (occ > 0)
                // occにはαスピンの電子のみ
                if (occ > 0.0 && occ <= static_cast<double>(2 * l + 1)) {
                    auto const rnl = m_eigProb[l].GetEigFun(n, r); // R_{n, \ell}(r)
                    // occが2 * l + 1まではαスピンの電子数はoccと等しい
                    rhoL += occ * rnl * rnl;
                }
                // occにβスピンの電子が含まれる
                else {
                    auto const rnl = m_eigProb[l].GetEigFun(n, r); // R_{n, \ell}(r)
                    // αスピンの電子数は2 * l + 1と等しい
                    rhoL += static_cast<double>(2 * l + 1) * rnl * rnl;
                }
            }

            // Sum up all constituents
            rho += rhoL;
        }

        return rho;
    }

    template <>
    double KohnSham<util::Spin::Beta>::Get(double r) const
    {
        auto const rc = m_db->GetDouble("Atom_Rc");
        auto const Lmax = m_stateSet->GetLmax();

        if (r >= rc) {
            return 0.0;
        }

        auto rho = 0.0;
        // For all quantum angular menetum numbers
        for (auto l = 0U; l < Lmax; l++)
        {
            auto rhoL = 0.0;

            // For all states for fixed "L"
            for (auto n = 0U; n < m_eigNo[l]; n++)
            {
                auto const occ = m_stateSet->Occ(l, n);
                //if (occ > 0)
                // occが2 * l + 1より大きければ、occにβスピンの電子が含まれる
                if (occ > static_cast<double>(2 * l + 1)) {
                    auto const rnl = m_eigProb[l].GetEigFun(n, r); // R_{n, \ell}(r)
                    // βスピンの電子数は、occから2 * l + 1分差し引いたもの
                    rhoL += (occ - static_cast<double>(2 * l + 1)) * rnl * rnl;
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


        // For each angular quantum number
        for (l = 0; l < Lmax; l++)
        {
            m_db->GetPath(path);
            path += (boost::format(".L=%lu") % static_cast<unsigned long>(l)).str();
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
