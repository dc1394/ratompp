#include "stdafx.h"
#include "kohnsham.h"
#include <type_traits>

// May 25th, 2014 Added by dc1394
namespace ks {
    //
    // Constructor
    //
    template <util::Spin Spin>
    KohnSham<Spin>::KohnSham(const ParamDb* db, StateSet* stateSet) : m_db(db), m_stateSet(stateSet)
    {
    }

    //
    // Destructor
    //
    template <util::Spin Spin>
    KohnSham<Spin>::~KohnSham(void)
    {
    }

    //
    // Solver configurations.
    // Required parameters are properly set
    //
    //void KohnSham::Config(util::Fun1D* pot)
    template <util::Spin Spin>
    void KohnSham<Spin>::Config(std::shared_ptr<const util::Fun1D> const & pot)
    {
        const double rc = m_db->GetDouble("Atom_Rc");
        const size_t eigNode = m_db->GetSize_t("Solver_EigNode");
        const size_t eigDeg = m_db->GetSize_t("Solver_EigDeg");

        const double gamma = 0.5; // Paramter for radial Kohn-Sham equation
        const size_t Lmax = m_stateSet->GetLmax();
        size_t l;

        // May 25th, 2014 Modified by dc1394
        //m_radPot.m_pot = pot;
        m_radPot.getPot() = pot;

        m_eigProb.resize(Lmax);
        for (l = 0; l < Lmax; l++)
        {
            m_eigProb[l].Define(gamma, &m_radPot);
            m_eigProb[l].GenMeshLin(0, rc, eigNode, eigDeg);
        }

        m_eigNo.resize(Lmax);
        for (l = 0; l < Lmax; l++)
            m_eigNo[l] = m_stateSet->GetNmax(l);
    }

    //
    // Solves linear eqigenvalue problem
    //
    template <util::Spin Spin>
    void KohnSham<Spin>::Solve(void)
    {
        const double abstol = m_db->GetDouble("Solver_EigAbsTol");
        const double absMaxCoef = m_db->GetDouble("Solver_EigAbsMaxCoef");
        const bool adapt = m_db->GetBool("Solver_EigAdapt");

        double eigVal;

        for (size_t l = 0; l < m_eigProb.size(); l++)
        {
            m_radPot.m_l = l;

            if (adapt)
                m_eigProb[l].SolveAdapt(m_eigNo[l], abstol, absMaxCoef);
            else
                m_eigProb[l].Solve(m_eigNo[l], abstol);

            // Sets eigenvalues of states
            for (size_t n = 0; n < m_eigNo[l]; n++)
            {
                eigVal = m_eigProb[l].GetEigVal(n);
                m_stateSet->SetEigVal(l, n, eigVal);
            }
        }
    }

    // May 23rd, 2014 Modified by dc1394
    //
    // Returns value of electron density for radius "r"
    //
    template <util::Spin Spin>
    double KohnSham<Spin>::Get(double r) const
    {
        // �X�s���ɉ������֐����Ăяo��
        return Get<Spin>(r);

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

    template <util::Spin Spin>
    template <util::Spin S, alpha_enabler<S> T>
    double KohnSham<Spin>::Get(double r) const
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
                // occ�ɂ̓��X�s���̓d�q�̂�
                if (occ > 0.0 && occ <= static_cast<double>(Lmax)) {
                    rnl = m_eigProb[l].GetEigFun(n, r); // R_{n, \ell}(r)
                    // Lmax�܂ł̓��X�s���̓d�q����occ�Ɠ�����
                    rhoL += occ * rnl * rnl;
                }
                // occ�Ƀ��X�s���̓d�q���܂܂��
                else if (occ > static_cast<double>(Lmax)) {
                    rnl = m_eigProb[l].GetEigFun(n, r); // R_{n, \ell}(r)
                    // ���X�s���̓d�q����Lmax�Ɠ�����
                    rhoL += static_cast<double>(Lmax) * rnl * rnl;
                }
            }

            // Sum up all constituents
            rho += rhoL;
        }

        return rho;
    }

    template <util::Spin Spin>
    template <util::Spin S, beta_enabler<S> T>
    double KohnSham<Spin>::Get(double r) const
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
                // occ��Lmax���傫����΁Aocc�Ƀ��X�s���̓d�q���܂܂��
                if (occ > static_cast<double>(Lmax)) {
                    rnl = m_eigProb[l].GetEigFun(n, r); // R_{n, \ell}(r)
                    // ���X�s���̓d�q���́Aocc����Lmax����������������
                    rhoL += (occ - static_cast<double>(Lmax)) * rnl * rnl;
                }
            }

            // Sum up all constituents
            rho += rhoL;
        }

        return rho;
    }

    //
    // Writes eigenfunctions into file.
    // Lobato coeffictienst are written as well.
    //
    template <util::Spin Spin>
    void KohnSham<Spin>::WriteEigen(void) const
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
    template double KohnSham<util::Spin::Alpha>::Get<util::Spin::Alpha>(double r) const;
    template double KohnSham<util::Spin::Beta>::Get<util::Spin::Beta>(double r) const;
}
