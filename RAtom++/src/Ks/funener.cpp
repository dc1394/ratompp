#include "stdafx.h"
#include "funener.h"
#include <utility>

namespace ks {
    //
    // Constructor
    //
    FunEner::FunEner(std::pair<const std::shared_ptr<const Pot<util::IsSpin::Alpha>>,
                     const std::shared_ptr<const Pot<util::IsSpin::Beta>>> pot, size_t type)
        : m_pot(std::move(pot)), m_type(type)
    {
    }

    //
    // Destructor
    //
    FunEner::~FunEner(void)
    {
    }

    //
    // Returns for radius "r" values of integrand for evaluation of terms of total energy
    //
    double FunEner::Get(double r) const
    {
        double v;

        assert(r > 0);
        switch (m_type)
        {
        case 0: v = GetTotal(r);	break;
        case 1: v = GetNucleus(r);	break;
        case 2: v = GetHartree(r);	break;
        case 3: v = GetExch(r);		break;
        case 4: v = GetCorr(r);		break;
        case 5: v = GetKinetic(r);	break;
        default: assert(0); v = 0;	break;
        }
        return v;
    }


    //
    // Returns for radius "r" values of integrand for evaluation of total energy
    //
    double FunEner::GetTotal(double r) const
    {
        // αスピンとβスピンの密度の合計
        const double rho = m_pot.first->GetRho<util::IsSpin::Alpha>(r) +
                           m_pot.second->GetRho<util::IsSpin::Beta>(r);
        //const double rhoT = RhoTilde(r, rho);
        double ex, co, vh;

        // March 7th, 2014	Modified by dc1394 
        //ex = m_pot->m_exch->EdiffV(rhoT, 0);
        ex = m_pot.first->m_exch->EdiffV(r) + m_pot.second->m_exch->EdiffV(r);
        //co = m_pot->m_corr->EdiffV(rhoT, 0);
        co = m_pot.first->m_corr->EdiffV(r) + m_pot.second->m_corr->EdiffV(r);
        vh = m_pot.first->Vh(r);

        if ((vh - m_pot.second->Vh(r)) > 1.0E-12)
            throw std::runtime_error("Debug Error!");

        return (ex + co - 0.5 * vh) * rho;
    }

    //
    // Integrand for energy "Nucleus"
    //
    double FunEner::GetNucleus(double r) const
    {
        // αスピンとβスピンの密度の合計
        const double rho = m_pot.first->GetRho<util::IsSpin::Alpha>(r) +
            m_pot.second->GetRho<util::IsSpin::Beta>(r);

        const double vn = m_pot.first->Vn(r);
        if ((vn - m_pot.second->Vn(r)) > 1.0E-12)
            throw std::runtime_error("Debug Error!");

        return vn * rho;
    }

    //
    // Integrand for energy "Hartree"
    //
    double FunEner::GetHartree(double r) const
    {
        // αスピンとβスピンの密度の合計
        const double rho = m_pot.first->GetRho<util::IsSpin::Alpha>(r) +
            m_pot.second->GetRho<util::IsSpin::Beta>(r);
        const double vh = m_pot.first->Vh(r);

        if ((vh - m_pot.second->Vh(r)) > 1.0E-12)
            throw std::runtime_error("Debug Error!");
        return 0.5 * vh * rho;
    }

    //
    // Integrand for energy "Exchenge"
    //
    double FunEner::GetExch(double r) const
    {
        // αスピンとβスピンの密度の合計
        const double rho = m_pot.first->GetRho<util::IsSpin::Alpha>(r) +
            m_pot.second->GetRho<util::IsSpin::Beta>(r);
        //const double rhoT = RhoTilde(r, rho);

        const double ex = m_pot.first->Ex(r);
        if ((ex - m_pot.second->Ex(r)) > 1.0E-12)
            throw std::runtime_error("Debug Error!");

        // March 7th, 2014	Modified by dc1394 
        //return m_pot->Ex(rhoT, 0) * rho;
        return ex * rho;
    }

    //
    // Integrand for energy "Correlation"
    //
    double FunEner::GetCorr(double r) const
    {
        // αスピンとβスピンの密度の合計
        const double rho = m_pot.first->GetRho<util::IsSpin::Alpha>(r) +
            m_pot.second->GetRho<util::IsSpin::Beta>(r);

        const double ec = m_pot.first->Ec(r);
        if ((ec - m_pot.second->Ec(r)) > 1.0E-12)
            throw std::runtime_error("Debug Error!");

        //const double rhoT = RhoTilde(r, rho);
        // March 7th, 2014	Modified by dc1394 
        //return m_pot->Ec(rhoT, 0) * rho;
        return ec * rho;
    }

    //
    // Integrand for energy "Kinetic"
    //
    double FunEner::GetKinetic(double r) const
    {
        // αスピン密度
        const double rhoalpha = m_pot.first->GetRho<util::IsSpin::Alpha>(r);
        // βスピン密度
        const double rhobeta = m_pot.second->GetRho<util::IsSpin::Beta>(r);
        // αスピンとβスピンの密度の合計
        const double rho = rhoalpha + rhobeta;
        //const double rhoT = RhoTilde(r, rho);

        //return (m_pot->Vx(rhoT, 0) + m_pot->Vc(rhoT, 0) + m_pot->Vn(r) + m_pot->Vh(r)) * rho;

        // March 7th, 2014	Modified by dc1394 
        return (m_pot.first->Vx(r) + m_pot.first->Vc(r) +
                m_pot.second->Vx(r) + m_pot.second->Vc(r) +
                m_pot.first->Vn(r) + m_pot.first->Vh(r)) * rho;
    }

    //
    // Returns value of function \tilde{\rho}(r)
    //
    /*double FunEner::RhoTilde(double r, double rho) const
    {
        return rho / (M_4PI * r * r);
    }*/
}
