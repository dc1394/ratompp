#include "stdafx.h"
#include "funener.h"
#include <utility>  // for std::move

namespace ks {
    //
    // Constructor
    //
    FunEner::FunEner(std::pair< std::shared_ptr< Pot<util::Spin::Alpha> const> const,
        std::shared_ptr< Pot<util::Spin::Beta> const> const> const & pot, std::size_t type)
        : m_pot(std::move(pot)), m_type(type)
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
        // αスピン密度
        auto const rhoalpha = m_pot.first->GetRho(r);
        // βスピン密度
        auto const rhobeta = m_pot.second->GetRho(r);
        //const double rhoT = RhoTilde(r, rho);

        // March 7th, 2014	Modified by dc1394 
        //ex = m_pot->m_exch->EdiffV(rhoT, 0);
        auto const exalpha = m_pot.first->Exch()->EdiffV(r);
        auto const exbeta = m_pot.second->Exch()->EdiffV(r);
        //co = m_pot->m_corr->EdiffV(rhoT, 0);
        auto const coalpha = m_pot.first->Corr()->EdiffV(r);
        auto const cobeta = m_pot.second->Corr()->EdiffV(r);
        auto const vh = m_pot.first->Vh(r);

        return (exalpha + coalpha - 0.5 * vh) * rhoalpha + (exbeta + cobeta - 0.5 * vh) * rhobeta;
    }

    //
    // Integrand for energy "Nucleus"
    //
    double FunEner::GetNucleus(double r) const
    {
        // αスピンとβスピンの密度の合計
        auto const rho = m_pot.first->GetRho(r) + m_pot.second->GetRho(r);

        auto const vn = m_pot.first->Vn(r);

        return vn * rho;
    }

    //
    // Integrand for energy "Hartree"
    //
    double FunEner::GetHartree(double r) const
    {
        // αスピンとβスピンの密度の合計
        auto const rho = m_pot.first->GetRho(r) + m_pot.second->GetRho(r);
        auto const vh = m_pot.first->Vh(r);

        return 0.5 * vh * rho;
    }

    //
    // Integrand for energy "Exchenge"
    //
    double FunEner::GetExch(double r) const
    {
        // αスピン密度
        auto const rhoalpha = m_pot.first->GetRho(r);
        // βスピン密度
        auto const rhobeta = m_pot.second->GetRho(r);
        //const double rhoT = RhoTilde(r, rho);

        auto const exalpha = m_pot.first->Ex(r);
        auto const exbeta = m_pot.second->Ex(r);

        // March 7th, 2014	Modified by dc1394 
        //return m_pot->Ex(rhoT, 0) * rho;
        return exalpha * rhoalpha + exbeta * rhobeta;
    }

    //
    // Integrand for energy "Correlation"
    //
    double FunEner::GetCorr(double r) const
    {
        // αスピン密度
        auto const rhoalpha = m_pot.first->GetRho(r);
        // βスピン密度
        auto const rhobeta = m_pot.second->GetRho(r);
        //const double rhoT = RhoTilde(r, rho);

        auto const ecalpha = m_pot.first->Ec(r);
        auto const ecbeta = m_pot.second->Ec(r);

        //const double rhoT = RhoTilde(r, rho);
        // March 7th, 2014	Modified by dc1394 
        //return m_pot->Ec(rhoT, 0) * rho;
        return ecalpha * rhoalpha + ecbeta * rhobeta;
    }

    //
    // Integrand for energy "Kinetic"
    //
    double FunEner::GetKinetic(double r) const
    {
        // αスピン密度
        auto const rhoalpha = m_pot.first->GetRho(r);
        // βスピン密度
        auto const rhobeta = m_pot.second->GetRho(r);
        // αスピンとβスピンの密度の合計
        auto const rho = rhoalpha + rhobeta;
        //const double rhoT = RhoTilde(r, rho);

        //return (m_pot->Vx(rhoT, 0) + m_pot->Vc(rhoT, 0) + m_pot->Vn(r) + m_pot->Vh(r)) * rho;

        return (m_pot.first->Vx(r) + m_pot.first->Vc(r)) * rhoalpha +
               (m_pot.second->Vx(r) + m_pot.second->Vc(r)) * rhobeta +
               (m_pot.first->Vn(r) + m_pot.first->Vh(r)) * rho;
    }
}
