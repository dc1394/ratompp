#ifndef __RATOM_FUNENER_H__
#define __RATOM_FUNENER_H__


/** \brief Helper function used to evaluated terms of total energy
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/
// March 16th, 2014 Modified by dc1394

#include "pot.h"

namespace ks {
    class FunEner final : public util::Fun1D
    {
    public:
        FunEner(std::pair<std::shared_ptr<Pot<util::Spin::Alpha>>,
                std::shared_ptr<Pot<util::Spin::Beta>>> const & pot, size_t type);
        virtual ~FunEner() = default;

        virtual double Get(double r) const;

    private:
        double GetTotal(double r) const;
        double GetNucleus(double r) const;
        double GetHartree(double r) const;
        double GetExch(double r) const;
        double GetCorr(double r) const;
        double GetKinetic(double r) const;

        double RhoTilde(double r, double rho) const;

    private:
        // Interaction potencial
        const std::pair<const std::shared_ptr<const Pot<util::Spin::Alpha>>,
                        const std::shared_ptr<const Pot<util::Spin::Beta>>> m_pot;

        // Flag defining term of evaluated total energy
        size_t m_type;
    };
}

#endif

