#ifndef __RATOM_ENERGY_H__
#define __RATOM_ENERGY_H__


/** \brief Evaluates of terms of total energy
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/


#include "pot.h"
#include "paramdb.h"
#include "stateset.h"

namespace ks {
    class Energy final
    {
    public:
        Energy(std::pair< std::shared_ptr<Pot<util::Spin::Alpha> const> const,
               std::shared_ptr<Pot<util::Spin::Beta> const> const> && pot,
               std::shared_ptr<StateSet const> && ss,
               std::shared_ptr<ParamDb const> const & db);
        ~Energy() = default;

        void SetNode(std::vector<double> const & node);
        void WriteEnergy(FILE* out);

    private:
        double Total() const;
        double Nucleus() const;
        double Hartree() const;
        double Exch() const;
        double Corr() const;
        double Kinetic() const;

        double Integ(size_t type) const;

    private:
        // Gaussa quadratures
        std::shared_ptr<Int1DGauss const> const m_gauss;

        // Coordinaes of nodes for approximation of electron density
        std::vector<double> m_node;

        // DFT potencial
        std::pair<std::shared_ptr<Pot<util::Spin::Alpha> const> const,
                  std::shared_ptr<Pot<util::Spin::Beta> const> const> const m_pot;

        // Integration is on the interval [0, rc]
        double const m_rc;

        // Atomic eigenstates
        std::shared_ptr<StateSet const> const m_ss;
    };
}

#endif

