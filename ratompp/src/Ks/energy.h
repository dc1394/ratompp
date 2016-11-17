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
        Energy(std::pair<std::shared_ptr<const Pot<util::Spin::Alpha>>,
               std::shared_ptr<const Pot<util::Spin::Beta>>> const & pot,
               std::shared_ptr<const StateSet> const & ss,
               std::shared_ptr<const ParamDb> const & db);
        virtual ~Energy() = default;

        void SetNode(const std::vector<double>& node);
        void WriteEnergy(FILE* out);

    private:
        double Total();
        double Nucleus();
        double Hartree();
        double Exch();
        double Corr();
        double Kinetic();

        double Integ(size_t type);

    private:
        // DFT potencial
        std::pair<std::shared_ptr<const Pot<util::Spin::Alpha>>,
                  std::shared_ptr<const Pot<util::Spin::Beta>>> const m_pot;

        // Atomic eigenstates
        std::shared_ptr<const StateSet> const m_ss;

        // Integration is on the interval [0, rc]
        double m_rc;

        // Gaussa quadratures
        std::shared_ptr<Int1DGauss> m_gauss;

        // Coordinaes of nodes for approximation of electron density
        std::vector<double> m_node;
    };
}

#endif
