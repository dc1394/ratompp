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
        Energy(std::pair<const std::shared_ptr<const Pot<util::Spin::Alpha>>,
               const std::shared_ptr<const Pot<util::Spin::Beta>>> pot, const StateSet* ss, const ParamDb* db);
        virtual ~Energy(void);

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
        const std::pair<const std::shared_ptr<const Pot<util::Spin::Alpha>>,
                        const std::shared_ptr<const Pot<util::Spin::Beta>>> m_pot;

        // Atomic eigenstates
        const StateSet* m_ss;

        // Integration is on the interval [0, rc]
        double m_rc;

        // Gaussa quadratures
        Int1DGauss *m_gauss;

        // Coordinaes of nodes for approximation of electron density
        std::vector<double> m_node;
    };
}

#endif

