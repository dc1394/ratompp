#include "stdafx.h"
#include "energy.h"
#include "funener.h"
#include <utility>

namespace ks {
    //
    // Construktor
    //
    Energy::Energy(std::pair<const std::shared_ptr<const Pot<util::Spin::Alpha>>,
                   const std::shared_ptr<const Pot<util::Spin::Beta>>> pot,
                   const StateSet* ss, const ParamDb* db) :
        m_pot(std::move(pot)), m_ss(ss)
    {
        m_rc = atof(db->Get("Atom_Rc"));
        m_gauss = new Int1DGauss(3 * db->GetLong("Rho_Deg"));
    }

    //
    // Destructor
    //
    Energy::~Energy(void)
    {
        delete m_gauss;
    }

    //
    // Sets array of nodes needed for approximation of electron density
    // This array is used for numerical integration
    //
    void Energy::SetNode(const std::vector<double>& node)
    {
        m_node = node;

# ifdef _DEBUG
        printf("\nRho mesh:\n");
        for (size_t i = 0; i < m_node.size(); ++i)
            printf("i=%d  node=%f\n", i, m_node[i]);
# endif
    }


    //
    // Returns total energy of system
    //
    double Energy::Total()
    {
        return m_ss->EigenEnerg() + Integ(0);
    }

    //
    // Returns interaction energy  of electrons with atomic core.
    //
    double Energy::Nucleus()
    {
        return Integ(1);
    }

    //
    // Returns Hartree energy
    //
    double Energy::Hartree()
    {
        return Integ(2);
    }

    //
    // Returns exchange energy
    //
    double Energy::Exch()
    {
        return Integ(3);
    }

    //
    // Returns correlation energy
    //
    double Energy::Corr()
    {
        return Integ(4);
    }

    //
    // Returns kinetic energy
    //
    double Energy::Kinetic()
    {
        return m_ss->EigenEnerg() - Integ(5);
    }


    //
    // Returns integral of type "type"
    //
    double Energy::Integ(size_t type)
    {
        // const double absTol = 1E-14;
        FunEner ener(m_pot, type);
        double val = 0;

        for (size_t i = 0; i < m_node.size() - 1; ++i)
        {
            val += m_gauss->Calc(ener, m_node[i], m_node[i + 1]);
            // val += m_gauss->Adapt(ener, m_node[i], m_node[i + 1], absTol);
        }
        return val;
    }

    //
    // Writes information about energy int file
    //
    void Energy::WriteEnergy(FILE* out)
    {
        double v, ex, ec;

        fprintf(out, "\n\n");
        fprintf(out, "==========================================================\n");
        fprintf(out, "     E N E R G Y \n");
        fprintf(out, "----------------------------------------------------------\n");

        v = Total();
        fprintf(out, "\t Etot   = %15.7lf Ha =  %13.4lf eV\n", v, v * M_EV);

        v = Kinetic();
        fprintf(out, "\t Ekin   = %15.7lf Ha =  %13.4lf eV\n", v, v * M_EV);

        v = Hartree();
        fprintf(out, "\t Ecoul  = %15.7lf Ha =  %13.4lf eV\n", v, v * M_EV);

        v = Nucleus();
        fprintf(out, "\t Eenucl = %15.7lf Ha =  %13.4lf eV\n", v, v * M_EV);

        ex = Exch();
        fprintf(out, "\t Eexch  = %15.7lf Ha =  %13.4lf eV\n", ex, ex * M_EV);

        ec = Corr();
        fprintf(out, "\t Ecorr  = %15.7lf Ha =  %13.4lf eV\n", ec, ec * M_EV);

        fprintf(out, "\t Exc    = %15.7lf Ha =  %13.4lf eV\n", ex + ec, (ex + ec) * M_EV);

        fprintf(out, "==========================================================\n");
    }
}
