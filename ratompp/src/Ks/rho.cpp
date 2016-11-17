#include "stdafx.h"
#include "rho.h"

// March 18th, 2014 Added by dc1394
#include <utility>

namespace ks {
    //
    // Constructor
    //
    Rho::Rho(std::shared_ptr<const ParamDb> const & db) : m_db(db)
    {
        auto const rhoDeg = m_db->GetSize_t("Rho_Deg");

        m_gauss = std::make_unique<Int1DGauss>(2 * rhoDeg);
    }


    //
    // Initialization of electron density
    //
    template <>
    void Rho::Init<util::Spin::Alpha>(void)
    {
        const bool def = m_db->GetBool("Rho0_Default");
        const double M = m_db->GetDouble("Atom_Proton");
        // const double N = m_db->GetDouble("Atom_Electron");
        double c, alpha;

        if (def)
        {
            double w;
            if (M < 50)
                w = M;
            else
                w = 50;
            // Default inicjalization
            c = w * w * w * w / 16;
            alpha = 0.5 * w;
            c *= (M / w);
        }
        else
        {
            // User defined initialization
            c = util::HelpFun::M_4PI * m_db->GetDouble("Rho0_c");
            alpha = m_db->GetDouble("Rho0_Alpha");
        }

        //std::shared_ptr<const RhoInit> pf(std::make_shared<const RhoInit>(c, alpha));
        Calc(std::make_shared<const RhoInit>(c, alpha));

        //
        // After initialization integral of electron densities shuld be equal to number of electrons
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n");
        //	printf("+  Number of electrons = %.2lf\n", N);
        printf("+  Number of protons (electrons) = %.2lf\n", M);
        printf("+  Applied 'Rho0' gives %.6lf electrons\n", Integ());
        printf("+++++++++++++++++++++++++++++++++++++++++++++++++++\n\n");
    }

    template <>
    void Rho::Init<util::Spin::Beta>(void)
    {
        const bool def = m_db->GetBool("Rho0_Default");
        const double M = m_db->GetDouble("Atom_Proton");
        // const double N = m_db->GetDouble("Atom_Electron");
        double c, alpha;

        if (def)
        {
            double w;
            if (M < 50)
                w = M;
            else
                w = 50;
            // Default inicjalization
            c = w * w * w * w / 16;
            alpha = 0.5 * w;
            c *= (M / w);
        }
        else
        {
            // User defined initialization
            c = util::HelpFun::M_4PI * m_db->GetDouble("Rho0_c");
            alpha = m_db->GetDouble("Rho0_Alpha");
        }

        //std::shared_ptr<const RhoInit> pf(std::make_shared<const RhoInit>(c, alpha));
        Calc(std::make_shared<const RhoInit>(c, alpha));
    }

    //
    // returns integral
    //
    // \int \rho(r) d \vec{r} = \int_0^{\infty} \rho(r) dr
    //
    double Rho::Integ(void) const
    {
        std::vector<double> node(GetNode());

        double val = 0.0;
        auto const loop = node.size() - 1;
        for (auto i = 0U; i < loop; ++i)
            val += m_gauss->Calc(*this, node[i], node[i + 1]);

        return val;
    }

    //
    // Calculates approximation of electron density based on function "f"
    //
    void Rho::Calc(std::shared_ptr<const util::Fun1D> && f)
    {
        const double rc = m_db->GetDouble("Atom_Rc");
        const size_t rhoDeg = m_db->GetSize_t("Rho_Deg");
        const double rhoDelta = m_db->GetDouble("Rho_Delta");

        m_approx.Define(rhoDeg, std::move(f));
        m_approx.SolveAdapt(0, rc, rhoDelta);

# ifdef _DEBUG
        printf("\n\nIntegrated density gives %.6lf electrons\n\n", Integ());
# endif
    }

    //
    // Returns electron density for radius "r"
    //
    double Rho::Get(double r) const
    {
        // March 19th, 2014 Modified by dc1394
        //# ifdef _DEBUG
        //	const double rc = m_db->GetDouble("Atom_Rc");
        //	assert(r <= rc);
        //# endif

        const double v = m_approx.Get(r);

        // Sometimes, because of approximation, the approximated value is less than zero,
        // what has no phisical meaning
        return std::max(v, 0.0);
    }


    //
    // returns value of helper function $\tilde{\rho}(r)$
    //
    double Rho::GetRhoTilde(double r) const
    {
        const double rho = Get(r);

        if (r > 0)
            return rho / (util::HelpFun::M_4PI * r * r);

        // Linear extrapolation for (r == 0)
        const double eps = 1E-12;
        double v1, v2;

        v1 = GetRhoTilde(eps);
        v2 = GetRhoTilde(2 * eps);
        return 2 * v1 - v2;
    }

    //
    // Returns nodes used for perfoming the approximation
    //
    std::vector<double> Rho::GetNode() const
    {
        return m_approx.GetNode();
    }
    
    //
    // Writes calculated electron density for each state int file
    //
    void Rho::Write(void) const
    {
        const size_t outRhoNode = m_db->GetSize_t("Out_RhoNode");
        FILE* out;
        double dr, r;

        {
            // Nodes coordinates and approximation coefficients for electron density
            out = m_db->OpenFile("rho.app", "wt");
            fprintf(out, "#\n");
            fprintf(out, "# Orders of intervals is random,\n");
            fprintf(out, "# You must sord on column 'r_i', in order to get ordered intervals.\n");
            fprintf(out, "#\n");
            m_approx.WriteCoef(out);
            fclose(out);
        }

    {
        // Storing values of electron densities
        std::vector<double> node(GetNode());
        out = m_db->OpenFile("rho", "wt");
        //fprintf(out, "%16s \t %16s \t %16s \n", "R", "Rho", "RhoTilde");
        for (auto i = 0U; i < node.size() - 1; ++i)
        {
            dr = (node[i + 1] - node[i]) / outRhoNode;
            r = node[i];
            for (auto k = 0U; k < outRhoNode; ++k)
            {
                fprintf(out, "%.15f,%.15f,%.15f\n", r, Get(r), GetRhoTilde(r));
                r += dr;
            }
        }
        fclose(out);
    }
    }

    // April 3rd, 2014 Added by dc1394
    // Returns value of electron density (1st derivative) for radius "r"
    //
    double Rho::GetDeriv(double r) const
    {
        return m_approx.GetDeriv(r);
    }



    // April 3rd, 2014 Added by dc1394
    // Returns value of electron density (2nd derivative) for radius "r"
    //
    double Rho::Get2ndDeriv(double r) const
    {
        return m_approx.Get2ndDeriv(r);
    }
}
