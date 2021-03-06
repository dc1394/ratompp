#include "stdafx.h"
#include "rho.h"
#include <iostream>         // for std::cout
#include <boost/format.hpp> // for boost::format

namespace ks {
    //
    // Constructor
    //
    Rho::Rho(std::shared_ptr<ParamDb const> && db) : m_db(std::move(db))
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
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n";
        //	printf("+  Number of electrons = %.2lf\n", N);
        std::cout << boost::format("+  Number of protons (electrons) = %.2lf\n") % M;
        std::cout << boost::format("+  Applied 'Rho0' gives %.6lf electrons\n") % Integ();
        std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
    }

    template <>
    void Rho::Init<util::Spin::Beta>()
    {
        Calc(std::make_shared<const RhoInit>(0.0, 0.0));
    }

    //
    // returns integral
    //
    // \int \rho(r) d \vec{r} = \int_0^{\infty} \rho(r) dr
    //
    double Rho::Integ() const
    {
        std::vector<double> node(GetNode());

        auto val = 0.0;
        auto const loop = node.size() - 1;
        for (std::size_t i = 0UL; i < loop; ++i)
            val += m_gauss->Calc(*this, node[i], node[i + 1]);

        return val;
    }

    //
    // Calculates approximation of electron density based on function "f"
    //
    void Rho::Calc(std::shared_ptr<util::Fun1D const> && f)
    {
        auto const rc = m_db->GetDouble("Atom_Rc");
        auto const rhoDeg = m_db->GetSize_t("Rho_Deg");
        auto const rhoDelta = m_db->GetDouble("Rho_Delta");

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
    void Rho::Write() const
    {
        auto const outRhoNode = m_db->GetSize_t("Out_RhoNode");
        
        {
            // Nodes coordinates and approximation coefficients for electron density
            auto out = m_db->OpenFile("rho.app", "wt");
            std::fprintf(out.get(), "#\n");
            std::fprintf(out.get(), "# Orders of intervals is random,\n");
            std::fprintf(out.get(), "# You must sord on column 'r_i', in order to get ordered intervals.\n");
            std::fprintf(out.get(), "#\n");
            m_approx.WriteCoef(out.get());
        }

        {
            // Storing values of electron densities
            std::vector<double> node(GetNode());
            auto out = m_db->OpenFile("rho", "wt");
            //fprintf(out, "%16s \t %16s \t %16s \n", "R", "Rho", "RhoTilde");
            for (std::size_t i = 0UL; i < node.size() - 1; ++i)
            {
                auto const dr = (node[i + 1] - node[i]) / outRhoNode;
                auto r = node[i];
                for (auto k = 0U; k < outRhoNode; ++k)
                {
                    std::fprintf(out.get(), "%.15f,%.15f,%.15f\n", r, Get(r), GetRhoTilde(r));
                    r += dr;
                }
            }
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
