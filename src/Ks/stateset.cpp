#include "stdafx.h"
#include "stateset.h"
#include <boost/algorithm/string/classification.hpp>    // for boost::is_any_of
#include <boost/algorithm/string/split.hpp>             // for boost::algorithm::split

namespace ks {
    //
    // The ground electronic configurations of the neutral atoms
    // http://physics.nist.gov/PhysRefData/DFTdata/configuration.html
    //
#define He "1s2"
#define Ne He" 2s2 2p6"
#define Ar Ne" 3s2 3p6"
#define Kr Ar" 3d10 4s2 4p6"
#define Xe Kr" 4d10 5s2 5p6"
#define Rn Xe" 4f14 5d10 6s2 6p6"
//
    const char* StateSet::m_atom[93] =
    {
        "EMPTY",
        "1s1",					// 01, H
        He,						// 02, He
        He" 2s1",				// 03, Li
        He" 2s2",				// 04, Be
        He" 2s2 2p1",			// 05, B
        He" 2s2 2p2",			// 06, C
        He" 2s2 2p3",			// 07, N
        He" 2s2 2p4",			// 08, O
        He" 2s2 2p5",			// 09, F

        Ne,						// 10, Ne
        Ne" 3s1",				// 11, Na
        Ne" 3s2",				// 12, Mg
        Ne" 3s2 3p1",			// 13, Al
        Ne" 3s2 3p2",			// 14, Si
        Ne" 3s2 3p3",			// 15, P
        Ne" 3s2 3p4",			// 16, S
        Ne" 3s2 3p5",			// 17, Cl
        Ar,						// 18, Ar
        Ar" 4s1",				// 19, K

        Ar" 4s2",				// 20, Ca
        Ar" 3d1 4s2",			// 21, Sc
        Ar" 3d2 4s2",			// 22, Ti
        Ar" 3d3 4s2",			// 23, V
        Ar" 3d5 4s1",			// 24, Cr
        Ar" 3d5 4s2",			// 25, Mn
        Ar" 3d6 4s2",			// 26, Fe
        Ar" 3d7 4s2",			// 27, Co
        Ar" 3d8 4s2",			// 28, Ni
        Ar" 3d10 4s1",			// 29, Cu

        Ar" 3d10 4s2",			// 30, Zn
        Ar" 3d10 4s2 4p1",		// 31, Ga
        Ar" 3d10 4s2 4p2",		// 32, Ge
        Ar" 3d10 4s2 4p3",		// 33, As
        Ar" 3d10 4s2 4p4",		// 34, Se
        Ar" 3d10 4s2 4p5",		// 35, Br
        Kr,						// 36, Kr
        Kr" 5s1",				// 37, Rb
        Kr" 5s2",				// 38, Sr
        Kr" 4d1 5s2",			// 39, Y

        Kr" 4d2 5s2",			// 40, Zr
        Kr" 4d4 5s1",			// 41, Nb
        Kr" 4d5 5s1",			// 42, Mo
        Kr" 4d5 5s2",			// 43, Tc
        Kr" 4d7 5s1",			// 44, Ru
        Kr" 4d8 5s1",			// 45, Rh
        Kr" 4d10",				// 46, Pd
        Kr" 4d10 5s1",			// 47, Ag
        Kr" 4d10 5s2",			// 48, Cd
        Kr" 4d10 5s2 5p1",		// 49, In

        Kr" 4d10 5s2 5p2",		// 50, Sn
        Kr" 4d10 5s2 5p3",		// 51, Sb
        Kr" 4d10 5s2 5p4",		// 52, Te
        Kr" 4d10 5s2 5p5",		// 53, I
        Xe,						// 54, Xe
        Xe" 6s1",				// 55, Cs
        Xe" 6s2",				// 56, Ba
        Xe" 5d1 6s2",			// 57, La
        Xe" 4f1 5d1 6s2",		// 58, Ce
        Xe" 4f3 6s2",			// 59, Pr

        Xe" 4f4 6s2",			// 60, Nd
        Xe" 4f5 6s2",			// 61, Pm
        Xe" 4f6 6s2",			// 62, Sm
        Xe" 4f7 6s2",			// 63, Eu
        Xe" 4f7 5d1 6s2",		// 64, Gd
        Xe" 4f9 6s2",			// 65, Tb
        Xe" 4f10 6s2",			// 66, Dy
        Xe" 4f11 6s2",			// 67, Ho
        Xe" 4f12 6s2",			// 68, Er
        Xe" 4f13 6s2",			// 69, Tm

        Xe" 4f14 6s2",			// 70, Yb
        Xe" 4f14 5d1 6s2",		// 71, Lu
        Xe" 4f14 5d2 6s2",		// 72, Hf
        Xe" 4f14 5d3 6s2",		// 73, Ta
        Xe" 4f14 5d4 6s2",		// 74, W
        Xe" 4f14 5d5 6s2",		// 75, Re
        Xe" 4f14 5d6 6s2",		// 76, Os
        Xe" 4f14 5d7 6s2",		// 77, Ir
        Xe" 4f14 5d9 6s1",		// 78, Pt
        Xe" 4f14 5d10 6s1",		// 79, Au

        Xe" 4f14 5d10 6s2",		// 80, Hg
        Xe" 4f14 5d10 6s2 6p1",	// 81, Tl
        Xe" 4f14 5d10 6s2 6p2",	// 82, Pb
        Xe" 4f14 5d10 6s2 6p3",	// 83, Bi
        Xe" 4f14 5d10 6s2 6p4",	// 84, Po
        Xe" 4f14 5d10 6s2 6p5",	// 85, At
        Rn,						// 86, Rn
        Rn" 7s1",				// 87, Fr
        Rn" 7s2",				// 88, Ra
        Rn" 6d1 7s2",			// 89, Ac

        Rn" 6d2 7s2",			// 90, Th
        Rn" 5f2 6d1 7s2",		// 91, Pa
        Rn" 5f3 6d1 7s2"		// 92, U
    };
#undef He
#undef Ne
#undef Ar
#undef Kr
#undef Xe
#undef Rn

    //
    // Constructor
    //
    StateSet::StateSet(std::shared_ptr<ParamDb> const & db) : m_db(db)
    {
        auto const proton = m_db->GetSize_t("Atom_Proton");
        auto const delim = " ";

        assert(proton <= 92);
        std::string config(m_atom[proton]);

        std::vector<std::string> tokens;
        boost::algorithm::split(tokens, config, boost::is_any_of(delim)); // ƒJƒ“ƒ}‚Å•ªŠ„
        for (auto && str : tokens) {
            push_back(State(str));
        }
    }

    //
    // Destructor
    //
    StateSet::~StateSet(void)
    {
    }

    //
    // Returns referenc for eigenstated defined by pair of quantum numbers(l, n)
    //
    const State& StateSet::Find(size_t l, size_t n) const
    {
        for (size_t i = 0; i < size(); ++i)
        {
            const State& s = (*this)[i];
            if ((s.L() == l) && (s.N() == n))
                return s;
        }

        assert(0); // Ther is no state (n, l)
        return front();
    }

    //
    // Returns occupation factor for state (l, n)
    //
    double StateSet::Occ(size_t l, size_t n) const
    {
        return Find(l, n).Occ();
    }

    //
    // Returns name of state (l, n)
    //
    std::string StateSet::Name(size_t l, size_t n) const
    {
        return Find(l, n).Name();
    }


    //
    // Returns maximal angular quantum number
    //
    size_t StateSet::GetLmax(void) const
    {
        size_t lMax = 0;

        for (size_t i = 0; i < size(); ++i)
        {
            const State& s = (*this)[i];
            if (s.L() > lMax)
                lMax = s.L();
        }
        return lMax + 1;
    }


    //
    // Returns maximum main quantum number of considered states
    // for given angular momentum number "l".
    //
    size_t StateSet::GetNmax(size_t l) const
    {
        size_t nMax = 0;

        for (size_t i = 0; i < size(); ++i)
        {
            const State& s = (*this)[i];
            if (s.L() != l)
                continue;
            if (s.N() > nMax)
                nMax = s.N();
        }
        return nMax + 1;
    }

    //
    // Returns energy of eigenstates including occupation factor and degenercy.
    //
    double StateSet::EigenEnerg(void) const
    {
        double ret = 0;

        for (size_t i = 0; i < size(); ++i)
        {
            const State& s = (*this)[i];
            ret += s.EigVal() * s.Occ();
        }

        return ret;
    }

    //
    // Returns sum of eigenvalues of states.
    // Non occupied einegstates are counted as well.
    //
    double StateSet::EigenSum(void) const
    {
        double ret = 0;

        for (size_t i = 0; i < size(); ++i)
        {
            const State& s = (*this)[i];
            ret += s.EigVal();
        }

        return ret;
    }


    //
    // Insters particular state (l, n, eigVal, occ) into set of states.
    //
    void StateSet::SetEigVal(size_t l, size_t n, double eigVal)
    {
        for (size_t i = 0; i < size(); ++i)
        {
            State& s = (*this)[i];
            if ((s.L() == l) && (s.N() == n))
            {
                s.EigVal(eigVal);
                return;
            }
        }

        // The state must be!
        assert(0);
    }

    //
    // Writes information about all states into file "out"
    //
    void StateSet::WriteSates(FILE* out) const
    {
        fprintf(out, "\n\n");
        fprintf(out, "===================================================================\n");
        fprintf(out, "     E I G E N V A L U E S\n");
        // fprintf(out, "+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +\n");
        fprintf(out, "-------------------------------------------------------------------\n");
        fprintf(out, "%16s %19s %20s\n", "State", "Value [Ha]", "Value [eV]");
        fprintf(out, "-------------------------------------------------------------------\n");
        for (size_t i = 0; i < size(); ++i)
        {
            const State& s = (*this)[i];

            fprintf(out, "   (n=%lu, L=%lu) %-5s %15.7lf Ha = %15.7lf eV\n",
                static_cast<unsigned long>(s.N()),
                static_cast<unsigned long>(s.L()),
                s.Name().c_str(),
                s.EigVal(),
                s.EigVal() * M_EV);
        }
        fprintf(out, "===================================================================\n");
    }
}
