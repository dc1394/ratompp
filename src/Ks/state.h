#ifndef __RATOM_STATE_H__
#define __RATOM_STATE_H__

/** \brief State of the atom
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/

namespace ks {
    class State final
    {
     public:
        State() = default;
        explicit State(std::string const & name);
        ~State() = default;

        std::string Name() const { return m_name; };
        std::size_t Occ() const { return m_occ; }
        std::size_t N() const { return m_n; }
        std::size_t L() const { return m_l; }
        double EigVal() const { return m_eigVal; }

        void EigVal(double val) { m_eigVal = val; }

        bool operator< (const State& s) const;

    private:
        std::size_t GetL(char ll) const;

    private:
        // Main quantum number
        std::size_t m_n = 0;

        // Angular quantum number
        std::size_t m_l = 0;

        // Eigenvalue (energy of state)
        double m_eigVal = 0.0;

        // Number of electrons on the state
        std::size_t m_occ = 0;
        
        // Name of the state
        std::string m_name;
    };
}

#endif

