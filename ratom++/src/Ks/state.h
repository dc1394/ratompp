#ifndef __RATOM_STATE_H__
#define __RATOM_STATE_H__

#include "../Util/spin.h"
#include <utility>

/** \brief State of the atom
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/

namespace ks {
    class State final
    {
        State(State const &) = delete;
        State & operator=(State const &) = delete;

    public:
        State(void);
        //explicit State(const char* name);
        explicit State(std::string const & name);
        ~State(void);

        std::string Name(void) const { return m_name; };
        size_t Occ(void) const { return m_occ; }
        size_t N(void) const { return m_n; }
        size_t L(void) const { return m_l; }
        double EigVal(void) const { return m_eigVal; }

        void EigVal(double val) { m_eigVal = val; }

        bool operator< (const State& s) const;

    private:
        size_t GetL(char ll) const;

    private:
        // Main quantum number
        size_t m_n;

        // Angular quantum number
        size_t m_l;

        // Eigenvalue (energy of state)
        double m_eigVal;

        // Number of electrons on the state
        //size_t m_occ;
        std::pair<std::size_t, std::size_t> m_occ;
        
        // Name of the state
        std::string m_name;
    };
}

#endif

