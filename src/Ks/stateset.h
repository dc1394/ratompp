#ifndef __RATOM_STATESET_H__
#define __RATOM_STATESET_H__


/** \brief Stores infomation about eigenfunctions of atom
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/

#include "state.h"
#include "paramdb.h"
#include <memory>   // for std::shared_ptr

namespace ks {
    class StateSet : private std::vector<State>
    {
    public:
        StateSet(std::shared_ptr<ParamDb const> && db);
        ~StateSet() = default;

        std::size_t GetLmax(void) const;
        std::size_t GetNmax(std::size_t l) const;

        double Occ(std::size_t l, std::size_t n) const;
        std::string Name(std::size_t l, std::size_t n) const;

        double EigenEnerg(void) const;
        double EigenSum(void) const;

        void SetEigVal(std::size_t l, std::size_t n, double eigVal);

        void WriteSates(FILE* out) const;

    private:
        std::size_t Key(std::size_t l, std::size_t n) const;
        const State& Find(std::size_t l, std::size_t n) const;

    private:
        // Database of parameters
        std::shared_ptr<ParamDb const> const m_db;

        // The ground electronic configurations of the neutral elements
        static const char* m_atom[93];
    };
}

#endif

