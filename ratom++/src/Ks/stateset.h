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
        StateSet(std::shared_ptr<ParamDb> const & db);
        ~StateSet(void);

        size_t GetLmax(void) const;
        size_t GetNmax(size_t l) const;

        double Occ(size_t l, size_t n) const;
        std::string Name(size_t l, size_t n) const;

        double EigenEnerg(void) const;
        double EigenSum(void) const;

        void SetEigVal(size_t l, size_t n, double eigVal);

        void WriteSates(FILE* out) const;

    private:
        size_t Key(size_t l, size_t n) const;
        const State& Find(size_t l, size_t n) const;

    private:
        // Database of parameters
        std::shared_ptr<ParamDb> const m_db;

        // The ground electronic configurations of the neutral elements
        static const char* m_atom[93];
    };
}

#endif

