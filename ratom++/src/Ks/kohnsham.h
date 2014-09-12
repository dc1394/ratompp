#ifndef __RATOM_KOHNSHAM_H__
#define __RATOM_KOHNSHAM_H__


/** \brief Linear Kohn-Sham equation
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/
// March 16th, 2014 Modified by dc1394

#include "paramdb.h"
#include "stateset.h"

#include "../Util/isspin.h"
#include <type_traits>

extern void* enabler;
// May 26th, 2014 Modified by dc1394
namespace ks {
    template <util::IsSpin T, typename X = void>
    using alpha_enabler = typename std::enable_if<T == util::IsSpin::Alpha, X>::type*&;
    template <util::IsSpin T, typename X = void>
    using beta_enabler = typename std::enable_if<T == util::IsSpin::Beta, X>::type*&;

    template <util::IsSpin Spin>
    class KohnSham : public util::Fun1D
    {
    private:
        // May 25th, 2014 Added by dc1394
        // Interaction potential used in radial Kohn-Sham equation
        class PotRad : public util::Fun1D
        {
            // May 25th, 2014 Added by dc1394
            std::shared_ptr<const util::Fun1D> m_pot;

        public:
            PotRad(void) : /*m_pot(NULL),*/ m_l(0) {}
            PotRad(std::shared_ptr<const util::Fun1D> const & pot) : m_pot(pot), m_l(0) {}
            virtual ~PotRad() { }

            virtual double Get(double r) const
            {
                assert(r > 0);
                return m_pot->Get(r) + m_l * (m_l + 1) / (2 * r * r);
            }
            std::shared_ptr<const util::Fun1D> & getPot()
            {
                return m_pot;
            }
        public:
            // May 25th, 2014 Modified by dc1394
            // Interaction potential
            //const util::Fun1D* m_pot;

            // Angular quantum number
            size_t m_l;
        };


    public:
        KohnSham(const ParamDb* db, StateSet* m_stateSet);
        ~KohnSham(void);

        // May 25th, 2014 Modified by dc1394
        //void Config(util::Fun1D* pot);
        void Config(const std::shared_ptr<const util::Fun1D> & pot);
        void Solve(void);
        double Get(double r) const;
        template <util::IsSpin S, alpha_enabler<S> T = enabler>
        double Get(double r) const;
        template <util::IsSpin S, beta_enabler<S> T = enabler>
        double Get(double r) const;
        void WriteEigen(void) const;

    private:
        void CalcOcc();

    private:
        // One solver for each angular quantum number L
        std::vector<EigProb> m_eigProb;

        // Effective interaction potential used in radial Kohn-Sham equation
        PotRad m_radPot;

        // Number of calculated eigenvalues for each angular quantum number L
        std::vector<size_t> m_eigNo;

        // database of parameters
        const ParamDb* m_db;

        // Set of eigenstates
        StateSet* m_stateSet;
    };
}
#endif

