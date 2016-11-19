#ifndef __RATOM_RHO_H__
#define __RATOM_RHO_H__


/** \brief Electron density of atom. It is represented as a piecewise polynomial function.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/
// March 16th, 2014 Modified by dc1394

#include "../Util/spin.h"
#include "paramdb.h"
#include <memory>       // for std::shared_ptr
#include <boost/mpl/int.hpp>

namespace ks {
    class Rho : public util::Fun1D
    {
    public:
        Rho(std::shared_ptr<const ParamDb> const & db);
        virtual ~Rho() = default;

        virtual double Get(double r) const;
        double GetRhoTilde(double r) const;
        // April 3rd, 2014 Added by dc1394
        virtual double GetDeriv(double r) const;
        virtual double Get2ndDeriv(double r) const;

        void Calc(std::shared_ptr<const util::Fun1D> && f);

        template <util::Spin S>
        void Init();
                
        std::vector<double> GetNode() const;
        void Write(void) const;

    private:
        double Integ() const;

    private:
        //! A private member variable.
        /*!
            Function approximation
        */
        fem1d::Approx m_approx;

        //! A private member variable.
        /*!
            Database of parameters
        */
        std::shared_ptr<const ParamDb> const m_db;

        //! A private member variable.
        /*!
            Number of calculated eigenvalues for each angular quantum number L
        */
        std::vector<size_t> m_eigNo;

        // Gauss quadratures
        std::unique_ptr<Int1DGauss> m_gauss;

    private:
        class RhoInit : public util::Fun1D
        {
            RhoInit(RhoInit const &) = delete;
            RhoInit& operator=(RhoInit const &) = delete;

        public:
            RhoInit(double c, double alpha) : m_c(c), m_alpha(alpha) { }
            virtual ~RhoInit() { }
            // èâä˙ñßìx
            virtual double Get(double r) const
            {
                return r * r * m_c * exp(-m_alpha * r);
            }
        public:
            double m_c, m_alpha;
        };
    };
}

#endif

