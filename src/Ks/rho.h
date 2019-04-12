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

namespace ks {
    class Rho final : public util::Fun1D
    {
    public:
        Rho(std::shared_ptr<const ParamDb> && db);
        ~Rho() override = default;

        double Get(double r) const override;
        double GetRhoTilde(double r) const;
        double GetDeriv(double r) const override;
        double Get2ndDeriv(double r) const override;

        void Calc(std::shared_ptr<const util::Fun1D> && f);

        template <util::Spin S>
        void Init();
                
        std::vector<double> GetNode() const;
        void Write() const;

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
        std::unique_ptr<const Int1DGauss> m_gauss;

    private:
        class RhoInit final : public util::Fun1D
        {
            RhoInit(RhoInit const &) = delete;
            RhoInit& operator=(RhoInit const &) = delete;

        public:
            RhoInit(double c, double alpha) : m_c(c), m_alpha(alpha) { }
            ~RhoInit() override = default;
            // 初期密度
            double Get(double r) const override
            {
                return r * r * m_c * exp(-m_alpha * r);
            }
        public:
            double m_c, m_alpha;
        };
    };
}

#endif
