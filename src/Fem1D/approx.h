#ifndef __RATOM_APPROX_H__
#define __RATOM_APPROX_H__


/** \brief Solves approximation problem on interval [a, b] for function $f$ by Lobatto polynomials.
*
* \author Zbigniew Romanowski [romz]
*
*/
// March 16th, 2014 Modified by dc1394

#include "mesh.h"
#include "heap.h"
#include "heapelt.h"
#include <memory>       // for std::shared_ptr, std::unique_ptr 

namespace fem1d {
    class Approx final : public util::Fun1D, private Mesh, private Lobatto
    {
        // Auxilary class
        class FunF2 final : public util::Fun1D
        {
        public:
            FunF2(std::shared_ptr<util::Fun1D const> const & f, const Element& e, double fa, double fb) : m_f(f), m_e(e), m_fa(fa), m_fb(fb) { }
            ~FunF2() override = default;
            double Get(double r) const override
            {
                auto const s = m_e.Xinv(r);
                auto const b0 = 0.5 * (1 - s), b1 = 0.5 * (1 + s);
                auto const fTilde = m_f->Get(r) - m_fa * b0 - m_fb * b1;
                return fTilde * fTilde;
            }
        public:
            std::shared_ptr<util::Fun1D const> const m_f;
            Element const & m_e;
            double m_fa, m_fb;
        };

    public:
        Approx() = default;
        Approx(size_t M, std::shared_ptr<util::Fun1D const> && f);
        ~Approx() override = default;

        void Define(size_t M, std::shared_ptr<util::Fun1D const> && f);
        double Get(double x) const override;

        double GetDeriv(double r) const override;
        double Get2ndDeriv(double r) const override;

        void SolveAdapt(double a, double b, double delta);

        std::vector<double> GetNode() const;

        void WriteCoef(FILE* out) const;

    private:

        double Solve(double a, double b);
        void SolveAdaptOld(double a, double b, double delta, double h);


        double FindB(double a, double b, double delta, double h);
        double CalcB(const Element& elt, size_t i, double fa, double fb) const;
        double Ftilde(double x, double s, double fa, double fb) const;
        double IntegF2(const Element& elt, double fa, double fb) const;
        double RunBisect(double a, double ra, double rb, double delta);
        void SetMesh(const std::vector<double>& x);

        void SetCoef(double a, double b, std::vector<double>& c);

    private:

        // Right hand size of system of equations
        std::unique_ptr<Vec> m_b;

        // Searched approximation coefficients for one element
        std::unique_ptr<Vec> m_c;

        // Searched approximation coefficients for all elements
        std::vector<double> m_cAll;

        // Matrix of system of equations
        std::unique_ptr<ClpMtxBand> m_K;

        // Approximation degree
        std::size_t m_M = 0;

        // Aproksymowana funkcja
        std::shared_ptr<util::Fun1D const> m_f;

        // Heap
        Heap<HeapElt> m_heap;
    };
}

#endif

