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
#include <memory>

namespace fem1d {
    class Approx : public util::Fun1D, private Mesh, private Lobatto
    {
    private:
        // Auxilary class
        class FunF2 : public util::Fun1D
        {
        public:
            FunF2(std::shared_ptr<const util::Fun1D> const & f, const Element& e, double fa, double fb) : m_f(f), m_e(e), m_fa(fa), m_fb(fb) { }
            virtual ~FunF2() { }
            virtual double Get(double r) const
            {
                double s = m_e.Xinv(r);
                double b0 = 0.5 * (1 - s), b1 = 0.5 * (1 + s);
                double fTilde = m_f->Get(r) - m_fa * b0 - m_fb * b1;
                return fTilde * fTilde;
            }
        public:
            const std::shared_ptr<const util::Fun1D> m_f;
            const Element& m_e;
            double m_fa, m_fb;
        };

    public:
        Approx();
        Approx(size_t M, std::shared_ptr<const util::Fun1D> f);
        virtual ~Approx(void);

        void Define(size_t M, std::shared_ptr<const util::Fun1D> f);
        virtual double Get(double x) const;

        // April 3rd, 2014 Added by dc1394
        virtual double GetDeriv(double r) const;
        virtual double Get2ndDeriv(double r) const;

        void SolveAdapt(double a, double b, double delta);

        void GetNode(std::vector<double>& node) const;

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
        Vec* m_b;

        // Searched approximation coefficients for one element
        Vec* m_c;

        // Searched approximation coefficients for all elements
        std::vector<double> m_cAll;

        // Matrix of system of equations
        ClpMtxBand* m_K;

        // Approximation degree
        size_t m_M;

        // Aproksymowana funkcja
        std::shared_ptr<const util::Fun1D> m_f;

        // Heap
        Heap<HeapElt> m_heap;
    };
}

#endif

