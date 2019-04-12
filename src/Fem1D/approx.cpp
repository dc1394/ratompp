#include "stdafx.h"
#include "approx.h"
#include <utility>  // for std::move

namespace fem1d {
    //
    // Constructor
    //
    Approx::Approx()
    {
        m_M = 0;
    }

    //
    // Constructor
    //
    Approx::Approx(size_t M, std::shared_ptr<const util::Fun1D> && f)
    {
        Define(M, std::move(f));
    }
    
    //
    // Defines approximation problem
    //
    void Approx::Define(size_t M, std::shared_ptr<const util::Fun1D> && f)
    {
        const size_t DIAG = 2; // Number of super diagonals

        m_f = std::move(f);

        assert(M >= 2);

        // Basis functions wit indexes $i=0$, $i=1$ are skipped
        // Only basis functions with indexes $2, 3, ... , M$ are considered.
        M = M - 1;
        m_M = M;

        m_b.reset();
        m_b = std::make_unique<Vec>(M);

        m_c.reset();
        m_c = std::make_unique<Vec>(M);

        m_K.reset();
        m_K = std::make_unique<ClpMtxBand>(M, DIAG, 0);

        size_t i, j;
        for (i = 0; i < M; i++)
        {
            const size_t jMax = std::min(i + DIAG, M - 1);
            for (j = i; j <= jMax; j++)
                m_K->Set(i, j) = Memi(i + 2, j + 2); // In order to keep consistency "2" must be added
        }
    }



    //
    // Solves approximation problem on interval [a,b].
    // Returns approximation error.
    //
    double Approx::Solve(double a, double b)
    {
        assert(b > a);
        size_t i;
        Element elt;
        const double fa = m_f->Get(a);
        const double fb = m_f->Get(b);

        // Define the element. Order "0" is not used!
        elt.Set(a, b, 0);

        for (i = 0; i < m_M; i++)
            m_b->Set(i) = CalcB(elt, i + 2, fa, fb);

        m_K->SolveSymPos(m_b, m_c);

        //	m_K->Write("K.mtx");
        //	m_b->Write("b.vec");
        //	m_c->Write("c.vec");

        double kij, sum = 0;
        for (i = 0; i < m_M; i++)
        {
            for (size_t j = 0; j < m_M; j++)
            {
                // Upper triangular matrix is stored only
                if (j >= i)
                    kij = m_K->Get(i, j);
                else
                    kij = m_K->Get(j, i);
                sum += m_c->Get(i) * m_c->Get(j) * kij;
            }
        }

        return ::sqrt(::fabs(IntegF2(elt, fa, fb) - elt.Jac() * sum));
    }

    //
    // Returns approximated value for "a <= x <= b".
    //
    double Approx::Get(double x) const
    {
        /*
            // The equation must be solved!
            assert(m_c);

            assert(IsInRange(x));
            const size_t n = FindElt(x);
            const Element& e = Elt(n);

            // s - wspolrzedna lokalna na lemencie "e"
            const double s = min(max(e.Xinv(x), -1), 1); // MIN, MAX - To avoid the rounding errors

            // Sum over all basis function with support on the element $e$
            double val = 0;

            const size_t m = n * (m_M + 2); // Numer wspolczynnika
            for(size_t i = 0; i <= m_M; i++)
            {
            val += m_cAll[m + i] * Basis(i, s);
            }

            //	{
            //	// Wartosci na koncach elementu "e"
            //	const double fa = m_f->Get(e.X(-1)); // lewy koniec
            //	const double fb = m_f->Get(e.X( 1)); // prawy koniec
            //	val += fa * Basis(0, e.Xinv(x)) + fb * Basis(1, e.Xinv(x));
            //	}
            //

            return val;
            */

        // Sum over all basis function with support on the element $e$




        for (size_t i = 0; i < m_heap.Size(); i++)
        {
            const HeapElt& e = m_heap[i];
            if (e.m_left <= x && x <= e.m_right)
            {
                double val = 0;
                // s - local variable for element "e"
                double s = e.Xinv(x);
                s = std::min(std::max(s, -1.0), 1.0); // MIN, MAX - To avoid the rounding errors

                for (size_t j = 0; j < m_M + 2; j++)
                    val += e.m_coef[j] * Basis(j, s);

                return val;
            }
        }

        assert(0);
        return 0;

    }

    // March 16th, 2014 Added by dc1394
    // Returns approximated value (1st derivative) for "a <= x <= b".
    //
    double Approx::GetDeriv(double x) const
    {
        for (size_t i = 0; i < m_heap.Size(); i++)
        {
            const HeapElt& e = m_heap[i];

            if (e.m_left <= x && x <= e.m_right)
            {
                double val = 0;
                // s - local variable for element "e"
                double s = e.Xinv(x);
                s = std::min(std::max(s, -1.0), 1.0); // MIN, MAX - To avoid the rounding errors

                for (size_t j = 0; j < m_M + 2; j++)
                    val += e.m_coef[j] * BasisD1(j, s);

                return val / e.Getc2();
            }
        }

        assert(0);
        return 0;

    }

    // March 16th, 2014 Added by dc1394
    // Returns approximated value (2nd derivative) for "a <= x <= b". 
    //
    double Approx::Get2ndDeriv(double x) const
    {
        for (size_t i = 0; i < m_heap.Size(); i++)
        {
            const HeapElt& e = m_heap[i];

            if (e.m_left <= x && x <= e.m_right)
            {
                double val = 0;
                // s - local variable for element "e"
                double s = e.Xinv(x);
                s = std::min(std::max(s, -1.0), 1.0); // MIN, MAX - To avoid the rounding errors

                for (size_t j = 0; j < m_M + 2; j++)
                    val += e.m_coef[j] * BasisD2(j, s);

                return val / (e.Getc2() * e.Getc2());
            }
        }

        assert(0);
        return 0;

    }

    //
    // Returns value $b[i]$. Gauss quadratures are applied.
    // 
    double Approx::CalcB(const Element& elt, size_t i, double fa, double fb) const
    {
        double w, s, b = 0;

        assert(m_xGauss.size() == m_wGauss.size());

        for (size_t n = 0; n < m_xGauss.size(); n++)
        {
            s = m_xGauss[n];
            w = m_wGauss[n];
            b += w * Basis(i, s) * Ftilde(elt.X(s), s, fa, fb);
        }
        return b;
    }

    //
    // Returns value of function $\tilde{f}$ at $x$
    // Local variable "s" corresponds to global variable "x"
    //
    double Approx::Ftilde(double x, double s, double fa, double fb) const
    {
        assert(m_f);

        return m_f->Get(x) - fa * Basis(0, s) - fb * Basis(1, s);
    }

    //
    // Returns value of integral for function f^2.
    // Gauss quadratures applied.
    // 
    double Approx::IntegF2(const Element& elt, double fa, double fb) const
    {

        double v, w, s, res = 0;

        assert(m_xGauss.size() == m_wGauss.size());

        for (size_t n = 0; n < m_xGauss.size(); n++)
        {
            s = m_xGauss[n];
            w = m_wGauss[n];
            v = Ftilde(elt.X(s), s, fa, fb);
            res += w * v * v;
        }
        return elt.Jac() * res;

        // double ww1 = elt.Jac() * res;

        /*
        Int1DGauss integ(10);
        FunF2 f2(m_f, elt, fa, fb);

        double ww2 = elt.Jac() * integ.Adapt(f2, XFront(), XBack(), 1E-6);
        return ww2;
        */
    }

    //
    // Zwraca wartosc wezla $r_0$, dla ktorego aproksymacja na przedziale $[a, r_0]$ daje blad mniejszy niz $delta$.
    // Argument $h$ oznacza wartosc pierwszego kroku podczas szukania przedzialu dla algorytmu bisekcji
    //
    double Approx::FindB(double a, double b, double delta, double h)
    {
        double ra, rb, eps = 0;

        assert(delta > 0);

        ra = 0;
        rb = 0;
        while (true)
        {
            rb += h;
            eps = Solve(a, a + rb);
            if (eps > delta)
                break;

            if (a + rb > b)
                return b;

            ra = rb;
            // h = 2 * h;
        }
        return RunBisect(a, ra, rb, delta);
    }

    //
    // Metoda bisekcji znajduje punkt $r_0$ na przedziale [r_a, r_b]
    //
    double Approx::RunBisect(double a, double ra, double rb, double delta)
    {
        const double prec = 1E-4;
        double x = 0.0, fa, fb, fx;

        assert(rb > ra);

        if (ra > 0)
            fa = Solve(a, a + ra) - delta;
        else
            fa = -delta;

        fb = Solve(a, a + rb) - delta;

        assert(fa * fb < 0);

        while (rb - ra > prec)
        {
            x = 0.5 * (ra + rb);
            fx = Solve(a, a + x) - delta;
            if (fa * fx <= 0)
            {
                rb = x;
                fb = fx;
            }
            else
            {
                ra = x;
                fa = fx;
            }
        }
        return a + x;
    }

    //
    // Znajduje wspolrzedne wezlow a = r_0 < r_1 < .. r_N = b, tkaich, ze aproksymacja 
    // na kazdym odcinku [r_i, r_{i+1}] ma mniejszy blad niz delta
    //
    void Approx::SolveAdaptOld(double a, double b, double delta, double h)
    {
        double r0;
        size_t i;
        std::vector<double> r;

        m_cAll.resize(0);
        r.push_back(a);
        while (true)
        {
            r0 = FindB(a, b, delta, h);

            // Wspolczynniki dla \psi_0 i dla \psi_1
            m_cAll.push_back(m_f->Get(a));
            m_cAll.push_back(m_f->Get(r0));

            for (i = 0; i < m_M; i++) // Pozostale wspolczynniki
                m_cAll.push_back(m_c->Get(i));

            if (r0 >= b)
            {
                r.push_back(b);
                break;
            }
            else
            {
                r.push_back(r0);
                a = r0;
            }

        }
        SetMesh(r);
    }


    void Approx::SetMesh(const std::vector<double>& x)
    {
        std::vector<size_t> degree(x.size() - 1, m_M);

        Mesh::Set(x, degree);
        Mesh::CreateCnnt(BndrType_Dir, BndrType_Dir);
    }


    //
    // Znajduje wspolrzedne wezlow a = r_0 < r_1 < .. r_N = b, tkaich, ze aproksymacja 
    // na kazdym odcinku [r_i, r_{i+1}] ma mniejszy blad niz delta
    //
    void Approx::SolveAdapt(double a, double b, double maxDelta)
    {
        double x, delta;
        // size_t i;
        // std::vector<double> r, c(m_M + 2);
        std::vector<double> c(m_M + 2);
        HeapElt e;

        // m_cAll.resize(0);
        //	r.push_back(a);
        //	r.push_back(b);

        m_heap.Clear();
        m_heap.Reserve(50);

        delta = Solve(a, b);
        SetCoef(a, b, c);
        m_heap.Push(HeapElt(a, b, delta, c));

        size_t ii = 0;
        while (true)
        {
            if (m_heap.Top().m_delta < maxDelta && ii > 10)
                break;

            ii++;
            m_heap.Pop(e);

            a = e.m_left;
            b = e.m_right;
            x = 0.5 * (a + b);

            delta = Solve(a, x);
            SetCoef(a, x, c);
            m_heap.Push(HeapElt(a, x, delta, c));

            delta = Solve(x, b);
            SetCoef(x, b, c);
            m_heap.Push(HeapElt(x, b, delta, c));


            //		r.push_back(x);

            //		// Wspolczynniki dla \psi_0 i dla \psi_1
            //		m_cAll.push_back(m_f->Get(a));
            //		m_cAll.push_back(m_f->Get(x));
            //
            //		for(i = 0; i < m_M; i++) // Pozostale wspolczynniki
            //			m_cAll.push_back(m_c->Get(i));
        }
        //	std::sort(r.begin(), r.end());
        //	SetMesh(r);

        // printf("   ii=%d\n", ii);


        //	for(size_t i = 0; i < m_heap.Size(); i++)
        //	{
        //		const HeapElt& e = m_heap[i];
        //		printf("i=%2d  Left=%6.3f  Right=%6.3f   Delta =%6.3E\n", i, e.m_left, e.m_right, e.m_delta);
        //	}
        //	printf("\n\n");

# ifdef _DEBUG
        printf("Heap size = %d\n", m_heap.Size());
# endif

    }


    void Approx::SetCoef(double a, double b, std::vector<double>& c)
    {
        // Wspolczynniki dla \psi_0 i dla \psi_1
        c[0] = m_f->Get(a);
        c[1] = m_f->Get(b);

        for (size_t i = 0; i < m_M; i++) // Pozostale wspolczynniki
            c[i + 2] = m_c->Get(i);
    }


    std::vector<double> Approx::GetNode() const
    {
        std::vector<double> node;
        node.reserve(m_heap.Size());

        for (size_t i = 0; i < m_heap.Size(); i++)
        {
            const HeapElt& e = m_heap[i];
            node.push_back(e.m_left);
            node.push_back(e.m_right);
        }

        std::sort(node.begin(), node.end());
        std::vector<double>::iterator newEnd;
        newEnd = std::unique(node.begin(), node.end());
        node.erase(newEnd, node.end());

        return node;
    }


    void Approx::WriteCoef(FILE* out) const
    {
        size_t k, i;

        fprintf(out, "%4s \t %16s \t %16s", "i", "r_i", "r_{i+1}");
        for (k = 0; k <= m_M; ++k)
            fprintf(out, " \t              c_%lu", static_cast<unsigned long>(k));
        fprintf(out, "\n");

        for (i = 0; i < m_heap.Size(); ++i)
        {
            const HeapElt& e = m_heap[i];

            fprintf(out, "%4lu \t %16.6E \t %16.6E", static_cast<unsigned long>(i), e.m_left, e.m_right);
            for (k = 0; k <= m_M; ++k)
                fprintf(out, " \t %16.6E", e.m_coef[k]);
            fprintf(out, "\n");
        }

    }
}
