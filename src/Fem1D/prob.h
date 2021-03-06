#ifndef __RATOM_PROB_H__
#define __RATOM_PROB_H__


/** \brief Represents function from R^3 into R.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/
/** \brief Represents problem solved by Fem1D library.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/


#include "mesh.h"
#include <memory>   // for std::shared_ptr

class Prob : protected Mesh, protected Lobatto
{
public:
	Prob() = default;
	Prob(Bndr left, Bndr right, double gamma, std::shared_ptr<util::Fun1D const> const & g, std::shared_ptr<util::Fun1D const> const & f);
	virtual ~Prob() = default;

	void DefineProb(Bndr left, Bndr right, double gamma, std::shared_ptr<util::Fun1D const> const & g, std::shared_ptr<util::Fun1D const> const & f);

	void GenMeshLin(double a, double b, std::size_t nodeNo, std::size_t degree);
	void GenLinExp(double a, double b, int m, std::size_t degree);
	void SetMesh(const std::vector<double>& x, const std::vector<std::size_t>& degree);
	void AddToMesh(const std::vector<std::size_t>& eltToSplit);

	std::size_t Dim() const;

	
protected:
	double CalcO(const Element& e, std::size_t i, std::size_t j) const;
	double CalcB(const Element& e, std::size_t i) const;
	double CalcS(const Element& e, std::size_t i, std::size_t j) const;


protected:
        // Left boundary condition
	Bndr m_left;

        // Right boundary condition
	Bndr m_right;

private:
	// Right hand side, function f(x). USED only for ODE
	std::shared_ptr<util::Fun1D const> m_f;

	// Function $g(x)$
	std::shared_ptr<util::Fun1D const> m_g;

	// Constant \gamma 
	double m_gamma = 0.0;
};

#endif

