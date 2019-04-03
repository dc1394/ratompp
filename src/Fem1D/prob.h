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
	Prob();
	Prob(Bndr left, Bndr right, double gamma, std::shared_ptr<util::Fun1D> const & g, std::shared_ptr<util::Fun1D> const & f);
	virtual ~Prob() = default;

	void DefineProb(Bndr left, Bndr right, double gamma, std::shared_ptr<util::Fun1D> const & g, std::shared_ptr<util::Fun1D> const & f);

	void GenMeshLin(double a, double b, size_t nodeNo, size_t degree);
	void GenLinExp(double a, double b, int m, size_t degree);
	void SetMesh(const std::vector<double>& x, const std::vector<size_t>& degree);
	void AddToMesh(const std::vector<size_t>& eltToSplit);

	size_t Dim() const;

	
protected:
	double CalcO(const Element& e, size_t i, size_t j) const;
	double CalcB(const Element& e, size_t i) const;
	double CalcS(const Element& e, size_t i, size_t j) const;


protected:
        // Left boundary condition
	Bndr m_left;

        // Right boundary condition
	Bndr m_right;

private:
	// Right hand side, function f(x). USED only for ODE
	std::shared_ptr<util::Fun1D> m_f;

	// Function $g(x)$
	std::shared_ptr<util::Fun1D> m_g;

	// Constant \gamma 
	double m_gamma;
};

#endif

