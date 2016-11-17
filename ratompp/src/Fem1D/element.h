#ifndef __RATOM_ELEMENT_H__
#define __RATOM_ELEMENT_H__



/** \brief One dimensional element. 
*
* \author Zbigniew Romanowski [romz]
*
*/

class Element
{
public:
	Element(void);
	~Element(void);

	size_t P() const { return m_dof.size() - 1; }
	double X(double s) const { return m_c1 + s * m_c2; }
	double Xinv(double x) const { return (x - m_c1) / m_c2; } // X^{-1}

	size_t DofNo() const { return m_dof.size(); }

	void Set(double x0, double x1, size_t p);

//	bool IsNeumann() const;

	void Write(FILE* out) const;

        // Jacobian
	double Jac() const { return m_c2; }

	size_t PsiId(size_t i) const;

public:
        // DOF - DEGREE OF FREEDOM
        // The length of this vector is (p + 1), where "p" is the maximal degree of applied Lobatto functions
	std::vector<int> m_dof;

private:
        // (x[m+1] + x[m]) / 2
	double m_c1;

        // Jacobian: (x[m+1] - x[m]) / 2
	double m_c2;
};

#endif

