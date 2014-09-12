#ifndef __RATOM_BNDR_H__
#define __RATOM_BNDR_H__


/** \brief Boundary conditions
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/



//
// Vertex boundary conditions
//
enum BndrType
{
	BndrType_Emp = 0, // Empty
	BndrType_Dir, // Dirichlet
	BndrType_Neu // Neumann
};


class Bndr
{
public:
	Bndr(void);
	Bndr(BndrType type, double val);
	~Bndr(void);

public:
        // Type of boundary condition
	BndrType m_type;

        // Value of boundary contition
	double m_val;
};

#endif

