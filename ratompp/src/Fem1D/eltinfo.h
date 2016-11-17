#ifndef __RATOM_ELTINFO_H__
#define __RATOM_ELTINFO_H__


/** \brief Information about element. Used by adaptive algorithm.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/


class EltInfo
{
public:
	EltInfo(void);
	~EltInfo(void);


        // Needed for sorting
	bool operator< (const EltInfo& e) const
	{
		return (m_eltId < e.m_eltId);
	}


        // Needed for duplicate elimination
	bool operator== (const EltInfo& e) const
	{
		return (e.m_eltId == m_eltId);
	}


public:
        // The largest from the smallest coefficients for eigenvalue "m_eigVal"
	double m_maxMinCoef;

        // Element ID
	size_t m_eltId;

        // Eigenvalue ID
	size_t m_eigVal;
};

#endif

