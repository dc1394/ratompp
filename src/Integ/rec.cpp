#include "stdafx.h"
#include "rec.h"


Rec::Rec(void)
{
}

Rec::~Rec(void)
{
}

//!
//! Zmniejsz prostkota wzdluz osi X. Brana jest lewa polowa.
//!
void Rec::HalfLeft()
{
	m_x1 = 0.5 * (m_x1 + m_x0);
}

//!
//! Zmniejsz prostkota wzdluz osi X. Brana jest prawa polowa.
//!
void Rec::HalfRight()
{
	m_x0 = 0.5 * (m_x1 + m_x0);
}

//!
//! Zmniejsz prostkota wzdluz osi Y. Brana jest gorna polowa.
//!
void Rec::HalfTop()
{
	m_y0 = 0.5 * (m_y1 + m_y0);
}


//!
//! Zmniejsz prostkota wzdluz osi Y. Brana jest dolna polowa.
//!
void Rec::HalfBottom()
{
	m_y1 = 0.5 * (m_y1 + m_y0);
}

//!
//! Operator porwnania
//!
bool Rec::operator<(const Rec& rec) const
{
	return (m_errAbs < rec.m_errAbs);
}

/*
//!
//! Tworzy cwiartke. Numery cwiartek podane sa na rysunku
//! 
//! *----------*----------* <-- y1
//! |          |          |
//! |    4     |    1     |
//! |          |          |
//! *----------*----------*
//! |          |          |
//! |    3     |    2     |
//! |          |          |
//! |          |          |
//! *----------*----------* <-- y0
//! ^                     ^
//! x0                    x1
void Rec::Quad(int i)
{
const double xMid = 0.5 * (m_x0 + m_x1);
const double yMid = 0.5 * (m_y0 + m_y1);

	if(i == 1)
	{
		m_x0 = xMid;
		m_y0 = yMid;
	}
	else if(i == 2)
	{
		m_x0 = xMid;
		m_y1 = yMid;
	}
	else if(i == 3)
	{
		m_x1 = xMid;
		m_y1 = yMid;
	}
	else if(i == 4)
	{
		m_x1 = xMid;
		m_y0 = yMid;
	}
	else
		assert(0);
}
*/

