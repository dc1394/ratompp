#include "stdafx.h"
#include "brick.h"


Brick::Brick(void)
{
}

Brick::~Brick(void)
{
}

//!
//! Zmniejsz prostkota wzdluz osi X (bierze gore).
//!
void Brick::HalfXup(void)
{
	m_x1 = 0.5 * (m_x1 + m_x0);
}

//!
//! Zmniejsz prostkota wzdluz osi X (bierze dol).
//!
void Brick::HalfXdown(void)
{
	m_x0 = 0.5 * (m_x1 + m_x0);
}



//!
//! Zmniejsz prostkota wzdluz osi Y (bierze gore).
//!
void Brick::HalfYup(void)
{
	m_y1 = 0.5 * (m_y1 + m_y0);
}

//!
//! Zmniejsz prostkota wzdluz osi Y (bierze dol).
//!
void Brick::HalfYdown(void)
{
	m_y0 = 0.5 * (m_y1 + m_y0);
}

//!
//! Zmniejsz prostkota wzdluz osi Z (bierze gore).
//!
void Brick::HalfZup(void)
{
	m_z1 = 0.5 * (m_z1 + m_z0);
}

//!
//! Zmniejsz prostkota wzdluz osi Z (bierze dol).
//!
void Brick::HalfZdown(void)
{
	m_z0 = 0.5 * (m_z1 + m_z0);
}



//!
//! Operator porwnania
//!
bool Brick::operator<(const Brick& rec) const
{
	return (m_errAbs < rec.m_errAbs);
	// return (m_errRel < rec.m_errRel);
}


