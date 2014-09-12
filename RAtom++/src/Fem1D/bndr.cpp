#include "stdafx.h"
#include "bndr.h"

//
// Constructor
//
Bndr::Bndr(void) : m_type(BndrType_Emp), m_val(0)
{
}

//
// Constructor
//
Bndr::Bndr(BndrType type, double val) : m_type(type), m_val(val)
{
}

//
// Destructor
//
Bndr::~Bndr(void)
{
}
