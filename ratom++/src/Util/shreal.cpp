#include "stdafx.h"
#include "shreal.h"
#include "helpfun.h"

double ShReal::m_clm[16];
double ShReal::m_slm[16];

//
// Constructor
//
ShReal::ShReal(void)
{
static bool init = false;

	if(!init)
	{	
		// l = 0
		m_clm[0] = 1. / sqrt(4 * util::HelpFun::M_PI);
	
		// l = 1
		m_clm[1] = sqrt(3. / (4 * util::HelpFun::M_PI));
		m_clm[2] = sqrt(3. / (4 * util::HelpFun::M_PI));
		m_clm[3] = sqrt(3. / (4 * util::HelpFun::M_PI));
	
		// l = 2
		m_clm[4] = sqrt(15. / ( 4 * util::HelpFun::M_PI));
		m_clm[5] = sqrt(15. / ( 4 * util::HelpFun::M_PI));
		m_clm[6] = sqrt( 5. / (16 * util::HelpFun::M_PI));
		m_clm[7] = sqrt(15. / ( 4 * util::HelpFun::M_PI));
		m_clm[8] = sqrt(15. / (16 * util::HelpFun::M_PI));
	
		// l = 3
		m_clm[ 9] = sqrt( 35. / (32 * util::HelpFun::M_PI));
		m_clm[10] = sqrt(105. / ( 4 * util::HelpFun::M_PI));
		m_clm[11] = sqrt( 21. / (32 * util::HelpFun::M_PI));
		m_clm[12] = sqrt(  7. / (16 * util::HelpFun::M_PI));
		m_clm[13] = sqrt( 21. / (32 * util::HelpFun::M_PI));
		m_clm[14] = sqrt(105. / (16 * util::HelpFun::M_PI));
		m_clm[15] = sqrt( 35. / (32 * util::HelpFun::M_PI));

		// -------------------------------------------------------------

		// l = 0
		m_slm[0] = 1. / sqrt(4 * util::HelpFun::M_PI);
	
		// l = 1
		m_slm[1] = sqrt(3. / (8 * util::HelpFun::M_PI));
		m_slm[2] = sqrt(3. / (4 * util::HelpFun::M_PI));
		m_slm[3] = sqrt(3. / (8 * util::HelpFun::M_PI));
	
		// l = 2
		m_slm[4] = sqrt(15. / (32 * util::HelpFun::M_PI));
		m_slm[5] = sqrt(15. / (32 * util::HelpFun::M_PI));
		m_slm[6] = sqrt( 5. / (64 * util::HelpFun::M_PI));
		m_slm[7] = sqrt(15. / (32 * util::HelpFun::M_PI));
		m_slm[8] = sqrt(15. / (32 * util::HelpFun::M_PI));
	
		// l = 3
		m_slm[ 9] = sqrt( 35. / (  64 * util::HelpFun::M_PI));
		m_slm[10] = sqrt(105. / (  32 * util::HelpFun::M_PI));
		m_slm[11] = sqrt( 21. / (1024 * util::HelpFun::M_PI));
		m_slm[12] = sqrt(  7. / ( 256 * util::HelpFun::M_PI));
		m_slm[13] = sqrt( 21. / (1024 * util::HelpFun::M_PI));
		m_slm[14] = sqrt(105. / (  32 * util::HelpFun::M_PI));
		m_slm[15] = sqrt( 35. / (  64 * util::HelpFun::M_PI));


		init = true;
	}
}

//
// Destructor
//
ShReal::~ShReal(void)
{
}

//
// Zwraca wartosc rzeczywistej harmoniki sferycznej w punkcie (x, y, z) w karezjanskim ukladzie wspolrzednych.
// Musi zachodzic zaleznosc r = sqrt(x*x + y*y + z*z)
// Promien "r" jest przekazywany aby uniknac wielokrotnego obliczania pierwiastka
//
double ShReal::Get(int l, int m, double x, double y, double z, double r) const
{
const int idx = Idx(l, m);
//const double r2 = r * r, r3 = r * r * r;
double v;

	assert(l >= 0 && l <= 3);
	assert(labs(m) <= l);

	if(r < FLT_EPSILON)
		return m_clm[0];

	switch(idx)
	{
		// l = 0
		case 0: v = 1; break;

		// l = 1
		case 1: v = y / r; break; // m = -1
		case 2: v = z / r; break; // m =  0
		case 3: v = x / r; break; // m =  1

		// l = 2
		case 4: v = x * y             / (r*r); break; // m = -2
		case 5: v = y * z             / (r*r); break; // m = -1
		case 6: v = (3 * z * z - r*r) / (r*r); break; // m =  0
		case 7: v = x * z             / (r*r); break; // m =  1
		case 8: v = (x - y) * (x + y) / (r*r); break; // m =  2

		// l = 3
		case  9: v = y * (3 * x * x - y * y)     / (r*r*r); break; // m = -3
		case 10: v = x * y * z                   / (r*r*r); break; // m = -2
		case 11: v = y * (5 * z * z - r * r)     / (r*r*r); break; // m = -1
		case 12: v = z * (5 * z * z - 3 * r * r) / (r*r*r); break; // m =  0
		case 13: v = x * (5 * z * z - r * r)     / (r*r*r); break; // m =  1
		case 14: v = z * (x - y) * (x + y)       / (r*r*r); break; // m =  2
		case 15: v = x * (x * x - 3 * y * y)     / (r*r*r); break; // m =  3

		default: assert(0); v = 0; break;

	}
	return m_clm[idx] * v;
}


//
// Zwraca wartosc wspolczynnika unormowania
//
double ShReal::Clm(int l, int m) const
{
	return m_clm[Idx(l, m)];
}


//
// Zwraca unikalny indeks odpowiadajacy parze (l, m)
//
int ShReal::Idx(int l, int m) const
{
	// Latwo sprawdzic, ze zachodzi
	// (l,  m) -> idx
	// (0,  0) -> 0
	// (1, -1) -> 1
	// (1,  0) -> 2
	// (1,  1) -> 3
	// (2, -2) -> 4
	// (2, -1) -> 5
	// (2,  0) -> 6
	// (2,  1) -> 7
	// (2,  2) -> 8
	return l * l + (l + m);
}


//
// Zwraca wartosc rzeczywistej harmoniki sferycznej w punkcie (theta, phi) w sferycznym ukladzie wspolrzednych.
//
double ShReal::Get(int l, int m, double theta, double phi) const
{
	return Plm(l, m, theta) * Qm(m, phi);
}


//
// Zwraca czesc zalezna od kata \theta w sferycznym ukladzie wspolrzednych
//
double ShReal::Plm(int l, int m, double theta) const
{
const int idx = Idx(l, m);
double v;

	assert(l >= 0 && l <= 3);
	assert(labs(m) <= l);

	switch(idx)
	{
		// l = 0
		case 0: v = 1; break;

		// l = 1
		case 1: v = sin(theta); break; // m = -1
		case 2: v = cos(theta); break; // m =  0
		case 3: v = sin(theta); break; // m =  1

		// l = 2
		case 4: v = sin(theta); v = v * v;	break; // m = -2
		case 5: v = sin(2 * theta);			break; // m = -1
		case 6: v = 1 + 3 * cos(2 * theta); break; // m =  0
		case 7: v = sin(2 * theta);			break; // m =  1
		case 8: v = sin(theta); v = v * v;  break; // m =  2


		// l = 3
		case 15: // m = +3
		case  9: // m = -3
			v = sin(theta); 
			v = v * v * v; 
		break; 

		case 14: // m = +2
		case 10: // m = -2
			v = sin(theta); 
			v = cos(theta) * v * v;
		break; 

		case 13: // m = +1
		case 11: // m = -1
			v = sin(theta) + 5 * sin(3 * theta); 
		break; 

		case 12: // m =  0
			v = 3 * cos(theta) + 5 * cos(3 * theta);
		break; 

		default: assert(0); v = 0; break;

	}
	return m_slm[idx] * v;

}

//
// Zwraca czesc zalezna od kata \theta w sferycznym ukladzie wspolrzednych
//
double ShReal::Qm(int m, double phi) const
{
static const double c = sqrt(2.);

	if(m > 0)
		return c * sin(m * phi);

	if(m < 0)
		return c * cos(-m * phi);

	return 1;
}

