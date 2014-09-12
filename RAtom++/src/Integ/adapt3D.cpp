#include "stdafx.h"
#include "adapt3D.h"


Adapt3D::Adapt3D(double epsAbs, size_t initSize) : m_epsAbs(epsAbs), m_pool(initSize), m_set(initSize)
{
	WeightAbsci();
}

Adapt3D::~Adapt3D(void)
{
}

//!
//! Oblicza wagi i wspolrzedne kwadratur
//!
void Adapt3D::WeightAbsci()
{
const int n = 3;
const double twondm = pow(2., n);

	m_q5 = sqrt(9. / 19.);
	m_q4 = sqrt(9. / 10.);
	m_q3 = m_q4;
	m_q2 = sqrt(9. / 70.);
	// m_q1 - nie ma takiej zmiennej	
	

	m_w7[0] = twondm * (12824 - 9120 * n + 400 * n * n) / 19683.;
	m_w7[1] = twondm * (980. / 6561.);
	m_w7[2] = twondm * (1820 - 400 * n) / 19683.;
	m_w7[3] = twondm * (200. / 19683.);
	m_w7[4] = 6859. / 19683.;

	m_w5[0] = twondm * (729 - 950 * n + 50 * n * n) / 729.;
	m_w5[1] = twondm * (245. / 486.);
	m_w5[2] = twondm * (265 - 100 * n) / 1458.;
	m_w5[3] = twondm * (25. / 729.);
}



//!
//! Zwraca wartosc calki z funkcji $f$ na prostopadloscianie $[x0, x1] \times [y0, y1] \times [z0, z1]$.
//!
double Adapt3D::Calc(const Fun3D& f, double x0, double x1, double y0, double y1, double z0, double z1)
{
Brick *br1, *br2;

	Clear();

	// Pierwsze wywolanie kwadratury
	br1 = NewBrick();
	br1->m_x0 = x0;
	br1->m_x1 = x1;
	br1->m_y0 = y0;
	br1->m_y1 = y1;
	br1->m_z0 = z0;
	br1->m_z1 = z1;

	Rule(f, br1);
	m_set.Add(br1); // Dodaje brick do kolejki priorytetowej

	size_t ii = 0;
	while(1)
	{
		if(ii > 10) // Musi byc wykonanych co najniej 10 podzialow
		{
			if(m_set.ErrAbs() < m_epsAbs)
				break;
		}
		ii++;

		br1 = m_set.Top(); // Pobieram brick z kolejki o najwiekszym bledzie
		br2 = NewBrick();
		if(br2 == NULL)
		{
			fprintf(stderr, "ARAPT WARNING. Max size = %lu of priority queue exceeded!\n\n", static_cast<unsigned long>(m_set.MaxSize()));
			break;
		}

		*br2 = *br1; // Kopiuje obeikty

		

		// Podzial na dwa prostopadlosciany
		if(br1->m_divAx == 0) // podzial wzdluz osi X
		{
			br1->HalfXup(); 
			br2->HalfXdown(); 
		}
		else if(br1->m_divAx == 1) // podzial wzdluz osi Y
		{
			br1->HalfYup();
			br2->HalfYdown();
		}
		else if(br1->m_divAx == 2) // podzial wzdluz osi Z
		{
			br1->HalfZup();
			br2->HalfZdown();
		}
		else
			assert(0);

		// Policzenie kwadratur
		Rule(f, br1);
		Rule(f, br2);

		// Dodanie do kolejki dwoch brikow
		m_set.Add(br1);
		m_set.Add(br2);

		// printf("Adapt3D::Quad = %18.8E\n", m_set.Quad());
	}

	return m_set.Quad();
}

/*
void Adapt3D::RuleEx(const Fun3D& fun, Brick* br)
{
	br->m_quad = Rule(fun, br->m_x0, br->m_x1, br->m_y0, br->m_y1, br->m_z0, br->m_z1, 
		&(br->m_errAbs), &(br->m_divAx));
}
*/

//!
//! Olicza kwadrature piatego i siodmego rzedu. Kwadratura obliczana jest z funkcji $f$
//! na prostopadloscianie $P = [x0, x1] \times [y0, y1] \times [z0, z1]$.
//! Prostopadloscian ma powierzchnie rownolegle 
//! do plaszczyzn ukladu wspolrzednych.
//!
void Adapt3D::Rule(const Fun3D& fun, Brick* br) 
{
const double x0 = br->m_x0, x1 = br->m_x1, y0 = br->m_y0, y1 = br->m_y1, z0 = br->m_z0, z1 = br->m_z1;
const double ratio = (m_q2 * m_q2) / (m_q4 * m_q4);
const int ndim = 3;
const double width[3] = {0.5 * (x1 - x0), 0.5 * (y1 - y0), 0.5 * (z1 - z0)};
const double center[3] = {x0 + width[0], y0 + width[1], z0 + width[2]};

double p[3], m[3], z[3], sum[5], f[4], df[2], difmax, dif;


// ******  COMPUTE FUN(0,0,0)
	sum[0] = fun.Get(center[0], center[1], center[2]);




// *****  COMPUTE SYMMETRIC SUMS OF FUN(LAMDA2,0,0) AND FUN(LAMDA3,0,0), AND MAXIMUM FOURTH DIFFERENCE 

    difmax = -1.;
	sum[1] = 0.;
    sum[2] = 0.;

	z[0] = center[0];
	z[1] = center[1];
	z[2] = center[2];

	for(int j = 0; j < ndim; j++)
	{
		z[j] = center[j] + m_q2 * width[j];
		f[0] = fun.Get(z[0], z[1], z[2]);

		z[j] = center[j] - m_q2 * width[j];
		f[1] = fun.Get(z[0], z[1], z[2]);

		z[j] = center[j] - m_q3 * width[j];
		f[2] = fun.Get(z[0], z[1], z[2]);

		z[j] = center[j] + m_q3 * width[j];
		f[3] = fun.Get(z[0], z[1], z[2]);

		sum[1] += (f[0] + f[1]);
		sum[2] += (f[2] + f[3]);

		df[0] = f[0] + f[1] - 2 * sum[0];
		df[1] = f[2] + f[3] - 2 * sum[0];

		dif = fabs(df[0] - ratio * df[1]);
		if(dif > difmax) // Poszukiwanie najwiekszej roznicy dzielonej
		{
			difmax = dif;
			br->m_divAx = j;
		}
		z[j] = center[j]; // Ustalenie wspolrzednej punktu
	}



// *****  COMPUTE SYMMETRIC SUM OF FUNCTN(LAMDA4,LAMDA4,0)

	p[0] = center[0] + m_q4 * width[0];
	m[0] = center[0] - m_q4 * width[0];
	p[1] = center[1] + m_q4 * width[1];
	m[1] = center[1] - m_q4 * width[1];
	p[2] = center[2] + m_q4 * width[2];
	m[2] = center[2] - m_q4 * width[2];

	{
	const double zz = center[2];

	sum[3]  = fun.Get(p[0], p[1], zz);
	sum[3] += fun.Get(p[0], m[1], zz);
	sum[3] += fun.Get(m[0], p[1], zz);
	sum[3] += fun.Get(m[0], m[1], zz);
	}

	{
	const double zz = center[1];
	
	sum[3] += fun.Get(p[0], zz, p[2]);
	sum[3] += fun.Get(p[0], zz, m[2]);
	sum[3] += fun.Get(m[0], zz, p[2]);
	sum[3] += fun.Get(m[0], zz, m[2]);
	}

	{
	const double zz = center[0];
	
	sum[3] += fun.Get(zz, p[1], p[2]);
	sum[3] += fun.Get(zz, p[1], m[2]);
	sum[3] += fun.Get(zz, m[1], p[2]);
	sum[3] += fun.Get(zz, m[1], m[2]);
	}



// *****  COMPUTE SYMMETRIC SUM OF FUNCTN(LAMDA5,LAMDA5,LAMDA5)


	p[0] = center[0] + m_q5 * width[0];
	m[0] = center[0] - m_q5 * width[0];

	p[1] = center[1] + m_q5 * width[1];
	m[1] = center[1] - m_q5 * width[1];

	p[2] = center[2] + m_q5 * width[2];
	m[2] = center[2] - m_q5 * width[2];

	sum[4]  = fun.Get(p[0], p[1], p[2]);
	sum[4] += fun.Get(p[0], p[1], m[2]);
	sum[4] += fun.Get(p[0], m[1], p[2]);
	sum[4] += fun.Get(p[0], m[1], m[2]);

	sum[4] += fun.Get(m[0], p[1], p[2]);
	sum[4] += fun.Get(m[0], p[1], m[2]);
	sum[4] += fun.Get(m[0], m[1], p[2]);
	sum[4] += fun.Get(m[0], m[1], m[2]);


// *****  COMPUTE FIFTH AND SEVENTH DEGREE RULES AND ERROR

	const double rgnvol = width[0] * width[1] * width[2];
	double v5, v7;

    v5 = rgnvol * (m_w5[0] * sum[0] + m_w5[1] * sum[1] + m_w5[2] * sum[2] + m_w5[3] * sum[3]);
    v7 = rgnvol * (m_w7[0] * sum[0] + m_w7[1] * sum[1] + m_w7[2] * sum[2] + m_w7[3] * sum[3] + m_w7[4] * sum[4]);
    br->m_errAbs = fabs(v7 - v5);
	// br->errRel = fabs((v7 - v5) / v7);

	// Zostalow wykonanych 33 obliczen funkcji
	m_funCall += 33;

	// Zapisanie wartosci kwadratury
	br->m_quad = v7;
}




//
// Zwraca koncowa liczba prostopadloscianow
//
size_t Adapt3D::BrickNo() const
{
	return m_set.EltNo();
}

//
// Zwraca calkowita liczba wywolan funkcji podcalkowej
//
size_t Adapt3D::FunCall() const
{
	return m_funCall;
}

Brick* Adapt3D::NewBrick()
{
	return m_pool.New();
}

//
// Przygotowanie do obliczen
//
void Adapt3D::Clear()
{
	m_funCall = 0;
	m_set.Clear();
	m_pool.Clear();
}


