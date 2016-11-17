#include "stdafx.h"
#include "adapt2D.h"


Adapt2D::Adapt2D(double epsAbs, size_t initSize, Adapt2DRule rule) : m_funCall(0), m_epsAbs(epsAbs), 
	m_set(initSize), m_pool(initSize), m_rule(rule)
{
	WeightAbsci();
}

Adapt2D::~Adapt2D(void)
{
}

//!
//! Oblicza wagi i wspolrzedne kwadratur
//!
void Adapt2D::WeightAbsci()
{
}



//!
//! Zwraca wartosc calki z funkcji $f$ na prostokacie $[x0, x1] \times [y0, y1]$.
//!
double Adapt2D::Calc(const Fun2D& f, double x0, double x1, double y0, double y1)
{
Rec *r1, *r2;

	Clear();

	// Pierwsze wywolanie kwadratury
	r1 = NewRec();
	r1->m_x0 = x0;
	r1->m_x1 = x1;
	r1->m_y0 = y0;
	r1->m_y1 = y1;

	Rule(f, r1);
	m_set.Add(r1); // Dodaje prostokat do kolejki priorytetowej

	size_t ii = 0;
	while(1)
	{
		if(ii > 10) // Musi byc wykonanych co najniej 10 podzialow
		{
			if(m_set.ErrAbs() < m_epsAbs)
				break;
		}
		ii++;

		r1 = m_set.Top(); // Pobieram brick z kolejki o najwiekszym bledzie
		r2 = NewRec();
		if(r2 == NULL)
		{
			fprintf(stderr, "ARAPT WARNING. Max size = %lu of priority queue exceeded!\n\n", static_cast<unsigned long>(m_set.MaxSize()));
			break;
		}

		*r2 = *r1; // Kopiuje obeikty

		

		// Podzial na dwa prostokaty
		if(r1->m_divAx == 0) // podzial wzdluz osi X
		{
			r1->HalfLeft(); 
			r2->HalfRight(); 
		}
		else // podzial wzdluz osi Y
		{
			r1->HalfTop();
			r2->HalfBottom();
		}

		// Policzenie kwadratur
		Rule(f, r1);
		Rule(f, r2);

		// Dodanie do kolejki dwoch brikow
		m_set.Add(r1);
		m_set.Add(r2);

		// printf("Adapt2D::Quad = %18.8E\n", m_set.Quad());
	}

	return m_set.Quad();
}


//!
//! Olicza kwadrature piatego i siodmego rzedu. Kwadratura obliczana jest z funkcji $f$
//! na prostopadloscianie $P = [x0, x1] \times [y0, y1] \times [z0, z1]$.
//! Prostopadloscian ma powierzchnie rownolegle 
//! do plaszczyzn ukladu wspolrzednych.
//!
void Adapt2D::Rule(const Fun2D& fun, Rec* rec) 
{
	if(m_rule == Adapt2DRule_57)
		return Rule57(fun, rec);

	if(m_rule == Adapt2DRule_79)
		return Rule79(fun, rec);

	if(m_rule == Adapt2DRule_911)
		return Rule911(fun, rec);

	if(m_rule == Adapt2DRule_913)
		return Rule913(fun, rec);

	assert(0);
}



//
// Zwraca koncowa liczba prostopadloscianow
//
size_t Adapt2D::RecNo() const
{
	return m_set.EltNo();
}

//
// Zwraca calkowita liczba wywolan funkcji podcalkowej
//
size_t Adapt2D::FunCall() const
{
	return m_funCall;
}

Rec* Adapt2D::NewRec()
{
	return m_pool.New();
}

//
// Przygotowanie do obliczen
//
void Adapt2D::Clear()
{
	m_funCall = 0;
	m_set.Clear();
	m_pool.Clear();
}


//!
//! Olicza kwadrature piategi i siodmeg rzedu. Kwadratura obliczana jest z funkcji $f$
//! na prostokacie $P = [x0, x1] \times [y0, y1]$. Prostokat ma boki rownolegle do osi ukladu wspolrzednych.
//!
void Adapt2D::Rule57(const Fun2D& fun, Rec* rec) 
{
const int ndim = 2;

const double x0 = rec->m_x0, x1 = rec->m_x1, y0 = rec->m_y0, y1 = rec->m_y1;

const double w7[5] = {
	-0.77549154092363968907,
     0.59746989788142051517,
	 0.20728547477518670934,
     0.040644210740232688106,
	 0.34847330183407000965};

const double w5[4] = {
	-5.32784636488340192044,
	 2.01646090534979423868,
	 0.17832647462277091907,
	 0.13717421124828532236};
    
const double lam[4] = { 0.0, sqrt(9. / 70.), sqrt(9. / 10.), sqrt(9. / 19.)};

const double ratio = (lam[1] * lam[1]) / (lam[2] * lam[2]);
const double width[2] = {0.5 * (x1 - x0), 0.5 * (y1 - y0)};
const double center[2] = {x0 + width[0], y0 + width[1]};

// double z[5], width[2], center[2], sum[5], f[4], df[2], difmax, dif, x[2], y[2];
double sum[5], f[4], df[2], x[2], y[2];


// ***** COMPUTE $f(0, 0)$
	sum[0] = fun.Get(center[0], center[1]);

// ***** COMPUTE SYMMETRIC SUMS OF $f(\lambda_1, 0)$ AND $f(\lambda_2, 0)$, AND MAXIMUM FOURTH DIFFERENCE 


	double difmax, dif, z[5];
    difmax = -1.;
	sum[1] = 0.;
    sum[2] = 0.;
	z[0] = center[0];
	z[1] = center[1];

	for(int j = 0; j < ndim; j++)
	{
		z[j] = center[j] - lam[1] * width[j];
		f[0] = fun.Get(z[0], z[1]);

		z[j] = center[j] + lam[1] * width[j];
		f[1] = fun.Get(z[0], z[1]);

		z[j] = center[j] - lam[2] * width[j];
		f[2] = fun.Get(z[0], z[1]);

		z[j] = center[j] + lam[2] * width[j];
		f[3] = fun.Get(z[0], z[1]);

		sum[1] += (f[0] + f[1]);
		sum[2] += (f[2] + f[3]);

		df[0] = f[0] + f[1] - 2 * sum[0];
		df[1] = f[2] + f[3] - 2 * sum[0];

		dif = fabs(df[0] - ratio * df[1]);
		if(dif > difmax) // Poszukiwanie najwiekszej roznicy dzielonej
		{
			difmax = dif;
			rec->m_divAx = j;
		}
		z[j] = center[j]; // Ustalenie wspolrzednej punktu
	}
/*
	// \lambda_1 and \lambbda_2
	for(int i = 1; i <= 2; i++)
	{
		x[0] = center[i - 1] + lam[i] * width[i - 1];
		x[1] = center[i - 1] - lam[i] * width[i - 1];

		y[0] = center[i - 1] + lam[i] * width[i - 1];
		y[1] = center[i - 1] - lam[i] * width[i - 1];

		f[0] = fun.Get(x[0], center[1]);
		f[1] = fun.Get(x[1], center[1]);
		f[2] = fun.Get(center[0], y[0]);
		f[3] = fun.Get(center[0], y[1]);

		sum[i] = f[0] + f[1] + f[2] + f[3];

		// Obliczenie roznicy dzielonej
		const double u1 = f[0] + f[1] - 2 * sum[0];
		const double u2 = f[2] + f[3] - 2 * sum[0];
		df[i - 1] = fabs(u1 - ratio * u2);
	}

	// Poszukiwanie najwiekszej roznicy dzielonej
	if(df[0] > df[1])
		*divAx = 1;
	else
		*divAx = 0;
*/

// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_2, \lambda_2)$

	x[0] = center[0] + lam[2] * width[0];
	x[1] = center[0] - lam[2] * width[0];

	y[0] = center[1] + lam[2] * width[1];
	y[1] = center[1] - lam[2] * width[1];

	sum[3]  = fun.Get(x[0], y[0]);
	sum[3] += fun.Get(x[0], y[1]);
	sum[3] += fun.Get(x[1], y[0]);
	sum[3] += fun.Get(x[1], y[1]);


// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_3, \lambda_3)$


	x[0] = center[0] + lam[3] * width[0];
	x[1] = center[0] - lam[3] * width[0];

	y[0] = center[1] + lam[3] * width[1];
	y[1] = center[1] - lam[3] * width[1];

	sum[4]  = fun.Get(x[0], y[0]);
	sum[4] += fun.Get(x[0], y[1]);
	sum[4] += fun.Get(x[1], y[0]);
	sum[4] += fun.Get(x[1], y[1]);


// *****  COMPUTE FIFTH AND SEVENTH DEGREE RULES AND ERROR

	double rgnvol = width[0] * width[1];
	double v5, v7;

    v5 = rgnvol * (w5[0] * sum[0] + w5[1] * sum[1] + w5[2] * sum[2] + w5[3] * sum[3]);
    v7 = rgnvol * (w7[0] * sum[0] + w7[1] * sum[1] + w7[2] * sum[2] + w7[3] * sum[3] + w7[4] * sum[4]);
    rec->m_errAbs = fabs(v7 - v5);

	// Zostalow wykonanych 17 obliczen funkcji
	m_funCall += 17;

	// Zapisanie wartosci kwadratury
	rec->m_quad = v7;
}



//!
//! Olicza kwadrature piategi i siodmeg rzedu. Kwadratura obliczana jest z funkcji $f$
//! na prostokacie $P = [x0, x1] \times [y0, y1]$. Prostokat ma boki rownolegle do osi ukladu wspolrzednych.
//!
void Adapt2D::Rule79(const Fun2D& fun, Rec* rec) 
{
const int ndim = 2;

const double x0 = rec->m_x0, x1 = rec->m_x1, y0 = rec->m_y0, y1 = rec->m_y1;

const double w9[6] = {
	0.32363456790123456790,
	0.13478507238752090312,
	0.27228653255075070182,
	0.22908540422399111713,
	0.056134348862428635955,
	0.11340000000000000000};
	
const double w7[6] = {
	0.67592092205970042525,
	0.23092842785903867626,
	0.,
	0.43953907332966785983,
	0.082373073956971141166,
	0.039089597169698608216};
	    
const double lam[3] = {0., 0.90617984593866399280, 0.53846931010568309104};

const double ratio = (lam[1] * lam[1]) / (lam[2] * lam[2]);
const double width[2] = {0.5 * (x1 - x0), 0.5 * (y1 - y0)};
const double center[2] = {x0 + width[0], y0 + width[1]};
double sum[6], f[4], df[2], x[2], y[2];


// ***** COMPUTE $f(0, 0)$
	sum[0] = fun.Get(center[0], center[1]);

// ***** COMPUTE SYMMETRIC SUMS OF $f(\lambda_1, 0)$ AND $f(\lambda_2, 0)$, AND MAXIMUM FOURTH DIFFERENCE 

	double difmax, dif, z[5];
    difmax = -1.;
	sum[1] = 0.;
    sum[2] = 0.;
	z[0] = center[0];
	z[1] = center[1];

	for(int j = 0; j < ndim; j++)
	{
		z[j] = center[j] - lam[1] * width[j];
		f[0] = fun.Get(z[0], z[1]);

		z[j] = center[j] + lam[1] * width[j];
		f[1] = fun.Get(z[0], z[1]);

		z[j] = center[j] - lam[2] * width[j];
		f[2] = fun.Get(z[0], z[1]);

		z[j] = center[j] + lam[2] * width[j];
		f[3] = fun.Get(z[0], z[1]);

		sum[1] += (f[0] + f[1]);
		sum[2] += (f[2] + f[3]);

		df[0] = f[0] + f[1] - 2 * sum[0];
		df[1] = f[2] + f[3] - 2 * sum[0];

		dif = fabs(df[0] - ratio * df[1]);
		if(dif > difmax) // Poszukiwanie najwiekszej roznicy dzielonej
		{
			difmax = dif;
			rec->m_divAx = j;
		}
		z[j] = center[j]; // Ustalenie wspolrzednej punktu
	}

// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_2, \lambda_2)$
	x[0] = center[0] + lam[2] * width[0];
	x[1] = center[0] - lam[2] * width[0];

	y[0] = center[1] + lam[2] * width[1];
	y[1] = center[1] - lam[2] * width[1];

	sum[3]  = fun.Get(x[0], y[0]);
	sum[3] += fun.Get(x[0], y[1]);
	sum[3] += fun.Get(x[1], y[0]);
	sum[3] += fun.Get(x[1], y[1]);



// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_1, \lambda_1)$
	x[0] = center[0] + lam[1] * width[0];
	x[1] = center[0] - lam[1] * width[0];

	y[0] = center[1] + lam[1] * width[1];
	y[1] = center[1] - lam[1] * width[1];

	sum[4]  = fun.Get(x[0], y[0]);
	sum[4] += fun.Get(x[0], y[1]);
	sum[4] += fun.Get(x[1], y[0]);
	sum[4] += fun.Get(x[1], y[1]);


// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_1, \lambda_2)$
	x[0] = center[0] + lam[1] * width[0];
	x[1] = center[0] - lam[1] * width[0];
	y[0] = center[1] + lam[2] * width[1];
	y[1] = center[1] - lam[2] * width[1];

	sum[5]  = fun.Get(x[0], y[0]);
	sum[5] += fun.Get(x[0], y[1]);
	sum[5] += fun.Get(x[1], y[0]);
	sum[5] += fun.Get(x[1], y[1]);

	x[0] = center[0] + lam[2] * width[0];
	x[1] = center[0] - lam[2] * width[0];
	y[0] = center[1] + lam[1] * width[1];
	y[1] = center[1] - lam[1] * width[1];

	sum[5] += fun.Get(x[0], y[0]);
	sum[5] += fun.Get(x[0], y[1]);
	sum[5] += fun.Get(x[1], y[0]);
	sum[5] += fun.Get(x[1], y[1]);



// *****  COMPUTE SEVENTH AND NINTH DEGREE RULES AND ERROR

	double rgnvol = width[0] * width[1];
	double v7 = 0, v9 = 0;

	for(int i = 0; i < 6; i++)
	{
		v7 += w7[i] * sum[i];
		v9 += w9[i] * sum[i];
	}
	v7 *= rgnvol;
	v9 *= rgnvol;

    rec->m_errAbs = fabs(v9 - v7);

	// Zostalow wykonanych 25 obliczen funkcji
	m_funCall += 25;

	// Zapisanie wartosci kwadratury
	rec->m_quad = v9;
}


void Adapt2D::Rule911(const Fun2D& fun, Rec* rec) 
{
double v11, v9;

	v11 = Rule11(fun, rec);
	v9  = Rule9 (fun, rec);

	rec->m_errAbs = fabs(v11 - v9);

	m_funCall += (20 + 28);
	rec->m_quad = v11;
}

void Adapt2D::Rule913(const Fun2D& fun, Rec* rec) 
{
double v13, v9;

	v13 = Rule13(fun, rec);
	v9  = Rule9 (fun, rec);

	rec->m_errAbs = fabs(v13 - v9);

	m_funCall += (37 + 20);
	rec->m_quad = v13;
}





//!
//! Olicza kwadrature piategi i siodmeg rzedu. Kwadratura obliczana jest z funkcji $f$
//! na prostokacie $P = [x0, x1] \times [y0, y1]$. Prostokat ma boki rownolegle do osi ukladu wspolrzednych.
//!
double Adapt2D::Rule11(const Fun2D& fun, Rec* rec) const 
{
const int ndim = 2;

const double x0 = rec->m_x0, x1 = rec->m_x1, y0 = rec->m_y0, y1 = rec->m_y1;

const double w11[6] = {
	0.0176679598882646,
    0.2322248008989674,
	0.0715516745178401,
    0.2192868905662522,
    0.2965842326220580,
	0.0813422207533089};

const double lam11[7] = {
	0.8989737240828844,
	0.7632367891419969,
	0.8949648832822285,
	0.6322452037101434,
	0.2797353125538562,
	0.9602661668053869,
	0.4347413023856830};

const double ratio = (lam11[0] * lam11[0]) / (lam11[1] * lam11[1]);
const double width[2] = {0.5 * (x1 - x0), 0.5 * (y1 - y0)};
const double center[2] = {x0 + width[0], y0 + width[1]};
double sum[6], f[4], df[2], x[2], y[2], f0;


// ***** COMPUTE $f(0, 0)$
	f0 = fun.Get(center[0], center[1]);

// ***** COMPUTE SYMMETRIC SUMS OF $f(\lambda_0, 0)$ AND $f(\lambda_1, 0)$, AND MAXIMUM FOURTH DIFFERENCE 

	double difmax, dif, z[5];
    difmax = -1.;
	sum[0] = 0.;
    sum[1] = 0.;
	z[0] = center[0];
	z[1] = center[1];

	for(int j = 0; j < ndim; j++)
	{
		z[j] = center[j] - lam11[0] * width[j];
		f[0] = fun.Get(z[0], z[1]);

		z[j] = center[j] + lam11[0] * width[j];
		f[1] = fun.Get(z[0], z[1]);

		z[j] = center[j] - lam11[1] * width[j];
		f[2] = fun.Get(z[0], z[1]);

		z[j] = center[j] + lam11[1] * width[j];
		f[3] = fun.Get(z[0], z[1]);

		sum[0] += (f[0] + f[1]);
		sum[1] += (f[2] + f[3]);

		df[0] = f[0] + f[1] - 2 * f0;
		df[1] = f[2] + f[3] - 2 * f0;

		dif = fabs(df[0] - ratio * df[1]);
		if(dif > difmax) // Poszukiwanie najwiekszej roznicy dzielonej
		{
			difmax = dif;
			rec->m_divAx = j;
		}
		z[j] = center[j]; // Ustalenie wspolrzednej punktu
	}

// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_2, \lambda_2)$
	x[0] = center[0] + lam11[2] * width[0];
	x[1] = center[0] - lam11[2] * width[0];

	y[0] = center[1] + lam11[2] * width[1];
	y[1] = center[1] - lam11[2] * width[1];

	sum[2]  = fun.Get(x[0], y[0]);
	sum[2] += fun.Get(x[0], y[1]);
	sum[2] += fun.Get(x[1], y[0]);
	sum[2] += fun.Get(x[1], y[1]);



// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_3, \lambda_3)$
	x[0] = center[0] + lam11[3] * width[0];
	x[1] = center[0] - lam11[3] * width[0];

	y[0] = center[1] + lam11[3] * width[1];
	y[1] = center[1] - lam11[3] * width[1];

	sum[3]  = fun.Get(x[0], y[0]);
	sum[3] += fun.Get(x[0], y[1]);
	sum[3] += fun.Get(x[1], y[0]);
	sum[3] += fun.Get(x[1], y[1]);


// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_4, \lambda_4)$
	x[0] = center[0] + lam11[4] * width[0];
	x[1] = center[0] - lam11[4] * width[0];

	y[0] = center[1] + lam11[4] * width[1];
	y[1] = center[1] - lam11[4] * width[1];

	sum[4]  = fun.Get(x[0], y[0]);
	sum[4] += fun.Get(x[0], y[1]);
	sum[4] += fun.Get(x[1], y[0]);
	sum[4] += fun.Get(x[1], y[1]);


// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_5, \lambda_6)$
	x[0] = center[0] + lam11[5] * width[0];
	x[1] = center[0] - lam11[5] * width[0];
	y[0] = center[1] + lam11[6] * width[1];
	y[1] = center[1] - lam11[6] * width[1];

	sum[5]  = fun.Get(x[0], y[0]);
	sum[5] += fun.Get(x[0], y[1]);
	sum[5] += fun.Get(x[1], y[0]);
	sum[5] += fun.Get(x[1], y[1]);

	x[0] = center[0] + lam11[6] * width[0];
	x[1] = center[0] - lam11[6] * width[0];
	y[0] = center[1] + lam11[5] * width[1];
	y[1] = center[1] - lam11[5] * width[1];

	sum[5] += fun.Get(x[0], y[0]);
	sum[5] += fun.Get(x[0], y[1]);
	sum[5] += fun.Get(x[1], y[0]);
	sum[5] += fun.Get(x[1], y[1]);



// *****  COMPUTE ELEVENTH RULE
	double v11 = 0;
	for(int i = 0; i < 6; i++)
		v11 += w11[i] * sum[i];

	return v11 * width[0] * width[1];
}

//!
//! Olicza kwadrature piategi i siodmeg rzedu. Kwadratura obliczana jest z funkcji $f$
//! na prostokacie $P = [x0, x1] \times [y0, y1]$. Prostokat ma boki rownolegle do osi ukladu wspolrzednych.
//!
double Adapt2D::Rule9(const Fun2D& fun, Rec* rec) const
{
const int ndim = 2;

const double x0 = rec->m_x0, x1 = rec->m_x1, y0 = rec->m_y0, y1 = rec->m_y1;

const double w9[4] = {
	0.0716134247098111,
    0.4540903525515453,
	0.0427846154667780,
	0.2157558036359328};

const double lam9[5] = {
	0.9845398119422523,
    0.4888863428423724,
	0.9395672874215217,
	0.8367103250239890,
	0.5073767736746132};

const double width[2] = {0.5 * (x1 - x0), 0.5 * (y1 - y0)};
const double center[2] = {x0 + width[0], y0 + width[1]};
double sum[4], f[4], x[2], y[2];


// ***** COMPUTE SYMMETRIC SUMS OF $f(\lambda_0, 0)$ AND $f(\lambda_1, 0)$
	double z[5];
	sum[0] = 0.;
    sum[1] = 0.;
	z[0] = center[0];
	z[1] = center[1];

	for(int j = 0; j < ndim; j++)
	{
		z[j] = center[j] - lam9[0] * width[j];
		f[0] = fun.Get(z[0], z[1]);

		z[j] = center[j] + lam9[0] * width[j];
		f[1] = fun.Get(z[0], z[1]);

		z[j] = center[j] - lam9[1] * width[j];
		f[2] = fun.Get(z[0], z[1]);

		z[j] = center[j] + lam9[1] * width[j];
		f[3] = fun.Get(z[0], z[1]);

		sum[0] += (f[0] + f[1]);
		sum[1] += (f[2] + f[3]);

		z[j] = center[j]; // Ustalenie wspolrzednej punktu
	}

// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_2, \lambda_2)$
	x[0] = center[0] + lam9[2] * width[0];
	x[1] = center[0] - lam9[2] * width[0];

	y[0] = center[1] + lam9[2] * width[1];
	y[1] = center[1] - lam9[2] * width[1];

	sum[2]  = fun.Get(x[0], y[0]);
	sum[2] += fun.Get(x[0], y[1]);
	sum[2] += fun.Get(x[1], y[0]);
	sum[2] += fun.Get(x[1], y[1]);


// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_3, \lambda_4)$
	x[0] = center[0] + lam9[3] * width[0];
	x[1] = center[0] - lam9[3] * width[0];
	y[0] = center[1] + lam9[4] * width[1];
	y[1] = center[1] - lam9[4] * width[1];

	sum[3]  = fun.Get(x[0], y[0]);
	sum[3] += fun.Get(x[0], y[1]);
	sum[3] += fun.Get(x[1], y[0]);
	sum[3] += fun.Get(x[1], y[1]);

	x[0] = center[0] + lam9[4] * width[0];
	x[1] = center[0] - lam9[4] * width[0];
	y[0] = center[1] + lam9[3] * width[1];
	y[1] = center[1] - lam9[3] * width[1];

	sum[3] += fun.Get(x[0], y[0]);
	sum[3] += fun.Get(x[0], y[1]);
	sum[3] += fun.Get(x[1], y[0]);
	sum[3] += fun.Get(x[1], y[1]);



// *****  COMPUTE NINTH RULE
	double v9 = 0;
	for(int i = 0; i < 4; i++)
		v9 += w9[i] * sum[i];

	return v9 * width[0] * width[1];
}


//!
//! Olicza kwadrature piategi i siodmeg rzedu. Kwadratura obliczana jest z funkcji $f$
//! na prostokacie $P = [x0, x1] \times [y0, y1]$. Prostokat ma boki rownolegle do osi ukladu wspolrzednych.
//!
double Adapt2D::Rule13(const Fun2D& fun, Rec* rec) const 
{
const int ndim = 2;

const double x0 = rec->m_x0, x1 = rec->m_x1, y0 = rec->m_y0, y1 = rec->m_y1;

const double w13[8] = {
	0.2995235559387052,
    0.0331100668669073,
	0.1802214941550577,
    0.0391672789603492,
    0.1387748348777338,
	0.2268881207335663,
	0.0365739576550950,
	0.1169047000557597};

const double lam13[9] = {
	0.9909890363004588,
	0.6283940712305196,
	0.9194861553393097,
	0.6973201917871096,
	0.3805687186904865,
	0.9708504361720127,
	0.6390348393207334,
	0.8623637916722844,
	0.3162277660168378};


const double ratio = (lam13[0] * lam13[0]) / (lam13[1] * lam13[1]);
const double width[2] = {0.5 * (x1 - x0), 0.5 * (y1 - y0)};
const double center[2] = {x0 + width[0], y0 + width[1]};
double sum[8], f[4], df[2], x[2], y[2];


// ***** COMPUTE $f(0, 0)$
	sum[0] = fun.Get(center[0], center[1]);

// ***** COMPUTE SYMMETRIC SUMS OF $f(\lambda_0, 0)$ AND $f(\lambda_1, 0)$, AND MAXIMUM FOURTH DIFFERENCE 

	double difmax, dif, z[5];
    difmax = -1.;
	sum[1] = 0.;
    sum[2] = 0.;
	z[0] = center[0];
	z[1] = center[1];

	for(int j = 0; j < ndim; j++)
	{
		z[j] = center[j] - lam13[0] * width[j];
		f[0] = fun.Get(z[0], z[1]);

		z[j] = center[j] + lam13[0] * width[j];
		f[1] = fun.Get(z[0], z[1]);

		z[j] = center[j] - lam13[1] * width[j];
		f[2] = fun.Get(z[0], z[1]);

		z[j] = center[j] + lam13[1] * width[j];
		f[3] = fun.Get(z[0], z[1]);

		sum[1] += (f[0] + f[1]);
		sum[2] += (f[2] + f[3]);

		df[0] = f[0] + f[1] - 2 * sum[0];
		df[1] = f[2] + f[3] - 2 * sum[0];

		dif = fabs(df[0] - ratio * df[1]);
		if(dif > difmax) // Poszukiwanie najwiekszej roznicy dzielonej
		{
			difmax = dif;
			rec->m_divAx = j;
		}
		z[j] = center[j]; // Ustalenie wspolrzednej punktu
	}

// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_2, \lambda_2)$
	x[0] = center[0] + lam13[2] * width[0];
	x[1] = center[0] - lam13[2] * width[0];

	y[0] = center[1] + lam13[2] * width[1];
	y[1] = center[1] - lam13[2] * width[1];

	sum[3]  = fun.Get(x[0], y[0]);
	sum[3] += fun.Get(x[0], y[1]);
	sum[3] += fun.Get(x[1], y[0]);
	sum[3] += fun.Get(x[1], y[1]);



// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_3, \lambda_3)$
	x[0] = center[0] + lam13[3] * width[0];
	x[1] = center[0] - lam13[3] * width[0];

	y[0] = center[1] + lam13[3] * width[1];
	y[1] = center[1] - lam13[3] * width[1];

	sum[4]  = fun.Get(x[0], y[0]);
	sum[4] += fun.Get(x[0], y[1]);
	sum[4] += fun.Get(x[1], y[0]);
	sum[4] += fun.Get(x[1], y[1]);


// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_4, \lambda_4)$
	x[0] = center[0] + lam13[4] * width[0];
	x[1] = center[0] - lam13[4] * width[0];

	y[0] = center[1] + lam13[4] * width[1];
	y[1] = center[1] - lam13[4] * width[1];

	sum[5]  = fun.Get(x[0], y[0]);
	sum[5] += fun.Get(x[0], y[1]);
	sum[5] += fun.Get(x[1], y[0]);
	sum[5] += fun.Get(x[1], y[1]);


// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_5, \lambda_6)$
	x[0] = center[0] + lam13[5] * width[0];
	x[1] = center[0] - lam13[5] * width[0];
	y[0] = center[1] + lam13[6] * width[1];
	y[1] = center[1] - lam13[6] * width[1];

	sum[6]  = fun.Get(x[0], y[0]);
	sum[6] += fun.Get(x[0], y[1]);
	sum[6] += fun.Get(x[1], y[0]);
	sum[6] += fun.Get(x[1], y[1]);

	x[0] = center[0] + lam13[6] * width[0];
	x[1] = center[0] - lam13[6] * width[0];
	y[0] = center[1] + lam13[5] * width[1];
	y[1] = center[1] - lam13[5] * width[1];

	sum[6] += fun.Get(x[0], y[0]);
	sum[6] += fun.Get(x[0], y[1]);
	sum[6] += fun.Get(x[1], y[0]);
	sum[6] += fun.Get(x[1], y[1]);


// *****  COMPUTE SYMMETRIC SUM OF $f(\lambda_7, \lambda_8)$
	x[0] = center[0] + lam13[7] * width[0];
	x[1] = center[0] - lam13[7] * width[0];
	y[0] = center[1] + lam13[8] * width[1];
	y[1] = center[1] - lam13[8] * width[1];

	sum[7]  = fun.Get(x[0], y[0]);
	sum[7] += fun.Get(x[0], y[1]);
	sum[7] += fun.Get(x[1], y[0]);
	sum[7] += fun.Get(x[1], y[1]);

	x[0] = center[0] + lam13[8] * width[0];
	x[1] = center[0] - lam13[8] * width[0];
	y[0] = center[1] + lam13[7] * width[1];
	y[1] = center[1] - lam13[7] * width[1];

	sum[7] += fun.Get(x[0], y[0]);
	sum[7] += fun.Get(x[0], y[1]);
	sum[7] += fun.Get(x[1], y[0]);
	sum[7] += fun.Get(x[1], y[1]);


// *****  COMPUTE ELEVENTH RULE
	double v13 = 0;
	for(int i = 0; i < 8; i++)
		v13 += w13[i] * sum[i];

	return v13 * width[0] * width[1];
}
