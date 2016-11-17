#ifndef __RATOM_ELSET_H__
#define __RATOM_ELSET_H__



/**
*	\brief Zbior elementow wykorzystywanych podczas calkowania adaptacyjnego.
*		
*	\author Zbigniew Romanowski, ROMZ
*
*	\version 29-Sty-2009 [romz]
*
*/

#include "priorityqueue.h"

template <typename T>
class EltSet
{
public:
	EltSet(size_t initSize);
	~EltSet(void);

	double Quad(void) const;
	size_t EltNo(void) const;

	void Clear(void);
	void Add(T* br);

	T* Top(void);

	double ErrAbs(void) const;

	size_t MaxSize(void) const { return m_pq.MaxSize(); }

private:

	// Suma bledow (BEZWZGLEDNYCH) ze wszystkich prostopadloscianow
	double m_errAbs;

	// Suma kwadratur ze wszystkich prostopadloscianow
	double m_quad;

	// Kolejka priorytetowa zawiera wskazniki
 	PriorityQueue<T> m_pq;
};

template <typename T>
EltSet<T>::EltSet(size_t initSize) : m_errAbs(0), m_quad(0), m_pq(initSize)
{
}

template <typename T>
EltSet<T>::~EltSet(void)
{
}

//!
//! Zwraca wartosc bledy BEZWZGLEDNEGO
//!
template <typename T>
double EltSet<T>::ErrAbs(void) const
{
	return m_errAbs;
}

//!
//! Zwraca wartosc kwadratury
//!
template <typename T>
double EltSet<T>::Quad(void) const
{
	return m_quad;
}

//!
//! Zwraca liczbe prostokatow w zbiorze
//!
template <typename T>
size_t EltSet<T>::EltNo(void) const
{
	return m_pq.EltNo();
}


//!
//! Dodaje prostopadloscian do kolejki
//!
template <typename T>
void EltSet<T>::Add(T* br)
{
	m_errAbs += br->m_errAbs;
	m_quad += br->m_quad;

	m_pq.Add(br);
}

//!
//! Usuwa wszystkie elementy ze zbioru
//!
template <typename T>
void EltSet<T>::Clear(void)
{
	m_pq.Clear();

	m_errAbs = 0;
	m_quad = 0;
}


//!
//! Zwraca wskaznik do prostopadloscian o najwiekszym bledzie BEZWZGLEDNYM. 
//! Zwracany prostopadloscian jest usuwany jest z kolejki priorytetowej.
//!
template <typename T>
T* EltSet<T>::Top(void)
{
	assert(!m_pq.Empty());

	T* br = m_pq.Top(); // Wskaznik do prostopadloscianu o najwiekszym bledzie

	// Usuniecie wkladu prostopadloscianu "br"
	m_errAbs -= br->m_errAbs;
	m_quad -= br->m_quad;

	// Zwracam wskaznik prostopadloscianu "br"
	return br;
}


#endif

