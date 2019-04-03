#ifndef __RATOM_RIORITYQUEUE_H__
#define __RATOM_RIORITYQUEUE_H__



//
// Zbigniew Romanowski: 06 Sie 2008
//
// Implementacja kolejki priorytetowej opartej na kopcu. 
// Elementy kolejki przechowywane sa w tablicy. Tym samym zuzywana jest tylko pamiec na przechowywanie elementow kolejki
// i nie ma dodatkowych narzutow zwiazanych z przechowywaniem wskaznikow (tak jak w drzewie).
//
// Opis tego algorytmu znajduje sie w ksiazkach:
// 1. R. Sedgewick "Algorytmy w C++", 1999, wydawnictwo RM
// 2. T. H. Cormen, Ch. E. Leiserson, R. L. Rivest "Wprowadzenie do algorytmow", 1998, WNT
//
// Zakladam, ze:
//		1. Typ T ma funkcje Pri() const, ktora zwraca priorytet elementu
//		2. Typ zwracany przez funkcje Pri() mozna porownywac operatorem mniejszosci "<"
//

// #include "stdafx.h"

template <typename T>
class PriorityQueue
{ 
public:
	PriorityQueue(size_t initSize);
	~PriorityQueue(void);

	bool Empty(void) const;
	void Add(T* t);
	T* Top(void);
	void Clear(void);
	size_t EltNo(void) const;

	size_t MaxSize(void) const { return m_tab.size(); }

private:
	void FixUp(size_t k);
	void FixDown(size_t k);
	void Exch(size_t n, size_t k);

private:

	// Tablica zawierajca wskazniki na elementy kolejki.
	//	Wykorzystywane sa elementy zaczynajace sie od indeksu jeden!
	std::vector<T*> m_tab;

	// Liczba elementow w kolejce
	size_t m_eltNo;

};

//
// Konstruktor
// initSize - poczatkowy rozmiar kolejki
//
template <typename T>
PriorityQueue<T>::PriorityQueue(size_t initSize) : m_tab(initSize), m_eltNo(0)
{
	// printf("PriorityQueue::size = %d\n", m_tab.size());
	assert(initSize > 1);
}


template <typename T>
PriorityQueue<T>::~PriorityQueue(void)
{
}

//
// Zamienia elementy w tablicy. Zamienia element "n" na element "k"
//
template <typename T>
void PriorityQueue<T>::Exch(size_t n, size_t k)
{
//	assert(n != k);

	T* tmp = m_tab[n];
	m_tab[n] = m_tab[k];
	m_tab[k] = tmp;
}

//
// Ukopcowanie zstepujace (z gory do dolu)
//
template <typename T>
void PriorityQueue<T>::FixDown(size_t k)
{
const size_t n = m_eltNo;
size_t j;

	while(2 * k <= n)
	{
		j = 2 * k;
		if(j < n && m_tab[j]->Pri() < m_tab[j + 1]->Pri())
			j++;
		if(!(m_tab[k]->Pri() < m_tab[j]->Pri()))
			break;
		Exch(k, j);
		k = j;
	}

}

//
// Ukopcowanie wstepujace (z dolu do gory)
//
template <typename T>
void PriorityQueue<T>::FixUp(size_t k)
{
	while(k > 1 && m_tab[k / 2]->Pri() < m_tab[k]->Pri())
	{
		Exch(k, k / 2);
		k = k / 2;
	}
}

//
// Zwraca "true", jezeli kolejka jest pusta
//
template <typename T>
bool PriorityQueue<T>::Empty(void) const
{
	return (m_eltNo == 0);
}

//
// Dodanie elementu "t" do kolejki
//
template <typename T>
void PriorityQueue<T>::Add(T* t)
{
	// Musi byc odpowiedni rozmiar tablicy
	assert(m_eltNo < m_tab.size());

	m_eltNo++;
	m_tab[m_eltNo] = t;
	FixUp(m_eltNo);
}

//
// Zwraca element o najwiekszej wartosci jednoczescie usuwa go z kolejki
//
template <typename T>
T* PriorityQueue<T>::Top(void)
{
T* ret = m_tab[1]; // Ten element ma najwiekszy priorytet

	assert(m_eltNo > 0); // Kolejka nie moze byc pusta
	
	Exch(1, m_eltNo);
	m_eltNo--;
	FixDown(1);

	return ret;
}

//
// Usuwa wszystkie elementy z kolejki
//
template <typename T>
void PriorityQueue<T>::Clear(void)
{
	// m_tab.resize(1);
	m_eltNo = 0;
}

//
// Zwraca liczbe elementow w kolejce
//
template <typename T>
size_t PriorityQueue<T>::EltNo(void) const
{
	return m_eltNo;
}





#endif // PRIORITYQUEUE_H
