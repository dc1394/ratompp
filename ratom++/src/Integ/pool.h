#ifndef __RATOM_POOL_H__
#define __RATOM_POOL_H__


//
// Zbigniew Romanowski: 12 Sty 2009
//
// Prosty sposob przydzielania pamieci bazujacy na wektorze
// Elementy sa wykorzystywane jednokrotnie
// Dedykowane dla kolejki priorytetowej i algorytmu calkowania adaptacyjnego


template <typename T>
class Pool 
{
public:
	Pool(size_t initSize) : m_tab(initSize), m_eltId(-1) { }
	~Pool(void) { }
	
	T* New();
	void Clear() { m_eltId = -1; }

private:
	// Tablica pamieci
	std::vector<T> m_tab;

	// Identyfikator w tablicy
	size_t m_eltId;
};


template <typename T>
T* Pool<T>::New()
{
	m_eltId++; 
	if(m_tab.size() - 1 > m_eltId)
		return &(m_tab[m_eltId]); 

	return NULL;
}




#endif 

