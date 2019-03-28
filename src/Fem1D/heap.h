#ifndef __RATOM_HEAP_H__
#define __RATOM_HEAP_H__


/** \brief Heap. Used by adaptive algorithm.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/


template <typename T>
class Heap : private std::vector<T>
{
public:
	Heap(void) { };
	~Heap(void) { };

	void Push(const T& e)
	{
		push_back(e);
		push_heap(this->begin(), this->end());
	}

	void Pop(T& e)
	{
		pop_heap(this->begin(), this->end());
		e = this->back();
		this->pop_back();
	}

	const T& Top() const
	{
		return this->front();
	}

	bool Empty(void) const
	{
		return this->empty();
	}

	size_t Size(void) const
	{
		return this->size();
	}

	void Clear()
	{
		this->clear();
	}

	const T& operator[](size_t i) const
	{
		return std::vector<T>::operator[](i);
	}

        void Reserve(size_t s)
        {
            this->reserve(s);
        }
};

#endif

