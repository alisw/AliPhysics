#ifndef CUSTOMQUEUE_H
#define CUSTOMQUEUE_H

#include <deque>

template <class T, int depth = 10>
class CustomQueue
{
public:
	CustomQueue() : fCollection(){};

	void Fill(T &object)
	{
		if (IsFull())
		{
			fCollection.pop_front();
		}
		fCollection.push_back(object);
	};

	T& GetElement(int index)
	{
		if (index < 0 || index >=  GetSize())
		{
			::Fatal("CustomQueue::GetElement", "index out of range");
		} else {
			return fCollection[index];
		}
	}

	int GetSize() { return (int)fCollection.size(); }

	int GetDepth() { return depth; }

	bool IsFull() { return (int)fCollection.size() == depth; }

private:
	std::deque<T> fCollection;
};

#endif
