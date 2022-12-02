#ifndef CUSTOMQUEUE_H
#define CUSTOMQUEUE_H

#include "TError.h"
#include <deque>

template <class T>
class CustomQueue
{
public:
  CustomQueue() : fCollection(), fDepth(5){};

  CustomQueue(unsigned int depth) : fCollection()
  {
    fDepth = depth;
  }

  void Fill(T &object)
  {
    if (IsFull())
    {
      fCollection.pop_front();
    }
    fCollection.push_back(object);
  };

  T &GetElement(int index)
  {
    if (index < 0 || index >= GetSize())
    {
      ::Fatal("CustomQueue::GetElement", "index out of range");
      return fCollection[0];
    }
    else
    {
      return fCollection[index];
    }
  }

  int GetSize() { return (int)fCollection.size(); }

  void SetDepth(unsigned int depth)
  {
    fDepth = depth;
    fCollection.clear();
  }

  int GetDepth() { return (int)fDepth; }

  bool IsFull() { return fCollection.size() == fDepth; }

  bool IsEmpty() { return (int)fCollection.size() == 0; }

private:
  std::deque<T> fCollection;
  unsigned int fDepth;

  ClassDefNV(CustomQueue, 1);
};

#endif
