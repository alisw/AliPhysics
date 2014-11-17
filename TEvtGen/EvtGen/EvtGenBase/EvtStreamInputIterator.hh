/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtStreamInputIterator.hh,v 1.2 2009-03-16 16:41:09 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Adapters are used to convert various types of input streams 
// into an iteratable interface.

#ifndef EVT_STREAM_INPUT_ITERATOR_HH
#define EVT_STREAM_INPUT_ITERATOR_HH

#include "EvtGenBase/EvtStreamAdapter.hh"

#include <iterator>
#include <cstddef>

using std::input_iterator_tag;

template <class Point>
class EvtStreamInputIterator {  
public:

  typedef input_iterator_tag iterator_category;
  typedef Point              value_type;
  typedef ptrdiff_t          difference_type;
  typedef const Point*       pointer;
  typedef const Point&       reference;

  EvtStreamInputIterator() 
    : _counter(0)
  {}
  
  EvtStreamInputIterator(const EvtStreamInputIterator& other) 
    : _counter(other._counter ? other._counter->clone() : 0),
      _currentValue(other._currentValue)
  {} 

  EvtStreamInputIterator(EvtStreamAdapter<Point>& counter)
    : _counter(counter.clone())
  {
    _currentValue = _counter->currentValue(); 
  }
  
  ~EvtStreamInputIterator()
  {
    if(_counter) delete _counter;
  }
  
  reference operator*() const 
  {
    return _currentValue; 
  }
  
  EvtStreamInputIterator& operator++() 
  {
    _read();
    return *this;
  }
  
  EvtStreamInputIterator operator++(int)
  {
    EvtStreamInputIterator tmp = *this;
    _read();
    return tmp;
  }
  
  bool operator==(const EvtStreamInputIterator& other) const 
  {
    // Equality is only defined for two past the end iterators
    return (pastEnd() && other.pastEnd());
  }
  
protected:

  EvtStreamAdapter<Point>* _counter;
  value_type _currentValue;
  
  bool pastEnd() const 
  {
    bool ret = true;
    if(_counter) ret = _counter->pastEnd();
    return ret;
  }

  // Advances the iterator

  void _read() {
    
    _counter->advance();
    _currentValue = _counter->currentValue();
  }
};


// For adaptable generators these shorthand functions can be used 
// to construct iterators.

template <class Generator>
EvtStreamInputIterator<typename Generator::result_type> iter(Generator gen, int N = 0)
{
  typedef typename Generator::result_type Point;
  EvtGenStreamAdapter<Point,Generator> counter(gen,N);
  return EvtStreamInputIterator<Point>(counter);
}


#endif



