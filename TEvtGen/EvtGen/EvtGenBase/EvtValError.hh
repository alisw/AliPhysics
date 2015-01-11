/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtValError.hh,v 1.2 2009-03-16 16:34:38 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Value and its associated error. E.g. this could be interval size and 
// the error associated with numerical integration.

#ifndef EVT_VAL_ERROR_HH
#define EVT_VAL_ERROR_HH

#include <iostream>
#include <assert.h>
#include <math.h>

class EvtValError {

public:

  EvtValError();
  EvtValError(double val);
  EvtValError(double val, double err);
  EvtValError(const EvtValError& other);
  ~EvtValError();

  inline int valueKnown() const { return _valKnown; }
  inline double value() const { assert(_valKnown); return _val; }
  inline int errorKnown() const { return _errKnown; }
  inline double error() const { assert(_errKnown); return _err; }

  double prec() const;
  void operator=(const EvtValError& other);
  void operator*=(const EvtValError& other);
  void operator/=(const EvtValError& other);
  void operator+=(const EvtValError& other);
  void operator*=(double c);

  void print(std::ostream&) const;

private:

  int    _valKnown;
  double _val;
  int    _errKnown;
  double _err;

};


EvtValError operator*(const EvtValError& x1, const EvtValError& x2);
EvtValError operator/(const EvtValError& x1, const EvtValError& x2);
EvtValError operator+(const EvtValError& x1, const EvtValError& x2);
EvtValError operator*(const EvtValError& x,double c);
EvtValError operator*(double c,const EvtValError& x);

std::ostream& operator<<(std::ostream&, const EvtValError&);

// Perform an accept/reject fraction count

template <class InputIterator, class Predicate>
EvtValError accept_reject(InputIterator it, InputIterator end, Predicate pred)
{
  int itsTried = 0;
  int itsPassed = 0;
  while(it != end) {
    
    itsTried++;
    if(pred(*it++)) itsPassed++;
  }

  return EvtValError(((double) itsPassed)/((double) itsTried),sqrt(itsPassed)/((double) itsTried));
}

#endif
