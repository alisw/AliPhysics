#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtValError.cpp,v 1.3 2009-03-16 15:39:28 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

#include <assert.h>
#include <math.h>
#include <iostream>
#include "EvtGenBase/EvtValError.hh"
using std::endl;
using std::ostream;

EvtValError::EvtValError() 
  : _valKnown(0), _val(0.), _errKnown(0), _err(0.)
{}

EvtValError::EvtValError(double val)
  : _valKnown(1), _val(val), _errKnown(0), _err(0.)
{}

EvtValError::EvtValError(double val, double err)
  : _valKnown(1), _val(val), _errKnown(1), _err(err)
{}
  
EvtValError::EvtValError(const EvtValError& other) 
  : _valKnown(other._valKnown), _val(other._val), 
  _errKnown(other._errKnown), _err(other._err)
{}

EvtValError::~EvtValError()
{}

double EvtValError::prec() const 
{ 
  assert(_valKnown && _errKnown); 
  return ( _val != 0) ? _err/_val : 0; 
}

void EvtValError::operator=(const EvtValError& other)
{
  _valKnown = other._valKnown;
  _val = other._val;
  _errKnown = other._errKnown;
  _err = other._err;
}

void EvtValError::operator*=(const EvtValError& other)
{
  assert(_valKnown && other._valKnown);

  // Relative errors add in quadrature
  if(_errKnown && other._errKnown)
    _err = _val * other._val * sqrt(prec()*prec() + other.prec() * other.prec());
  else _errKnown = 0;
  
  // Modify the value  
  _val *= other._val;
}

void EvtValError::operator/=(const EvtValError& other)
{
  assert(_valKnown && other._valKnown && other._val != 0.);

  // Relative errors add in quadrature
  if(_errKnown && other._errKnown)
    _err = _val/other._val * sqrt(prec()*prec() + other.prec() * other.prec());
  else _errKnown = 0;
  
  // Modify the value  
  _val /= other._val;
}


void EvtValError::print(ostream& os) const
{
  if(_valKnown) os << _val;
  else os << "Undef";
  os << " +/- ";
  if(_errKnown) os << _err;
  else os << "Undef";
  os << endl;
}


void EvtValError::operator+=(const EvtValError& other)
{
  assert(_valKnown); assert(other._valKnown);
  _val += other._val;
  
    // add errors in quadrature
  
  if(_errKnown && other._errKnown) {

    _err = sqrt(_err*_err + other._err*other._err);
  }
  else {
    
      _errKnown = 0;
  }
}

void EvtValError::operator*=(double c) {
  
  assert(_valKnown);
  _val *= c;
  if(_errKnown) _err*=c;
}


EvtValError operator*(const EvtValError& x1, const EvtValError& x2)
{
  EvtValError ret(x1);
  ret *= x2;
  return ret;
}

EvtValError operator/(const EvtValError& x1, const EvtValError& x2)
{
  EvtValError ret(x1);
  ret /= x2;
  return ret;
}


EvtValError operator+(const EvtValError& x1, const EvtValError& x2)
{
  EvtValError ret(x1);
  ret += x2;
  return ret;
}


EvtValError operator*(const EvtValError& x,double c) 
{
  EvtValError ret(x);
  ret*=c;
  return ret;
}


EvtValError operator*(double c,const EvtValError& x) 
{
  EvtValError ret(x);
  ret*=c;
  return ret;
}


ostream& operator<<(ostream& os, const EvtValError& other)
{
  other.print(os);
  return os;
}


