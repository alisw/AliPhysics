#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtTwoBodyKine.cpp,v 1.3 2009-03-16 15:37:54 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

#include <iostream>
#include <assert.h>
#include <math.h>
#include "EvtGenBase/EvtTwoBodyKine.hh"
#include "EvtGenBase/EvtReport.hh"
using std::endl;
using std::ostream;


EvtTwoBodyKine::EvtTwoBodyKine()
  : _mA(0.), _mB(0.), _mAB(0.)
{}

EvtTwoBodyKine::EvtTwoBodyKine(double mA, double mB, double mAB)
  : _mA(mA), _mB(mB), _mAB(mAB)
{
  if(mAB < mA + mB) {

    report(Severity::Info,"EvtGen") << mAB << " < " << mA << " + " << mB << endl;
    assert(0);
  }
}

EvtTwoBodyKine::EvtTwoBodyKine(const EvtTwoBodyKine& other)
  : _mA(other._mA), _mB(other._mB), _mAB(other._mAB)
{}

EvtTwoBodyKine::~EvtTwoBodyKine()
{}


double EvtTwoBodyKine::m(Index i) const
{
  double ret = _mAB;
  if(A == i) ret = _mA;
  else
    if(B == i) ret = _mB;
  
  return ret;
}


double EvtTwoBodyKine::p(Index i) const
{ 
  double p0 = 0.;

  if(i == AB) {

    double x = _mAB*_mAB - _mA*_mA - _mB*_mB;
    double y = 2*_mA*_mB;
    p0 = sqrt(x*x - y*y)/2./_mAB;
  }
  else 
    if(i == A) {

      double x = _mA*_mA - _mAB*_mAB - _mB*_mB;
      double y = 2*_mAB*_mB;
      p0 = sqrt(x*x - y*y)/2./_mA;
    }
    else {

      double x = _mB*_mB - _mAB*_mAB - _mA*_mA;
      double y = 2*_mAB*_mA;
      p0 = sqrt(x*x - y*y)/2./_mB;
    }

  return p0;
}


double EvtTwoBodyKine::e(Index i, Index j) const
{
  double ret = m(i);
  if(i != j) {

    double pD = p(j);
    ret = sqrt(ret*ret + pD*pD);
  }
  return ret;
}


void EvtTwoBodyKine::print(ostream& os) const
{
  os << " mA = " << _mA << endl;
  os << " mB = " << _mB << endl;
  os << "mAB = " << _mAB << endl;
}


ostream& operator<<(ostream& os, const EvtTwoBodyKine& p)
{
  p.print(os);
  return os;
}
