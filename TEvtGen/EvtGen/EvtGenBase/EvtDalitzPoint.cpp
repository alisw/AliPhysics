#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtDalitzPoint.cpp,v 1.3 2009-03-16 15:53:27 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

#include "EvtGenBase/EvtPatches.hh"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "EvtGenBase/EvtDalitzPoint.hh"
using namespace EvtCyclic3;

EvtDalitzPoint::EvtDalitzPoint() 
  : _mA(-1.), _mB(-1.), _mC(-1.), _qAB(-1.), _qBC(-1.), _qCA(-1.)
{}

EvtDalitzPoint::EvtDalitzPoint(double mA, double mB, double mC, double qAB, double qBC, double qCA)
  : _mA(mA), _mB(mB), _mC(mC), _qAB(qAB), _qBC(qBC), _qCA(qCA)
{}

// Constructor from Zemach coordinates

EvtDalitzPoint::EvtDalitzPoint(double mA, double mB, double mC, 
			       EvtCyclic3::Pair i, 
			       double qres, double qhel, double qsum)
  : _mA(mA), _mB(mB), _mC(mC)
{
  double qi = qres + qsum/3.;
  double qj = -qres/2. + qhel + qsum/3.;
  double qk = -qres/2. - qhel + qsum/3.;
  
  if(i == AB) { _qAB = qi; _qBC = qj; _qCA = qk; }
  else if(i == BC) { _qAB = qk; _qBC = qi; _qCA = qj; }
  else if(i == CA) { _qAB = qj; _qBC = qk; _qCA = qi; }
}

EvtDalitzPoint::EvtDalitzPoint(const EvtDalitzPlot& dp, const EvtDalitzCoord& x)
  : _mA(dp.m(A)), _mB(dp.m(B)), _mC(dp.m(C))
{
  if(x.pair1() == AB) _qAB = x.q1();
  else
    if(x.pair2() == AB) _qAB = x.q2();
    else _qAB = dp.sum() - x.q1() - x.q2();

  if(x.pair1() == BC) _qBC = x.q1();
  else
    if(x.pair2() == BC) _qBC = x.q2();
    else _qBC = dp.sum() - x.q1() - x.q2();

  if(x.pair1() == CA) _qCA = x.q1();
  else
    if(x.pair2() == CA) _qCA = x.q2();
    else _qCA = dp.sum() - x.q1() - x.q2();
			
}

EvtDalitzPoint::EvtDalitzPoint(const EvtDalitzPoint& other)
  : _mA(other._mA), _mB(other._mB), _mC(other._mC),
    _qAB(other._qAB), _qBC(other._qBC), _qCA(other._qCA)
{}

EvtDalitzPoint::~EvtDalitzPoint()
{}

double EvtDalitzPoint::q(EvtCyclic3::Pair i) const
{
  double ret = _qAB;
  if(BC == i) ret = _qBC;
  else
    if(CA == i) ret = _qCA;
  
  return ret;
}

double EvtDalitzPoint::m(EvtCyclic3::Index i) const
{
  double ret = _mA;
  if(B == i) ret = _mB;
  else
    if(C == i) ret = _mC;
  
  return ret;
}

// Zemach variables

double EvtDalitzPoint::qres(EvtCyclic3::Pair i) const
{
  return (2.*q(i) - q(EvtCyclic3::prev(i)) - q(EvtCyclic3::next(i)))/3.;
}
double EvtDalitzPoint::qhel(EvtCyclic3::Pair i) const
{
  Pair j = next(i);
  Pair k = prev(i);
  return (q(j) - q(k))/2.;
}
double EvtDalitzPoint::qsum() const
{
  return _qAB + _qBC + _qCA;
}


double EvtDalitzPoint::qMin(EvtCyclic3::Pair i, EvtCyclic3::Pair j) const
{
  EvtDalitzPlot dp = getDalitzPlot();
  return dp.qMin(i,j,q(j));
}

double EvtDalitzPoint::qMax(EvtCyclic3::Pair i, EvtCyclic3::Pair j) const
{
  EvtDalitzPlot dp = getDalitzPlot();
  return dp.qMax(i,j,q(j));
}
  
double EvtDalitzPoint::pp(EvtCyclic3::Index i, EvtCyclic3::Index j) const 
{
  if(i == j) return m(i)*m(i);  
  else return (q(combine(i,j)) - m(i)*m(i) - m(j)*m(j))/2.;
}

double EvtDalitzPoint::e(EvtCyclic3::Index i, EvtCyclic3::Pair j) const
{ 
  EvtDalitzPlot dp = getDalitzPlot();
  return dp.e(i,j,q(j)); 
}

double EvtDalitzPoint::p(EvtCyclic3::Index i, EvtCyclic3::Pair j) const
{ 
  EvtDalitzPlot dp = getDalitzPlot();
  return dp.p(i,j,q(j)); 
}

double EvtDalitzPoint::cosTh(EvtCyclic3::Pair pairAng, EvtCyclic3::Pair pairRes) const
{
  EvtDalitzPlot dp = getDalitzPlot();
  return dp.cosTh(pairAng,q(pairAng),pairRes,q(pairRes));
}

  
EvtDalitzCoord EvtDalitzPoint::getDalitzPoint(EvtCyclic3::Pair i, EvtCyclic3::Pair j) const
{
  return EvtDalitzCoord(i,q(i),j,q(j));
}


EvtDalitzPlot EvtDalitzPoint::getDalitzPlot() const
{
  return EvtDalitzPlot(_mA,_mB,_mC,bigM());
}

bool EvtDalitzPoint::isValid() const
{ 
  // Check masses

  double M = bigM();
  if(_mA < 0 || _mB < 0 || _mC < 0 || M <= 0) return false;
  if(M < _mA + _mB + _mC) return false;

  // Check that first coordinate is within absolute limits
 
  bool inside = false; 
  EvtDalitzPlot dp = getDalitzPlot();

  if(dp.qAbsMin(AB) <= _qAB && _qAB <= dp.qAbsMax(AB)) 
    if(qMin(BC,AB) <= _qBC && _qBC <= qMax(BC,AB))
      inside = true;
  
  return inside;
};

double EvtDalitzPoint::bigM() const
{
  return sqrt(_qAB+_qBC+_qCA - _mA*_mA - _mB*_mB - _mC*_mC);
}


void EvtDalitzPoint::print() const
{
  getDalitzPlot().print();
  printf("%f %f %f\n",_qAB,_qBC,_qCA);
}


