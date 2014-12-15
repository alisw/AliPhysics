//-----------------------------------------------------------------------
// File and Version Information: 
//      $Id: EvtDalitzPlot.cpp,v 1.3 2009-03-16 15:53:27 robbep Exp $
// 
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 1998 Caltech, UCSB
//
// Module creator:
//      Alexei Dvoretskii, Caltech, 2001-2002.
//-----------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

// Global 3-body Dalitz decay kinematics as defined by the mass
// of the mother and the daughters. Spins are not considered.

#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtDalitzPlot.hh"
#include "EvtGenBase/EvtTwoBodyVertex.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtDecayMode.hh"

using namespace EvtCyclic3;


EvtDalitzPlot::EvtDalitzPlot()
  : _mA(0.), _mB(0.), _mC(0.), _bigM(0.),
    _ldel(0.), _rdel(0.)
{}


EvtDalitzPlot::EvtDalitzPlot(double mA, double mB, double mC, double bigM,
			     double ldel, double rdel)
  : _mA(mA), _mB(mB), _mC(mC), _bigM(bigM),
    _ldel(ldel), _rdel(rdel) 
{
  sanityCheck();
}

EvtDalitzPlot::EvtDalitzPlot(const EvtDecayMode& mode, double ldel, double rdel )
{
  _mA = EvtPDL::getMeanMass(EvtPDL::getId(mode.dau(A)));
  _mB = EvtPDL::getMeanMass(EvtPDL::getId(mode.dau(B)));
  _mC = EvtPDL::getMeanMass(EvtPDL::getId(mode.dau(C)));
  _bigM = EvtPDL::getMeanMass(EvtPDL::getId(mode.mother()));

  _ldel = ldel;
  _rdel = rdel;
  
  sanityCheck();
}


EvtDalitzPlot::EvtDalitzPlot(const EvtDalitzPlot& other) 
  : _mA(other._mA), _mB(other._mB), _mC(other._mC), _bigM(other._bigM),
    _ldel(other._ldel), _rdel(other._rdel)
{}


EvtDalitzPlot::~EvtDalitzPlot()
{}


bool EvtDalitzPlot::operator==(const EvtDalitzPlot& other) const
{
  bool ret = false;
  if(_mA == other._mA && 
     _mB == other._mB &&
     _mC == other._mC &&
     _bigM == other._bigM) ret = true;

  return ret;
}


const EvtDalitzPlot* EvtDalitzPlot::clone() const
{
  return new EvtDalitzPlot(*this);
}


void EvtDalitzPlot::sanityCheck() const
{
  if(_mA < 0 || _mB < 0 || _mC < 0 || _bigM <= 0 || _bigM - _mA - _mB - _mC < 0.) {
    
    printf("Invalid Dalitz plot %f %f %f %f\n",_mA,_mB,_mC,_bigM);
    assert(0);
  }
  assert(_ldel <= 0.);
  assert(_rdel >= 0.);
}


double EvtDalitzPlot::m(Index i) const {

  double m = _mA;
  if(i == B) m = _mB;
  else
    if(i == C) m = _mC;

  return m;
}


double EvtDalitzPlot::sum() const
{
  return _mA*_mA + _mB*_mB + _mC*_mC + _bigM*_bigM;
}


double EvtDalitzPlot::qAbsMin(Pair i) const
{
  Index j = first(i);
  Index k = second(i);

  return (m(j) + m(k))*(m(j) + m(k));
}


double EvtDalitzPlot::qAbsMax(Pair i) const
{
  Index j = other(i);
  return (_bigM-m(j))*(_bigM-m(j));
}


double EvtDalitzPlot::qResAbsMin(EvtCyclic3::Pair i) const
{
  return qAbsMin(i) - sum()/3.;
}

double EvtDalitzPlot::qResAbsMax(EvtCyclic3::Pair i) const
{
  return  qAbsMax(i) - sum()/3.;
}

double EvtDalitzPlot::qHelAbsMin(EvtCyclic3::Pair i) const
{
  Pair j = next(i);
  Pair k = prev(i);
  return  (qAbsMin(j) - qAbsMax(k))/2.;
}

double EvtDalitzPlot::qHelAbsMax(EvtCyclic3::Pair i) const
{
  Pair j = next(i);
  Pair k = prev(i);
  return  (qAbsMax(j) - qAbsMin(k))/2.;
}


double EvtDalitzPlot::mAbsMin(Pair i) const
{
  return sqrt(qAbsMin(i));
}


double EvtDalitzPlot::mAbsMax(Pair i) const
{
  return sqrt(qAbsMax(i));
}


// parallel

double EvtDalitzPlot::qMin(Pair i, Pair j, double q) const
{
  if(i == j) return q;

  else {

    // Particle pair j defines the rest-frame
    // 0 - particle common to r.f. and angle calculations
    // 1 - particle belonging to r.f. but not angle
    // 2 - particle not belonging to r.f.

    Index k0 = common(i,j);
    Index k2 = other(j);
    Index k1 = other(k0,k2);

    // Energy, momentum of particle common to rest-frame and angle
    EvtTwoBodyKine jpair(m(k0),m(k1),sqrt(q)); 
    double pk = jpair.p();
    double ek = jpair.e(EvtTwoBodyKine::A,EvtTwoBodyKine::AB);


    // Energy and momentum of the other particle
    EvtTwoBodyKine mother(sqrt(q),m(k2),bigM());
    double ej = mother.e(EvtTwoBodyKine::B,EvtTwoBodyKine::A);
    double pj = mother.p(EvtTwoBodyKine::A);


    // See PDG 34.4.3.1
    return (ek+ej)*(ek+ej) - (pk+pj)*(pk+pj);
  }
}


// antiparallel

double EvtDalitzPlot::qMax(Pair i, Pair j, double q) const
{

  if(i == j) return q;
  else {

    // Particle pair j defines the rest-frame
    // 0 - particle common to r.f. and angle calculations
    // 1 - particle belonging to r.f. but not angle
    // 2 - particle not belonging to r.f.
    
    Index k0 = common(i,j);
    Index k2 = other(j);
    Index k1 = other(k0,k2); 

    // Energy, momentum of particle common to rest-frame and angle
    EvtTwoBodyKine jpair(m(k0),m(k1),sqrt(q)); 
    double ek = jpair.e(EvtTwoBodyKine::A,EvtTwoBodyKine::AB);
    double pk = jpair.p();

    // Energy and momentum of the other particle
    EvtTwoBodyKine mother(sqrt(q),m(k2),bigM());
    double ej = mother.e(EvtTwoBodyKine::B,EvtTwoBodyKine::A);
    double pj = mother.p(EvtTwoBodyKine::A);

    
    // See PDG 34.4.3.1
    return (ek+ej)*(ek+ej) - (pk-pj)*(pk-pj);
  }
}


double EvtDalitzPlot::getArea(int N, Pair i, Pair j) const
{
  // Trapezoidal integral over qi. qj can be calculated.
  // The first and the last point are zero, so they are not counted

  double dh = (qAbsMax(i) - qAbsMin(i))/((double) N);
  double sum = 0;

  int ii;
  for(ii=1;ii<N;ii++) {

    double x = qAbsMin(i) + ii*dh;
    double dy = qMax(j,i,x) - qMin(j,i,x);
    sum += dy;
  }

  return sum * dh;
}


double EvtDalitzPlot::cosTh(EvtCyclic3::Pair i1, double q1, EvtCyclic3::Pair i2, double q2) const
{
  if(i1 == i2) return 1.;
  
  double qmax = qMax(i1,i2,q2);
  double qmin = qMin(i1,i2,q2);
  
  double cos = (qmax + qmin - 2*q1)/(qmax - qmin);
  
  return cos;
}


double EvtDalitzPlot::e(Index i, Pair j, double q) const
{
  if(i == other(j)) {
 
    // i does not belong to pair j

    return (bigM()*bigM()-q-m(i)*m(i))/2/sqrt(q);
  }
  else {
    
    // i and k make pair j

    Index k;
    if(first(j) == i) k = second(j);
    else k = first(j); 

    double e = (q + m(i)*m(i) - m(k)*m(k))/2/sqrt(q);	       
    return e;
  }
}


double EvtDalitzPlot::p(Index i, Pair j, double q) const
{
  double en = e(i,j,q);
  double p2 = en*en - m(i)*m(i);
  
  if(p2 < 0) {
    printf("Bad value of p2 %f %d %d %f %f\n",p2,i,j,en,m(i));
    assert(0);
  }

  return sqrt(p2);  
}


double EvtDalitzPlot::q(EvtCyclic3::Pair i1, double cosTh, EvtCyclic3::Pair i2, double q2) const
{
  if(i1 == i2) return q2;

  EvtCyclic3::Index f = first(i1);
  EvtCyclic3::Index s = second(i1);
  return m(f)*m(f) + m(s)*m(s) + 2*e(f,i2,q2)*e(s,i2,q2) - 2*p(f,i2,q2)*p(s,i2,q2)*cosTh;
}


double EvtDalitzPlot::jacobian(EvtCyclic3::Pair i, double q) const
{
  return 2*p(first(i),i,q)*p(other(i),i,q);  // J(BC) = 2pA*pB = 2pA*pC
}


EvtTwoBodyVertex EvtDalitzPlot::vD(Pair iRes, double m0, int L) const
{
  return EvtTwoBodyVertex(m(first(iRes)),
			  m(second(iRes)),m0,L);
}


EvtTwoBodyVertex EvtDalitzPlot::vB(Pair iRes, double m0, int L) const
{
  return EvtTwoBodyVertex(m0,m(other(iRes)),bigM(),L);
}


void EvtDalitzPlot::print() const
{
  // masses
  printf("Mass  M    %f\n",bigM());
  printf("Mass mA    %f\n",_mA);
  printf("Mass mB    %f\n",_mB);
  printf("Mass mC    %f\n",_mC);
  // limits
  printf("Limits qAB %f : %f\n",qAbsMin(AB),qAbsMax(AB));
  printf("Limits qBC %f : %f\n",qAbsMin(BC),qAbsMax(BC));
  printf("Limits qCA %f : %f\n",qAbsMin(CA),qAbsMax(CA));
  printf("Sum q       %f\n",sum());
  printf("Limit qsum  %f : %f\n",qSumMin(),qSumMax());
}


