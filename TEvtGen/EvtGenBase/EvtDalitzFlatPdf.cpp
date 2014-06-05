#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtDalitzFlatPdf.cpp,v 1.3 2009-03-16 15:53:27 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtDalitzFlatPdf.hh"

EvtDalitzFlatPdf::EvtDalitzFlatPdf(const EvtDalitzPlot& dp)
  : EvtPdf<EvtDalitzPoint>(), _dp(dp)
{}

EvtDalitzFlatPdf::EvtDalitzFlatPdf(const EvtDalitzFlatPdf& other)
  : EvtPdf<EvtDalitzPoint>(other), _dp(other._dp)
{}

EvtDalitzFlatPdf::~EvtDalitzFlatPdf()
{}

EvtPdf<EvtDalitzPoint>* EvtDalitzFlatPdf::clone() const
{
  return new EvtDalitzFlatPdf(*this);
}

double EvtDalitzFlatPdf::pdf(const EvtDalitzPoint&) const
{
  return 1.;
}
  
EvtValError EvtDalitzFlatPdf::compute_integral(int N) const
{
  return EvtValError(_dp.getArea(N),0.);
}

EvtDalitzPoint EvtDalitzFlatPdf::randomPoint()
{
  // To obtain a uniform distribution generate 
  // in terms of q's. Generate in a box that circumscribes the 
  // Dalitz plot. Accept points inside. If there are two 
  // many unsuccessful attempts it's a hint that the Dalitz plot
  // area is tiny compared to the box. It's a pathological
  // case. Abort.
  
  EvtCyclic3::Pair pair1 = EvtCyclic3::BC;
  EvtCyclic3::Pair pair2 = EvtCyclic3::CA;
  
  int n = 0;
  int maxTries = 1000;
  while(n++ < maxTries) {
    
    double q1 = EvtRandom::Flat(_dp.qAbsMin(pair1),_dp.qAbsMax(pair2));
    double q2 = EvtRandom::Flat(_dp.qAbsMin(pair2),_dp.qAbsMax(pair2)); 
    
    EvtDalitzCoord point(pair1,q1,pair2,q2);
    EvtDalitzPoint x(_dp,point);
    
    if(x.isValid()) return x;
  }
  
  printf("No point generated for dalitz plot after %d tries\n",maxTries);
  return EvtDalitzPoint();
}
