#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtDalitzResPdf.cpp,v 1.3 2009-03-16 15:54:07 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

#include <stdio.h>
#include <math.h>
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtDalitzResPdf.hh"
#include "EvtGenBase/EvtDalitzCoord.hh"
#include "EvtGenBase/EvtRandom.hh"
using namespace EvtCyclic3;

EvtDalitzResPdf::EvtDalitzResPdf(const EvtDalitzPlot& dp, 
				 double _m0, double _g0, EvtCyclic3::Pair pair)
  : EvtPdf<EvtDalitzPoint>(), 
  _dp(dp), _m0(_m0), _g0(_g0), _pair(pair)
{}


EvtDalitzResPdf::EvtDalitzResPdf(const EvtDalitzResPdf& other)
  : EvtPdf<EvtDalitzPoint>(other), 
  _dp(other._dp),_m0(other._m0), _g0(other._g0), _pair(other._pair)
{}


EvtDalitzResPdf::~EvtDalitzResPdf()
{}

EvtValError EvtDalitzResPdf::compute_integral(int N) const
{
  assert(N != 0);
  
  EvtCyclic3::Pair i = _pair;
  EvtCyclic3::Pair j = EvtCyclic3::next(i);

  // Trapezoidal integral

  double dh = (_dp.qAbsMax(j) - _dp.qAbsMin(j))/((double) N);
  double sum = 0;
  
  int ii;
  for(ii=1;ii<N;ii++) {
    
    double x = _dp.qAbsMin(j) + ii*dh;
    double min = (_dp.qMin(i,j,x) - _m0*_m0)/_m0/_g0;
    double max = (_dp.qMax(i,j,x) - _m0*_m0)/_m0/_g0;
    double itg = 1/EvtConst::pi*(atan(max) - atan(min));
    sum += itg;
  }
  EvtValError ret(sum*dh,0.); 
  
  return ret;
}


EvtDalitzPoint EvtDalitzResPdf::randomPoint()
{
  // Random point generation must be done in a box encompassing the 
  // Dalitz plot


  EvtCyclic3::Pair i = _pair;
  EvtCyclic3::Pair j = EvtCyclic3::next(i);  
  double min = 1/EvtConst::pi*atan((_dp.qAbsMin(i) - _m0*_m0)/_m0/_g0);
  double max = 1/EvtConst::pi*atan((_dp.qAbsMax(i) - _m0*_m0)/_m0/_g0);

  int n = 0;
  while(n++ < 1000) {

    double qj = EvtRandom::Flat(_dp.qAbsMin(j),_dp.qAbsMax(j));
    double r = EvtRandom::Flat(min,max);
    double qi = tan(EvtConst::pi*r)*_g0*_m0 + _m0*_m0;
    EvtDalitzCoord x(i,qi,j,qj);
    EvtDalitzPoint ret(_dp,x);
    if(ret.isValid()) return ret;
  }
  
  // All generated points turned out to be outside of the Dalitz plot
  // (in the outer box)
  
  printf("No point generated for dalitz plot after 1000 tries\n");
  return EvtDalitzPoint(0.,0.,0.,0.,0.,0.);
}


double EvtDalitzResPdf::pdf(const EvtDalitzPoint& x) const
{
  EvtCyclic3::Pair i = _pair;
  double dq = x.q(i) - _m0*_m0;
  return 1/EvtConst::pi*_g0*_m0/(dq*dq + _g0*_g0*_m0*_m0);
}


double EvtDalitzResPdf::pdfMaxValue() const
{
  return 1/(EvtConst::pi*_g0*_m0);
}

