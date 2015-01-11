#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtBreitWignerPdf.cpp,v 1.3 2009-03-16 15:55:55 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Breit-Wigner shape PDF. If the width is zero it degenerates into a delta 
// function. The integral and its inverse can be still evaluated. 

#include <assert.h>
#include <stdio.h>
#include <math.h>
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtBreitWignerPdf.hh"
#include "EvtGenBase/EvtConst.hh"

EvtBreitWignerPdf::EvtBreitWignerPdf(double min, double max, double m0, double g0)
  : EvtIntegPdf1D(min,max), _m0(m0), _g0(g0)
{}


EvtBreitWignerPdf::EvtBreitWignerPdf(const EvtBreitWignerPdf& other)
  : EvtIntegPdf1D(other), _m0(other._m0), _g0(other._g0)
{}


EvtBreitWignerPdf::~EvtBreitWignerPdf()
{}


double EvtBreitWignerPdf::pdf(const EvtPoint1D& x) const
{
  double m = x.value();
  if((0 == (m - _m0)) && (0. == _g0)) {

    printf("Delta function Breit-Wigner\n");
    assert(0);
  }
  
  double ret = _g0/EvtConst::twoPi/((m-_m0)*(m-_m0)+_g0*_g0/4);
  
  return ret;
}


double EvtBreitWignerPdf::pdfIntegral(double m) const
{
  double itg = 0;
  if(_g0 == 0) {

    if(m > _m0) itg = 1.;
    else
      if(m < _m0) itg = 0.;
      else
	itg = 0.5;
  }
  else itg = atan((m-_m0)/(_g0/2.))/EvtConst::pi + 0.5; 

  return itg;
}


double EvtBreitWignerPdf::pdfIntegralInverse(double x) const
{
  if(x < 0 || x > 1) {
    
    printf("Invalid integral value %f\n",x);
    assert(0);
  }

  double m = _m0;
  if(_g0 != 0) m = _m0 + (_g0/2.)*tan(EvtConst::pi*(x-0.5));
 
  return m;
}




