#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtIntegPdf1D.cpp,v 1.3 2009-03-16 15:48:09 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

#include <assert.h>
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtMacros.hh"
#include "EvtGenBase/EvtIntegPdf1D.hh"

EvtIntegPdf1D::EvtIntegPdf1D(double min, double max)
  : EvtPdf<EvtPoint1D>(), _min(min), _max(max)
{
  assert(min <= max);
}

EvtIntegPdf1D::EvtIntegPdf1D(const EvtIntegPdf1D& other)
  : EvtPdf<EvtPoint1D>(other), _min(other._min), _max(other._max)
{}

EvtIntegPdf1D::~EvtIntegPdf1D() 
{}

EvtValError EvtIntegPdf1D::compute_integral() const
{
  double x1 = pdfIntegral(_min); 
  double x2 = pdfIntegral(_max);
  return EvtValError(x2-x1,0.);
}


EvtPoint1D EvtIntegPdf1D::randomPoint()
{
  double itgmin = pdfIntegral(_min);
  double itgmax = pdfIntegral(_max);
  double itgrnd = EvtRandom::Flat(itgmin,itgmax);

  return EvtPoint1D(_min,_max,pdfIntegralInverse(itgrnd));
}


