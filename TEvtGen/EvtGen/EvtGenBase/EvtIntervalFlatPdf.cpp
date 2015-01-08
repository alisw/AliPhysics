#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtIntervalFlatPdf.cpp,v 1.3 2009-03-16 15:51:08 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

#include "EvtGenBase/EvtPatches.hh"
#include <assert.h>
#include "EvtGenBase/EvtIntervalFlatPdf.hh"
#include "EvtGenBase/EvtRandom.hh"

EvtIntervalFlatPdf::EvtIntervalFlatPdf(double min, double max)
  : EvtPdf<EvtPoint1D>(), _min(min), _max(max)
{
  assert(max >= min);
}

EvtIntervalFlatPdf::EvtIntervalFlatPdf(const EvtIntervalFlatPdf& other)
  : EvtPdf<EvtPoint1D>(other), _min(other._min), _max(other._max)
{}

EvtIntervalFlatPdf::~EvtIntervalFlatPdf()
{}

EvtPdf<EvtPoint1D>* EvtIntervalFlatPdf::clone() const
{
  return new EvtIntervalFlatPdf(*this);
}

double EvtIntervalFlatPdf::pdf(const EvtPoint1D&) const
{
  return 1.;
}
  
EvtValError EvtIntervalFlatPdf::compute_integral() const
{
  return EvtValError(_max-_min,0.);
}
  
EvtPoint1D EvtIntervalFlatPdf::randomPoint()
{
  return EvtPoint1D(_min,_max,EvtRandom::Flat(_min,_max));
}
