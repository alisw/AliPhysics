#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtPoint1D.cpp,v 1.3 2009-03-16 15:44:41 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

// Point on a finite 1-D interval. isValid shows whether for a given specification,
// the coordinate _value is inside the interval defined by _min, _max.

#include <stdio.h>
#include "EvtGenBase/EvtPoint1D.hh"

EvtPoint1D::EvtPoint1D()
  : _min(0.), _max(-1.), _value(0.), _valid(false)
{}

EvtPoint1D::EvtPoint1D(double value)
  : _min(0.), _max(-1.), _value(value), _valid(true)
{}

EvtPoint1D::EvtPoint1D(double min, double max, double value)
  : _min(min), _max(max), _value(value), _valid((_min <= _value && _value <= _max) ? true : false)
{} 
  
EvtPoint1D::~EvtPoint1D()
{}

void EvtPoint1D::print() const
{
  printf("%f (%f : %f)\n",_value,_min,_max);
}

