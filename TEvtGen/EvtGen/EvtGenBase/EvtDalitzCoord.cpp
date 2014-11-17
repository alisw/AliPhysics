#include "EvtGenBase/EvtPatches.hh"
/*******************************************************************************
 * Project: BaBar detector at the SLAC PEP-II B-factory
 * Package: EvtGenBase
 *    File: $Id: EvtDalitzCoord.cpp,v 1.3 2009-03-16 15:55:13 robbep Exp $
 *  Author: Alexei Dvoretskii, dvoretsk@slac.stanford.edu, 2001-2002
 *
 * Copyright (C) 2002 Caltech
 *******************************************************************************/

#include <assert.h>
#include <iostream>
#include "EvtGenBase/EvtDalitzCoord.hh"
using std::endl;
using std::ostream;
using EvtCyclic3::Pair;


// For coordinates it's good to alway have a
// default ctor. Initialize to something invalid.

EvtDalitzCoord::EvtDalitzCoord()
  : _i1(EvtCyclic3::AB), _i2(EvtCyclic3::BC), _q1(-1.), _q2(-1.)
{}

EvtDalitzCoord::EvtDalitzCoord(const EvtDalitzCoord& other)
  : _i1(other._i1), _i2(other._i2), _q1(other._q1), _q2(other._q2)
{}


EvtDalitzCoord::EvtDalitzCoord(Pair i1, double q1, Pair i2, double q2)
  : _i1(i1), _i2(i2),_q1(q1),_q2(q2)
{} 


EvtDalitzCoord::~EvtDalitzCoord()
{}


bool EvtDalitzCoord::operator==(const EvtDalitzCoord& other) const
{
  return (_i1 == other._i1 && _i2 == other._i2 && 
	  _q1 == other._q1 && _q2 == other._q2);
}

void EvtDalitzCoord::print(ostream& os) const
{
  os << _i1 << " " << _q1 << endl;
  os << _i2 << " " << _q2 << endl;
}


ostream& operator<<(ostream& os,const EvtDalitzCoord& p)
{
  p.print(os);
  return os;
}

