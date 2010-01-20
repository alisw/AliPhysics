//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtComplex.cc
//
// Description: EvtComlex.cc
//
// Modification history:
//
//    RYD  December 5, 1998            Created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <math.h>
#include "EvtGenBase/EvtComplex.hh"
using std::ostream;


ostream& operator<<(ostream& s, const EvtComplex& c){

  s<<"("<<c._rpart<<","<<c._ipart<<")";
  return s;

}

EvtComplex& EvtComplex::operator*=(EvtComplex c){

  double r=_rpart*c._rpart-_ipart*c._ipart;
  double i=_rpart*c._ipart+_ipart*c._rpart;
  
  _rpart=r;
  _ipart=i;

  return *this;

}

EvtComplex& EvtComplex::operator/=(EvtComplex c){

  double inv=1.0/(c._rpart*c._rpart+c._ipart*c._ipart);

  double r=inv*(_rpart*c._rpart+_ipart*c._ipart);
  double i=inv*(_rpart*c._ipart-_ipart*c._rpart);

  _rpart=r;
  _ipart=i;

  return *this;

}




