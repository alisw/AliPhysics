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
// Module: EvtId.cc
//
// Description: Class for particle Id used in EvtGen.
//
// Modification history:
//
//    RYD     May 26, 1998        Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtPDL.hh"
#include <iostream>
using std::ostream;


ostream& operator<<(ostream& s, const EvtId& id){

  s<<"(Id="<<id._id<<" Alias="<<id._alias<<")";

  return s;

}


int EvtId::isConjugate(const EvtId & id) const {
  return EvtPDL::getStdHep(*this) == - EvtPDL::getStdHep(id);
}

std::string EvtId::getName() const {

  return EvtPDL::name(*this);

}
