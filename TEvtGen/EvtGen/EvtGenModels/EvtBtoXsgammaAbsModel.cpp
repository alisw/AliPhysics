//--------------------------------------------------------------------------
//
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Module: EvtBtoXsgammaAbsModel.cc
//
// Description:
//      B->Xs gamma model base class.
//
// Modification history:
//
//    Jane Tinslay            March 21, 2000      Module Created
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"

#include <stdlib.h>
#include "EvtGenModels/EvtBtoXsgammaAbsModel.hh"

EvtBtoXsgammaAbsModel::~EvtBtoXsgammaAbsModel() {}

void EvtBtoXsgammaAbsModel::init(int, double*) {

  //This default version of init does nothing;
  //A specialized version of this function can be
  //supplied for each decay model to do initialization.

  return;

}
