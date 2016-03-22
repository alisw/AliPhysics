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
// Module: EvtGenBase/EvtSpinType.cc
//
// Description: Class for enumarating the different types of
//              particles and the number of states they have.
//
// Modification history:
//
//    RYD     Jan 26,  2006         Module created
//
//------------------------------------------------------------------------


#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSpinType.hh"


int EvtSpinType::getSpin2(spintype stype){
  
  switch (stype){
  case SCALAR: case STRING:
    return 0;
  case DIRAC: case NEUTRINO:
    return 1;
  case VECTOR: case PHOTON: 
    return 2;
  case RARITASCHWINGER:
    return 3;
  case TENSOR:
    return 4;
  case SPIN5HALF:
    return 5;
  case SPIN3:
    return 6;
  case SPIN7HALF:
    return 7;
  case SPIN4:
    return 8;
  default:
    report(Severity::Error,"EvtGen")<<"Unknown spintype in EvtSpinType!"<<std::endl;
    return 0;
  }
  
}



int EvtSpinType::getSpinStates(spintype stype){

  switch (stype){
  case SCALAR: case STRING: case NEUTRINO:
    return 1;
  case DIRAC: case PHOTON:
    return 2;
  case VECTOR: 
    return 3;
  case RARITASCHWINGER:
    return 4;
  case TENSOR:
    return 5;
  case SPIN5HALF:
    return 6;
  case SPIN3:
    return 7;
  case SPIN7HALF:
    return 8;
  case SPIN4:
    return 9;
  default:
    report(Severity::Error,"EvtGen")<<"Unknown spintype in EvtSpinType!"<<std::endl;
    return 0;
  }

}
