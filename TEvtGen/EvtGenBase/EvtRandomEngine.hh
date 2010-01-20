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
// Module: EvtGen/EvtRandom.hh
//
// Description:Class to generate random numbers. Single member
//             function random is expected to return a random
//             number in the range ]0..1[.
//
// Modification history:
//
//    RYD     December 25, 1999         Module created
//    RYD     October 2, 2006           Converted to a pure interface class 
//
//------------------------------------------------------------------------

#ifndef EVTRANDOMENGINE_HH
#define EVTRANDOMENGINE_HH

#include <TRandom.h>

class EvtRandomEngine{

public:

  virtual ~EvtRandomEngine() {};

  virtual double random()=0;

private:

};

//define class for generating random numbers (I put this in place of the commented part)
class EvtNUMRandomEngine:public EvtRandomEngine{
public:
  double random(){return  gRandom->Rndm();}
};


//old lines: those are in the macro testEvtGen.cc    
/*
//define class for generating random numbers
class EvtCLHEPRandomEngine:public EvtRandomEngine{
public:
  double random();
};

#ifdef NOCLHEPNAMESPACE
namespace CLHEP {
typedef HepJamesRandom HepJamesRandom ;
}
#endif

double EvtCLHEPRandomEngine::random(){
  static CLHEP::HepJamesRandom randengine;
  return randengine.flat();
}
*/

#endif


