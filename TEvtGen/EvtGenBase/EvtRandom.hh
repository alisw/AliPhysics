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
// Description:Source of random numbers for EvtGen code
//
// Modification history:
//
//    RYD     March 24, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTRANDOM_HH
#define EVTRANDOM_HH


class EvtRandomEngine;

class EvtRandom{

public:

  static double Flat();
  static double Flat(double max);
  static double Flat(double min, double max);

  //generate unit Gaussian
  static double Gaussian();

  static double random();
  
  //This class does not take ownership of the random engine;
  //the caller needs to make sure that the engine is not
  //destroyed.
  static void setRandomEngine(EvtRandomEngine* randomEngine);

private:

  static EvtRandomEngine* _randomEngine;

};

#endif

