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
// Module: EvtGen/EvtYmSToYnSpipiCLEO.hh
//
// Description: This model is based on matrix element method used by
//              CLEO in Phys.Rev.D76:072001,2007 to model the dipion mass
//              and helicity angle distribution in the decays Y(mS) -> pi pi Y(nS),
//              where m,n are integers and m>n and m<4.
//              This model has two parameters, Re(B/A) and Im(B/A), which
//              are coefficients of the matrix element's terms determined by
//              the CLEO fits.
//
// Example:
//
// Decay  Upsilon(3S)
//  1.0000    Upsilon      pi+  pi-     YMSTOYNSPIPICLEO -2.523 1.189;
// Enddecay
// Decay  Upsilon(3S)
//  1.0000    Upsilon(2S)  pi+  pi-     YMSTOYNSPIPICLEO -0.395 0.001;
// Enddecay
// Decay  Upsilon(2S)
//  1.0000    Upsilon      pi+  pi-     YMSTOYNSPIPICLEO -0.753 0.000;
// Enddecay
//
//   --> the order of parameters is: Re(B/A) Im(B/A)
//
// Modification history:
//
//    SEKULA  Jan. 28, 2008         Module created
//
//------------------------------------------------------------------------

#ifndef EVTYMSTOYNSPIPICLEO_HH
#define EVTYMSTOYNSPIPICLEO_HH

// #include "EvtGenBase/EvtDecayProb.hh"
#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtYmSToYnSpipiCLEO:public  EvtDecayAmp  {
  //EvtDecayProb  {

public:

  EvtYmSToYnSpipiCLEO() {}
  virtual ~EvtYmSToYnSpipiCLEO();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p); 
  void init();
  void initProbMax();

};

#endif

