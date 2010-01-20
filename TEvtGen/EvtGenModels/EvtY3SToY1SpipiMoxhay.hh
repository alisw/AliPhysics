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
// Module: EvtGen/EvtVVpipiMoxhay.hh
//
// Description: This model is based on the proposal by Tuan and Lipkin
//              (Phys.Lett.B206:349-353,1988) and the subsequent model
//              by Moxhay (Phys.Rev.D39:3497,1989) for the dipion spectrum
//              in Y(3S) -> pi+ pi- Y(1S). Please Note: in Moxhay's paper,
//              he wrote the fitted value of the parameter Im(B)/A as
//              -0.2983. However, using his quoted value leads to the wrong
//              spectrum. Changing the sign of his quoted Im(B)/A fixes the
//              shape and reproduces his result. Therefore, please pass
//              Im(B)/A = 0.2983 and Re(B)/A = 0.2196 to get the correct shape
//              based on his fit to the CLEO data.
//
// Example:
//
// Decay  Upsilon(3S)
//  1.0000    Upsilon  pi+  pi-     Y3STOY1SPIPIMOXHAY 0.2196 0.2983;
// Enddecay
//
//   --> the order of parameters is: Re(B)/A Im(B)/A
//
//
// Modification history:
//
//    SEKULA  November 02, 2007         Module created
//
//------------------------------------------------------------------------

#ifndef EVTY3STOY1SPIPIMOXHAY_HH
#define EVTY3STOY1SPIPIMOXHAY_HH

#include "EvtGenBase/EvtDecayProb.hh"

class EvtParticle;

class EvtY3SToY1SpipiMoxhay:public  EvtDecayProb  {

public:

  EvtY3SToY1SpipiMoxhay() {}
  virtual ~EvtY3SToY1SpipiMoxhay();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p); 
  void init();
  void initProbMax();

};

#endif

