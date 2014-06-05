//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2000      Caltech
//
// Module: 
// This is the base class for form factors in b->sll transitions.
//
// Description:
//
// Modification history:
//
//    RYD     January 5, 2000         Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTOSLLFF_HH
#define EVTBTOSLLFF_HH

class EvtId;

class EvtbTosllFF{

 public:

  virtual ~EvtbTosllFF( ) { } ;

  virtual void getScalarFF(EvtId /*parent*/, EvtId /*daught*/,double /*t*/, 
                           double /*mass*/, double& /*fp*/,double& /*f0*/,
			   double& /*ft*/) {return;}
  virtual void getVectorFF(EvtId /*parent*/, EvtId /*daught*/,double /*t*/, 
			   double /*mass*/, double& /*a1*/,double& /*a2*/,
			   double& /*a0*/, double& /*v*/,double& /*t1*/, 
			   double& /*t2*/, double& /*t3*/ ) {return;}

};

#endif
