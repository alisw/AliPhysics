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
// Module: 
// Description: Form factors for b->sll according to Ali '02 et al.
//              PRD 66 34002
//
// Modification history:
//
//    Ryd     March 30, 2003         Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTOSLLALIFF_HH
#define EVTBTOSLLALIFF_HH

#include "EvtGenModels/EvtbTosllFF.hh"

class EvtId;


class EvtbTosllAliFF : public EvtbTosllFF {

public:

  EvtbTosllAliFF();

  void getScalarFF(EvtId parent, EvtId daught,double t, double mass, 
		   double& fp,double& f0,double& ft);
  void getVectorFF(EvtId parent, EvtId daught,double t, double mass, 
		   double& a1,double& a2,double& a0, double& v,
		   double& t1, double& t2, double& t3 );



private:

};

#endif

