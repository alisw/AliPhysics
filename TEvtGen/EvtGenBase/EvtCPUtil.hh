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
// Module: EvtGen/EvtCPUtil.hh
//
// Description:Class to hold CP physics utilities.
//
// Modification history:
//
//    RYD     March 24, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTCPUTIL_HH
#define EVTCPUTIL_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtPatches.hh"
class EvtParticle;
class EvtId;

class EvtCPUtil{

public:

  static void fractB0CP(EvtComplex Af, EvtComplex Abarf, 
			double deltam, double beta, double &fract);

  static void fractB0nonCP(EvtComplex Af, EvtComplex Abarf, 
			   EvtComplex Afbar, EvtComplex Abarfbar, 
			   double deltam, double beta, int flip, 
			   double &fract);

  static void OtherB(EvtParticle *p, double &t, EvtId &otherb);

  static void OtherB(EvtParticle *p, double &t, EvtId &otherb, double probB0);

  //id is the produced particle
  //t returns the lifetime of the particle
  //and mix will be 1 if it mixed otherwise 0
  static void incoherentMix(const EvtId id, double &t, int &mix);



};


#endif

