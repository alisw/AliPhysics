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
// Module: EvtGen/EvtRadCorr.hh
//
// Description:FSR interface 
//
// Modification history:
//
//    Lange   April 27, 2002 Created
//
//------------------------------------------------------------------------

#ifndef EVTRADCORR_HH
#define EVTRADCORR_HH


class EvtAbsRadCorr;
class EvtParticle;

class EvtRadCorr{

public:
  EvtRadCorr();
  ~EvtRadCorr();

  static void doRadCorr(EvtParticle *p);
  
  //This class does not take ownership of the fsr engine;
  //the caller needs to make sure that the engine is not
  //destroyed.
  static void setRadCorrEngine(EvtAbsRadCorr* fsrEngine);
  static bool alwaysRadCorr();
  static bool neverRadCorr();
  static void setAlwaysRadCorr();
  static void setNeverRadCorr();
  static void setNormalRadCorr();

private:

  static EvtAbsRadCorr* _fsrEngine;
  static bool _alwaysRadCorr;
  static bool _neverRadCorr;
};

#endif

