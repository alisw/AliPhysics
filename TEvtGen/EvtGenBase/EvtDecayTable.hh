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
// Module: EvtGen/EvtDecayTable.hh
//
// Description: Class to read in and handle the decays available
//              to EvtGen for each particle, and the model to be
//              used for each one.
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------

#ifndef EVTDECAYTABLE_HH
#define EVTDECAYTABLE_HH

#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtParticleDecayList.hh"
#include <vector>

class EvtId;


typedef EvtDecayBase* EvtDecayBasePtr;

class EvtDecayTable{

public:

  static int getNMode(int ipar);

  static EvtDecayBase* getDecay(int ipar, int imode);

  static void readDecayFile(const std::string dec_name, bool verbose=true);

  static int findChannel(EvtId parent,std::string model,int ndaug, 
			 EvtId *daugs,
			 int narg, std::string *args);
  
  static int inChannelList(EvtId parent, int ndaug, EvtId *daugs);

  static EvtDecayBase *getDecayFunc(EvtParticle *p);

  static void printSummary();

  static void checkConj();

  static const std::vector<EvtParticleDecayList> & decaytable() {return _decaytable;};

private:

  static std::vector<EvtParticleDecayList> _decaytable;

};

#endif


