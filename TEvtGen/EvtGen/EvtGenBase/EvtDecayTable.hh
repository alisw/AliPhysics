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

class EvtDecayTable{

public:

  static EvtDecayTable* getInstance();

  int getNMode(int ipar);

  EvtDecayBase* getDecay(int ipar, int imode);

  void readDecayFile(const std::string dec_name, bool verbose=true);
  void readXMLDecayFile(const std::string dec_name, bool verbose=true);

  bool stringToBoolean(std::string valStr);
  void checkParticle(std::string particle);

  int findChannel(EvtId parent,std::string model,int ndaug, 
		  EvtId *daugs,
		  int narg, std::string *args);
  
  int inChannelList(EvtId parent, int ndaug, EvtId *daugs);

  EvtDecayBase *getDecayFunc(EvtParticle *p);

  void printSummary();

  void checkConj();

  std::vector<EvtParticleDecayList> getDecayTable() {return _decaytable;};

  EvtDecayBase* findDecayModel(int aliasInt, int modeInt);
  EvtDecayBase* findDecayModel(EvtId id, int modeInt);

  bool hasPythia(int aliasInt);
  bool hasPythia(EvtId id);

  int getNModes(int aliasInt);
  int getNModes(EvtId id);

  std::vector<std::string> splitString(std::string& theString, 
				       std::string& splitter);

protected:  

  EvtDecayTable();
  ~EvtDecayTable();

private:

  std::vector<EvtParticleDecayList> _decaytable;

  EvtDecayTable(const EvtDecayTable&) {};
  //EvtDecayTable& operator=(const EvtDecayTable&) {};

};

#endif


