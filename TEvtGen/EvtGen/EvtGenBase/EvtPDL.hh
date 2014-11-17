//--------------------------------------------------------------------------
// 
// Environment: 
// This software is part of the EvtGen package developed jointly 
// for the BaBar and CLEO collaborations.  If you use all or part 
// of it, please give an appropriate acknowledgement.
// 
// Copyright Information: See EvtGen/COPYRIGHT 
// Copyright (C) 1998 Caltech, UCSB 
// 
// Module: EvtGen/EvtPDL.hh 
// 
// Description:Class to keep track of particle properties.  
// 
// Modification history: 
//
// DJL/RYD September 25, 1996 Module created 
//
//------------------------------------------------------------------------

#ifndef EVTPDL_HH
#define EVTPDL_HH

#include "EvtGenBase/EvtPartProp.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtStringHash.hh"
#include <vector>
#include <map>

const int SPIN_NAME_LENGTH = 100;

class EvtPDL {

public:

  EvtPDL();  

  ~EvtPDL();

  void read(const char* fname);
  void readPDT(const std::string fname);

  
  static double getMeanMass(EvtId i );
  static double getMass(EvtId i );
  static double getRandMass(EvtId i, EvtId *parId, int nDaug, EvtId *dauId, EvtId *othDaugId,double maxMass, double *dauMasses );
  static double getMassProb(EvtId i, double mass, double massPar, int nDaug, double *massDau);

  static double getMaxMass(EvtId i );
  static double getMinMass(EvtId i );
  //the number we got from PDT
  static double getMaxRange(EvtId i );
  static double getWidth(EvtId i );
  static double getctau(EvtId i );
  static int getStdHep(EvtId id );
  static int getLundKC(EvtId id );

  // Function to retrieve EvtId from PythiaID  
  static EvtId evtIdFromLundKC(int pythiaId );
  static EvtId evtIdFromStdHep(int stdhep );
  static EvtId chargeConj(EvtId id );
  static int chg3(EvtId i );
  static EvtSpinType::spintype getSpinType(EvtId i );
  static EvtId getId(const std::string& name );
  static std::string name(EvtId i);
  static void alias(EvtId num,const std::string& newname);
  static void aliasChgConj(EvtId a,EvtId abar);
  static size_t entries();
  static EvtId getEntry(int i);
  static void reSetMass(EvtId i, double mass);
  static void reSetWidth(EvtId i, double width);
  static void reSetMassMin(EvtId i, double mass);
  static void reSetMassMax(EvtId i,double mass);
  static void reSetBlatt(EvtId i,double blatt);
  static void reSetBlattBirth(EvtId i,double blatt);
  static void includeBirthFactor(EvtId i,bool yesno);
  static void includeDecayFactor(EvtId i,bool yesno);
  static void changeLS(EvtId i, std::string &newLS );
  static void setPWForDecay(EvtId i, int spin, EvtId d1, EvtId d2);
  static void setPWForBirthL(EvtId i, int spin, EvtId par, EvtId othD);
private:

  void setUpConstsPdt();

  static unsigned int    _firstAlias;
  static int    _nentries;

  static std::vector<EvtPartProp>& partlist() {
    static std::vector<EvtPartProp> s_partlist;
    return s_partlist;
  }

  static std::map<std::string, int> _particleNameLookup;
  
}; // EvtPDL.h

#endif


