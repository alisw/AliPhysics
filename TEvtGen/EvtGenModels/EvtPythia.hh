//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See BelEvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: BelEvtGen/EvtJetSet.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//    RS          October 28, 2002        copied from JETSET module
//
//------------------------------------------------------------------------

#ifndef EVTPYTHIA_HH
#define EVTPYTHIA_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"
#include "EvtGenBase/EvtParticle.hh"
#include <string>

#include <iosfwd>

typedef EvtDecayBase* EvtDecayBasePtr;

class EvtPythia:public  EvtDecayIncoherent  {

public:

  EvtPythia();
  virtual ~EvtPythia();
  
  std::string getName();
  EvtDecayBase* clone();
  void decay(EvtParticle *p); 
  
  std::string commandName();
  void command(std::string cmd);
  
  void init();
  
  void initProbMax();
  
  //initialize jetset; sets up decay table and
  //paramters. Static so it can be invoked from
  //from EvtJscont.
  static void pythiaInit(int f);
  static void pythiacont(double *,int *, int *,
			 double *,double *,double *,double *);
  
private:
  
  void store(EvtDecayBase* jsdecay);
  void fixPolarizations(EvtParticle* p);
  static void MakePythiaFile(char* fname);
  static void WritePythiaParticle(std::ofstream &outdec,EvtId ipar,EvtId iparname,int &first);
  static void WritePythiaEntryHeader(std::ofstream &outdec, int lundkc,
				     EvtId evtnum,std::string name,
				     int chg, int cchg, int spin2,double mass,
				     double width, double maxwidth,double ctau,
				     int stable,double rawbrfrsum);
  static bool diquark(int);
  static double NominalMass(int);
  static int njetsetdecays;
  static EvtDecayBasePtr* jetsetdecays;
  static int ntable;
  
  static int ncommand;
  static int lcommand;
  static std::string* commands;
};

#endif




