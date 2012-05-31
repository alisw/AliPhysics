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
// Module: EvtGen/EvtSecondary.hh
//
// Description:Class store decays of secondary particles
//
// Modification history:
//
//    RYD     March. 12, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSECONDARY_HH
#define EVTSECONDARY_HH


const int EVTSECONDARYLENGTH =100;

class EvtParticle;
#include <iosfwd>

class EvtSecondary {

public:

  EvtSecondary(){}
  ~EvtSecondary(){}

  void init();

  int getStdHepIndex(int i) {return _stdhepindex[i];}
  int getD1(int i) {return _id1[i];}
  int getD2(int i) {return _id2[i];}
  int getD3(int i) {return _id3[i];}
  
  int getNPart();
  void createSecondary(int stdhepindex,EvtParticle* prnt);

  friend std::ostream& operator<<(std::ostream& s, const EvtSecondary& secondary);  

private:

  int _npart;
  int _stdhepindex[EVTSECONDARYLENGTH];
  int _id1[EVTSECONDARYLENGTH];
  int _id2[EVTSECONDARYLENGTH];
  int _id3[EVTSECONDARYLENGTH];
  
}; 

#endif

