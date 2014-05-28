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
// Module: EvtGen/EvtParticleDecayList.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTPARTICLEDECAYLIST_HH
#define EVTPARTICLEDECAYLIST_HH

#include "EvtGenBase/EvtParticleDecay.hh"

typedef EvtParticleDecay* EvtParticleDecayPtr;

class EvtParticleDecayList{

public:

  EvtParticleDecayList(){ 
   _decaylist=0;
    _nmode=0;
    _rawbrfrsum=0;
  }

  EvtParticleDecayList(const EvtParticleDecayList &o);

  ~EvtParticleDecayList();

  EvtParticleDecayList& operator=(const EvtParticleDecayList &o);

  int getNMode() const {return _nmode;}

  void setNMode(int nmode);

  EvtDecayBase* getDecayModel(EvtParticle *p);
  EvtDecayBase* getDecayModel(int imode);

  EvtParticleDecay& getDecay(int nchannel) const;

  double getRawBrfrSum() {return _rawbrfrsum;}
  void setRawBrfrSum(double rawbrfrsum) {_rawbrfrsum=rawbrfrsum;}
  
  void makeChargeConj(EvtParticleDecayList* conjDecayList);

  void removeDecay();

  void alocateDecay(int nmode){
    _decaylist= new EvtParticleDecayPtr[nmode];
  }

  void removeMode(EvtDecayBase* decay);

  void addMode(EvtDecayBase* decay,double brfr,double massmin);
  void finalize();

  void printSummary();

  bool isJetSet() const ;

private:

  EvtParticleDecayPtr* _decaylist;

  double _rawbrfrsum;
  int _nmode;

};

#endif

