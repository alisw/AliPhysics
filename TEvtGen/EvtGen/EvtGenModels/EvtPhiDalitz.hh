
#ifndef EVTPHIDALITZ_HH
#define EVTPHIDALITZ_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtPhiDalitz:public  EvtDecayAmp  {

public:

  EvtPhiDalitz() {}
  virtual ~EvtPhiDalitz();

  std::string getName();
  EvtDecayBase* clone();

  void init();

  void decay(EvtParticle *p); 

private:
  double _mRho;
  double _gRho;
  double _aD;
  double _phiD;
  double _aOmega;
  double _phiOmega;
  int _locPip;
  int _locPim;
  int _locPi0;

};

#endif
