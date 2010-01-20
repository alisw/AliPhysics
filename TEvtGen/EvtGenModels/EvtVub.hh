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
// Module: EvtGen/EvtVub.hh
//
// Description:
// Class to generate inclusive B to X_u l nu decays according to various
// decay models. Implemtented are ACCM, parton-model and a QCD model.
//
// Modification history:
//
//   Sven Menke     January 17, 2001         Module created
//
//------------------------------------------------------------------------

#ifndef EVTVUB_HH
#define EVTVUB_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"

#include <vector>

class EvtParticle;
class EvtVubdGamma;

class EvtVub:public  EvtDecayIncoherent  {

public:
  
  EvtVub() : _masses(0) , _weights(0) , _dGamma(0) {}
  virtual ~EvtVub();

  std::string getName();

  EvtDecayBase* clone();

  void initProbMax();

  void init();

  void decay(EvtParticle *p); 

private:

  double _mb;     // the b-quark pole mass in GeV (try 4.65 to 4.9)
  double _a;      // Parameter for the Fermi Motion (1.29 is good)
  double _alphas; // Strong Coupling at m_b (around 0.24)
  double _dGMax;  // max dGamma*p2 value;
  int    _nbins;
  int    _storeQplus;
  double * _masses;
  double * _weights;

  EvtVubdGamma *_dGamma; // calculates the decay rate
  double findPFermi();
  std::vector<double> _pf;
};

#endif

