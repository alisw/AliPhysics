//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Module: EvtGenModels/EvtBcToNPi.hh
//
// Description: General decay model for Bc -> V + npi and Bc -> P + npi
//
// Modification history:
//
//    A.Berezhnoy, A.Likhoded, A.Luchinsky  April 2011   Module created
//
//------------------------------------------------------------------------

#ifndef EvtBcToNPi_HH
#define EvtBcToNPi_HH

#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <string>

class EvtBcToNPi: public EvtDecayAmp {

public:
  
  EvtBcToNPi(bool printAuthorInfo = false);
  virtual ~EvtBcToNPi();

  std::string getName();

  EvtDecayBase* clone();

  void initProbMax();

  void init();

  void decay(EvtParticle *p); 

protected:
  
  int nCall;
  double maxAmp2;

  // Bc form factors
  double _maxProb;
  double FA0_N, FA0_c1, FA0_c2;
  double FAm_N, FAm_c1, FAm_c2;
  double FAp_N, FAp_c1, FAp_c2;
  double FV_N, FV_c1, FV_c2;
  
  double Fp_N, Fp_c1, Fp_c2;
  double Fm_N, Fm_c1, Fm_c2;
  
  // W -> pi... form factors
  double _beta;
  double  _mRho;
  double  _gammaRho;
  double  _mRhopr;
  double  _gammaRhopr;
  double _mA1;
  double _gammaA1;

  double _ee(double M, double m1, double m2);
  double _pp(double M, double m1, double m2);
  EvtComplex Fpi( EvtVector4R q1, EvtVector4R q2);
  double pi3G(double m2,int dupD);

private:

  void printAuthorInfo();

};

#endif

