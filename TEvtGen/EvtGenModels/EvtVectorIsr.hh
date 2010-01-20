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
// Module: EvtGen/EvtVectorIsr2.hh
//
// Description: 
//   This is a special decay model to generate e+e- -> phi gamma + soft gammas
//   using soft collinear ISR calculation from AfkQed
//   This is implemented as a decay of the VPHO.
//
// Modification history:
//
//    Joe Izen        Oct, 2005             Soft Colinear Photons (secondary ISR) ported from AfkQed
//    Joe Izen        Dec  16, 2002         Fix cos_theta distribution - prevents boom at cos_theta=+/-1 
//    RYD/Adriano     June 16, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTVECTORISR_HH
#define EVTVECTORISR_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"

class EvtParticle;


class EvtVectorIsr:public  EvtDecayIncoherent  {

public:

  EvtVectorIsr() {}
  virtual ~EvtVectorIsr();


  std::string getName();

  EvtDecayBase* clone();

  void decay(EvtParticle *p); 

  void init();

  void initProbMax();

  double ckhrad1(double xx, double a, double b);
  
  void ckhrad(const double& e_beam,const double& q2_min,double& e01,double& e02,double& f);


private:

  double csfrmn,csbkmn;
  double fmax;
  bool firstorder;
};

#endif











