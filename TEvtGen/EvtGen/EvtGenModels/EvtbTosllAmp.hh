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
// Module: EvtGen/EvtbTosllAmp.hh
//
// Description:
//
// Modification history:
//
//    RYD     January 30 2000         Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTOSLLAMP_HH
#define EVTBTOSLLAMP_HH

class EvtAmp;
class EvtId;
class EvtbTosllFF;
class EvtParticle;
class EvtComplex;

class EvtbTosllAmp{

 public:
  virtual ~EvtbTosllAmp() { } ;

  //Daughters are initialized and have been added to the parent.
  //No need to carry around the daughters seperately!

  virtual void CalcAmp( EvtParticle *parent, EvtAmp& amp,
			EvtbTosllFF *formFactors )=0;

  double CalcMaxProb( EvtId parent, EvtId meson, EvtId lepton,
		      EvtId nudaug, EvtbTosllFF *formFactors,
		      double& poleSize);

  EvtComplex GetC7Eff(double q2, bool nnlo=true);
  EvtComplex GetC9Eff(double q2, bool nnlo=true, bool btod=false);
  EvtComplex GetC10Eff(double q2, bool nnlo=true);

  double dGdsProb(double mb, double ms, double ml,
                                  double s);

  double dGdsdupProb(double mb, double ms, double ml,
                                     double s,  double u);

};

#endif



