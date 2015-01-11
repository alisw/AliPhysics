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
// Module: EvtGen/EvtSemiLeptonicBaryonAmp.hh
//
// Description:
//
// Modification history:
//
//    Lange Oct 20, 2004 Created
//
//------------------------------------------------------------------------

#ifndef EVTSEMILEPTONICBARYONAMP_HH
#define EVTSEMILEPTONICBARYONAMP_HH

#include "EvtGenBase/EvtSemiLeptonicAmp.hh"

class EvtParticle;
class EvtAmp;
class EvtSemiLeptonicFF;

class EvtVector4C;
class EvtVector4R;
class EvtDiracSpinor;
class EvtRaritaSchwinger;

class EvtSemiLeptonicBaryonAmp:public EvtSemiLeptonicAmp {

 public:

  virtual ~EvtSemiLeptonicBaryonAmp();

  //Daughters are initialized and have been added to the parent.
  //No need to carry around the daughters seperately!
  void CalcAmp( EvtParticle *parent,EvtAmp& amp,
		EvtSemiLeptonicFF *FormFactors );

  void CalcAmp( EvtParticle *parent, EvtAmp& amp,
		EvtSemiLeptonicFF *FormFactors,
		EvtComplex r00, EvtComplex r01,
		EvtComplex r10, EvtComplex r11 );
  
  double CalcMaxProb( EvtId parent, EvtId meson, EvtId lepton,
                      EvtId nudaug, EvtSemiLeptonicFF *FormFactors,
                      EvtComplex r00, EvtComplex r01,
                      EvtComplex r10, EvtComplex r11);


 private:

  EvtVector4C EvtBaryonVACurrent( const EvtDiracSpinor& Bf,
				  const EvtDiracSpinor& Bi, 
				  EvtVector4R parent, 
				  EvtVector4R daught, 
				  const double *ff, int pflag);

  EvtVector4C EvtBaryonVARaritaCurrent( const EvtRaritaSchwinger& Bf_vect,
					const EvtDiracSpinor& Bi, 
					EvtVector4R parent, 
					EvtVector4R daught, 
					const double *ff, int pflag);

};

#endif


