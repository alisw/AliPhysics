//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
//
// Module: EvtGen/EvtBtoXsllUtil.hh
//
// Description:
// Class to generate inclusive non-resonant B -> Xs l+ l- decays.
//
// Modification history:
//
//    Stephane Willocq    Jan 19, 2001   Module created
//    Stephane Willocq    Nov  6, 2003   Update Wilson Coeffs
//
//------------------------------------------------------------------------

#ifndef EVTBTOXSLLUTIL_HH
#define EVTBTOXSLLUTIL_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDecayIncoherent.hh"
#include "EvtGenBase/EvtParticle.hh"

class EvtComplex;

class EvtBtoXsllUtil{

public:

  EvtComplex GetC7Eff0(double sh, bool nnlo=true);
  EvtComplex GetC7Eff1(double sh, double mb, bool nnlo=true);
  EvtComplex GetC9Eff0(double sh, double mb, bool nnlo=true, bool btod=false);
  EvtComplex GetC9Eff1(double sh, double mb, bool nnlo=true, bool btod=false);
  EvtComplex GetC10Eff(double sh, bool nnlo=true);

  double dGdsProb(double mb, double ms, double ml,
                  double s);

  double dGdsdupProb(double mb, double ms, double ml,
                     double s,  double u);
  
  double FermiMomentum( double pf );
  
  double FermiMomentumProb( double pb, double pf );

};

#endif

