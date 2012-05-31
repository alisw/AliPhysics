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
// Module: EvtGen/EvtDecayIncoherent.hh
//
// Description: Base class for models that calculate
//              decay kinematics and do not do any accept/reject.
//              Useful e.g. for interface to other generators
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EvtDecayIncoherent_HH
#define EvtDecayIncoherent_HH

#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtParticle.hh"



class EvtDecayIncoherent : public EvtDecayBase{

public:

  void makeDecay(EvtParticle* p, bool recursive=true);

  virtual ~EvtDecayIncoherent() {}

  void setDaughterSpinDensity(int daughter)
  { spinDensitySet[daughter]=1; return;}

  int isDaughterSpinDensitySet(int daughter) 
  {return spinDensitySet[daughter];}

private:

  int spinDensitySet[MAX_DAUG];

};




#endif



