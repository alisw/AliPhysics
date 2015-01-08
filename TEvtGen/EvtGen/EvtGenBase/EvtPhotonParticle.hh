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
// Module: EvtGen/EvtPhotonParticle.hh
//
// Description:Class to describe photons
//
// Modification history:
//
//    DJL/RYD     Sept. 25, 1996         Module created
//
//------------------------------------------------------------------------

#ifndef EVTPHOTONPARTICLE_HH
#define EVTPHOTONPARTICLE_HH

#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtParticle.hh"
class EvtId;

//Class to handle massless spin 1 particles.

class EvtPhotonParticle: public EvtParticle {

public:

  EvtPhotonParticle();
  virtual ~EvtPhotonParticle();

  void init(EvtId part_n,double e,double px,double py,double pz);
  void init(EvtId part_n,const EvtVector4R& p4);

  //Return polarization vectors
  EvtVector4C epsParentPhoton(int i); 
  EvtVector4C epsPhoton(int i); 


  EvtSpinDensity rotateToHelicityBasis() const;
  EvtSpinDensity rotateToHelicityBasis(double alpha,
				       double beta,
				       double gamma) const;

private:

  EvtVector4C eps1,eps2;
  int _evalBasis;

  EvtPhotonParticle(const EvtPhotonParticle& photon);
  EvtPhotonParticle& operator=(const EvtPhotonParticle& photon);

};

#endif

