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
// Module: EvtGen/EvtScalarParticle.hh
//
// Description:Class to describe all spin 0 particles.
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSCALARPARTICLE_HH
#define EVTSCALARPARTICLE_HH

#include "EvtGenBase/EvtParticle.hh"
class EvtId;


class EvtScalarParticle: public EvtParticle {

public:

  EvtScalarParticle() {}
  virtual ~EvtScalarParticle();

  void init(EvtId part_n,double e,double px,double py,double pz);
  void init(EvtId part_n,const EvtVector4R& p);

  EvtSpinDensity rotateToHelicityBasis() const;
  EvtSpinDensity rotateToHelicityBasis(double alpha,
				       double beta,
				       double gamma) const;
   
private:

  EvtScalarParticle(const EvtScalarParticle& scalar);
  EvtScalarParticle& operator=(const EvtScalarParticle& scalar);

};

#endif

