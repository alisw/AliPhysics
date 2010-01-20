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
// Module: EvtGen/EvtDiracParticle.hh
//
// Description:EvtDiracParticle particles i.e. spin 1/2 particles.
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------

#ifndef EVTDIRACPARTICLE_HH
#define EVTDIRACPARTICLE_HH

#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtParticle.hh"

class EvtId;
class EvtVector4R;

class EvtDiracParticle:public EvtParticle {

public:

  
  EvtDiracParticle();
  virtual ~EvtDiracParticle();
  void init(EvtId part_n,const EvtVector4R& p4);
  void init(EvtId part_n,const EvtVector4R& p4,
	    const EvtDiracSpinor &,const EvtDiracSpinor &,
	    const EvtDiracSpinor &,const EvtDiracSpinor &);
  EvtDiracSpinor spParent(int i) const {return _spinorParent[i];}
  EvtDiracSpinor sp(int i) const {return _spinorRest[i];}  
  EvtSpinDensity rotateToHelicityBasis() const;
  EvtSpinDensity rotateToHelicityBasis(double alpha,
				       double beta,
				       double gamma) const;

private:

  EvtDiracSpinor _spinorRest[2];
  EvtDiracSpinor _spinorParent[2];
  EvtDiracParticle(const EvtDiracParticle& d);
  EvtDiracParticle& operator=(const EvtDiracParticle& d);

};
#endif

