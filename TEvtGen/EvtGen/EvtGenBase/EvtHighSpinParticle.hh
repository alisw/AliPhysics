//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2000      Caltech, UCSB
//
// Module: EvtGen/EvtHighSpinParticle.hh
//
// Description:Class to describe particle with spin>2.
//
// Modification history:
//
//    RYD     August 8, 2000         Module created
//
//------------------------------------------------------------------------

#ifndef EVTHIGHSPINPARTICLE_HH
#define EVTHIGHSPINPARTICLE_HH

#include "EvtGenBase/EvtParticle.hh"

class EvtId;

class EvtHighSpinParticle: public EvtParticle {

public:

  EvtHighSpinParticle() {}
  virtual ~EvtHighSpinParticle();

  void init(EvtId id,const EvtVector4R& p);

  EvtSpinDensity rotateToHelicityBasis() const;
  EvtSpinDensity rotateToHelicityBasis(double alpha,
				       double beta,
				       double gamma) const;

   
private:

  EvtHighSpinParticle(const EvtHighSpinParticle& highSpin);
  EvtHighSpinParticle& operator=(const EvtHighSpinParticle& highSpin);

};

#endif





