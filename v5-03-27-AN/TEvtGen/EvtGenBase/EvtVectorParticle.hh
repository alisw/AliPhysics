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
// Module: EvtGen/EvtVectorParticle.hh
//
// Description: Class to describe vector particles.
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------

#ifndef EVTVECTORPARTICLE_HH
#define EVTVECTORPARTICLE_HH

#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtParticle.hh"

class EvtId;

class EvtVectorParticle: public EvtParticle {

public:

  EvtVectorParticle() {}
  virtual ~EvtVectorParticle();

  void init(EvtId part_n,double e,double px,double py,double pz);
  void init(EvtId part_n,const EvtVector4R& p);
  void init(EvtId part_n,const EvtVector4R& p,
	    const EvtVector4C&,const EvtVector4C&,const EvtVector4C&);
  EvtVector4C epsParent(int i) const   {return boostTo(_eps[i],this->getP4());}
  EvtVector4C eps(int i) const {return _eps[i];} 
  EvtSpinDensity rotateToHelicityBasis() const;
  EvtSpinDensity rotateToHelicityBasis(double alpha,
				       double beta,
				       double gamma) const;
 
private:
  
  EvtVector4C _eps[3];

  EvtVectorParticle(const EvtVectorParticle& vector);
  EvtVectorParticle& operator=(const EvtVectorParticle& vector);
  
};

#endif

