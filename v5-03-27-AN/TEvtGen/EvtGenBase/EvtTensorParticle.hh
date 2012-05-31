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
// Module: EvtGen/EvtTensorParticle.hh
//
// Description: Class to describe tensor ( spin 2 ) particles.
//
// Modification history:
//
//    DJL/RYD     Sept. 25, 1996       Module created
//
//------------------------------------------------------------------------

#ifndef EVTTENSORPARTICLE_HH
#define EVTTENSORPARTICLE_HH

#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtParticle.hh"

class EvtTensorParticle: public EvtParticle {
  
public:

  EvtTensorParticle() {}
  virtual ~EvtTensorParticle();

  void init(EvtId part_n,double e,double px,double py,double pz);
  void init(EvtId part_n,const EvtVector4R& p4);
  void init(EvtId part_n,const EvtVector4R& p4,
	    const EvtTensor4C&,const EvtTensor4C&,const EvtTensor4C&,
	    const EvtTensor4C&,const EvtTensor4C&);
  //Returns polarization tensors.
  EvtTensor4C epsTensorParent(int i) const; 
  EvtTensor4C epsTensor(int i) const; 

  EvtSpinDensity rotateToHelicityBasis() const;
  EvtSpinDensity rotateToHelicityBasis(double alpha,
				       double beta,
				       double gamma) const;

  
private:
  
  EvtTensor4C eps[5];//eps1,eps2,eps3,eps4,eps5; 

  EvtTensorParticle(const EvtTensorParticle& tensor);  
  EvtTensorParticle& operator=(const EvtTensorParticle& tensor);  

};

#endif

