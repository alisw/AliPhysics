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
// Module: EvtGen/EvtDiracSpinor.hh
//
// Description:Class to manipulate dirac spinors
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------

#ifndef EVTDIRACSPINOR_HH
#define EVTDIRACSPINOR_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVector3R.hh"

class EvtTensor4C;
class EvtVector4C;
class EvtDiracSpinor ;


class EvtDiracSpinor {

  friend EvtDiracSpinor rotateEuler(const EvtDiracSpinor& sp,
				 double alpha,double beta,double gamma);
  friend EvtDiracSpinor boostTo(const EvtDiracSpinor& sp,
			     const EvtVector4R p4);
  friend EvtDiracSpinor boostTo(const EvtDiracSpinor& sp,
			     const EvtVector3R boost);  
  friend EvtVector4C EvtLeptonVACurrent
	(const EvtDiracSpinor& d,const EvtDiracSpinor& dp);
  friend EvtVector4C EvtLeptonVCurrent 
	(const EvtDiracSpinor& d,const EvtDiracSpinor& dp);
  friend EvtVector4C EvtLeptonACurrent 
	(const EvtDiracSpinor& d,const EvtDiracSpinor& dp);
  friend EvtComplex  EvtLeptonSCurrent 
	(const EvtDiracSpinor& d,const EvtDiracSpinor& dp);
  friend EvtComplex  EvtLeptonPCurrent 
	(const EvtDiracSpinor& d,const EvtDiracSpinor& dp);
  friend EvtTensor4C  EvtLeptonTCurrent 
	(const EvtDiracSpinor& d,const EvtDiracSpinor& dp);
  friend EvtDiracSpinor operator+(const EvtDiracSpinor& u1,
				  const EvtDiracSpinor& u2); 
  friend EvtDiracSpinor operator-(const EvtDiracSpinor& u1,
				  const EvtDiracSpinor& u2); 
  friend EvtDiracSpinor operator*(const EvtComplex& c,
				  const EvtDiracSpinor& d);

  friend EvtComplex operator*(const EvtDiracSpinor& d ,
                              const EvtDiracSpinor& dp ) ;
 
  friend std::ostream& operator<<(std::ostream& s, const EvtDiracSpinor& c);  

public:

  inline EvtDiracSpinor();
  EvtDiracSpinor(const EvtComplex& sp0,const EvtComplex& sp1,
		 const EvtComplex& sp2,const EvtComplex& sp3);
  virtual ~EvtDiracSpinor();
  inline EvtDiracSpinor(const EvtDiracSpinor& dspinor);
  inline EvtDiracSpinor& operator=(const EvtDiracSpinor& dspinor);

  inline EvtDiracSpinor& operator+=(const EvtDiracSpinor& u2);
  inline EvtDiracSpinor& operator-=(const EvtDiracSpinor& u2);

  void set(const EvtComplex& sp0,const EvtComplex& sp1,
	   const EvtComplex& sp2,const EvtComplex& sp3);
  void set_spinor(int i,const EvtComplex& sp);
  const EvtComplex& get_spinor(int i) const; 
  EvtDiracSpinor conj() const;
  void applyRotateEuler(double alpha,double beta,double gamma);
  void applyBoostTo(const EvtVector4R& p4);
  void applyBoostTo(const EvtVector3R& boost);
  EvtDiracSpinor adjoint() const;
  
private:

  EvtComplex spinor[4];

};

EvtDiracSpinor::EvtDiracSpinor(){

  spinor[0]=EvtComplex(); spinor[1]=EvtComplex();
  spinor[2]=EvtComplex(); spinor[3]=EvtComplex();

}

EvtDiracSpinor::EvtDiracSpinor(const EvtDiracSpinor& dspinor){

  spinor[0]=dspinor.spinor[0];
  spinor[1]=dspinor.spinor[1];
  spinor[2]=dspinor.spinor[2];
  spinor[3]=dspinor.spinor[3];

}

EvtDiracSpinor& EvtDiracSpinor::operator=(const EvtDiracSpinor& dspinor){

  spinor[0]=dspinor.spinor[0];
  spinor[1]=dspinor.spinor[1];
  spinor[2]=dspinor.spinor[2];
  spinor[3]=dspinor.spinor[3];

  return *this;

}

inline EvtDiracSpinor& EvtDiracSpinor::operator+=(const EvtDiracSpinor& u2){

  spinor[0]+=u2.spinor[0];
  spinor[1]+=u2.spinor[1];
  spinor[2]+=u2.spinor[2];
  spinor[3]+=u2.spinor[3];
  
  return *this; 
}

inline EvtDiracSpinor operator+(const EvtDiracSpinor& u1,
				const EvtDiracSpinor& u2){
  
  return EvtDiracSpinor(u1)+=u2;

}

inline EvtDiracSpinor& EvtDiracSpinor::operator-=(const EvtDiracSpinor& u2){

  spinor[0]-=u2.spinor[0];
  spinor[1]-=u2.spinor[1];
  spinor[2]-=u2.spinor[2];
  spinor[3]-=u2.spinor[3];
  
  return *this; 
}

inline EvtDiracSpinor operator-(const EvtDiracSpinor& u1,
				const EvtDiracSpinor& u2){
  
  return EvtDiracSpinor(u1)-=u2;

}

#endif


