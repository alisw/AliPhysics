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
// Module: EvtGen/EvtRaritaSchwing.hh
//
// Description:Class to handle spin 3/2
//
// Modification history:
//
//    RYD     April 23, 2000         Module created
//
//------------------------------------------------------------------------

#ifndef EVTRARITASCHWINGER_HH
#define EVTRARITASCHWINGER_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"

class EvtRaritaSchwinger ;
EvtRaritaSchwinger rotateEuler(const EvtRaritaSchwinger& rs,
			       double alpha,double beta,double gamma);
EvtRaritaSchwinger boostTo(const EvtRaritaSchwinger& rs,
			   const EvtVector4R p4);
EvtRaritaSchwinger boostTo(const EvtRaritaSchwinger& rs,
			   const EvtVector3R boost);
EvtRaritaSchwinger dirProd(EvtVector4R v,EvtDiracSpinor u);
EvtRaritaSchwinger dirProd(EvtVector4C v,EvtDiracSpinor u);
EvtRaritaSchwinger operator+(const EvtRaritaSchwinger& u1,
			     const EvtRaritaSchwinger& u2); 
EvtRaritaSchwinger operator-(const EvtRaritaSchwinger& u1,
			     const EvtRaritaSchwinger& u2); 
EvtComplex operator*(const EvtRaritaSchwinger& u1,
		     const EvtRaritaSchwinger& u2);

class EvtRaritaSchwinger{

  friend EvtRaritaSchwinger rotateEuler(const EvtRaritaSchwinger& rs,
					double alpha,double beta,double gamma);
  friend EvtRaritaSchwinger boostTo(const EvtRaritaSchwinger& rs,
				    const EvtVector4R p4);
  friend EvtRaritaSchwinger boostTo(const EvtRaritaSchwinger& rs,
				    const EvtVector3R boost);

  friend EvtRaritaSchwinger dirProd(EvtVector4R v,EvtDiracSpinor u);
  friend EvtRaritaSchwinger dirProd(EvtVector4C v,EvtDiracSpinor u);

  friend EvtRaritaSchwinger operator+(const EvtRaritaSchwinger& u1,
				  const EvtRaritaSchwinger& u2); 
  friend EvtRaritaSchwinger operator-(const EvtRaritaSchwinger& u1,
				  const EvtRaritaSchwinger& u2); 

  friend EvtComplex operator*(const EvtRaritaSchwinger& u1,
				  const EvtRaritaSchwinger& u2); 

public:

  inline EvtRaritaSchwinger();
  virtual ~EvtRaritaSchwinger();
  inline EvtRaritaSchwinger(const EvtRaritaSchwinger& rs);
  inline EvtRaritaSchwinger& operator=(const EvtRaritaSchwinger& rs);

  void set(int i,int j,const EvtComplex& sp);

  void applyRotateEuler(double alpha,double beta,double gamma);
  void applyBoostTo(const EvtVector4R p4);
  void applyBoostTo(const EvtVector3R boost);
 
  EvtRaritaSchwinger& operator+=(const EvtRaritaSchwinger& u2);
  EvtRaritaSchwinger& operator-=(const EvtRaritaSchwinger& u2);

  EvtComplex get(int i,int j) const; 
  friend std::ostream& operator<<(std::ostream& s, const EvtRaritaSchwinger& rs); 

  EvtVector4C getVector(int i) const;
  EvtDiracSpinor getSpinor(int i) const;

  void setVector(int i,const EvtVector4C& v);
  void setSpinor(int i,const EvtDiracSpinor& sp);


  
private:
  
  //First index in spinor index, second is Lorentz index.
  EvtComplex _rs[4][4];

};

EvtRaritaSchwinger::EvtRaritaSchwinger(){

  int i,j;
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      _rs[i][j]=0.0;
    }
  }

}

EvtRaritaSchwinger::EvtRaritaSchwinger(const EvtRaritaSchwinger& rs){

  int i,j;
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      _rs[i][j]=rs._rs[i][j];
    }
  }

}

EvtRaritaSchwinger& EvtRaritaSchwinger::operator=(const EvtRaritaSchwinger& rs){

  int i,j;
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      _rs[i][j]=rs._rs[i][j];
    }
  }

  return *this;

}

#endif


