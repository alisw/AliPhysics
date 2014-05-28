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
// Module: EvtTensor4C.cc
//
// Description: Implementation of tensor particles.
//
// Modification history:
//
//    DJL/RYD   September 25,1996           Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <math.h>
#include <assert.h>
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
using std::endl;
using std::ostream;



EvtTensor4C::EvtTensor4C( const EvtTensor4C& t1 ) {

  int i,j;
  
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      t[i][j] = t1.t[i][j];
    }
  }

}

EvtTensor4C::~EvtTensor4C() { }

const EvtTensor4C& EvtTensor4C::g(){

  static EvtTensor4C g_metric(1.0,-1.0,-1.0,-1.0);

  return g_metric;

}

EvtTensor4C& EvtTensor4C::operator=(const EvtTensor4C& t1) {
  int i,j;
  
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      t[i][j] = t1.t[i][j];
    }
  }
  return *this;
}

EvtTensor4C EvtTensor4C::conj() const {
  EvtTensor4C temp;
  
  int i,j;
  
  for(i=0;i<4;i++) {
    for(j=0;j<4;j++) {
      temp.set(j,i,::conj(t[i][j]));
    }
  }
  return temp;
}


EvtTensor4C rotateEuler(const EvtTensor4C& rs,
			double alpha,double beta,double gamma){

  EvtTensor4C tmp(rs);
  tmp.applyRotateEuler(alpha,beta,gamma);
  return tmp;

}

EvtTensor4C boostTo(const EvtTensor4C& rs,
		    const EvtVector4R p4){

  EvtTensor4C tmp(rs);
  tmp.applyBoostTo(p4);
  return tmp;

}

EvtTensor4C boostTo(const EvtTensor4C& rs,
		    const EvtVector3R boost){

  EvtTensor4C tmp(rs);
  tmp.applyBoostTo(boost);
  return tmp;

}

void EvtTensor4C::applyBoostTo(const EvtVector4R& p4){

  double e=p4.get(0);

  EvtVector3R boost(p4.get(1)/e,p4.get(2)/e,p4.get(3)/e);

  applyBoostTo(boost);

  return;

}


void EvtTensor4C::applyBoostTo(const EvtVector3R& boost){

  double bx,by,bz,gamma,b2;
  double lambda[4][4];
  EvtComplex tt[4][4];

  bx=boost.get(0);
  by=boost.get(1);
  bz=boost.get(2);

  double bxx=bx*bx;
  double byy=by*by;
  double bzz=bz*bz;

  b2=bxx+byy+bzz;


  if (b2==0.0){
    return;
  }

  assert(b2<1.0);

  gamma=1.0/sqrt(1-b2);


  int i,j,k;
  
  
  if (b2==0.0){
    return ;
  }
  
  lambda[0][0]=gamma;
  lambda[0][1]=gamma*bx;
  lambda[1][0]=gamma*bx;
  lambda[0][2]=gamma*by;
  lambda[2][0]=gamma*by;
  lambda[0][3]=gamma*bz;
  lambda[3][0]=gamma*bz;

  lambda[1][1]=1.0+(gamma-1.0)*bx*bx/b2;
  lambda[2][2]=1.0+(gamma-1.0)*by*by/b2;
  lambda[3][3]=1.0+(gamma-1.0)*bz*bz/b2;
  
  lambda[1][2]=(gamma-1.0)*bx*by/b2;
  lambda[2][1]=(gamma-1.0)*bx*by/b2;
  
  lambda[1][3]=(gamma-1.0)*bx*bz/b2;
  lambda[3][1]=(gamma-1.0)*bx*bz/b2;
  
  lambda[3][2]=(gamma-1.0)*bz*by/b2;
  lambda[2][3]=(gamma-1.0)*bz*by/b2;
  
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      tt[i][j] = EvtComplex(0.0);
      for(k=0;k<4;k++){
        tt[i][j]=tt[i][j]+lambda[j][k]*t[i][k]; 
      }
    }
  }
  
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      t[i][j] = EvtComplex(0.0);
      for(k=0;k<4;k++){
        t[i][j]=t[i][j]+lambda[i][k]*tt[k][j]; 
      }
    }
  }
  
}

void EvtTensor4C::zero(){
  int i,j;
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      t[i][j]=EvtComplex(0.0,0.0);
    }
  }
}



ostream& operator<<(ostream& s,const EvtTensor4C& t){

  int i,j;
  s<< endl;
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      s << t.t[i][j];
    }
   s << endl;
  }
  return s;
}

void EvtTensor4C::setdiag(double g00, double g11, double g22, double g33){
  t[0][0]=EvtComplex(g00);
  t[1][1]=EvtComplex(g11);
  t[2][2]=EvtComplex(g22);
  t[3][3]=EvtComplex(g33);
  t[0][1] = EvtComplex(0.0);
  t[0][2] = EvtComplex(0.0);
  t[0][3] = EvtComplex(0.0);
  t[1][0] = EvtComplex(0.0);
  t[1][2] = EvtComplex(0.0);
  t[1][3] = EvtComplex(0.0);
  t[2][0] = EvtComplex(0.0);
  t[2][1] = EvtComplex(0.0);
  t[2][3] = EvtComplex(0.0);
  t[3][0] = EvtComplex(0.0);
  t[3][1] = EvtComplex(0.0);
  t[3][2] = EvtComplex(0.0);
}


EvtTensor4C& EvtTensor4C::operator+=(const EvtTensor4C& t2){
  
  int i,j;
  
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      t[i][j]+=t2.get(i,j);
    }
  }
  return *this;
}

EvtTensor4C& EvtTensor4C::operator-=(const EvtTensor4C& t2){

  int i,j;
  
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      t[i][j]-=t2.get(i,j);
    }
  }
  return *this;
}


EvtTensor4C& EvtTensor4C::operator*=(const EvtComplex& c) {
  int i,j;
  
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      t[i][j]*=c;
    }
  }
  return *this;
}


EvtTensor4C operator*(const EvtTensor4C& t1,const EvtComplex& c){

  return EvtTensor4C(t1)*=c;

}

EvtTensor4C operator*(const EvtComplex& c,const EvtTensor4C& t1){

  return EvtTensor4C(t1)*=c;

}


EvtTensor4C& EvtTensor4C::operator*=(double d) {
  int i,j;
  
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      t[i][j]*=EvtComplex(d,0.0);
    }
  }
  return *this;
}


EvtTensor4C operator*(const EvtTensor4C& t1, double d){

  return EvtTensor4C(t1)*=EvtComplex(d,0.0);

}

EvtTensor4C operator*(double d, const EvtTensor4C& t1){

  return EvtTensor4C(t1)*=EvtComplex(d,0.0);

}

EvtComplex cont(const EvtTensor4C& t1,const EvtTensor4C& t2){

  EvtComplex sum(0.0,0.0);
  int i,j;

  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      if ((i==0&&j!=0) || (j==0&&i!=0)) {
	sum -= t1.t[i][j]*t2.t[i][j];
      } else {
	sum += t1.t[i][j]*t2.t[i][j]; 
      }
    }
  }

  return sum;
}


EvtTensor4C EvtGenFunctions::directProd(const EvtVector4C& c1,
                                        const EvtVector4C& c2){ 
  EvtTensor4C temp;
  int i,j;
  
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      temp.set(i,j,c1.get(i)*c2.get(j));
    }
  }
  return temp;
}


EvtTensor4C EvtGenFunctions::directProd(const EvtVector4C& c1,
                                        const EvtVector4R& c2){ 
  EvtTensor4C temp;
  int i,j;
  
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      temp.set(i,j,c1.get(i)*c2.get(j));
    }
  }
  return temp;
}


EvtTensor4C EvtGenFunctions::directProd(const EvtVector4R& c1,
                                        const EvtVector4R& c2){ 

  EvtTensor4C temp;
  int i,j;
  
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      temp.t[i][j]=EvtComplex(c1.get(i)*c2.get(j),0.0);
    }
  }
  return temp;
}

EvtTensor4C& EvtTensor4C::addDirProd(const EvtVector4R& p1,const EvtVector4R& p2){ 

  int i,j;
  
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      t[i][j]+=p1.get(i)*p2.get(j);
    }
  }
  return *this;
}


EvtTensor4C dual(const EvtTensor4C& t2){ 
  
  EvtTensor4C temp;
  
  temp.set(0,0,EvtComplex(0.0,0.0));
  temp.set(1,1,EvtComplex(0.0,0.0));
  temp.set(2,2,EvtComplex(0.0,0.0));
  temp.set(3,3,EvtComplex(0.0,0.0));
  
  temp.set(0,1,t2.get(3,2)-t2.get(2,3));
  temp.set(0,2,-t2.get(3,1)+t2.get(1,3));
  temp.set(0,3,t2.get(2,1)-t2.get(1,2));
  
  temp.set(1,2,-t2.get(3,0)+t2.get(0,3));
  temp.set(1,3,t2.get(2,0)-t2.get(0,2));
  
  temp.set(2,3,-t2.get(1,0)+t2.get(0,1));
  
  temp.set(1,0,-temp.get(0,1));
  temp.set(2,0,-temp.get(0,2));
  temp.set(3,0,-temp.get(0,3));
  
  temp.set(2,1,-temp.get(1,2));
  temp.set(3,1,-temp.get(1,3));

  temp.set(3,2,-temp.get(2,3));
  
  return temp;
  
}


EvtTensor4C conj(const EvtTensor4C& t2) { 
  EvtTensor4C temp;
  
  int i,j;

  for(i=0;i<4;i++){ 
    for(j=0;j<4;j++){ 
      temp.set(i,j,::conj((t2.get(i,j))));
    }
  }
  
  return temp;
}


EvtTensor4C cont22(const EvtTensor4C& t1,const EvtTensor4C& t2){ 
  EvtTensor4C temp;

  int i,j;
  EvtComplex c;
  
  for(i=0;i<4;i++){ 
    for(j=0;j<4;j++){ 
      c=t1.get(i,0)*t2.get(j,0)-t1.get(i,1)*t2.get(j,1)
	-t1.get(i,2)*t2.get(j,2)-t1.get(i,3)*t2.get(j,3);
      temp.set(i,j,c);
    }
  }
  
  return temp;
}

EvtTensor4C cont11(const EvtTensor4C& t1,const EvtTensor4C& t2){ 
  EvtTensor4C temp;
  
  int i,j;
  EvtComplex c;
  
  for(i=0;i<4;i++){ 
    for(j=0;j<4;j++){ 
        c=t1.get(0,i)*t2.get(0,j)-t1.get(1,i)*t2.get(1,j)
	  -t1.get(2,i)*t2.get(2,j)-t1.get(3,i)*t2.get(3,j);
	temp.set(i,j,c);
    }
  }
  
  return temp;
}


EvtVector4C EvtTensor4C::cont1(const EvtVector4C& v4) const {
  EvtVector4C temp;
  
  int i;
  
  for(i=0;i<4;i++){
    temp.set(i,t[0][i]*v4.get(0)-t[1][i]*v4.get(1)
	     -t[2][i]*v4.get(2)-t[3][i]*v4.get(3));
  }
  
  return temp;
} 

EvtVector4C EvtTensor4C::cont2(const EvtVector4C& v4) const {
  EvtVector4C temp;

  int i;
  
  for(i=0;i<4;i++){
    temp.set(i,t[i][0]*v4.get(0)-t[i][1]*v4.get(1)
	     -t[i][2]*v4.get(2)-t[i][3]*v4.get(3));
  }
  
  return temp;
} 


EvtVector4C EvtTensor4C::cont1(const EvtVector4R& v4) const {
  EvtVector4C temp;
  
  int i;
  
  for(i=0;i<4;i++){
    temp.set(i,t[0][i]*v4.get(0)-t[1][i]*v4.get(1)
	     -t[2][i]*v4.get(2)-t[3][i]*v4.get(3));
  }

  return temp;
} 


EvtVector4C EvtTensor4C::cont2(const EvtVector4R& v4) const {
  EvtVector4C temp;
  
  int i;
  
  for(i=0;i<4;i++){
    temp.set(i,t[i][0]*v4.get(0)-t[i][1]*v4.get(1)
	     -t[i][2]*v4.get(2)-t[i][3]*v4.get(3));
  }
  
  return temp;
} 



void EvtTensor4C::applyRotateEuler(double phi,double theta,double ksi){

  EvtComplex tt[4][4];
  double sp,st,sk,cp,ct,ck;
  double lambda[4][4];

  sp=sin(phi);
  st=sin(theta);
  sk=sin(ksi);
  cp=cos(phi);
  ct=cos(theta);
  ck=cos(ksi);


  lambda[0][0]=1.0;
  lambda[0][1]=0.0;
  lambda[1][0]=0.0;
  lambda[0][2]=0.0;
  lambda[2][0]=0.0;
  lambda[0][3]=0.0;
  lambda[3][0]=0.0;

  lambda[1][1]= ck*ct*cp-sk*sp;
  lambda[1][2]=-sk*ct*cp-ck*sp;
  lambda[1][3]=st*cp;

  lambda[2][1]= ck*ct*sp+sk*cp;
  lambda[2][2]=-sk*ct*sp+ck*cp;
  lambda[2][3]=st*sp;

  lambda[3][1]=-ck*st;
  lambda[3][2]=sk*st;
  lambda[3][3]=ct;
  

  int i,j,k;

  
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      tt[i][j] = EvtComplex(0.0);
      for(k=0;k<4;k++){
        tt[i][j]+=lambda[j][k]*t[i][k]; 
      }
    }
  }
  
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      t[i][j] = EvtComplex(0.0);
      for(k=0;k<4;k++){
        t[i][j]+=lambda[i][k]*tt[k][j]; 
      }
    }
  }
  
}



