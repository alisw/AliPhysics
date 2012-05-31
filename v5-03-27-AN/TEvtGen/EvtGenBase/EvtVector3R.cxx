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
// Module: EvtVector3R.cc
//
// Description: Real implementation of 3-vectors
//
// Modification history:
//
//    RYD       September 5, 1997       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <math.h>
#include "EvtGenBase/EvtVector3R.hh"
using std::ostream;



EvtVector3R::~EvtVector3R(){}

EvtVector3R::EvtVector3R(){
  
  v[0]=v[1]=v[2]=0.0;
}


EvtVector3R::EvtVector3R(double x,double y, double z){
  
  v[0]=x; v[1]=y; v[2]=z;
}

EvtVector3R rotateEuler(const EvtVector3R& v,
			double alpha,double beta,double gamma){

  EvtVector3R tmp(v);
  tmp.applyRotateEuler(alpha,beta,gamma);
  return tmp;

}


void EvtVector3R::applyRotateEuler(double phi,double theta,double ksi){

  double temp[3];
  double sp,st,sk,cp,ct,ck;

  sp=sin(phi);
  st=sin(theta);
  sk=sin(ksi);
  cp=cos(phi);
  ct=cos(theta);
  ck=cos(ksi);

  temp[0]=( ck*ct*cp-sk*sp)*v[0]+( -sk*ct*cp-ck*sp)*v[1]+st*cp*v[2];
  temp[1]=( ck*ct*sp+sk*cp)*v[0]+(-sk*ct*sp+ck*cp)*v[1]+st*sp*v[2];
  temp[2]=-ck*st*v[0]+sk*st*v[1]+ct*v[2];


  v[0]=temp[0];
  v[1]=temp[1];
  v[2]=temp[2];
}

ostream& operator<<(ostream& s,const EvtVector3R& v){
 
  s<<"("<<v.v[0]<<","<<v.v[1]<<","<<v.v[2]<<")";

  return s;

}


EvtVector3R cross( const EvtVector3R& p1,const EvtVector3R& p2 ){

  //Calcs the cross product.  Added by djl on July 27, 1995.
  //Modified for real vectros by ryd Aug 28-96

  return EvtVector3R(p1.v[1]*p2.v[2] - p1.v[2]*p2.v[1],
		     p1.v[2]*p2.v[0] - p1.v[0]*p2.v[2],
		     p1.v[0]*p2.v[1] - p1.v[1]*p2.v[0]);

}

double EvtVector3R::d3mag() const

// returns the 3 momentum mag.
{
  double temp;

  temp = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
  temp = sqrt( temp );

  return temp;
} // r3mag

double EvtVector3R::dot ( const EvtVector3R& p2 ){

  double temp;

  temp = v[0]*p2.v[0];
  temp += v[0]*p2.v[0];
  temp += v[0]*p2.v[0];
 
  return temp;
} //dot





