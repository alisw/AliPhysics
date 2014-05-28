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
// Module: EvtdFunctionSingle.cc
//
// Description: Evaluates one Wigner d-Functions.
//
// Modification history:
//
//    fkw           February 2, 2001     changes to satisfy KCC
//    RYD            August 10, 2000         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <assert.h>
#include "EvtGenBase/EvtdFunctionSingle.hh"

EvtdFunctionSingle::EvtdFunctionSingle(){
  _j=0;
  _m1=0;
  _m2=0;
  _coef=0;
  _kmin=0;
  _kmax=0;
}


EvtdFunctionSingle::~EvtdFunctionSingle(){
  if (_coef!=0) delete [] _coef;
}


void EvtdFunctionSingle::init(int j,int m1,int m2){

  assert(abs(m2)>=abs(m1));
  assert(m2>=0);

  _j=j;
  _m1=m1;
  _m2=m2;

  _kmin=_m2-_m1;
  _kmax=_j-_m1;

  assert(_kmin<=_kmax);

  _coef=new double[(_kmax-_kmin)/2+1];

  int k;

  for(k=_kmin;k<=_kmax;k+=2){
    int sign=1;
    if ((k-_m2+_m1)%4!=0) sign=-sign;
    double tmp = fact((_j+_m2)/2)*fact((_j-_m2)/2)
                   *fact((_j+_m1)/2)*fact((_j-_m1)/2);
    _coef[(k-_kmin)/2]=sign*sqrt(tmp)/
      (fact((_j+_m2-k)/2)*fact(k/2)*fact((_j-_m1-k)/2)*fact((k-_m2+_m1)/2));

  }

}


double EvtdFunctionSingle::d(int j,int m1,int m2, double theta){

  assert(j==_j); _unused( j );
  assert(m1==_m1);
  assert(m2==_m2);

  double c2=cos(0.5*theta);
  double s2=sin(0.5*theta);

  double d=0.0;
  
  int k;
  for(k=_kmin;k<=_kmax;k+=2){
    d+=_coef[(k-_kmin)/2]*pow(c2,(2*_j-2*k+m2-m1)/2)*pow(s2,(2*k-m2+m1)/2);
  }

  return d;
  
}


int EvtdFunctionSingle::fact(int n){

  assert(n>=0);

  int f=1;

  int k;
  for(k=2;k<=n;k++) f*=k;

  return f;

}



