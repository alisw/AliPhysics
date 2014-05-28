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
// Module: EvtKine.cc
//
// Description: Evaluates the Wigner d-Functions.
//
// Modification history:
//
//    RYD            March 14, 1999         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <assert.h>
#include "EvtGenBase/EvtdFunction.hh"
#include "EvtGenBase/EvtdFunctionSingle.hh"


double EvtdFunction::d(int j,int m1,int m2, double theta){


  int m1p=m1;
  int m2p=m2;


  int signp=1;
  //make |m2p|>|m1p|
  if (abs(m2p)<abs(m1p)) {
    int tmp=m1p;
    m1p=m2p;
    m2p=tmp;
    if ((m1p-m2p)%4!=0) signp=-signp;
  } 

  //make m2p non-negative
  if (m2p<0) {
    m1p=-m1p;
    m2p=-m2p;
    if ((m1p-m2p)%4!=0) signp=-signp;
  }


  EvtdFunctionSingle df;

  df.init(j,m1p,m2p);

  double d=df.d(j,m1p,m2p,theta)*signp;

  return d;
  
}



