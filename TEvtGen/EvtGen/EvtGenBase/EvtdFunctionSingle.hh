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
// Module: EvtGen/EvtdFunction.hh
//
// Description:Evaluation of one Wigner d-Functions
//
// Modification history:
//
//    RYD      August 10, 2000         Module created
//
//------------------------------------------------------------------------

#ifndef EVTDFUNCTIONSINGLE_HH
#define EVTDFUNCTIONSINGLE_HH

class EvtdFunctionSingle{

public:

  EvtdFunctionSingle();
  ~EvtdFunctionSingle();

  void init(int j,int m1,int m2);

  double d(int j,int m1,int m2,double theta);
  
private:

  int fact(int n);
  
  int _j;
  int _m1;
  int _m2;

  double *_coef;


  int _kmin;
  int _kmax;

};

#endif

