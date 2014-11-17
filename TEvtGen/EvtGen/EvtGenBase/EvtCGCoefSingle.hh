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
// Module: EvtGen/EvtCGCoefSingle.hh
//
// Description:Evaluation of Clebsch-Gordon coef for a fixed j1 and j2.
//
// Modification history:
//
//    RYD      August 12, 2000         Module created
//
//------------------------------------------------------------------------

#ifndef EVTCGCOEFSINGLE_HH
#define EVTCGCOEFSINGLE_HH

#include <vector>


class EvtCGCoefSingle{

public:

  EvtCGCoefSingle(int j1, int j2) {init(j1,j2);}
  ~EvtCGCoefSingle();


  double coef(int J,int M,int j1,int j2,int m1,int m2);
  
private:

  void init(int j1,int j2);
  double& cg(int J,int M, int m1, int m2);

  int _j1;
  int _j2;

  int _Jmax;
  int _Jmin;

  std::vector<std::vector<std::vector<double> > > _table;

};

#endif


