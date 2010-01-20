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
// Module: EvtGen/EvtdFunction.hh
//
// Description:Evaluation of Wigner d-Functions
//
// Modification history:
//
//    RYD      March 14, 1999         Module created
//
//------------------------------------------------------------------------

#ifndef EVTDFUNCTION_HH
#define EVTDFUNCTION_HH

class EvtdFunction{

public:

 static double d(int j,int m1,int m2,double theta);

};

#endif

