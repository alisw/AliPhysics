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
// Module: EvtDecayParm.cc
//
// Description: Store decay parameters for one decay.
//
// Modification history:
//
//    RYD     April 5, 1997         Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctype.h>
#include "EvtGenBase/EvtDecayParm.hh"
#include <string>
using std::fstream;

void EvtDecayParm::init(fcnPtr pfcn, int ndaug, int *daugs, int narg,
	       double *args, std::string name) {

  int i;

  itsfcn=pfcn;
  itsndaug=ndaug;
  itsnarg=narg;
  
  itsdaugs=new int [itsndaug];
  for(i=0;i<itsndaug;i++){
    itsdaugs[i]=daugs[i];
  }
  itsargs=new double [itsnarg];
  for(i=0;i<itsnarg;i++){
    itsargs[i]=args[i];
  }
  modelname=name;
}

EvtDecayParm::EvtDecayParm() {

  itsfcn=0;
  itsndaug=0;
  itsnarg=0;
  itsdaugs=0;
  itsargs=0;

  modelname="**********";

}

EvtDecayParm::~EvtDecayParm() {

  if (itsdaugs!=0){
     delete [] itsdaugs;
  }

  if (itsargs!=0){
     delete [] itsargs;
  }

}

