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
// Module: EvtDDalitz.cc
//
// Description: Routine to handle three-body decays of D0/D0_bar or D+/D-
//
// Modification history:
//
//    NK     September 3, 1997         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtResonance.hh"
#include "EvtGenBase/EvtResonance2.hh"
#include "EvtGenModels/EvtDDalitz.hh"
#include <string>
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtFlatte.hh"
#include "EvtGenBase/EvtDecayTable.hh"
using std::endl;

EvtDDalitz::~EvtDDalitz() {}

std::string EvtDDalitz::getName(){
  
  return "D_DALITZ";     

}


EvtDecayBase* EvtDDalitz::clone(){

  return new EvtDDalitz;

}

void EvtDDalitz::init(){

  // check that there are 0 arguments

  static EvtId DM=EvtPDL::getId("D-");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId DSP=EvtPDL::getId("D_s+");
  static EvtId DSM=EvtPDL::getId("D_s-");
  static EvtId KM=EvtPDL::getId("K-");
  static EvtId KP=EvtPDL::getId("K+");
  static EvtId K0=EvtPDL::getId("K0");
  static EvtId KB=EvtPDL::getId("anti-K0");
  static EvtId KL=EvtPDL::getId("K_L0");
  static EvtId KS=EvtPDL::getId("K_S0");
  static EvtId PIM=EvtPDL::getId("pi-");
  static EvtId PIP=EvtPDL::getId("pi+");
  static EvtId PI0=EvtPDL::getId("pi0");

  static double MPI = EvtPDL::getMeanMass(PI0);
  static double MKP = EvtPDL::getMeanMass(KP);


  checkNArg(0);
  checkNDaug(3);

  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);
  checkSpinDaughter(2,EvtSpinType::SCALAR);

  EvtId parnum=getParentId();
  EvtId d1=getDaug(0);
  EvtId d2=getDaug(1);
  EvtId d3=getDaug(2);
  _flag=0;
  if ( parnum == D0 ) {
    //look for either a K- pi+ pi0 or K0bar pi+ pi-
    if ( d1==KM && d2==PIP && d3==PI0 ) { _flag=4; _d1=0; _d2=1; _d3=2;}
    if ( d1==KM && d3==PIP && d2==PI0 ) { _flag=4; _d1=0; _d2=2; _d3=1;}
    if ( d2==KM && d1==PIP && d3==PI0 ) { _flag=4; _d1=1; _d2=0; _d3=2;}
    if ( d2==KM && d3==PIP && d1==PI0 ) { _flag=4; _d1=1; _d2=2; _d3=0;}
    if ( d3==KM && d1==PIP && d2==PI0 ) { _flag=4; _d1=2; _d2=0; _d3=1;}
    if ( d3==KM && d2==PIP && d1==PI0 ) { _flag=4; _d1=2; _d2=1; _d3=0;}

    if ( d1==KB && d2==PIP && d3==PIM )  { _flag=3; _d1=0; _d2=2; _d3=1;}
    if ( d1==KB && d3==PIP && d2==PIM ) { _flag=3; _d1=0; _d2=1; _d3=2;}
    if ( d2==KB && d1==PIP && d3==PIM ) { _flag=3; _d1=1; _d2=2; _d3=0;}
    if ( d2==KB && d3==PIP && d1==PIM ) { _flag=3; _d1=1; _d2=0; _d3=2;}
    if ( d3==KB && d1==PIP && d2==PIM ) { _flag=3; _d1=2; _d2=1; _d3=0;}
    if ( d3==KB && d2==PIP && d1==PIM ) { _flag=3; _d1=2; _d2=0; _d3=1;}

    if ( d1==KL && d2==PIP && d3==PIM )  { _flag=3; _d1=0; _d2=2; _d3=1;}
    if ( d1==KL && d3==PIP && d2==PIM ) { _flag=3; _d1=0; _d2=1; _d3=2;}
    if ( d2==KL && d1==PIP && d3==PIM ) { _flag=3; _d1=1; _d2=2; _d3=0;}
    if ( d2==KL && d3==PIP && d1==PIM ) { _flag=3; _d1=1; _d2=0; _d3=2;}
    if ( d3==KL && d1==PIP && d2==PIM ) { _flag=3; _d1=2; _d2=1; _d3=0;}
    if ( d3==KL && d2==PIP && d1==PIM ) { _flag=3; _d1=2; _d2=0; _d3=1;}

    if ( d1==KS && d2==PIP && d3==PIM )  { _flag=3; _d1=0; _d2=2; _d3=1;}
    if ( d1==KS && d3==PIP && d2==PIM ) { _flag=3; _d1=0; _d2=1; _d3=2;}
    if ( d2==KS && d1==PIP && d3==PIM ) { _flag=3; _d1=1; _d2=2; _d3=0;}
    if ( d2==KS && d3==PIP && d1==PIM ) { _flag=3; _d1=1; _d2=0; _d3=2;}
    if ( d3==KS && d1==PIP && d2==PIM ) { _flag=3; _d1=2; _d2=1; _d3=0;}
    if ( d3==KS && d2==PIP && d1==PIM ) { _flag=3; _d1=2; _d2=0; _d3=1;}

    
    if ( d1==KS && d2==KP && d3==KM ) {_flag=5;_d1=0;_d2=2;_d3=1;}
    if ( d1==KS && d3==KP && d2==KM ) {_flag=5;_d1=0;_d2=1;_d3=2;}
    if ( d2==KS && d1==KP && d3==KM ) {_flag=5;_d1=1;_d2=2;_d3=0;}
    if ( d2==KS && d3==KP && d1==KM ) {_flag=5;_d1=1;_d2=0;_d3=2;}
    if ( d3==KS && d1==KP && d2==KM ) {_flag=5;_d1=2;_d2=1;_d3=0;}
    if ( d3==KS && d2==KP && d1==KM ) {_flag=5;_d1=2;_d2=0;_d3=1;}
    
    if ( d1==KL && d2==KP && d3==KM ) {_flag=5;_d1=0;_d2=2;_d3=1;}
    if ( d1==KL && d3==KP && d2==KM ) {_flag=5;_d1=0;_d2=1;_d3=2;}
    if ( d2==KL && d1==KP && d3==KM ) {_flag=5;_d1=1;_d2=2;_d3=0;}
    if ( d2==KL && d3==KP && d1==KM ) {_flag=5;_d1=1;_d2=0;_d3=2;}
    if ( d3==KL && d1==KP && d2==KM ) {_flag=5;_d1=2;_d2=1;_d3=0;}
    if ( d3==KL && d2==KP && d1==KM ) {_flag=5;_d1=2;_d2=0;_d3=1;}    
    
    if ( d1==K0 && d2==KP && d3==KM ) {_flag=5;_d1=0;_d2=2;_d3=1;}
    if ( d1==K0 && d3==KP && d2==KM ) {_flag=5;_d1=0;_d2=1;_d3=2;}
    if ( d2==K0 && d1==KP && d3==KM ) {_flag=5;_d1=1;_d2=2;_d3=0;}
    if ( d2==K0 && d3==KP && d1==KM ) {_flag=5;_d1=1;_d2=0;_d3=2;}
    if ( d3==K0 && d1==KP && d2==KM ) {_flag=5;_d1=2;_d2=1;_d3=0;}
    if ( d3==K0 && d2==KP && d1==KM ) {_flag=5;_d1=2;_d2=0;_d3=1;}   
  }
  if ( parnum == D0B ) {
    //look for either a K+ pi- pi0 or K0 pi+ pi-
    if ( d1==KP && d2==PIM && d3==PI0 )  { _flag=4; _d1=0; _d2=1; _d3=2;}
    if ( d1==KP && d3==PIM && d2==PI0 ) { _flag=4; _d1=0; _d2=2; _d3=1;}
    if ( d2==KP && d1==PIM && d3==PI0 ) { _flag=4; _d1=1; _d2=0; _d3=2;}
    if ( d2==KP && d3==PIM && d1==PI0 ) { _flag=4; _d1=1; _d2=2; _d3=0;}
    if ( d3==KP && d1==PIM && d2==PI0 ) { _flag=4; _d1=2; _d2=0; _d3=1;}
    if ( d3==KP && d2==PIM && d1==PI0 ) { _flag=4; _d1=2; _d2=1; _d3=0;}

    if ( d1==K0 && d2==PIP && d3==PIM )  { _flag=3; _d1=0; _d2=1; _d3=2;}
    if ( d1==K0 && d3==PIP && d2==PIM ) { _flag=3; _d1=0; _d2=2; _d3=1;}
    if ( d2==K0 && d1==PIP && d3==PIM ) { _flag=3; _d1=1; _d2=0; _d3=2;}
    if ( d2==K0 && d3==PIP && d1==PIM ) { _flag=3; _d1=1; _d2=2; _d3=0;}
    if ( d3==K0 && d1==PIP && d2==PIM ) { _flag=3; _d1=2; _d2=0; _d3=1;}
    if ( d3==K0 && d2==PIP && d1==PIM ) { _flag=3; _d1=2; _d2=1; _d3=0;}

    if ( d1==KL && d2==PIP && d3==PIM )  { _flag=3; _d1=0; _d2=1; _d3=2;}
    if ( d1==KL && d3==PIP && d2==PIM ) { _flag=3; _d1=0; _d2=2; _d3=1;}
    if ( d2==KL && d1==PIP && d3==PIM ) { _flag=3; _d1=1; _d2=0; _d3=2;}
    if ( d2==KL && d3==PIP && d1==PIM ) { _flag=3; _d1=1; _d2=2; _d3=0;}
    if ( d3==KL && d1==PIP && d2==PIM ) { _flag=3; _d1=2; _d2=0; _d3=1;}
    if ( d3==KL && d2==PIP && d1==PIM ) { _flag=3; _d1=2; _d2=1; _d3=0;}

    if ( d1==KS && d2==PIP && d3==PIM )  { _flag=3; _d1=0; _d2=1; _d3=2;}
    if ( d1==KS && d3==PIP && d2==PIM ) { _flag=3; _d1=0; _d2=2; _d3=1;}
    if ( d2==KS && d1==PIP && d3==PIM ) { _flag=3; _d1=1; _d2=0; _d3=2;}
    if ( d2==KS && d3==PIP && d1==PIM ) { _flag=3; _d1=1; _d2=2; _d3=0;}
    if ( d3==KS && d1==PIP && d2==PIM ) { _flag=3; _d1=2; _d2=0; _d3=1;}
    if ( d3==KS && d2==PIP && d1==PIM ) { _flag=3; _d1=2; _d2=1; _d3=0;}

    if ( d1==KS && d2==KP && d3==KM ) {_flag=5;_d1=0;_d2=1;_d3=2;}
    if ( d1==KS && d3==KP && d2==KM ) {_flag=5;_d1=0;_d2=2;_d3=1;}
    if ( d2==KS && d1==KP && d3==KM ) {_flag=5;_d1=1;_d2=0;_d3=2;}
    if ( d2==KS && d3==KP && d1==KM ) {_flag=5;_d1=1;_d2=2;_d3=0;}
    if ( d3==KS && d1==KP && d2==KM ) {_flag=5;_d1=2;_d2=0;_d3=1;}
    if ( d3==KS && d2==KP && d1==KM ) {_flag=5;_d1=2;_d2=1;_d3=0;}
    
    if ( d1==KL && d2==KP && d3==KM ) {_flag=5;_d1=0;_d2=1;_d3=2;}
    if ( d1==KL && d3==KP && d2==KM ) {_flag=5;_d1=0;_d2=2;_d3=1;}
    if ( d2==KL && d1==KP && d3==KM ) {_flag=5;_d1=1;_d2=0;_d3=2;}
    if ( d2==KL && d3==KP && d1==KM ) {_flag=5;_d1=1;_d2=2;_d3=0;}
    if ( d3==KL && d1==KP && d2==KM ) {_flag=5;_d1=2;_d2=0;_d3=1;}
    if ( d3==KL && d2==KP && d1==KM ) {_flag=5;_d1=2;_d2=1;_d3=0;}    
    
    if ( d1==K0 && d2==KP && d3==KM ) {_flag=5;_d1=0;_d2=1;_d3=2;}
    if ( d1==K0 && d3==KP && d2==KM ) {_flag=5;_d1=0;_d2=2;_d3=1;}
    if ( d2==K0 && d1==KP && d3==KM ) {_flag=5;_d1=1;_d2=0;_d3=2;}
    if ( d2==K0 && d3==KP && d1==KM ) {_flag=5;_d1=1;_d2=2;_d3=0;}
    if ( d3==K0 && d1==KP && d2==KM ) {_flag=5;_d1=2;_d2=0;_d3=1;}
    if ( d3==K0 && d2==KP && d1==KM ) {_flag=5;_d1=2;_d2=1;_d3=0;}
    
  }

  if ( parnum == DP ) {
    //look for K- pi+ pi+
    if ( d1==KB && d2==PIP && d3==PI0 )  { _flag=2; _d1=0; _d2=1; _d3=2;}
    if ( d1==KB && d3==PIP && d2==PI0 ) { _flag=2; _d1=0; _d2=2; _d3=1;}
    if ( d2==KB && d1==PIP && d3==PI0 ) { _flag=2; _d1=1; _d2=0; _d3=2;}
    if ( d2==KB && d3==PIP && d1==PI0 ) { _flag=2; _d1=1; _d2=2; _d3=0;}
    if ( d3==KB && d1==PIP && d2==PI0 ) { _flag=2; _d1=2; _d2=0; _d3=1;}
    if ( d3==KB && d2==PIP && d1==PI0 ) { _flag=2; _d1=2; _d2=1; _d3=0;}

    if ( d1==KL && d2==PIP && d3==PI0 )  { _flag=2; _d1=0; _d2=1; _d3=2;}
    if ( d1==KL && d3==PIP && d2==PI0 ) { _flag=2; _d1=0; _d2=2; _d3=1;}
    if ( d2==KL && d1==PIP && d3==PI0 ) { _flag=2; _d1=1; _d2=0; _d3=2;}
    if ( d2==KL && d3==PIP && d1==PI0 ) { _flag=2; _d1=1; _d2=2; _d3=0;}
    if ( d3==KL && d1==PIP && d2==PI0 ) { _flag=2; _d1=2; _d2=0; _d3=1;}
    if ( d3==KL && d2==PIP && d1==PI0 ) { _flag=2; _d1=2; _d2=1; _d3=0;}

    if ( d1==KS && d2==PIP && d3==PI0 )  { _flag=2; _d1=0; _d2=1; _d3=2;}
    if ( d1==KS && d3==PIP && d2==PI0 ) { _flag=2; _d1=0; _d2=2; _d3=1;}
    if ( d2==KS && d1==PIP && d3==PI0 ) { _flag=2; _d1=1; _d2=0; _d3=2;}
    if ( d2==KS && d3==PIP && d1==PI0 ) { _flag=2; _d1=1; _d2=2; _d3=0;}
    if ( d3==KS && d1==PIP && d2==PI0 ) { _flag=2; _d1=2; _d2=0; _d3=1;}
    if ( d3==KS && d2==PIP && d1==PI0 ) { _flag=2; _d1=2; _d2=1; _d3=0;}

    if ( d1==KM && d2==PIP && d3==PIP )  { _flag=1; _d1=0; _d2=1; _d3=2;}
    if ( d2==KM && d1==PIP && d3==PIP ) { _flag=1; _d1=1; _d2=0; _d3=2;}
    if ( d3==KM && d1==PIP && d2==PIP ) { _flag=1; _d1=2; _d2=0; _d3=1;}
  }

  if ( parnum == DM ) {
    //look for K- pi+ pi+
    if ( d1==K0 && d2==PIM && d3==PI0 )  { _flag=2; _d1=0; _d2=1; _d3=2;}
    if ( d1==K0 && d3==PIM && d2==PI0 ) { _flag=2; _d1=0; _d2=2; _d3=1;}
    if ( d2==K0 && d1==PIM && d3==PI0 ) { _flag=2; _d1=1; _d2=0; _d3=2;}
    if ( d2==K0 && d3==PIM && d1==PI0 ) { _flag=2; _d1=1; _d2=2; _d3=0;}
    if ( d3==K0 && d1==PIM && d2==PI0 ) { _flag=2; _d1=2; _d2=0; _d3=1;}
    if ( d3==K0 && d2==PIM && d1==PI0 ) { _flag=2; _d1=2; _d2=1; _d3=0;}

    if ( d1==KL && d2==PIM && d3==PI0 )  { _flag=2; _d1=0; _d2=1; _d3=2;}
    if ( d1==KL && d3==PIM && d2==PI0 ) { _flag=2; _d1=0; _d2=2; _d3=1;}
    if ( d2==KL && d1==PIM && d3==PI0 ) { _flag=2; _d1=1; _d2=0; _d3=2;}
    if ( d2==KL && d3==PIM && d1==PI0 ) { _flag=2; _d1=1; _d2=2; _d3=0;}
    if ( d3==KL && d1==PIM && d2==PI0 ) { _flag=2; _d1=2; _d2=0; _d3=1;}
    if ( d3==KL && d2==PIM && d1==PI0 ) { _flag=2; _d1=2; _d2=1; _d3=0;}

    if ( d1==KS && d2==PIM && d3==PI0 )  { _flag=2; _d1=0; _d2=1; _d3=2;}
    if ( d1==KS && d3==PIM && d2==PI0 ) { _flag=2; _d1=0; _d2=2; _d3=1;}
    if ( d2==KS && d1==PIM && d3==PI0 ) { _flag=2; _d1=1; _d2=0; _d3=2;}
    if ( d2==KS && d3==PIM && d1==PI0 ) { _flag=2; _d1=1; _d2=2; _d3=0;}
    if ( d3==KS && d1==PIM && d2==PI0 ) { _flag=2; _d1=2; _d2=0; _d3=1;}
    if ( d3==KS && d2==PIM && d1==PI0 ) { _flag=2; _d1=2; _d2=1; _d3=0;}

    if ( d1==KP && d2==PIM && d3==PIM )  { _flag=1; _d1=0; _d2=1; _d3=2;}
    if ( d2==KP && d1==PIM && d3==PIM ) { _flag=1; _d1=1; _d2=0; _d3=2;}
    if ( d3==KP && d1==PIM && d2==PIM ) { _flag=1; _d1=2; _d2=0; _d3=1;}
  }

  if ( parnum == DSP ) {
     if ( d1==KM && d2==KP && d3==PIP ) { _flag=6; _d1=0; _d2=1; _d3=2; }
     if ( d1==KM && d3==KP && d2==PIP ) { _flag=6; _d1=0; _d2=2; _d3=1; }
     if ( d2==KM && d1==KP && d3==PIP ) { _flag=6; _d1=1; _d2=0; _d3=2; }
     if ( d2==KM && d3==KP && d1==PIP ) { _flag=6; _d1=1; _d2=2; _d3=0; }
     if ( d3==KM && d1==KP && d2==PIP ) { _flag=6; _d1=2; _d2=0; _d3=1; }
     if ( d3==KM && d2==KP && d1==PIP ) { _flag=6; _d1=2; _d2=1; _d3=0; }
  }

  if ( parnum == DSM ) {
     if ( d1==KP && d2==KM && d3==PIM ) { _flag=6; _d1=0; _d2=1; _d3=2; }
     if ( d1==KP && d3==KM && d2==PIM ) { _flag=6; _d1=0; _d2=2; _d3=1; }
     if ( d2==KP && d1==KM && d3==PIM ) { _flag=6; _d1=1; _d2=0; _d3=2; }
     if ( d2==KP && d3==KM && d1==PIM ) { _flag=6; _d1=1; _d2=2; _d3=0; }
     if ( d3==KP && d1==KM && d2==PIM ) { _flag=6; _d1=2; _d2=0; _d3=1; }
     if ( d3==KP && d2==KM && d1==PIM ) { _flag=6; _d1=2; _d2=1; _d3=0; }
  }

  if ( _flag==6) {
    _kkpi_params.push_back(EvtFlatteParam(MPI, MPI, 0.406));
    _kkpi_params.push_back(EvtFlatteParam(MKP, MKP, 0.800));
  }

  if ( _flag==0) {
    report(ERROR,"EvtGen") << "EvtDDaltiz: Invalid mode."<<endl;
    assert(0);
  }
}

void EvtDDalitz::initProbMax() {

//probmax different for different modes!  

  if ( _flag==1 ) {setProbMax(9.6);}
  if ( _flag==2 ) {setProbMax(147.9);}
  if ( _flag==3 ) {setProbMax(5000.0);}
  if ( _flag==4 ) {setProbMax(3000.0);}
  if ( _flag==5 ) {setProbMax(10000000.0);}
  if ( _flag==6 ) {setProbMax(50000.0);}

}

void EvtDDalitz::decay( EvtParticle *p){

  static EvtId BP = EvtPDL::getId("B+");                                       
  static EvtId BM = EvtPDL::getId("B-");                                       
  static EvtId B0 = EvtPDL::getId("B0");                                       
  static EvtId B0B = EvtPDL::getId("anti-B0");         
  static EvtId DM=EvtPDL::getId("D-");
  static EvtId DP=EvtPDL::getId("D+");
  static EvtId D0=EvtPDL::getId("D0");
  static EvtId D0B=EvtPDL::getId("anti-D0");
  static EvtId KM=EvtPDL::getId("K-");
  static EvtId KP=EvtPDL::getId("K+");
  static EvtId K0=EvtPDL::getId("K0");
  static EvtId KB=EvtPDL::getId("anti-K0");
  static EvtId PIM=EvtPDL::getId("pi-");
  static EvtId PIP=EvtPDL::getId("pi+");
  static EvtId PI0=EvtPDL::getId("pi0");

  double oneby2 = 0.707106782;

  bool isBToDK=false; 
  if ( p -> getParent () ) {                                                   
    EvtId parId = p -> getParent()->getId ();                              
    if ( ( BP == parId ) || ( BM == parId ) || ( B0 == parId ) ||              
	               ( B0B == parId ) )
      if (EvtDecayTable::getDecayFunc(p->getParent())->getName() == "BTODDALITZCPK") isBToDK=true;   
  }                                                                            
  

//same structure for all of these decays

  p->initializePhaseSpace(getNDaug(),getDaugs());
  EvtVector4R moms1 = p->getDaug(_d1)->getP4();
  EvtVector4R moms2 = p->getDaug(_d2)->getP4();
  EvtVector4R moms3 = p->getDaug(_d3)->getP4();

  EvtVector4R p4_p;
  p4_p.set(p->mass(),0.0,0.0,0.0);

  EvtComplex amp(1.0,0.0);

//now determine which D and which decay

//data from Anjos et al, Phys.Rev.D 1993, v.48,num.1,p.56 (E691 resuls)
//for D+ -> K- pi+ pi+, and from Adler et al, Phys.Lett. B196 (1987), 107
//(Mark III results) for D+ -> K0bar pi+ pi0. 
  //CLEO results for D0->k-pi+pi0

  if ( _flag==1) {

//have a D+ -> K- pi+ pi+ decay, or charge conjugate
//Anjos etal e691 - Phys Rev D48, 56 (1993) 
    EvtResonance DplusRes11(p4_p,moms1,moms2,0.78,-60.0,0.0498,0.89610,1);
    EvtResonance DplusRes12(p4_p,moms3,moms1,0.78,-60.0,0.0498,0.89610,1);//K*(892)
    
    EvtResonance DplusRes21(p4_p,moms1,moms2,0.53,132.0,0.287,1.429,0);
    EvtResonance DplusRes22(p4_p,moms3,moms1,0.53,132.0,0.287,1.429,0);//K*(1430)
    
    EvtResonance DplusRes31(p4_p,moms1,moms2,0.47,-51.0,0.323,1.714,1);
    EvtResonance DplusRes32(p4_p,moms3,moms1,0.47,-51.0,0.323,1.714,1);//K*(1680)
    
    amp = amp + oneby2*(-DplusRes11.resAmpl()+DplusRes12.resAmpl()) + oneby2*(DplusRes21.resAmpl() + DplusRes22.resAmpl()) + oneby2*(-DplusRes31.resAmpl()+ DplusRes32.resAmpl());
    
 }
  
  if ( _flag==2) {

//have a D+ -> K0bar pi+ pi0 decay 
//adler etal MarkIII - Phys Lett B196, 107 (1987)    
// Results in this paper:
//   Kbar rho+    FitFraction = 68+/-8+/-12    Phase   0
//   Kbar* pi+                  19+/-6+/-6            43+/-23
//   nonres                     13+/-7+/-8           250+/-19   
// These numbers below seem not to be exactly the same
// the phases are equiv to -106=254 and 41
// 
    EvtResonance DplusKpipi0Res1(p4_p,moms2,moms3,1.00,0.00,0.1512,0.7699,1); //rho+  
    EvtResonance DplusKpipi0Res2(p4_p,moms3,moms1,0.8695,0.7191,0.0498,0.89159,1); //K*0
    
    amp = 0.9522*EvtComplex(cos(-1.8565),sin(-1.8565)) + 1.00*DplusKpipi0Res1.relBrWig(0) + 0.8695*EvtComplex(cos(0.7191),sin(0.7191))*DplusKpipi0Res2.relBrWig(1);
    
  }

  if(_flag==3) {
    // D0 -> K0 pi+ pi- + CC                                                                       
    // If it does not come from a B->DK, decay it as D0 or D0bar separatly                         
    // if p4_p is D0, moms1 is K0, moms2 is pi-, moms3 is pi+                                      
    // if p4_p is D0bar, moms1 is K0, moms2 is pi+, moms3 is pi-                                   

    if ( isBToDK ) {
      // Gamma angle in rad.                                                                       
      double gamma = EvtDecayTable::getDecayFunc( p->getParent() )
        -> getArg( 0 )  ;
      // Strong phase in rad.                                                                      
      double delta =  EvtDecayTable::getDecayFunc( p->getParent() )
        -> getArg( 1 )  ;
      // Ratio between B->D0K and B->D0barK                                                        
      double A     =  EvtDecayTable::getDecayFunc( p->getParent() )
        -> getArg( 2 )  ;

      EvtComplex Factor( fabs( A ) * cos ( delta ) ,
                         fabs( A ) * sin ( delta ) ) ;

      if ( ( p->getParent()->getId() == BP ) ||
           ( p->getParent()->getId() == B0 ) ) {
        // the ratio D/Dbar                                                                        
        Factor = Factor * EvtComplex( cos ( gamma ) , sin ( gamma ) ) ;
        if ( p->getId() == D0 ) {
          // the flavor of the particle has no meaning. But we need                                
          // it to know which daughter is pi+ or pi-                                               
          // M( B+ or B0 ) = f(Dbar) + factor * f(D)                                               
          // f(Dbar) = amplDtoK0PiPi(pD, K0, pi+, pi-)                                             
          // f(D)    = amplDtoK0PiPi(pD, K0, pi-, pi+)                                             
          // Then ...                                                
          amp = amplDtoK0PiPi( p4_p , moms1 , moms3 , moms2 ) +
            Factor * amplDtoK0PiPi( p4_p , moms1 , moms2 , moms3 ) ;
        }
        else {
          amp = amplDtoK0PiPi( p4_p , moms1 , moms2 , moms3 ) +
            Factor * amplDtoK0PiPi( p4_p , moms1 , moms3 , moms2 ) ;
        }
      }
      else if ( ( p->getParent() -> getId() == BM ) ||
                ( p->getParent() -> getId() == B0B ) ) {
        Factor = Factor * EvtComplex( cos ( gamma ) , - sin ( gamma ) ) ;
        // here M( B- or B0bar ) = f(D) + factor * f(Dbar) then ...                                
        if ( p->getId() == D0 ) {
          amp = amplDtoK0PiPi( p4_p , moms1 , moms2 , moms3 ) +
            Factor * amplDtoK0PiPi( p4_p , moms1 , moms3 , moms2 ) ;
        }
        else {
          amp = amplDtoK0PiPi( p4_p , moms1 , moms3 , moms2 ) +
            Factor * amplDtoK0PiPi( p4_p , moms1 , moms2 , moms3 ) ;
        }
      }
    }
    else {
      amp = amplDtoK0PiPi( p4_p , moms1 , moms2 , moms3 ) ;
    }
  }

  
  if(_flag==4) {

    EvtResonance2 DKpipi0Res1(p4_p,moms2,moms3,1.0  ,0.0   ,0.1507,0.770 ,1); //rho
    EvtResonance2 DKpipi0Res2(p4_p,moms1,moms2,0.39, -0.2  ,0.0505,0.8961,1); //k*0
    EvtResonance2 DKpipi0Res3(p4_p,moms1,moms3,0.44, 163.0 ,0.050 ,0.8915,1); //k*-
    
    EvtResonance2 DKpipi0Res4(p4_p,moms1,moms3,0.77 ,55.5  ,0.294 ,1.412 ,0); //k01430-
    EvtResonance2 DKpipi0Res5(p4_p,moms1,moms2,0.85 ,166.0 ,0.294 ,1.412 ,0); //k01430bar
    EvtResonance2 DKpipi0Res6(p4_p,moms2,moms3,2.5  ,171.0 ,0.240 ,1.700 ,1); //rho1700
    EvtResonance2 DKpipi0Res7(p4_p,moms1,moms3,2.5  ,103.0 ,0.322 ,1.717 ,1); //K*1680-
    
    
    
    double pi180inv = 1.0/EvtConst::radToDegrees;
    
    amp = EvtComplex(1.75*cos(31.2*pi180inv),1.75*sin(31.2*pi180inv)) 
      + DKpipi0Res1.resAmpl() + DKpipi0Res2.resAmpl() + DKpipi0Res3.resAmpl()
      + DKpipi0Res4.resAmpl() + DKpipi0Res5.resAmpl() 
      + DKpipi0Res6.resAmpl()
      + DKpipi0Res7.resAmpl();
    
  }
 
  if(_flag==5) {

    // D0 -> K0 K+ K- + CC                                                                         
    // If it does not come from a B->DK, decay it as D0 or D0bar separatly                         
    // if p4_p is D0, moms1 is K0, moms2 is pi-, moms3 is pi+                                      
    // if p4_p is D0bar, moms1 is K0, moms2 is pi+, moms3 is pi-                                   

    if ( isBToDK ){
      // Gamma angle in rad.                                                                       
      double gamma = EvtDecayTable::getDecayFunc( p->getParent() )
        -> getArg( 0 )  ;
      // Strong phase in rad.                                                                      
      double delta =  EvtDecayTable::getDecayFunc( p->getParent() )
        -> getArg( 1 )  ;
      // Ratio between B->D0K and B->D0barK                                                        
      double A     =  EvtDecayTable::getDecayFunc( p->getParent() )
        -> getArg( 2 )  ;

      EvtComplex Factor( fabs( A ) * cos ( delta ) ,
                         fabs( A ) * sin ( delta ) ) ;

      if ( ( p->getParent()->getId() == BP ) ||
           ( p->getParent()->getId() == B0 ) ) {
        // the ratio D/Dbar                                                                        
        Factor = Factor * EvtComplex( cos ( gamma ) , sin ( gamma ) ) ;
        if ( p->getId() == D0 ) {
          // the flavor of the particle has no meaning. But we need                                
          // it to know which daughter is pi+ or pi-                                               
          // M( B+ or B0 ) = f(Dbar) + factor * f(D)                                               
          // f(Dbar) = amplDtoK0PiPi(pD, K0, K+, K-)                                               
          // f(D)    = amplDtoK0PiPi(pD, K0, K-, K+)                                               
          // Then ...                                                                              
          amp = amplDtoK0KK( p4_p , moms1 , moms3 , moms2 ) +
            Factor * amplDtoK0KK( p4_p , moms1 , moms2 , moms3 ) ;
        }
        else {
          amp = amplDtoK0KK( p4_p , moms1 , moms2 , moms3 ) +
            Factor * amplDtoK0KK( p4_p , moms1 , moms3 , moms2 ) ;
        }
      }
      else if ( ( p->getParent() -> getId() == BM ) ||
                ( p->getParent() -> getId() == B0B ) ) {
        Factor = Factor * EvtComplex( cos ( gamma ) , - sin ( gamma ) ) ;
        // here M( B- or B0bar ) = f(D) + factor * f(Dbar) then ...                                
        if ( p->getId() == D0 ) {
          amp = amplDtoK0KK( p4_p , moms1 , moms2 , moms3 ) +
            Factor * amplDtoK0KK( p4_p , moms1 , moms3 , moms2 ) ;
        }
        else {
          amp = amplDtoK0KK( p4_p , moms1 , moms3 , moms2 ) +
            Factor * amplDtoK0KK( p4_p , moms1 , moms2 , moms3 ) ;
        }
      }
    }
    else {
      amp = amplDtoK0KK( p4_p , moms1 , moms2 , moms3 ) ;
    }
  }




  // Ds -> K K pi
  if(_flag==6) {
     EvtResonance2 DsKKpiRes1(p4_p, moms3, moms1, 1.0, 0.0, 0.0455, 0.8944, 1, true); // K*(892)
     EvtResonance2 DsKKpiRes2(p4_p, moms3, moms1, 1.48, 138., 0.290, 1.414, 0); // K*_0(1430)
     EvtFlatte     DsKKpiRes3(p4_p, moms1, moms2, 5.07, 156., 0.965, _kkpi_params); // f_0(980)
     EvtResonance2 DsKKpiRes4(p4_p, moms1, moms2, 1.15, -10., 0.00426, 1.019455, 1, true); // phi(1020)
     EvtResonance2 DsKKpiRes5(p4_p, moms1, moms2, 1.28, 53., 0.265, 1.350, 0); // f_0(1370)
     EvtResonance2 DsKKpiRes6(p4_p, moms1, moms2, 1.19, 87., 0.137, 1.724, 0); // f_0(1710)
     amp = DsKKpiRes1.resAmpl() + DsKKpiRes2.resAmpl() + DsKKpiRes3.resAmpl()
        + DsKKpiRes4.resAmpl() + DsKKpiRes5.resAmpl() + DsKKpiRes6.resAmpl();

  }


  vertex(amp);

  return ;
}

EvtComplex EvtDDalitz::amplDtoK0PiPi(EvtVector4R p4_p,  EvtVector4R moms1, 
                                     EvtVector4R moms2, EvtVector4R moms3) {

    //K*(892)-
    EvtResonance2 DK2piRes1(p4_p,moms1,moms2,1.418,-190.0,0.0508,0.89166,1);
    //K0*(1430)
    EvtResonance2 DK2piRes2(p4_p,moms1,moms2,1.818,-337.0,0.294 ,1.412  ,0);
    //K2*(1430)
    EvtResonance2 DK2piRes3(p4_p,moms1,moms2,0.909,  -5.0,0.0985,1.4256 ,2);
    //K*(1680)
    EvtResonance2 DK2piRes4(p4_p,moms1,moms2,5.091,-166.0,0.322 ,1.717  ,1);
    //DCS K*(892)
    EvtResonance2 DK2piRes5(p4_p,moms1,moms3,0.100, -19.0,0.0508,0.89166,1);
    
    //Rho
    EvtResonance2 DK2piRes6(p4_p,moms3,moms2,0.909,-340.0,0.1502,0.7693,1);
    //Omega
    EvtResonance2 DK2piRes7(p4_p,moms3,moms2,.0336,-226.0,0.00844,0.78257,1);
    //f0(980)
    EvtResonance2 DK2piRes8(p4_p,moms3,moms2,0.309,-152.0,0.05,0.977,0);
    //f0(1370)
    EvtResonance2 DK2piRes9(p4_p,moms3,moms2,1.636,-255.0,0.272,1.31,0);
    //f2(1270)
    EvtResonance2 DK2piRes10(p4_p,moms3,moms2,0.636,-32.0,0.1851,1.2754,2);
    
    return EvtComplex(1.0,0.0) + 
      DK2piRes1.resAmpl() + DK2piRes2.resAmpl() +
      DK2piRes3.resAmpl() + DK2piRes4.resAmpl() + 
      DK2piRes5.resAmpl() + DK2piRes6.resAmpl() + 
      DK2piRes7.resAmpl() + DK2piRes8.resAmpl() + 
      DK2piRes9.resAmpl() + DK2piRes10.resAmpl();
}

//
// BaBar decay amplitudes for D0->Ks K+ K-
//
// p4_p is D0
// moms1 is K0s
// moms2 is K+
// moms3 is K-
// Amplitudes and phases are taken from BaBar hep-ex/0207089
// with convention : Non Resonant = Amp 1. / Phase 0. 

EvtComplex EvtDDalitz::amplDtoK0KK(EvtVector4R p4_p,  EvtVector4R moms1, 
                                     EvtVector4R moms2, EvtVector4R moms3) {

    //phi
    EvtResonance DK0KKRes1( p4_p, moms2, moms3, 113.75, -40.0, 0.0043,
                            1.019456, 1 ) ;
    //a0(980)
    EvtResonance DK0KKRes2( p4_p, moms2, moms3, 152.25, 69.0, 0.1196 , 0.9847,
                            0 ) ;
    //f0(980)
    EvtResonance DK0KKRes3( p4_p, moms2, moms3, 30.5, -201.0, 0.05, 0.980 , 
                            0 ) ;
    //a0(980)+
    EvtResonance DK0KKRes4( p4_p, moms1, moms2, 85.75, -93.0, 0.1196 , 0.9847,
                            0 ) ;
    //a0(980)-
    EvtResonance DK0KKRes5( p4_p, moms3, moms1, 8. , -53.0 ,0.1196, 0.9847,
                            0 ) ;

    return EvtComplex(1.0,0.0) +
      DK0KKRes1.resAmpl() + DK0KKRes2.resAmpl() +
      DK0KKRes3.resAmpl() + DK0KKRes4.resAmpl() + 
      DK0KKRes5.resAmpl() ;

}
