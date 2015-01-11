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
// Module: EvtCBTo3piP00.cc
//
// Description: Routine to decay B+/-->pi0 pi0 pi+/-
//              and has CP violation.
//
// Modification history:
//
//    RYD,Versille     May 6, 1997         Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtCBTo3piP00.hh"
#include <string>

//Below you will have do modify the declaration to be appropriate
//for your new routine for the calculation of the amplitude

#ifdef WIN32
extern "C" {
  extern void EVT3PIONSP00(double *,int *,
			   double *,
			   double *,double *,
			   double *,double *,
			   double *,double *,double *,double *);
}
#else
extern "C" {
  extern void evt3pionsp00_(double *,int *,
			     double *,
			     double *,double *,
			     double *,double *,
			     double *,double *,double *,double *);
}
#endif

EvtCBTo3piP00::~EvtCBTo3piP00() {}

std::string EvtCBTo3piP00::getName(){

  return "CB3PI-P00";     

}


EvtDecayBase* EvtCBTo3piP00::clone(){

  return new EvtCBTo3piP00;

}

void EvtCBTo3piP00::init(){

  // check that there are 1 argument
  checkNArg(1);
  checkNDaug(3);

  checkSpinParent(EvtSpinType::SCALAR);
  
  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::SCALAR);
  checkSpinDaughter(2,EvtSpinType::SCALAR);

}



void EvtCBTo3piP00::initProbMax(){


  setProbMax(1.5);

}


void EvtCBTo3piP00::decay( EvtParticle *p ){

  //added by Lange Jan4,2000
  static EvtId BM=EvtPDL::getId("B-");
  static EvtId BP=EvtPDL::getId("B+");

  EvtParticle *pi1,*pi2,*pi3;

  p->makeDaughters(getNDaug(),getDaugs());
  pi1=p->getDaug(0);
  pi2=p->getDaug(1);
  pi3=p->getDaug(2);

  EvtVector4R p4[3];
  double alpha = getArg(0);
  int iset;
  static int first=1;

  if (first==1) {
    iset=10000;
    first=0;
  }
  else{
    iset=0;
  }

  double p4pi1[4],p4Gamma11[4],p4Gamma12[4];
  double p4Gamma21[4],p4Gamma22[4];

  double realA,imgA,realbarA,imgbarA;

  evt3pionsp00_(&alpha,&iset,
		 p4pi1,
		 p4Gamma11,p4Gamma12,
		 p4Gamma21,p4Gamma22,
		 &realA,&imgA,&realbarA,&imgbarA);

  p4[0].set(p4pi1[3],p4pi1[0],p4pi1[1],p4pi1[2]);
  p4[1].set(p4Gamma11[3]+p4Gamma12[3],
	    p4Gamma11[0]+p4Gamma12[0],
	    p4Gamma11[1]+p4Gamma12[1],
	    p4Gamma11[2]+p4Gamma12[2]);
  p4[2].set(p4Gamma21[3]+p4Gamma22[3],
	    p4Gamma21[0]+p4Gamma22[0],
	    p4Gamma21[1]+p4Gamma22[1],
	    p4Gamma21[2]+p4Gamma22[2]);

  pi1->init( getDaug(0), p4[0] );
  pi2->init( getDaug(1), p4[1] );
  pi3->init( getDaug(2), p4[2] );

  EvtComplex A(realA,imgA);
  EvtComplex Abar(realbarA, imgbarA);
   
  EvtComplex  amp;
  if(p->getId()==BP)
    {
      amp = A;
    }
  if(p->getId()==BM)
    {
      amp = Abar;
    }

  vertex(amp);

  return ;
}

