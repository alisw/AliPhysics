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
// Module: EvtJscontCDF.cc
//
// Description: Routine to generate e+e- --> q\barq  via Jetset
//
// Modification history:
//
//    PCK     August 4, 1997        Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <stdio.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenModels/EvtJscontCDF.hh"
#include "EvtGenModels/EvtJetSetCDF.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtReport.hh"
#include <string>
using namespace std;

extern "C" {
  extern void continuumcdf_(double *,int *,int *,int *,
		       double *,double *,double *,double *);
}

extern "C" {
  extern void lygive_(const char *cnfgstr,int length);
}


EvtJscontCDF::~EvtJscontCDF() {}

std::string EvtJscontCDF::getName(){

  return "JSCONTCDF";     

}

EvtDecayBase* EvtJscontCDF::clone(){

  return new EvtJscontCDF;

}

void EvtJscontCDF::init(){

  // check that there are 1 argument

  checkNArg(1,2);

}


void EvtJscontCDF::initProbMax(){

  noProbMax();

}


void EvtJscontCDF::decay( EvtParticle *p){

  EvtJetSetCDF::jetSetInit();
  static int first=1;

  if (first){
    first=0;

    float val=0.6;
    if ( getNArg()>1) {
      val=getArg(1);
    }
    char vak[20];
    sprintf(vak,"PARJ(13)=%f",val);
    std::string temp(vak);
    lygive_(temp.c_str(),strlen(temp.c_str()));
  }
  EvtVector4R p4[100];
  
  double energy=p->mass();

  int flavor;

  int i,more;
  int ndaugjs;
  int kf[100];
  EvtId id[100];
  int type[MAX_DAUG]; 

  flavor=(int)getArg(0);

  double px[100],py[100],pz[100],e[100];

  if ( p->getNDaug() != 0 ) { return;}
  do{

    continuumcdf_(&energy,&flavor,&ndaugjs,kf,px,py,pz,e);

    for(i=0;i<ndaugjs;i++){

      id[i]=EvtPDL::evtIdFromStdHep(kf[i]);

      type[i]=EvtPDL::getSpinType(id[i]);

      // have to protect against negative mass^2 for massless particles
      // i.e. neutrinos and photons.
      // this is uggly but I need to fix it right now....

      if (px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]>=e[i]*e[i]){

        e[i]=sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])+0.0000000000001;

      }

      p4[i].set(e[i],px[i],py[i],pz[i]);

    }

    int channel=EvtDecayTable::inChannelList(p->getId(),ndaugjs,id);

    more=((channel!=-1)&&(channel!=p->getChannel()));


  }while(more);

  p->makeDaughters(ndaugjs,id);

  for(i=0;i<ndaugjs;i++){
    p->getDaug(i)->init( id[i], p4[i] );
  }
  return ;
}




