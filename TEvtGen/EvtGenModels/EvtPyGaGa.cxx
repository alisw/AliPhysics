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
// Module: EvtPycont.cc
//
// Description: Routine to generate e+e- --> q\barq  via Jetset
//
// Modification history:
//
//    PCK     August 4, 1997        Module created
//    RS      October 28, 2002      copied from EvtJscont.cc
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenModels/EvtPyGaGa.hh"
#include "EvtGenModels/EvtPythia.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtReport.hh"
#include <string.h>
#include <iostream>

extern "C" {
  extern void pystat_(int &);
  extern struct
  {
    int dc[18];
  } decaych_;
}


EvtPyGaGa::~EvtPyGaGa()
{
  int i=1;
  pystat_(i);
}

std::string EvtPyGaGa::getName()
{
  return "PYGAGA";  
}

EvtDecayBase* EvtPyGaGa::clone()
{
  return new EvtPyGaGa;
}

void EvtPyGaGa::init()
{
  // check that there are 1 argument
  checkNArg(0);
  for( int i=0; i<18; i++)
    decaych_.dc[i]=0;
}

void EvtPyGaGa::initProbMax()
{
  noProbMax();
}

void EvtPyGaGa::decay( EvtParticle *p)
{
  EvtPythia::pythiaInit(1);
  EvtVector4R p4[100];
  
  double energy=p->mass();
  
  int i,more;
  int ndaugjs;
  int kf[100];
  EvtId id[100];
  int type[MAX_DAUG];
  
  double px[100],py[100],pz[100],e[100];
  
  if ( p->getNDaug() != 0 ) { return;}
  do{
    EvtPythia::pythiacont(&energy,&ndaugjs,kf,px,py,pz,e);
    
    for(i=0;i<ndaugjs;i++)
      {
	
	id[i]=EvtPDL::evtIdFromStdHep(kf[i]);
	
	type[i]=EvtPDL::getSpinType(id[i]);
	
	// have to protect against negative mass^2 for massless particles
	// i.e. neutrinos and photons.
	// this is uggly but I need to fix it right now....
	
	if (px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i]>=e[i]*e[i])
	  e[i]=sqrt(px[i]*px[i]+py[i]*py[i]+pz[i]*pz[i])+0.0000000000001;
	
	p4[i].set(e[i],px[i],py[i],pz[i]);
	
      }
    
    int channel=EvtDecayTable::inChannelList(p->getId(),ndaugjs,id);
    
    more=((channel!=-1)&&(channel!=p->getChannel()));
    
  }while(more);
  
  p->makeDaughters(ndaugjs,id);
  
  for(i=0;i<ndaugjs;i++)
    p->getDaug(i)->init( id[i], p4[i] );
  
  return ;
}

