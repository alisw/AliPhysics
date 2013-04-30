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
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenModels/EvtPycont.hh"
#include "EvtGenModels/EvtPythia.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtReport.hh"
#include <string>
#include <iostream>
using std::endl;

extern "C" {
  extern void pystat_(int &);
  extern struct
  {
    int dc[18];
  } decaych_;
}



EvtPycont::~EvtPycont() {
  //  int i=1;
  //  pystat_(i);
}

std::string EvtPycont::getName()
{
  return "PYCONT";  
}

EvtDecayBase* EvtPycont::clone()
{
  return new EvtPycont;
}

void EvtPycont::init()
{
  // check that there are 1 argument
  if ( getNArg() != 12 && getNArg() != 0 ) {
    report(ERROR,"EvtGen") << "EvtPYCONT expects "
			   << " 12 arguments (d u s c b t e nu_e mu nu_mu tau nu_tau) but found: "
			   << getNArg() <<endl;

  }
  checkNArg(0,12);

   for( int i=0; i<18; i++)
     decaych_.dc[i]=0;
   if ( getNArg() == 12 ) {
     decaych_.dc[0]=(int)getArg(0);
     decaych_.dc[1]=(int)getArg(1);
     decaych_.dc[2]=(int)getArg(2);
     decaych_.dc[3]=(int)getArg(3);
     decaych_.dc[4]=(int)getArg(4);
     decaych_.dc[5]=(int)getArg(5);
     decaych_.dc[10]=(int)getArg(6);
     decaych_.dc[11]=(int)getArg(7);
     decaych_.dc[12]=(int)getArg(8);
     decaych_.dc[13]=(int)getArg(9);
     decaych_.dc[14]=(int)getArg(10);
     decaych_.dc[15]=(int)getArg(11);
   }
   else{
     decaych_.dc[0]=1;
     decaych_.dc[1]=1;
     decaych_.dc[2]=1;
     decaych_.dc[3]=1;
   }

}

void EvtPycont::initProbMax()
{
  noProbMax();
}

void EvtPycont::decay( EvtParticle *p)
{
  EvtPythia::pythiaInit(0);
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

