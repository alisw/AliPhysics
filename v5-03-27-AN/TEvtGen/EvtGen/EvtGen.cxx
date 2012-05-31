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
// Module: EvtGen.cc
//
// Description: Main class to provide user interface to EvtGen.
//
// Modification history:
//
//    RYD     March 24, 1998        Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdio.h>
#include <fstream>
#include <math.h>
#include "EvtGenBase/EvtComplex.hh"
#include <stdlib.h>
#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVectorParticle.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtStdHep.hh"
#include "EvtGenBase/EvtSecondary.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtId.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtSimpleRandomEngine.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
//#include "CLHEP/Vector/LorentzVector.h"
#include "TLorentzVector.h"
#include "EvtGenModels/EvtModelReg.hh"
#include "EvtGenBase/EvtStatus.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtRadCorr.hh"
#include "EvtGenModels/EvtPHOTOS.hh"
using std::endl;
using std::fstream;
using std::ifstream;


extern "C" void begevtgenstore_(int *,int *,int *,int *,
				int *,int *,int *,int *,int *,
                                double *,double *,double *, 
                                double *,double *,double *, 
                                double *,double *,double *);

extern "C" {
extern void evtgen_(float svertex[3],float *e_cms,float *beta_zs,
                    int *mode);
}

EvtGen::~EvtGen(){

  //This is a bit uggly, should not do anything
  //in a destructor. This will fail if EvtGen is made a static
  //because then this destructor might be called _after_
  //the destructoin of objects that it depends on, e.g., EvtPDL.

  if (getenv("EVTINFO")){
    EvtDecayTable::printSummary();
  }

}

EvtGen::EvtGen(const char* const decayName,
	       const char* const pdtTableName,
	       EvtRandomEngine* randomEngine,
	       EvtAbsRadCorr* isrEngine,
	       const std::list<EvtDecayBase*>* extraModels){


  report(INFO,"EvtGen") << "Initializing EvtGen"<<endl;

  report(INFO,"EvtGen") << "Storing known decay models"<<endl;
  EvtModelReg dummy(extraModels);

  report(INFO,"EvtGen") << "Main decay file name  :"<<decayName<<endl;
  report(INFO,"EvtGen") << "PDT table file name   :"<<pdtTableName<<endl;
  
  report(INFO,"EvtGen") << "Initializing RadCorr=PHOTOS"<<endl;
  if (isrEngine==0){
    static EvtPHOTOS defaultRadCorrEngine;
    EvtRadCorr::setRadCorrEngine(&defaultRadCorrEngine);
  }
  else{
    EvtRadCorr::setRadCorrEngine(isrEngine);    
  }

  _pdl.readPDT(pdtTableName);

  EvtDecayTable::readDecayFile(decayName,false);

  if (randomEngine==0){
    static EvtSimpleRandomEngine defaultRandomEngine;
    EvtRandom::setRandomEngine((EvtRandomEngine*)&defaultRandomEngine);
    report(INFO,"EvtGen") <<"No random engine given in "
			  <<"EvtGen::EvtGen constructor, "
			  <<"will use default EvtSimpleRandomEngine."<<endl;
  }
  else{
    EvtRandom::setRandomEngine(randomEngine);    
  }

  report(INFO,"EvtGen") << "Done initializing EvtGen"<<endl;

}


void EvtGen::readUDecay(const char* const uDecayName){

  ifstream indec;

  if ( uDecayName[0] == 0) {
    report(INFO,"EvtGen") << "Is not reading a user decay file!"<<endl;
  }
  else{  
    indec.open(uDecayName);
    if (indec) {
      EvtDecayTable::readDecayFile(uDecayName,true);
    }    
    else{
      
      report(INFO,"EvtGen") << "Can not find UDECAY file '"
			    <<uDecayName<<"'.  Stopping"<<endl;
      ::abort();
    }
  }
  
}


void EvtGen::generateDecay(int stdhepid, 
			   EvtVector4R P, 
			   EvtVector4R D,
			   EvtStdHep *evtStdHep,
			   EvtSpinDensity *spinDensity ){

 
  EvtParticle *p;

  if(spinDensity==0){    
    p=EvtParticleFactory::particleFactory(EvtPDL::evtIdFromStdHep(stdhepid),P);
  }
  else{
    p=EvtParticleFactory::particleFactory(EvtPDL::evtIdFromStdHep(stdhepid),
					  P,*spinDensity);
  }

  generateDecay(p);
  //  p->Decay();

  evtStdHep->init();

  p->makeStdHep(*evtStdHep);
  
  evtStdHep->translate(D);
  
  p->deleteTree();


}

void EvtGen::generateDecay(EvtParticle *p){

  int times=0;
  do{
    times+=1;
    EvtStatus::initRejectFlag();

    p->decay();
    //ok then finish.
    if ( EvtStatus::getRejectFlag()==0 ) { 
      times=0;
    }
    else{   
      for (size_t ii=0;ii<p->getNDaug();ii++){
	EvtParticle *temp=p->getDaug(ii);
	temp->deleteTree();
      }
      p->resetFirstOrNot();
      p->resetNDaug();
      
    }

    if ( times==10000) {
      report(ERROR,"EvtGen") << "Your event has been rejected 10000 times!"<<endl;
      report(ERROR,"EvtGen") << "Will now abort."<<endl;
      ::abort();
      times=0;
    }
  } while (times);

}



//void EvtGen::generateEvent(int stdhepid,CLHEP::HepLorentzVector P,CLHEP::HepLorentzVector D){
void EvtGen::generateEvent(int stdhepid,TLorentzVector P,TLorentzVector D){
  EvtParticle *root_part;
  EvtVectorParticle *vector_part;
  
  vector_part=new EvtVectorParticle;
  EvtVector4R p_init;

  //p_init.set(P.t(),P.x(),P.y(),P.z());
  p_init.set(P.E(),P.Px(),P.Py(),P.Pz());
  vector_part->init(EvtPDL::evtIdFromStdHep(stdhepid),p_init);
  
  root_part=(EvtParticle *)vector_part;
  
  root_part->setVectorSpinDensity();      

  generateEvent(root_part,D);

  root_part->deleteTree();  

}

//void EvtGen::generateEvent(EvtParticle *root_part,CLHEP::HepLorentzVector D){
void EvtGen::generateEvent(EvtParticle *root_part,TLorentzVector D){
  int i;  
  
  static int nevent=0;
  
  nevent++;  

  static EvtStdHep evtstdhep;
  //  static EvtSecondary evtsecondary;

  int j;
  int istat;
  int partnum;
  double px,py,pz,e,m;
  double x,y,z,t;

  EvtVector4R p4,x4;
    
  generateDecay(root_part);
  //  root_part->Decay();
  
  int          npart=0;
  
  EvtId        list_of_stable[10];
  EvtParticle* stable_parent[10];
    
    
  list_of_stable[0]=EvtId(-1,-1);
  stable_parent[0]=0;

  evtstdhep.init();
  //  evtsecondary.init();
  // root_part->makeStdHep(evtstdhep,evtsecondary,list_of_stable);
  root_part->makeStdHep(evtstdhep);

  //report(INFO,"EvtGen") << evtstdhep;
  //report(INFO,"EvtGen") << evtsecondary;
  
  npart=evtstdhep.getNPart();

  for(i=0;i<evtstdhep.getNPart();i++){

    j=i+1;

    int jmotherfirst=evtstdhep.getFirstMother(i)+1;
    int jmotherlast=evtstdhep.getLastMother(i)+1;
    int jdaugfirst=evtstdhep.getFirstDaughter(i)+1;
    int jdauglast=evtstdhep.getLastDaughter(i)+1;

    partnum=evtstdhep.getStdHepID(i);

    istat=evtstdhep.getIStat(i);

    p4=evtstdhep.getP4(i);
    x4=evtstdhep.getX4(i);
      
    px=p4.get(1);
    py=p4.get(2);
    pz=p4.get(3);
    e=p4.get(0);
	  
    x=x4.get(1)+D.X();
    y=x4.get(2)+D.Y();
    z=x4.get(3)+D.Z();
    t=x4.get(0)+D.T();
      
    m=p4.mass();

    begevtgenstore_(&j,&nevent,&npart,&istat,
		    &partnum,&jmotherfirst,&jmotherlast,
		    &jdaugfirst,&jdauglast,
		    &px,&py,&pz,&e,
		    &m,&x,&y,&z,&t);
  }

}









