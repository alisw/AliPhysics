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
//    JBack   June 2011             Added HepMC event interface
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGen/EvtGen.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtRandomEngine.hh"
#include "EvtGenBase/EvtSimpleRandomEngine.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenModels/EvtModelReg.hh"
#include "EvtGenBase/EvtStatus.hh"
#include "EvtGenBase/EvtAbsRadCorr.hh"
#include "EvtGenBase/EvtRadCorr.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtHepMCEvent.hh"

#include "EvtGenModels/EvtNoRadCorr.hh"

#include <cstdlib>
#include <fstream>
#include <string>

using std::endl;
using std::fstream;
using std::ifstream;

EvtGen::~EvtGen(){

  //This is a bit ugly, should not do anything
  //in a destructor. This will fail if EvtGen is made a static
  //because then this destructor might be called _after_
  //the destruction of objects that it depends on, e.g., EvtPDL.

  if (getenv("EVTINFO")){
    EvtDecayTable::getInstance()->printSummary();
  }

}

EvtGen::EvtGen(const char* const decayName,
	       const char* const pdtTableName,
	       EvtRandomEngine* randomEngine,
	       EvtAbsRadCorr* isrEngine,
	       const std::list<EvtDecayBase*>* extraModels,
	       int mixingType, bool useXml){


  report(Severity::Info,"EvtGen") << "Initializing EvtGen"<<endl;

  if (randomEngine==0){
    static EvtSimpleRandomEngine defaultRandomEngine;
    EvtRandom::setRandomEngine(&defaultRandomEngine);
    report(Severity::Info,"EvtGen") <<"No random engine given in "
			  <<"EvtGen::EvtGen constructor, "
			  <<"will use default EvtSimpleRandomEngine."<<endl;
  }
  else{
    EvtRandom::setRandomEngine(randomEngine);    
  }

  report(Severity::Info,"EvtGen") << "Storing known decay models"<<endl;
  EvtModelReg dummy(extraModels);

  report(Severity::Info,"EvtGen") << "Main decay file name  :"<<decayName<<endl;
  report(Severity::Info,"EvtGen") << "PDT table file name   :"<<pdtTableName<<endl;
  
  _pdl.readPDT(pdtTableName);

  if(useXml) {
    EvtDecayTable::getInstance()->readXMLDecayFile(decayName,false);
  } else {
    EvtDecayTable::getInstance()->readDecayFile(decayName,false);
  }

  _mixingType = mixingType;
  report(Severity::Info,"EvtGen") << "Mixing type integer set to "<<_mixingType<<endl;
  EvtCPUtil::getInstance()->setMixingType(_mixingType);

  // Set the radiative correction engine

  if (isrEngine != 0) {

    EvtRadCorr::setRadCorrEngine(isrEngine);

  } else {

    // Owing to the pure abstract interface, we still need to define a concrete 
    // implementation of a radiative correction engine. Use one which does nothing.
    EvtAbsRadCorr* noRadCorr = new EvtNoRadCorr();
    EvtRadCorr::setRadCorrEngine(noRadCorr);

  }

  report(Severity::Info,"EvtGen") << "Done initializing EvtGen"<<endl;

}


void EvtGen::readUDecay(const char* const uDecayName, bool useXml){

  ifstream indec;

  if ( uDecayName[0] == 0) {
    report(Severity::Info,"EvtGen") << "Is not reading a user decay file!"<<endl;
  }
  else{  
    indec.open(uDecayName);
    if (indec) {
      if(useXml) {
        EvtDecayTable::getInstance()->readXMLDecayFile(uDecayName,true);
      } else {
        EvtDecayTable::getInstance()->readDecayFile(uDecayName,true);
      }
    }    
    else{
      
      report(Severity::Info,"EvtGen") << "Can not find UDECAY file '"
			    <<uDecayName<<"'.  Stopping"<<endl;
      ::abort();
    }
  }
  
}

EvtHepMCEvent* EvtGen::generateDecay(int PDGId, EvtVector4R refFrameP4,
				     EvtVector4R translation,
				     EvtSpinDensity* spinDensity) {

  EvtParticle* theParticle(0);

  if (spinDensity == 0 ){
    theParticle = EvtParticleFactory::particleFactory(EvtPDL::evtIdFromStdHep(PDGId),
						      refFrameP4);
  } else {
    theParticle = EvtParticleFactory::particleFactory(EvtPDL::evtIdFromStdHep(PDGId),
						      refFrameP4, *spinDensity);
  }

  generateDecay(theParticle);
  EvtHepMCEvent* hepMCEvent = new EvtHepMCEvent();
  hepMCEvent->constructEvent(theParticle, translation);

  theParticle->deleteTree();

  return hepMCEvent;

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
      report(Severity::Error,"EvtGen") << "Your event has been rejected 10000 times!"<<endl;
      report(Severity::Error,"EvtGen") << "Will now abort."<<endl;
      ::abort();
      times=0;
    }
  } while (times);

}
