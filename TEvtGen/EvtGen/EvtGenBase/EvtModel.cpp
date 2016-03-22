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
// Module: EvtModel.cc
//
// Description:
//
// Modification history:
//
//    RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtModel.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtParticleDecayList.hh"
#include "EvtGenBase/EvtParser.hh"
#include "EvtGenBase/EvtReport.hh"
#include <string>
using std::fstream;

EvtModel* EvtModel::_instance=0;

EvtModel::EvtModel() {

}

EvtDecayBase* EvtModel::getFcn(std::string model_name){

  EvtDecayBase *model=0;
  if ( _modelNameHash.find(model_name)!=_modelNameHash.end() ) {
    model=_modelNameHash[model_name];
  }

  if (model==0){
    report(Severity::Error,"EvtGen") << "Did not find the right model:"
			   <<model_name.c_str()<<"\n";
    return 0;
  }

  return model->clone();

}


void EvtModel::registerModel(EvtDecayBase* prototype){

  std::string modelName= prototype->getName();

  _modelNameHash[modelName]=prototype;

  std::string commandName=prototype->commandName();
  
  if (commandName!=""){

    _commandNameHash[commandName]=prototype;

  }

}

int EvtModel::isModel(std::string model_name){

  if ( _modelNameHash.find(model_name)!=_modelNameHash.end() ) {
    return 1;
  }
  return 0;
}


int EvtModel::isCommand(std::string cmd){

  if ( _commandNameHash.find(cmd)!=_commandNameHash.end() ) {
    return 1;
  }
  return 0;
}

void EvtModel::storeCommand(std::string cmd,std::string cnfgstr){

  EvtDecayBase *model=0;
  if ( _commandNameHash.find(cmd)!=_commandNameHash.end() ) {
    model=_commandNameHash[cmd];
  }

  assert(model!=0);

  model->command(cnfgstr);

}




