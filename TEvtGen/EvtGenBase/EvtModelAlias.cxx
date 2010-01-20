//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2002      Caltech, LLNL
//
// Module: EvtGen/EvtModelAlias.cc
//
// Description:Class to keep track of model aliases 
//             read in from the decay table
//
// Modification history:
//
//    Lange     January 19, 2002         Module created
//
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenBase/EvtModelAlias.hh"

EvtModelAlias::EvtModelAlias(std::string alias, std::string model, std::vector<std::string> args):

  _aliasName(alias)
  ,_model(model)
  ,_modelArgs(args)


{
}

EvtModelAlias::EvtModelAlias(const EvtModelAlias &copyMe) :

  _aliasName(copyMe._aliasName)
  ,_model(copyMe._model)
  ,_modelArgs(copyMe._modelArgs)

{

}

EvtModelAlias EvtModelAlias::operator=(const EvtModelAlias &copyMe) {

  _aliasName=copyMe._aliasName;
  _model=copyMe._model;
  _modelArgs = copyMe._modelArgs;

  return *this;
}

std::vector<std::string> EvtModelAlias::getArgList() {
  
  return _modelArgs;
}
