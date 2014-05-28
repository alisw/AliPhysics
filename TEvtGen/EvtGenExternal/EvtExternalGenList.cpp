//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2012     University of Warwick, UK
//
// Module: EvtExternalGenFactory
//
// Description: A factory type method to create engines for external physics
// generators like Pythia.
//
// Modification history:
//
//    John Back       Sept 2012           Module created
//
//------------------------------------------------------------------------------
//
#include <TSystem.h>
#include "EvtGenExternal/EvtExternalGenList.hh"

#include "EvtGenExternal/EvtExternalGenFactory.hh"
#include "EvtGenExternal/EvtPHOTOS.hh"
#include "EvtGenExternal/EvtPythia.hh"
#include "EvtGenExternal/EvtTauola.hh"

EvtExternalGenList::EvtExternalGenList(bool convertPythiaCodes, std::string pythiaXmlDir,
				       std::string photonType, bool useEvtGenRandom) {

  // Instantiate the external generator factory
  EvtExternalGenFactory* extFactory = EvtExternalGenFactory::getInstance();

  // Define the external generator "engines" here
  extFactory->definePhotosGenerator(photonType, useEvtGenRandom);

  if (pythiaXmlDir.size() < 1) {
    // If we have no string defined, check the value of the
    // PYTHIA8DATA environment variable which should be set to the 
    // xmldoc Pythia directory
    //char* pythiaDataDir = getenv("PYTHIA8DATA");
    char *pythiaDataDir = gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8175/xmldoc");
    if (pythiaDataDir != 0) {pythiaXmlDir = pythiaDataDir;}
  }

  extFactory->definePythiaGenerator(pythiaXmlDir, convertPythiaCodes,
				    useEvtGenRandom);

  extFactory->defineTauolaGenerator(useEvtGenRandom);  

}

EvtExternalGenList::~EvtExternalGenList() {
}

EvtAbsRadCorr* EvtExternalGenList::getPhotosModel() {

  // Define the Photos model, which uses the EvtPhotosEngine class.
  EvtPHOTOS* photosModel = new EvtPHOTOS();
  return photosModel;

}

std::list<EvtDecayBase*> EvtExternalGenList::getListOfModels() {

  // Create the Pythia and Tauola models, which use their own engine classes.
  EvtPythia* pythiaModel = new EvtPythia();
  EvtTauola* tauolaModel = new EvtTauola();

  std::list<EvtDecayBase*> extraModels;
  extraModels.push_back(pythiaModel);
  extraModels.push_back(tauolaModel);

  // Return the list of models
  return extraModels;

}
