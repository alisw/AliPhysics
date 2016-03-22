#ifdef EVTGEN_PYTHIA
//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2011      University of Warwick, UK
//
// Module: EvtPythiaEngine
//
// Description: Interface to the Pytha 8 external generator
//
// Modification history:
//
//    John Back       April 2011            Module created
//
//------------------------------------------------------------------------

#include "EvtGenExternal/EvtPythiaEngine.hh"

#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtDecayTable.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenBase/EvtParticleFactory.hh"
#include "EvtGenBase/EvtReport.hh"

#include "EvtGenBase/EvtExtGeneratorCommandsTable.hh"
#include "EvtGenExternal/EvtPythia6CommandConverter.hh"

#include "Pythia8/Event.h"


#include <iostream>
#include <sstream>

using std::endl;

EvtPythiaEngine::EvtPythiaEngine(std::string xmlDir, bool convertPhysCodes,
				 bool useEvtGenRandom) {

  // Create two Pythia generators. One will be for generic
  // Pythia decays in the decay.dec file. The other one will be to 
  // only decay aliased particles, which are in general "signal" 
  // decays different from those in the decay.dec file.
  // Even though it is not possible to have two different particle
  // versions in one Pythia generator, we can use two generators to 
  // get the required behaviour.

  report(Severity::Info,"EvtGen")<<"Creating generic Pythia generator"<<endl;
  _genericPythiaGen = new Pythia8::Pythia(xmlDir);
  _genericPartData = _genericPythiaGen->particleData;

  report(Severity::Info,"EvtGen")<<"Creating alias Pythia generator"<<endl;
  _aliasPythiaGen  = new Pythia8::Pythia(xmlDir);
  _aliasPartData = _aliasPythiaGen->particleData;

  _thePythiaGenerator = 0;
  _daugPDGVector.clear(); _daugP4Vector.clear();

  _convertPhysCodes = convertPhysCodes;

  // Specify if we are going to use the random number generator (engine)
  // from EvtGen for Pythia 8.
  _useEvtGenRandom = useEvtGenRandom;

  _evtgenRandom = new EvtPythiaRandom();

  _initialised = false;

}

EvtPythiaEngine::~EvtPythiaEngine() {

  delete _genericPythiaGen; _genericPythiaGen = 0;
  delete _aliasPythiaGen; _aliasPythiaGen = 0;

  delete _evtgenRandom; _evtgenRandom = 0;

  _thePythiaGenerator = 0;

  this->clearDaughterVectors();
  this->clearPythiaModeMap();
 
}

void EvtPythiaEngine::clearDaughterVectors() {
  _daugPDGVector.clear();
  _daugP4Vector.clear();
}

void EvtPythiaEngine::clearPythiaModeMap() {

  PythiaModeMap::iterator iter;
  for (iter = _pythiaModeMap.begin(); iter != _pythiaModeMap.end(); ++iter) {

    std::vector<int> modeVector = iter->second;
    modeVector.clear();

  }

  _pythiaModeMap.clear();

}

void EvtPythiaEngine::initialise() {

  if (_initialised) {return;}

  this->clearPythiaModeMap();

  this->updateParticleLists();

  // Hadron-level processes only (hadronized, string fragmentation and secondary decays).
  // We do not want to generate the full pp or e+e- event structure etc..
  _genericPythiaGen->readString("ProcessLevel:all = off");
  _aliasPythiaGen->readString("ProcessLevel:all = off");

  // Apply any other physics (or special particle) requirements/cuts etc..
  this->updatePhysicsParameters();

  // Set the random number generator
  if (_useEvtGenRandom == true) {

    _genericPythiaGen->setRndmEnginePtr(_evtgenRandom);
    _aliasPythiaGen->setRndmEnginePtr(_evtgenRandom);

  }

  _genericPythiaGen->init();
  _aliasPythiaGen->init();

  _initialised = true;

}

bool EvtPythiaEngine::doDecay(EvtParticle* theParticle) {

  // Store the mother particle within a Pythia8 Event object.
  // Then do the hadron level decays.
  // The EvtParticle must be a colour singlet (meson/baryon/lepton), i.e. not a gluon or quark

  // We delete any daughters the particle may have, since we are asking Pythia
  // to generate the decay anew. Also note that _any_ Pythia decay allowed for the particle
  // will be generated and not the specific Pythia decay mode that EvtGen has already
  // specified. This is necessary since we only want to initialise the Pythia decay table
  // once; all Pythia branching fractions for a given mother particle are renormalised to sum to 1.0.
  // In EvtGen decay.dec files, it may be the case that Pythia decays are only used
  // for some of the particle decays (i.e. Pythia BF sum < 1.0). As we loop over many events,
  // the total frequency for each Pythia decay mode will normalise correctly to what
  // we wanted via the specifications made to the decay.dec file, even though event-by-event
  // the EvtGen decay channel and the Pythia decay channel may be different.

  if (_initialised == false) {this->initialise();}
  
  if (theParticle == 0) {
    report(Severity::Info,"EvtGen")<<"Error in EvtPythiaEngine::doDecay. The mother particle is null. Not doing any Pythia decay."<<endl;
    return false;
  }

  // Delete EvtParticle daughters if they already exist
  if (theParticle->getNDaug() != 0) {
    bool keepChannel(false);
    theParticle->deleteDaughters(keepChannel);
  }

  EvtId particleId = theParticle->getId();
  int isAlias = particleId.isAlias();

  // Choose the generator depending if we have an aliased (parent) particle or not
  _thePythiaGenerator = _genericPythiaGen;
  if (isAlias == 1) {_thePythiaGenerator = _aliasPythiaGen;}

  //_thePythiaGenerator->settings.listChanged();

  // Need to use the reference to the Pythia8::Event object,
  // otherwise it will just return a new empty, default event object.
  Pythia8::Event& theEvent = _thePythiaGenerator->event;
  theEvent.reset();

  // Initialise the event to be the particle rest frame
  int PDGCode = EvtPDL::getStdHep(particleId);

  int status(1);
  int colour(0), anticolour(0);
  double px(0.0), py(0.0), pz(0.0);
  double m0 = theParticle->mass();
  double E = m0;

  theEvent.append(PDGCode, status, colour, anticolour, px, py, pz, E, m0);

  // Generate the Pythia event
  int iTrial(0);
  bool generatedEvent(false);
  for (iTrial = 0; iTrial < 10; iTrial++) {
    
    generatedEvent = _thePythiaGenerator->next();
    if (generatedEvent) {break;}
    
  }

  bool success(false);

  if (generatedEvent) {

    // Store the daughters for this particle from the Pythia decay tree.
    // This is a recursive function that will continue looping through
    // all available daughters until the first set of non-quark and non-gluon 
    // particles are encountered in the Pythia Event structure.

    // First, clear up the internal vectors storing the daughter
    // EvtId types and 4-momenta.
    this->clearDaughterVectors();

    // Now store the daughter info. Since this is a recursive function
    // to loop through the full Pythia decay tree, we do not want to create 
    // EvtParticles here but in the next step.
    this->storeDaughterInfo(theParticle, 1);

    // Now create the EvtParticle daughters of the (parent) particle.
    // We need to use the EvtParticle::makeDaughters function
    // owing to the way EvtParticle stores parent-daughter information.
    this->createDaughterEvtParticles(theParticle);

    //theParticle->printTree();
    //theEvent.list(true, true);

    success = true;

  } 

  return success;

}

void EvtPythiaEngine::storeDaughterInfo(EvtParticle* theParticle, int startInt) {
  
  Pythia8::Event& theEvent = _thePythiaGenerator->event;

  std::vector<int> daugList = theEvent.daughterList(startInt);

  std::vector<int>::iterator daugIter;
  for (daugIter = daugList.begin(); daugIter != daugList.end(); ++daugIter) {

    int daugInt = *daugIter;

    // Ask if the daughter is a quark or gluon. If so, recursively search again.
    Pythia8::Particle& daugParticle = theEvent[daugInt];

    if (daugParticle.isQuark() || daugParticle.isGluon()) {

      // Recursively search for correct daughter type
      this->storeDaughterInfo(theParticle, daugInt);

    } else {

      // We have a daughter that is not a quark nor gluon particle.
      // Make sure we are not double counting particles, since several quarks
      // and gluons make one particle.
      // Set the status flag for any "new" particle to say that we have stored it.
      // Use status flag = 1000 (within the user allowed range for Pythia codes).

      // Check that the status flag for the particle is not equal to 1000
      int statusCode = daugParticle.status();
      if (statusCode != 1000) {

	int daugPDGInt = daugParticle.id();

	double px = daugParticle.px();
	double py = daugParticle.py();
	double pz = daugParticle.pz();
	double E = daugParticle.e();
	EvtVector4R daughterP4(E, px, py, pz);

	// Now store the EvtId and 4-momentum in the internal vectors
	_daugPDGVector.push_back(daugPDGInt);
	_daugP4Vector.push_back(daughterP4);

	// Set the status flag for the Pythia particle to let us know
	// that we have already considered it to avoid double counting.
	daugParticle.status(1000);

      } // Status code != 1000

    }

  }

}

void EvtPythiaEngine::createDaughterEvtParticles(EvtParticle* theParent) {

  if (theParent == 0) {
    report(Severity::Info,"EvtGen")<<"Error in EvtPythiaEngine::createDaughterEvtParticles. The parent is null"<<endl;
    return;
  }

  // Get the list of Pythia decay modes defined for this particle id alias.
  // It would be easier to just use the decay channel number that Pythia chose to use 
  // for the particle decay, but this is not accessible from the Pythia interface at present.

  int nDaughters = _daugPDGVector.size();
  std::vector<EvtId> daugAliasIdVect(0);

  EvtId particleId = theParent->getId();
  // Check to see if we have an anti-particle. If we do, charge conjugate the particle id to get the
  // Pythia "alias" we can compare with the defined (particle) Pythia modes.
  int PDGId = EvtPDL::getStdHep(particleId);
  int aliasInt = particleId.getAlias();
  int pythiaAliasInt(aliasInt);

  if (PDGId < 0) {
    // We have an anti-particle.
    EvtId conjPartId = EvtPDL::chargeConj(particleId);
    pythiaAliasInt = conjPartId.getAlias();
  }

  std::vector<int> pythiaModes = _pythiaModeMap[pythiaAliasInt];

  // Loop over all available Pythia decay modes and find the channel that matches
  // the daughter ids. Set each daughter id to also use the alias integer.
  // This will then convert the Pythia generated channel to the EvtGen alias defined one.

  std::vector<int>::iterator modeIter;
  bool gotMode(false);

  for (modeIter = pythiaModes.begin(); modeIter != pythiaModes.end(); ++modeIter) {

    // Stop the loop if we have the right decay mode channel
    if (gotMode) {break;}

    int pythiaModeInt = *modeIter;

    EvtDecayBase* decayModel = EvtDecayTable::getInstance()->findDecayModel(aliasInt, pythiaModeInt);

    if (decayModel != 0) {

      int nModeDaug = decayModel->getNDaug();

      // We need to make sure that the number of daughters match
      if (nDaughters == nModeDaug) {

	int iModeDaug(0);
	for (iModeDaug = 0; iModeDaug < nModeDaug; iModeDaug++) {

	  EvtId daugId = decayModel->getDaug(iModeDaug);
	  int daugPDGId = EvtPDL::getStdHep(daugId);
	  // Pythia has used the right PDG codes for this decay mode, even for conjugate modes
	  int pythiaPDGId = _daugPDGVector[iModeDaug];

	  if (daugPDGId == pythiaPDGId) {
	    daugAliasIdVect.push_back(daugId);
	  }

	} // Loop over EvtGen mode daughters

	int daugAliasSize = daugAliasIdVect.size();
	if (daugAliasSize == nDaughters) {
	  // All daughter Id codes are accounted for. Set the flag to stop the loop.
	  gotMode = true;
	} else {
	  // We do not have the correct daughter ordering. Clear the id vector
	  // and try another mode.
	  daugAliasIdVect.clear();
	}

      } // Same number of daughters
    
    } // decayModel != 0

  } // Loop over available Pythia modes  

  if (gotMode == false) {

    // We did not find a match for the daughter aliases. Just use the normal PDG codes
    // from the Pythia decay result
    int iPyDaug(0);
    for (iPyDaug = 0; iPyDaug < nDaughters; iPyDaug++) {

      int daugPDGCode = _daugPDGVector[iPyDaug];
      EvtId daugPyId = EvtPDL::evtIdFromStdHep(daugPDGCode);
      daugAliasIdVect.push_back(daugPyId);

    }
  }

  // Make the EvtParticle daughters (with correct alias id's). Their 4-momenta are uninitialised.
  theParent->makeDaughters(nDaughters, daugAliasIdVect);

  // Now set the 4-momenta of the daughters.
  int iDaug(0);
  // Can use an iterator here, but we already had to use the vector size...
  for (iDaug = 0; iDaug < nDaughters; iDaug++) {

    EvtParticle* theDaughter = theParent->getDaug(iDaug);

    // Set the correct 4-momentum for each daughter particle.
    if (theDaughter != 0) {
      EvtId theDaugId = daugAliasIdVect[iDaug];
      const EvtVector4R theDaugP4 = _daugP4Vector[iDaug];
      theDaughter->init(theDaugId, theDaugP4);
    }
    
  }

}

void EvtPythiaEngine::updateParticleLists() {

  // Use the EvtGen decay table (decay/user.dec) to update particle entries 
  // for Pythia. Pythia 8 should use the latest PDG codes, so if the evt.pdl
  // file is up to date, just let Pythia 8 find the particle properties
  // knowing the PDG code integer. If we want to use evt.pdl for _all_
  // particle properties, then we need to make sure that this is up to date,
  // and modify the code in this class to read that data and use it...
  // Using the PDG code only also avoids the need to convert EvtGen particle names
  // to Pythia particle names.

  // Loop over all entries in the EvtPDL particle data table.
  // Aliases are added at the end with id numbers equal to the
  // original particle, but with alias integer = lastPDLEntry+1 etc..
  int iPDL;
  int nPDL = EvtPDL::entries();

  // Reset the _addedPDGCodes map that keeps track
  // of any new particles added to the Pythia input data stream
  _addedPDGCodes.clear();

  for (iPDL = 0; iPDL < nPDL; iPDL++) {

    EvtId particleId = EvtPDL::getEntry(iPDL);
    int aliasInt = particleId.getAlias();

    // Check which particles have a Pythia decay defined.
    // Get the list of all possible decays for the particle, using the alias integer.
    // If the particle is not actually an alias, aliasInt = idInt.

    // Should change isJetSet to isPythia eventually.
    bool hasPythiaDecays = EvtDecayTable::getInstance()->hasPythia(aliasInt);

    if (hasPythiaDecays) {

      int isAlias = particleId.isAlias();

      int PDGCode = EvtPDL::getStdHep(particleId);

      // Decide what generator to use depending on ether we have 
      // an alias particle or not
      _thePythiaGenerator = _genericPythiaGen;
      _theParticleData = _genericPartData;
      if (isAlias == 1) {
	_thePythiaGenerator = _aliasPythiaGen;
	_theParticleData = _aliasPartData;
      }

      // Find the Pythia particle name given the standard PDG code integer
      std::string dataName = _theParticleData.name(PDGCode);
      bool alreadyStored(false);
      if (_addedPDGCodes.find(abs(PDGCode)) != _addedPDGCodes.end()) {alreadyStored = true;}

      if (dataName == " " && alreadyStored == false) {

        // Particle and its antiparticle does not exist in the Pythia database.
	// Create a new particle, then create the new decay modes.
	this->createPythiaParticle(particleId, PDGCode);

      } else {
       
	// Particle exists in the Pythia database.
	// Could update mass/lifetime values here. For now just use Pythia defaults.

      }

      // For the particle, create the Pythia decay modes.
      // Update Pythia data tables.
      this->updatePythiaDecayTable(particleId, aliasInt, PDGCode);

    } // Loop over Pythia decays

  } // Loop over EvtPDL entries

  //report(Severity::Info,"EvtGen")<<"Writing out changed generic Pythia decay list"<<endl;
  //_genericPythiaGen->particleData.listChanged();

  //report(Severity::Info,"EvtGen")<<"Writing out changed alias Pythia decay list"<<endl;
  //_aliasPythiaGen->particleData.listChanged();

}

void EvtPythiaEngine::updatePythiaDecayTable(EvtId& particleId, int aliasInt, int PDGCode) {
  
  // Update the particle data table in Pythia.
  // The tables store information about the allowed decay modes
  // whre the PDGId for all particles must be positive; anti-particles are stored
  // with the corresponding particle entry.
  // Since we do not want to implement CP violation here, just use the same branching
  // fractions for particle and anti-particle modes.

  int nModes = EvtDecayTable::getInstance()->getNModes(aliasInt);
  int iMode(0);

  bool firstMode(true);

  // Only process positive PDG codes.
  if (PDGCode < 0) {return;} 

  // Keep track of which decay modes are Pythia decays for each aliasInt
  std::vector<int> pythiaModes(0);

  // Loop over the decay modes for this particle
  for (iMode = 0; iMode < nModes; iMode++) {
      
    EvtDecayBase* decayModel = EvtDecayTable::getInstance()->findDecayModel(aliasInt, iMode);

    if (decayModel != 0) {

      int nDaug = decayModel->getNDaug();

      // If the decay mode has no daughters, then that means that there will be 
      // no entries for any submode re-definitions for Pythia. 
      // This sometimes occurs for any mode using non-standard Pythia 6 codes.
      // Do not refine the decay mode, i.e. accept the Pythia 8 default (if it exists).
      if (nDaug > 0) {

	// Check to see if we have a Pythia decay mode
	std::string modelName = decayModel->getModelName();

	if (modelName == "PYTHIA") {

	  // Keep track which decay mode is a Pythia one. We need this in order to 
	  // reassign alias Id values for particles generated in the decay.
	  pythiaModes.push_back(iMode);

	  std::ostringstream oss;
	  oss.setf(std::ios::scientific);
	  // Write out the absolute value of the PDG code, since
	  // particles and anti-particles occupy the same part of the Pythia table.
	  oss << PDGCode;

	  if (firstMode) {
	    // Create a new channel
	    oss <<":oneChannel = ";
	    firstMode = false;
	  } else {
	    // Add the channel
	    oss <<":addChannel = ";
	  }

	  // Select all channels (particle and anti-particle).
	  // For CP violation, or different BFs for particle and anti-particle, 
	  // use options 2 or 3 (not here).
	  int onMode(1);
	  oss << onMode << " ";

	  double BF = decayModel->getBranchingFraction();
	  oss << BF << " ";
	  
	  // Need to convert the old Pythia physics mode integers with the new ones
	  // To do this, get the model argument and write a conversion method.
	  int modeInt = this->getModeInt(decayModel);
	  oss << modeInt;
	  
	  int iDaug(0);	
	  for (iDaug = 0; iDaug < nDaug; iDaug++) {
	    
	    EvtId daugId = decayModel->getDaug(iDaug);
	    int daugPDG = EvtPDL::getStdHep(daugId);
	    oss << " " << daugPDG;
	    
	  } // Daughter list

	  _thePythiaGenerator->readString(oss.str());
		
	} // is Pythia

      } else {

	report(Severity::Info,"EvtGen")<<"Warning in EvtPythiaEngine. Trying to redefine Pythia table for "
			     <<EvtPDL::name(particleId)<<" for a decay channel that has no daughters."
			     <<" Keeping Pythia default (if available)."<<endl;	  

      }
	
    } else {
	
      report(Severity::Info,"EvtGen")<<"Error in EvtPythiaEngine. decayModel is null for particle "
			   <<EvtPDL::name(particleId)<<" mode number "<<iMode<<endl;
	
    }

  } // Loop over modes

  _pythiaModeMap[aliasInt] = pythiaModes;

  // Now, renormalise the decay branching fractions to sum to 1.0
  std::ostringstream rescaleStr;
  rescaleStr.setf(std::ios::scientific);
  rescaleStr << PDGCode << ":rescaleBR = 1.0";
  
  _thePythiaGenerator->readString(rescaleStr.str());

}

int EvtPythiaEngine::getModeInt(EvtDecayBase* decayModel) {

  int tmpModeInt(0), modeInt(0);

  if (decayModel != 0) {

    int nVars = decayModel->getNArg();
    // Just read the first integer, which specifies the Pythia decay model. 
    // Ignore any other values.
    if (nVars > 0) {
      tmpModeInt = static_cast<int>(decayModel->getArg(0));
    }
  }

  if (_convertPhysCodes) {

    // Extra code to convert the old Pythia decay model integer MDME(ICC,2) to the new one.
    // This should be removed eventually after updating decay.dec files to use
    // the new convention.
    
    if (tmpModeInt == 0) {
      modeInt = 0; // phase-space
    } else if (tmpModeInt == 1) {
      modeInt = 1; // omega or phi -> 3pi
    } else if (tmpModeInt == 2) {
      modeInt = 11; // Dalitz decay
    } else if (tmpModeInt == 3) {
      modeInt = 2; // V -> PS PS
    } else if (tmpModeInt == 4) {
      modeInt = 92; // onium -> ggg or gg gamma
    } else if (tmpModeInt == 11) {
      modeInt = 42; // phase-space of hadrons from available quarks
    } else if (tmpModeInt == 12) {
      modeInt = 42; // phase-space for onia resonances
    } else if (tmpModeInt == 13) {
      modeInt = 43; // phase-space of at least 3 hadrons
    } else if (tmpModeInt == 14) {
      modeInt = 44; // phase-space of at least 4 hadrons
    } else if (tmpModeInt == 15) {
      modeInt = 45; // phase-space of at least 5 hadrons
    } else if (tmpModeInt >= 22 && tmpModeInt <= 30) {
      modeInt = tmpModeInt + 40; // phase space of hadrons with fixed multiplicity (modeInt - 60)
    } else if (tmpModeInt == 31) {
      modeInt = 42; // two or more quarks phase-space; one spectactor quark
    } else if (tmpModeInt == 32) {
      modeInt = 91; // qqbar or gg pair
    } else if (tmpModeInt == 33) {
      modeInt = 0; // triplet q X qbar, where X = gluon or colour singlet (superfluous, since g's are created anyway)
    } else if (tmpModeInt == 41) {
      modeInt = 21; // weak decay phase space, weighting nu_tau spectrum
    } else if (tmpModeInt == 42) {
      modeInt = 22; // weak decay V-A matrix element
    } else if (tmpModeInt == 43) {
      modeInt = 22; // weak decay V-A matrix element, quarks as jets (superfluous)
    } else if (tmpModeInt == 44) {
      modeInt = 22; // weak decay V-A matrix element, parton showers (superfluous)
    } else if (tmpModeInt == 48) {
      modeInt = 23; // weak decay V-A matrix element, at least 3 decay products
    } else if (tmpModeInt == 50) {
      modeInt = 0; // default behaviour
    } else if (tmpModeInt == 51) {
      modeInt = 0; // step threshold (channel switched off when mass daughters > mother mass
    } else if (tmpModeInt == 52 || tmpModeInt == 53) {
      modeInt = 0; // beta-factor threshold
    } else if (tmpModeInt == 84) {
      modeInt = 42; // unknown physics process - just use phase-space
    } else if (tmpModeInt == 101) {
      modeInt = 0; // continuation line
    } else if (tmpModeInt == 102) {
      modeInt = 0; // off mass shell particles.
    } else {
      report(Severity::Info,"EvtGen")<<"Pythia mode integer "<<tmpModeInt
			   <<" is not recognised. Using phase-space model"<<endl;
      modeInt = 0; // Use phase-space for anything else
    }

  } else {

    // No need to convert the physics mode integer code
    modeInt = tmpModeInt;

  }

  return modeInt;

}

void EvtPythiaEngine::createPythiaParticle(EvtId& particleId, int PDGCode) {

  // Use the EvtGen name, PDGId and other variables to define the new Pythia particle.
  EvtId antiPartId = EvtPDL::chargeConj(particleId);

  std::string aliasName = EvtPDL::name(particleId); // If not an alias, aliasName = normal name
  std::string antiName = EvtPDL::name(antiPartId);

  EvtSpinType::spintype spinType = EvtPDL::getSpinType(particleId);
  int spin = EvtSpinType::getSpin2(spinType);

  int charge = EvtPDL::chg3(particleId);

  // Must set the correct colour type manually here, since the evt.pdl file
  // does not store this information. This is required for quarks otherwise
  // Pythia cannot generate the decay properly.
  int PDGId = EvtPDL::getStdHep(particleId);
  int colour(0);
  if (PDGId == 21) {
    colour = 2; // gluons
  } else if (PDGId <= 8 && PDGId > 0) {
    colour = 1; // single quarks
  }

  double m0 = EvtPDL::getMeanMass(particleId);
  double mWidth = EvtPDL::getWidth(particleId);
  double mMin = EvtPDL::getMinMass(particleId);
  double mMax = EvtPDL::getMaxMass(particleId);

  double tau0 = EvtPDL::getctau(particleId);

  std::ostringstream oss;
  oss.setf(std::ios::scientific);
  int absPDGCode = abs(PDGCode);
  oss << absPDGCode << ":new = " << aliasName << " " << antiName << " "
      << spin << " " << charge << " " << colour << " " 
      << m0 << " " << mWidth << " " << mMin << " " << mMax << " "
      << tau0;
  
  // Pass this information to Pythia
  _thePythiaGenerator->readString(oss.str());

  // Also store the absolute value of the PDG entry
  // to keep track of which new particles have been added,
  // which also automatically includes the anti-particle.
  // We need to avoid creating new anti-particles when
  // they already exist when the particle was added.
  _addedPDGCodes[absPDGCode] = 1;

}

void EvtPythiaEngine::updatePhysicsParameters() {

  // Update any more Pythia physics (or special particle) requirements/cuts etc..
  // This should be used if any of the Pythia 6 parameters like JetSetPar MSTJ(i) = x
  // are needed. Such commands will need to be implemented using the new interface
  // pythiaGenerator->readString(cmd); Here cmd is a string telling Pythia 8
  // what physics parameters to change. This will need to be done for the generic and
  // alias generator pointers, as appropriate.

  // Set the multiplicity level for hadronic weak decays
  std::string multiWeakCut("ParticleDecays:multIncreaseWeak = 2.0");
  _genericPythiaGen->readString(multiWeakCut);
  _aliasPythiaGen->readString(multiWeakCut);

  // Set the multiplicity level for all other decays
  std::string multiCut("ParticleDecays:multIncrease = 4.5");
  _genericPythiaGen->readString(multiCut);
  _aliasPythiaGen->readString(multiCut);

  //Now read in any custom configuration entered in the XML
  GeneratorCommands commands = EvtExtGeneratorCommandsTable::getInstance()->getCommands("PYTHIA");
  GeneratorCommands::iterator it = commands.begin();

  for( ; it!=commands.end(); it++) {

    Command command = *it;
    std::vector<std::string> commandStrings;

    if(command["VERSION"] == "PYTHIA6") {
      report(Severity::Info,"EvtGen")<<"Converting Pythia 6 command: "<<command["MODULE"]<<"("<<command["PARAM"]<<")="<<command["VALUE"]<<"..."<<endl;
      commandStrings = convertPythia6Command(command);
    } else if(command["VERSION"] == "PYTHIA8") {
      commandStrings.push_back(command["MODULE"]+":"+command["PARAM"]+" = "+command["VALUE"]);
    } else {
      report(Severity::Error, "EvtGen") << "Pythia command received by EvtPythiaEngine has bad version:"<<endl;
      report(Severity::Error, "EvtGen") << "Received "<<command["VERSION"]<<" but expected PYTHIA6 or PYTHIA8."<<endl;
      report(Severity::Error, "EvtGen") << "The error is likely to be in EvtDecayTable.cpp"<<endl;
      report(Severity::Error, "EvtGen") << "EvtGen will now abort."<<endl;
      ::abort();
    }
    std::string generator = command["GENERATOR"];
    if(generator == "GENERIC" || generator == "Generic" || generator == "generic" ||
       generator == "BOTH" || generator == "Both" || generator == "both") {
      std::vector<std::string>::iterator it2 = commandStrings.begin();
      for( ; it2!=commandStrings.end(); it2++) {
        report(Severity::Info,"EvtGen")<<"Configuring generic Pythia generator: " << (*it2) << endl;
        _genericPythiaGen->readString(*it2);
      }
    }
    if(generator == "ALIAS" || generator == "Alias" || generator == "alias" ||
       generator == "BOTH" || generator == "Both" || generator == "both") {
      std::vector<std::string>::iterator it2 = commandStrings.begin();
      for( ; it2!=commandStrings.end(); it2++) {
        report(Severity::Info,"EvtGen")<<"Configuring alias Pythia generator: " << (*it2) << endl;
        _aliasPythiaGen->readString(*it2);
      }
    }
  }
}

#endif
