// DatabasePDG stores and handles PDG information
// The PDG particle definitions and decay channels are read
// in the begining from ASCII files
// PDG definitions loaded can be selected according to their
// mass and decay width

/*
  Copyright   : The FASTMC and SPHMC Collaboration
  Author      : Ionut Cristian Arsene 
  Affiliation : Oslo University, Norway & Institute for Space Sciences, Bucharest, Romania
  e-mail      : i.c.arsene@fys.uio.no
  Date        : 2007/05/30

  This class is using the particle and decay lists provided by the 
  THERMINATOR (Computer Physics Communications 174 669 (2006)) and
  SHARE (Computer Physics Communications 167 229 (2005)) collaborations.
*/


#ifndef DATABASEPDG_H
#include "DatabasePDG.h"
#endif

#include <iostream>
#include <cstring>
#include <fstream>
using std::cout;
using std::endl;
using std::strcmp;
using std::strcpy;
using std::ifstream;
using std::ios;

DatabasePDG::DatabasePDG():
  fNParticles(0),
  fUseCharmParticles(kTRUE),
  fMinimumWidth(0.0),
  fMaximumWidth(10.),
  fMinimumMass(0.0001),
  fMaximumMass(10.)
{
  // Default constructor, initialize members, set input files  
  
  strncpy(fParticleFilename, "particles.data", 256);
  strncpy(fDecayFilename, "tabledecay.txt", 256);
  for(Int_t i=0; i<kMaxParticles; i++) {
    fParticles[i] = new ParticlePDG();
    fStatus[i] = kFALSE;
  }
}

DatabasePDG::~DatabasePDG() {
  for(Int_t i=0; i<kMaxParticles; i++)
    if(fParticles[i]) 
      delete fParticles[i];
}

void DatabasePDG::SetParticleFilename(Char_t *filename) {
  if(strlen(filename)<255)
    strncpy(fParticleFilename, filename, 256);
}

void DatabasePDG::SetDecayFilename(Char_t *filename) {
  if(strlen(filename)<255)
    strncpy(fDecayFilename, filename, 256);
}

Bool_t DatabasePDG::LoadData() {
  return (LoadParticles() && LoadDecays());
}

Bool_t DatabasePDG::LoadParticles() {
  // Read particle definitions from the ascii file  

  ifstream particleFile;
  particleFile.open(fParticleFilename);
  if(!particleFile) {
    cout << "ERROR in DatabasePDG::LoadParticles() : The ASCII file containing the PDG particle list (\""
         << fParticleFilename << "\") was not found !! Exiting..." << endl;
    return kFALSE;
  }
  
  Char_t name[9];
  Double_t mass, width, spin, isospin, isospinZ, q, s, aq, as, c, ac;
  Int_t pdg;
  Int_t goodStatusParticles = 0;

  cout << "Info in DatabasePDG::LoadParticles() : Start loading particles with the following criteria:" << endl
       << "       Use particles containing charm quarks (1-yes;0-no) : " << fUseCharmParticles << endl
       << "       Mass range                                         : (" << fMinimumMass << "; " << fMaximumMass << ")" << endl
       << "       Width range                                        : (" << fMinimumWidth << "; " << fMaximumWidth << ")" << endl;
  
  particleFile.exceptions(ios::failbit);
  while(!particleFile.eof()) {
    try {
      particleFile >> name >> mass >> width >> spin >> isospin >> isospinZ >> q >> s >> aq >> as >> c >> ac >> pdg;
    }
    catch (ios::failure const &problem) {
      cout << problem.what() << endl;
      break;
    }
        
    fParticles[fNParticles]->SetName(name);
    fParticles[fNParticles]->SetPDG(pdg);
    fParticles[fNParticles]->SetMass(mass);
    fParticles[fNParticles]->SetWidth(width);
    fParticles[fNParticles]->SetSpin(spin);
    fParticles[fNParticles]->SetIsospin(isospin);
    fParticles[fNParticles]->SetIsospinZ(isospinZ);
    fParticles[fNParticles]->SetLightQNumber(q);
    fParticles[fNParticles]->SetStrangeQNumber(s);
    fParticles[fNParticles]->SetLightAQNumber(aq);
    fParticles[fNParticles]->SetStrangeAQNumber(as);
    fParticles[fNParticles]->SetCharmQNumber(c);
    fParticles[fNParticles]->SetCharmAQNumber(ac);
    
    fStatus[fNParticles] = kTRUE;
    // check if we want charmed particles
    if(!fUseCharmParticles && (c>0 || ac>0)) {
      fStatus[fNParticles] = kFALSE;
    }
    // check that the particle mass is inside accepted limits
    if(!(fMinimumMass<=mass && mass<=fMaximumMass)) {
      fStatus[fNParticles] = kFALSE;
    }
    // check that the particle width is inside accepted limits
    if(!(fMinimumWidth<=width && width<=fMaximumWidth)) {
      fStatus[fNParticles] = kFALSE;
    }
    if(fStatus[fNParticles]) goodStatusParticles++;
    fNParticles++;
  }
  particleFile.close();
  if(fNParticles==0) {
    cout << "Warning in DatabasePDG::LoadParticles(): No particles were found in the file specified!!" << endl;
    return kFALSE;
  }
  SortParticles();
  cout << "Info in DatabasePDG::LoadParticles(): Particle definitions found = " << fNParticles << endl
       << "                                      Good status particles      = " << goodStatusParticles << endl;
  return kTRUE;
}

Bool_t DatabasePDG::LoadDecays() {
  // Read the decay channel definitions from the ascii file

  ifstream decayFile;
  decayFile.open(fDecayFilename);
  if(!decayFile) {
    cout << "ERROR in DatabasePDG::LoadDecays() : The ASCII file containing the decays list (\""
         << fDecayFilename << "\") was not found !! Exiting..." << endl;
    return kFALSE;
  }
  
  Int_t motherPdg, daughterPdg[3];
  Double_t branching;
  
  decayFile.exceptions(ios::failbit);
  while(!decayFile.eof()) {
    motherPdg = 0;
    for(Int_t i=0; i<3; i++) daughterPdg[i] = 0;
    branching = -1.0;
    try {
      decayFile >> motherPdg;
      for(Int_t i=0; i<3; i++) 
        decayFile >> daughterPdg[i];
      decayFile >> branching;
    }
    catch (ios::failure const &problem) {
      cout << problem.what() << endl;
      break;
    }
    if((motherPdg!=0) && (daughterPdg[0]!=0) && (branching>=0)) {
      Int_t nDaughters = 0;
      for(Int_t i=0; i<3; i++)
        if(daughterPdg[i]!=0)
          nDaughters++;
      ParticlePDG* particle = GetPDGParticle(motherPdg);
      DecayChannel decay(motherPdg, branching, nDaughters, daughterPdg);
      particle->AddChannel(&decay);
    }
  }
  decayFile.close();
  Int_t nDecayChannels = 0;
  for(Int_t i=0; i<fNParticles; i++) {
    nDecayChannels += fParticles[i]->GetNDecayChannels();
  }
  cout << "Info in DatabasePDG::LoadDecays(): Number of decays found in the database is " << nDecayChannels << endl;
  return kTRUE;
}

ParticlePDG* DatabasePDG::GetPDGParticleByIndex(Int_t index) const {
  // Return a PDG particle definition based on its index in the particle list 

  if(index<0 || index>fNParticles) {
    cout << "Warning in DatabasePDG::GetPDGParticleByIndex(Int_t): Particle index is negative or too big !!" << endl
         << " It must be inside this range: [0, " << fNParticles-1 << "]" << endl
         << " Returning null pointer!!" << endl;
    return 0x0;
  }
  return fParticles[index];
}

Bool_t DatabasePDG::GetPDGParticleStatusByIndex(Int_t index) const {
  // Return the status of a PDG particle definition based on its index in the particle list
  // The status is kTRUE when a particle passed the mass, width and charm criteria
  // and kFALSE otherwise

  if(index<0 || index>fNParticles) {
    cout << "Warning in DatabasePDG::GetPDGParticleStatusByIndex(Int_t): Particle index is negative or too big !!" << endl
         << " It must be inside this range: [0, " << fNParticles-1 << "]" << endl
         << " Returning null pointer!!" << endl;
    return kFALSE;
  }
  return fStatus[index];
}

ParticlePDG* DatabasePDG::GetPDGParticle(Int_t pdg) const {
  // Return a PDG particle definition based on the PDG code (PYTHIA convention) 
  // If more than 1 definition with the asked PDG code is found then a warning is issued 

  Int_t nFindings = 0;
  Int_t firstTimeIndex = 0;
  for(Int_t i=0; i<fNParticles; i++) {
    if(pdg == fParticles[i]->GetPDG()) {
      if(nFindings == 0) firstTimeIndex = i;
      nFindings++;
    }
  }
  if(nFindings == 1) return fParticles[firstTimeIndex];
  if(nFindings == 0) {
    //cout << "Warning in DatabasePDG::GetPDGParticle(Int_t): The particle required with PDG = " << pdg
    //     << " was not found in the database!!" << endl;
    return 0x0;
  }
  if(nFindings >= 2) {
    cout << "Warning in DatabasePDG::GetPDGParticle(Int_t): The particle required with PDG = " << pdg
         << " was found with " << nFindings << " entries in the database. Check it out !!" << endl
	 << "Returning the first instance found" << endl;
    return fParticles[firstTimeIndex];
  }
  return 0x0;
}

Bool_t DatabasePDG::GetPDGParticleStatus(Int_t pdg) const {
  // Return the status of a PDG particle definition based on its PDG code (PYTHIA convention)
  // The status is kTRUE when a particle passed the mass, width and charm criteria
  // and kFALSE otherwise

  Int_t nFindings = 0;
  Int_t firstTimeIndex = 0;
  for(Int_t i=0; i<fNParticles; i++) {
    if(pdg == fParticles[i]->GetPDG()) {
      if(nFindings == 0) firstTimeIndex = i;
      nFindings++;
    }
  }
  if(nFindings == 1) return fStatus[firstTimeIndex];
  if(nFindings == 0) {
    //cout << "Warning in DatabasePDG::GetPDGParticle(Int_t): The particle required with PDG = " << pdg
    //     << " was not found in the database!!" << endl;
    return kFALSE;
  }
  if(nFindings >= 2) {
    cout << "Warning in DatabasePDG::GetPDGParticleStatus(Int_t): The particle status required for PDG = " << pdg
         << " was found with " << nFindings << " entries in the database. Check it out !!" << endl
	 << "Returning the status of first instance found" << endl;
    return fStatus[firstTimeIndex];
  }
  return kFALSE;
}

ParticlePDG* DatabasePDG::GetPDGParticle(Char_t* name) const {
  // Return a PDG particle definition based on its name
  // If more than 1 definition with the asked PDG code is found then a warning is issued

  Int_t nFindings = 0;
  Int_t firstTimeIndex = 0;
  for(Int_t i=0; i<fNParticles; i++) {
    if(!strcmp(name, fParticles[i]->GetName())) {
      if(nFindings == 0) firstTimeIndex = i;
      nFindings++;
    }
  }
  if(nFindings == 1) return fParticles[firstTimeIndex];
  if(nFindings == 0) {
    //cout << "Warning in DatabasePDG::GetPDGParticle(Char_t*): The particle required with name \"" << name
    //     << "\" was not found in the database!!" << endl;
    return 0x0;
  }
  if(nFindings >= 2) {
    cout << "Warning in DatabasePDG::GetPDGParticle(Char_t*): The particle required with name \"" << name
         << "\" was found with " << nFindings << " entries in the database. Check it out !!" << endl
	 << "Returning the first instance found" << endl;
    return fParticles[firstTimeIndex];
  }
  return 0x0;
}

Bool_t DatabasePDG::GetPDGParticleStatus(Char_t* name) const {
  // Return the status of a PDG particle definition based on its name
  // The status is kTRUE when a particle passed the mass, width and charm criteria
  // and kFALSE otherwise

  Int_t nFindings = 0;
  Int_t firstTimeIndex = 0;
  for(Int_t i=0; i<fNParticles; i++) {
    if(!strcmp(name, fParticles[i]->GetName())) {
      if(nFindings == 0) firstTimeIndex = i;
      nFindings++;
    }
  }
  if(nFindings == 1) return fStatus[firstTimeIndex];
  if(nFindings == 0) {
    //cout << "Warning in DatabasePDG::GetPDGParticle(Char_t*): The particle required with name \"" << name
    //     << "\" was not found in the database!!" << endl;
    return kFALSE;
  }
  if(nFindings >= 2) {
    cout << "Warning in DatabasePDG::GetPDGParticleStatus(Char_t*): The particle status required for name \"" << name
         << "\" was found with " << nFindings << " entries in the database. Check it out !!" << endl
	 << "Returning the first instance found" << endl;
    return fStatus[firstTimeIndex];
  }
  return kFALSE;
}

void DatabasePDG::DumpData(Bool_t dumpAll) const {
  // Printout all the information in the PDG database

  cout << "***********************************************************************************************************" << endl;
  cout << "Dumping all the information contained in the database..." << endl;
  Int_t nDecays = 0;
  Int_t nGoodStatusDecays = 0;
  Int_t nGoodStatusParticles = 0;
  for(Int_t currPart=0; currPart<fNParticles; currPart++) {
    nGoodStatusParticles += (fStatus[currPart] ? 1:0);
    nGoodStatusDecays += (fStatus[currPart] ? fParticles[currPart]->GetNDecayChannels() : 0);
    nDecays += fParticles[currPart]->GetNDecayChannels();
    if(!(dumpAll || (!dumpAll && fStatus[currPart]))) continue;
    cout << "###### Particle: " << fParticles[currPart]->GetName() << " with PDG code " << fParticles[currPart]->GetPDG() << endl;
    cout << "   status          = " << fStatus[currPart] << endl;
    cout << "   mass            = " << fParticles[currPart]->GetMass() << " GeV" << endl;
    cout << "   width           = " << fParticles[currPart]->GetWidth() << " GeV" << endl;
    cout << "   2*spin          = " << Int_t(2.*fParticles[currPart]->GetSpin()) << endl;
    cout << "   2*isospin       = " << Int_t(2.*fParticles[currPart]->GetIsospin()) << endl;
    cout << "   2*isospin3      = " << Int_t(2.*fParticles[currPart]->GetIsospinZ()) << endl;
    cout << "   u,d quarks      = " << Int_t(fParticles[currPart]->GetLightQNumber()) << endl;
    cout << "   s quarks        = " << Int_t(fParticles[currPart]->GetStrangeQNumber()) << endl;
    cout << "   c quarks        = " << Int_t(fParticles[currPart]->GetCharmQNumber()) << endl;
    cout << "   anti u,d quarks = " << Int_t(fParticles[currPart]->GetLightAQNumber()) << endl;
    cout << "   anti s quarks   = " << Int_t(fParticles[currPart]->GetStrangeAQNumber()) << endl;
    cout << "   anti c quarks   = " << Int_t(fParticles[currPart]->GetCharmAQNumber()) << endl;
    cout << "   baryon number   = " << Int_t(fParticles[currPart]->GetBaryonNumber()) << endl;
    cout << "   strangeness     = " << Int_t(fParticles[currPart]->GetStrangeness()) << endl;
    cout << "   charmness       = " << Int_t(fParticles[currPart]->GetCharmness()) << endl;
    cout << "   electric charge = " << Int_t(fParticles[currPart]->GetElectricCharge()) << endl;
    cout << "   full branching  = " << fParticles[currPart]->GetFullBranching() << endl;
    cout << "   decay modes     = " << fParticles[currPart]->GetNDecayChannels() << endl;
    for(Int_t currChannel=0; currChannel<fParticles[currPart]->GetNDecayChannels(); currChannel++) {
      cout << "   channel " << currChannel+1 << " with branching " << fParticles[currPart]->GetDecayChannel(currChannel)->GetBranching() << endl;
      cout << "   daughters PDG codes: ";
      Double_t daughtersMass = 0.0;
      for(Int_t currDaughter=0; currDaughter<fParticles[currPart]->GetDecayChannel(currChannel)->GetNDaughters(); currDaughter++) {
        cout << fParticles[currPart]->GetDecayChannel(currChannel)->GetDaughterPDG(currDaughter) << "\t";
	ParticlePDG *daughter = GetPDGParticle(fParticles[currPart]->GetDecayChannel(currChannel)->GetDaughterPDG(currDaughter));
        daughtersMass += daughter->GetMass();
      }
      cout << endl;
      cout << "   daughters sum mass = " << daughtersMass << endl;
    }
  }
  if(dumpAll) {
    cout << "Finished dumping information for " << fNParticles << " particles with " << nDecays << " decay channels in total." << endl;
    cout << "*************************************************************************************************************" << endl;
  }
  else {
    cout << "Finished dumping information for " << nGoodStatusParticles << "(" << fNParticles << ")" 
	 << " particles with " << nGoodStatusDecays << "(" << nDecays << ")" << " decay channels in total." << endl;
    cout << "*************************************************************************************************************" << endl;
  }
}

Int_t DatabasePDG::CheckImpossibleDecays(Bool_t dump) const {
  // Check the database for impossible decays

  Int_t nImpossibleDecays = 0;
  for(Int_t currPart=0; currPart<fNParticles; currPart++) {
    if(!fStatus[currPart]) continue;
    Int_t allChannels = fParticles[currPart]->GetNDecayChannels();
    Int_t allowedChannels = GetNAllowedChannels(fParticles[currPart], fParticles[currPart]->GetMass());
    if(dump) {
      cout << "Particle " << fParticles[currPart]->GetPDG() << " has " << allChannels << " decay channels specified in the database" << endl;
      cout << " Allowed channels assuming table mass = " << allowedChannels << endl;
    }
    if(dump && allChannels>0 && allowedChannels == 0) {
      cout << "**********************************************************************" << endl;
      cout << "       All channels for this particles are not allowed" << endl;
      cout << "**********************************************************************" << endl;
    }
    if(dump && fParticles[currPart]->GetWidth() > 0. && allChannels == 0) {
      cout << "**********************************************************************" << endl;
      cout << "    Particle has finite width but no decay channels specified" << endl;
      cout << "**********************************************************************" << endl;
    }
    for(Int_t currChannel=0; currChannel<fParticles[currPart]->GetNDecayChannels(); currChannel++) {
      Double_t motherMass = fParticles[currPart]->GetMass();
      Double_t daughtersSumMass = 0.;
      for(Int_t currDaughter=0; currDaughter<fParticles[currPart]->GetDecayChannel(currChannel)->GetNDaughters(); currDaughter++) {
        ParticlePDG *daughter = GetPDGParticle(fParticles[currPart]->GetDecayChannel(currChannel)->GetDaughterPDG(currDaughter));
        daughtersSumMass += daughter->GetMass();
      }
      if(daughtersSumMass >= motherMass) {
        nImpossibleDecays++;
        if(dump) {
          cout << "Imposible decay for particle " << fParticles[currPart]->GetPDG() << endl;
          cout << "  Channel: " << fParticles[currPart]->GetPDG() << " --> ";
          for(Int_t currDaughter=0; currDaughter<fParticles[currPart]->GetDecayChannel(currChannel)->GetNDaughters(); currDaughter++) {
            ParticlePDG *daughter = GetPDGParticle(fParticles[currPart]->GetDecayChannel(currChannel)->GetDaughterPDG(currDaughter));
            cout << daughter->GetPDG() << " ";
          }
          cout << endl;
          cout << "  Mother particle mass = " << motherMass << endl;
          cout << "  Daughters sum mass   = " << daughtersSumMass << endl;
        }
      }
    }
  }
  return nImpossibleDecays;
}

void DatabasePDG::SetUseCharmParticles(Bool_t flag) {
  // Switch on/off the use of charmed particles

  if(fNParticles>0) {
    fUseCharmParticles = flag;
    for(Int_t i=0; i<fNParticles; i++) {
      if(fParticles[i]->GetCharmQNumber()>0 || fParticles[i]->GetCharmAQNumber())		  
	fStatus[i] = flag;
    }
    SortParticles();
    return;
  }
  else
    fUseCharmParticles = flag;
  return;
}

void DatabasePDG::SetMinimumWidth(Double_t value) {
  // Set the minimum decay width for the particle definitions to be generated in the soft fireball
  if(fNParticles>0) {
    fMinimumWidth = value;
    for(Int_t i=0; i<fNParticles; i++) {
      if(fParticles[i]->GetWidth() < fMinimumWidth)		  
	fStatus[i] = kFALSE;
    }
    SortParticles();
    return;
  }
  else
    fMinimumWidth = value;
  return;
}

void DatabasePDG::SetMaximumWidth(Double_t value) {
  // Set the maximum decay width for the particle definitions to be generated in the soft fireball

  if(fNParticles>0) {
    fMaximumWidth = value;
    for(Int_t i=0; i<fNParticles; i++) {
      if(fParticles[i]->GetWidth() > fMaximumWidth)		  
	fStatus[i] = kFALSE;
    }
    SortParticles();
    return;
  }
  else
    fMaximumWidth = value;
  return;
}

void DatabasePDG::SetWidthRange(Double_t min, Double_t max) {
  // Set the decay width range for the particle definitions to be generated in the soft fireball

  if(fNParticles>0) {
    fMinimumWidth = min;
    fMaximumWidth = max;
    for(Int_t i=0; i<fNParticles; i++) {
      if((fParticles[i]->GetWidth()<fMinimumWidth) || (fParticles[i]->GetWidth()>fMaximumWidth))  
	fStatus[i] = kFALSE;
    }
    SortParticles();
    return;
  }
  else {
    fMinimumWidth = min;
    fMaximumWidth = max;
  }
  return;
}

void DatabasePDG::SetMinimumMass(Double_t value) {
  // Set the minimum mass for the particle definitions to be generated in the soft fireball

  if(fNParticles>0) {
    fMinimumMass = value;
    for(Int_t i=0; i<fNParticles; i++) {
      if(fParticles[i]->GetMass() < fMinimumMass)		  
	fStatus[i] = kFALSE;
    }
    SortParticles();
    return;
  }
  else
    fMinimumMass = value;
  return;
}

void DatabasePDG::SetMaximumMass(Double_t value) {
  // Set the maximum mass for the particle definitions to be generated in the soft fireball

  if(fNParticles>0) {
    fMaximumMass = value;
    for(Int_t i=0; i<fNParticles; i++) {
      if(fParticles[i]->GetMass() > fMaximumMass)		  
	fStatus[i] = kFALSE;
    }
    SortParticles();
    return;
  }
  else
    fMaximumMass = value;
  return;
}

void DatabasePDG::SetMassRange(Double_t min, Double_t max) {
  // Set the mass range for the particle definitions to be generated in the soft fireball

  if(fNParticles>0) {
    fMinimumMass = min;
    fMaximumMass = max;
    for(Int_t i=0; i<fNParticles; i++) {
      if((fParticles[i]->GetMass()<fMinimumMass) || (fParticles[i]->GetMass()>fMaximumMass))  
	fStatus[i] = kFALSE;
    }
    SortParticles();
    return;
  }
  else {
    fMinimumMass = min;
    fMaximumMass = max;
  }
  return;
}

void DatabasePDG::SortParticles() {
  // Sort the particle list so that those with kTRUE status will be always on top of the list

  if(fNParticles<2) {
    cout << "Warning in DatabasePDG::SortParticles() : No particles to sort. Load data first!!" << endl;
    return;
  }

  Int_t nGoodStatus = 0;
  for(Int_t i=0; i<fNParticles; i++)
    if(fStatus[i]) nGoodStatus++;
  if(nGoodStatus==fNParticles)    // if all particles have good status then there is nothing to do
    return;
  if(nGoodStatus==0)              // no good status particles, again nothing to do
    return;

  Int_t shifts = 1;
  while(shifts) {
    shifts = 0;
    for(Int_t i=0; i<fNParticles-1; i++) {
      if(!fStatus[i] && fStatus[i+1]) {   // switch if false status is imediately before a true status particle
	ParticlePDG *temporaryPointer = fParticles[i];
	fParticles[i] = fParticles[i+1];
	fParticles[i+1] = temporaryPointer;
	Bool_t temporaryStatus = fStatus[i];
	fStatus[i] = fStatus[i+1];
	fStatus[i+1] = temporaryStatus;
	shifts++;
      }
    }
  }
  return;
}

Int_t DatabasePDG::GetNParticles(Bool_t all) const {
  // Return the number of particle definitions in the database
  // If all is kTRUE then return number of all particle
  // If all is kFALSE then return the number of good status (kTRUE) particles

  if(all)
    return fNParticles;

  Int_t nGoodStatus = 0;
  for(Int_t i=0; i<fNParticles; i++)
    if(fStatus[i]) nGoodStatus++;
  return nGoodStatus;
}

void DatabasePDG::UseThisListOfParticles(Char_t *filename, Bool_t exclusive) {
  // Read a list of PDG codes from the file "filename" and mark them with good status (kTRUE)
  // while all the other will be marked kFALSE (only if exclusive = kTRUE)     

  if(fNParticles<1) {
    cout << "Error in DatabasePDG::UseThisListOfParticles(Char_t*, Bool_t) : You must load the data before calling this function!!" << endl;
    return;
  }

  ifstream listFile;
  listFile.open(filename);
  if(!listFile) {
    cout << "ERROR in DatabasePDG::UseThisListOfParticles(Char_t*, Bool_t) : The ASCII file containing the PDG codes list (\""
         << filename << "\") was not found !! Exiting..." << endl;
    return;
  }

  Bool_t flaggedIndexes[kMaxParticles];
  for(Int_t i=0; i<kMaxParticles; i++)
    flaggedIndexes[i] = kFALSE;
  Int_t pdg = 0;
  listFile.exceptions(ios::failbit);
  while(!listFile.eof()) {
    try {
      listFile >> pdg;
    }
    catch (ios::failure const &problem) {
      cout << problem.what() << endl;
      break;
    }
    Int_t found = 0;
    for(Int_t i=0; i<fNParticles; i++) {
      if(fParticles[i]->GetPDG()==pdg) {
	found++;
	flaggedIndexes[i] = kTRUE;
      }
    }
    if(!found) {
      cout << "Warning in DatabasePDG::UseThisListOfParticles(Char_t*, Bool_t) : The particle with PDG code "
	   << pdg << " was asked but not found in the database!!" << endl;
    }
    if(found>1) {
      cout << "Warning in DatabasePDG::UseThisListOfParticles(Char_t*, Bool_t) : The particle with PDG code "
	   << pdg << " was found more than once in the database!!" << endl;
    }
  }

  if(exclusive) {
    for(Int_t i=0; i<kMaxParticles; i++)
      fStatus[i] = flaggedIndexes[i];
  }
  else {
    for(Int_t i=0; i<kMaxParticles; i++)
      fStatus[i] = (fStatus[i] && flaggedIndexes[i]);
  }
  SortParticles();

  return;
}

Bool_t DatabasePDG::IsChannelAllowed(DecayChannel *channel, Double_t motherMass) const {
  // Check if the decay channel "channel" is allowed by using the mother particle mass "motherMass"

  Double_t daughtersSumMass = 0.0;
  for(Int_t i=0; i<channel->GetNDaughters(); i++)
    daughtersSumMass += GetPDGParticle(channel->GetDaughterPDG(i))->GetMass();
  if(daughtersSumMass<=motherMass)
    return kTRUE;
  return kFALSE;
}

Int_t DatabasePDG::GetNAllowedChannels(ParticlePDG *particle, Double_t motherMass) const {
  // Check how many decay channels are allowed for a given particle definition at a given mass

  Int_t nAllowedChannels = 0;
  for(Int_t i=0; i<particle->GetNDecayChannels(); i++) 
    nAllowedChannels += (IsChannelAllowed(particle->GetDecayChannel(i), motherMass) ? 1:0);
  
  return nAllowedChannels;
}
