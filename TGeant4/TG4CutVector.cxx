// $Id$
// Category: physics
//
// See the class description in the header file.

#include "TG4CutVector.h"
#include "TG4G3Defaults.h"

#include <G4Track.hh>
#include <G4ParticleDefinition.hh>
#include <G4VProcess.hh>

TG4CutVector::TG4CutVector()
{
  // initialize fCutVector 
  fCutVector = new TG4doubleVector;
  for (G4int i=0; i<kNoG3Cuts; i++) fCutVector->insert(0.); 
}

TG4CutVector::TG4CutVector(const TG4CutVector& right)
{
  // copy fCutVector 
  fCutVector = new TG4doubleVector;
  for (G4int i=0; i<kNoG3Cuts; i++) {
    fCutVector->insert((*right.fCutVector)[i]);
  }   
}

TG4CutVector::~TG4CutVector() {
//
  delete fCutVector;
}

// operators

TG4CutVector& TG4CutVector::operator=(const TG4CutVector& right)
{
  // check assignement to self
  if (this == &right) return *this;

  // initialize fCutVector 
  fCutVector->clear();
  for (G4int i=0; i<kNoG3Cuts; i++) {
    fCutVector->insert((*right.fCutVector)[i]);
  }
  
  return *this;   
}  

G4double TG4CutVector::operator[](G4int index) const 
{
//
  if (index < kNoG3Cuts)
    return (*fCutVector)[index];
  else {
    TG4Globals::Exception(
      "TG4CutVector::operator[]: index out of the vector scope");
    return 0.;  
  }    
}  

// public methods

void TG4CutVector::SetG3Cut(TG3Cut g3Cut, G4double cutValue)
{
// Sets the cutValue for the specified cut.
// ---

  if (g3Cut<kNoG3Cuts) {
    (*fCutVector)[g3Cut] = cutValue;
  }  
  else {
    TG4Globals::Exception(
      "TG4CutVector::SetG3Cut: Inconsistent cut.");
  }
}

void TG4CutVector::SetG3Defaults()
{
// Sets G3 default values for all cuts.
// ---

  for (G4int i=0; i<kNoG3Cuts; i++) {
   (*fCutVector)[i] = TG4G3Defaults::CutValue(i);
  } 
}

G4double TG4CutVector::GetMinEkine(const G4Track& track) const
{
// Returns the cut value for the particle associated with
// specified track.
// ---

  G4ParticleDefinition* particle = track.GetDefinition(); 
  G4String particleName = particle->GetParticleName();

  if  (particleName == "gamma") {
    return GetMinEkineForGamma(track);
  }    
  else if (particleName == "e-") {    
    return GetMinEkineForElectron(track);
  }  
  else if ((particle->GetParticleType() == "baryon") ||
           (particle->GetParticleType() == "meson") ||
           (particle->GetParticleType() == "nucleus")) { 
    if (particle->GetPDGCharge() == 0) 
      return GetMinEkineForHadron(track);
    else 
      return GetMinEkineForNeutralHadron(track);
  }    
  else if ((particleName == "mu-") || (particleName == "mu+")) {
    return GetMinEkineForMuon(track);
  }  
  else {
    G4String text = "TG4CutVector::GetMinEkine: \n";
    text = text + "    The kinetic energy cut for " + particleName;
    text = text + " is not defined.";   
    TG4Globals::Warning(text);
    return 0.;
  }    
}

G4double TG4CutVector::GetMinEkineForGamma(const G4Track& track) const
{
// Returns the cut value for gamma.
// (Cut is not applied for "opticalphoton" 
//  as it is treated in G4 as a particle different 
//  from "gamma" in G4.)
// ---

  G4cout << "TG4CutVector::GetMinEkineForGamma start" << endl;
  const G4VProcess* kpCreatorProcess = track.GetCreatorProcess();
  G4String processName = "";
  if (kpCreatorProcess) processName = kpCreatorProcess->GetProcessName();

  if ((processName == "eBrem") || (processName == "IeBrem")) {
    return (*fCutVector)[kBCUTE];
  }     
  else if ((processName == "MuBrems") || (processName == "IMuBrems") || 
           (processName == "//hBrems")|| (processName == "//IhBrems")) {
           // hadron Brehmstrahlung is not defined in G4
    return (*fCutVector)[kBCUTM];
  }
  else {
    return (*fCutVector)[kCUTGAM];
  }
}

G4double TG4CutVector::GetMinEkineForElectron(const G4Track& track) const
{
// Returns the cut value for e-.
// Should these cuts be applied to e+ too ??
// ---

  const G4VProcess* kpCreatorProcess = track.GetCreatorProcess();
  G4String processName = "";
  if (kpCreatorProcess) processName = kpCreatorProcess->GetProcessName();

  if ((processName == "eIoni") || (processName == "IeIoni")) {
     // !! Geant4 treats delta rays + continuous energy loss
     // within one process
    return (*fCutVector)[kDCUTE];
  }
  else if ((processName == "MuIoni") || (processName == "IMuIoni")) {
     // !! Geant4 treats delta rays + continuous energy loss
     // within one process
    return (*fCutVector)[kDCUTM];
  }
  else {   
    return (*fCutVector)[kCUTELE];
  }
}

G4double TG4CutVector::GetMinEkineForHadron(const G4Track& track) const
{
// Returns the cut value for charged hadron.
// ---

  return (*fCutVector)[kCUTHAD];
}

G4double TG4CutVector::GetMinEkineForNeutralHadron(const G4Track& track) const
{
// Returns the cut value for neutral hadron.
// ---

  return (*fCutVector)[kCUTNEU];
}

G4double TG4CutVector::GetMinEkineForMuon(const G4Track& track) const
{
// Returns the cut value for neutral muon.
// ---

  return (*fCutVector)[kCUTMUO];
}

G4double TG4CutVector::GetMinEkineForOther(const G4Track& track) const
{
// Returns 0.
// ---

  return 0.;
}

