// $Id$
// Category: physics
//
// See the class description in the header file.

#include "TG4PhysicsManager.h"
#include "TG4PhysicsList.h"
#include "TG4CutVector.h"
#include "TG4FlagVector.h"
#include "TG4G3Defaults.h"

#include <G4ParticleDefinition.hh>

#include <TDatabasePDG.h>

TG4PhysicsManager* TG4PhysicsManager::fgInstance = 0;

TG4PhysicsManager::TG4PhysicsManager()
  : fLock(false),
    fPhysicsList(0),
    fCutVector(0),
    fFlagVector(0) 
{ 
//
  if (fgInstance) {
    TG4Globals::Exception(
      "TG4PhysicsManager: attempt to create two instances of singleton.");
  }
      
  fgInstance = this;  
      
  // initialize fIsCutVector 
  fIsCutVector = new TG4boolVector;
  G4int i;
  //for (i=0; i<kNofParticlesWSP; i++) fIsCutVector->insert(false);  
  for (i=0; i<kNofParticlesWSP; i++) fIsCutVector->push_back(false);

  // initialize fIsFlagVector 
  fIsFlagVector = new TG4boolVector;
  //for (i=0; i<kNofParticlesWSP; i++) fIsFlagVector->insert(false);
  for (i=0; i<kNofParticlesWSP; i++) fIsFlagVector->push_back(false);

  // define fCutNameVector, fFlagNameVector
  FillG3CutNameVector();
  FillG3FlagNameVector();

  // fill process name map
  FillProcessMap();
}

TG4PhysicsManager::TG4PhysicsManager(const TG4PhysicsManager& right) {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4PhysicsManager singleton.");
}

TG4PhysicsManager::~TG4PhysicsManager() {
//
  delete fIsCutVector;
  delete fIsFlagVector;
}

// operators

TG4PhysicsManager& 
TG4PhysicsManager::operator=(const TG4PhysicsManager& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4PhysicsManager singleton.");
    
  return *this;  
}    
          
// private methods

void TG4PhysicsManager::LockException() const
{
// Gives exception in case of attempt to modified physics
// setup after physics manager was locked.
// ---

  G4String text = "TG4PhysicsManager: \n";
  text = text + "    It is too late to change physics setup. \n";
  text = text + "    PhysicsManager has been already locked.";
  TG4Globals::Exception(text);
}

void TG4PhysicsManager::FillG3CutNameVector()
{
// Defines fCutNameVector.
// ---

  fG3CutNameVector.insert("CUTGAM");
  fG3CutNameVector.insert("CUTELE");
  fG3CutNameVector.insert("CUTNEU");
  fG3CutNameVector.insert("CUTHAD");
  fG3CutNameVector.insert("CUTMUO");
  fG3CutNameVector.insert("BCUTE");
  fG3CutNameVector.insert("BCUTM"); 
  fG3CutNameVector.insert("DCUTE");
  fG3CutNameVector.insert("DCUTM");
  fG3CutNameVector.insert("PPCUTM");
}

void TG4PhysicsManager::FillG3FlagNameVector() 
{
// Defines fFlagNameVector.
// ---

  fG3FlagNameVector.insert("PAIR");
  fG3FlagNameVector.insert("COMP");
  fG3FlagNameVector.insert("PHOT");
  fG3FlagNameVector.insert("PFIS");
  fG3FlagNameVector.insert("DRAY");
  fG3FlagNameVector.insert("ANNI");
  fG3FlagNameVector.insert("BREM");
  fG3FlagNameVector.insert("HADR");
  fG3FlagNameVector.insert("MUNU");
  fG3FlagNameVector.insert("DCAY");
  fG3FlagNameVector.insert("LOSS");
  fG3FlagNameVector.insert("MULS");
}

void TG4PhysicsManager::FillProcessMap()
{
// Fills fProcessMap.
// The default G4 process names are used in the map.
// ---

  // multiple scattering
  fProcessMap.Add("msc",  kPMultipleScattering);
  fProcessMap.Add("Imsc", kPMultipleScattering);
    
  // continuous energy loss
  // !! including delta rays
  fProcessMap.Add("eIoni",  kPEnergyLoss);
  fProcessMap.Add("IeIoni", kPEnergyLoss);
  fProcessMap.Add("LowEnergyIoni", kPEnergyLoss);
  fProcessMap.Add("hIoni",  kPEnergyLoss);
  fProcessMap.Add("IhIoni", kPEnergyLoss);
  fProcessMap.Add("hLowEIoni", kPEnergyLoss);
  fProcessMap.Add("MuIoni", kPEnergyLoss);
  fProcessMap.Add("IMuIonisation", kPEnergyLoss);
  fProcessMap.Add("ionIoni",  kPEnergyLoss);
  fProcessMap.Add("ionLowEIoni",  kPEnergyLoss);
  fProcessMap.Add("PAIonisation",  kPEnergyLoss);
  
  // bending in mag. field
  // kPMagneticFieldL

  // particle decay
  fProcessMap.Add("Decay", kPDecay);
  
  // photon pair production or
  // muon direct pair production
  fProcessMap.Add("conv", kPPair);
  fProcessMap.Add("LowEnConversion", kPPair);
  fProcessMap.Add("MuPairProd", kPPair);
  fProcessMap.Add("IMuPairProduction", kPPair);

  // Compton scattering
  fProcessMap.Add("compt", kPCompton);
  fProcessMap.Add("LowEnCompton", kPCompton);
  fProcessMap.Add("polarCompt", kPCompton);

  // photoelectric effect
  fProcessMap.Add("phot", kPPhotoelectric);
  fProcessMap.Add("LowEnPhotoElec", kPPhotoelectric);

  // bremsstrahlung
  fProcessMap.Add("eBrem", kPBrem);
  fProcessMap.Add("IeBrem", kPBrem);
  fProcessMap.Add("MuBrem", kPBrem);
  fProcessMap.Add("IMuBremsstrahlung", kPBrem);
  fProcessMap.Add("LowEnBrem", kPBrem);

  // delta-ray production
  // kPDeltaRay
  // has to be distinguished from kPEnergyLoss on flight
  
  // positron annihilation
  fProcessMap.Add("annihil", kPAnnihilation);
  fProcessMap.Add("Iannihil", kPAnnihilation);

  // hadronic interaction
  // kPHadronic

  // nuclear evaporation
  // kPEvaporation
  
  // nuclear fission
  // kPNuclearFission

  // nuclear absorption
  fProcessMap.Add("PionMinusAbsorptionAtRest", kPNuclearAbsorption);
  fProcessMap.Add("PiMinusAbsorptionAtRest", kPNuclearAbsorption);
  fProcessMap.Add("KaonMinusAbsorption", kPNuclearAbsorption);         
  fProcessMap.Add("KaonMinusAbsorptionAtRest", kPNuclearAbsorption);         
  
  // antiproton annihilation
  fProcessMap.Add("AntiProtonAnnihilationAtRest", kPPbarAnnihilation);
  // fProcessMap.Add("AntiNeutronAnnihilationAtRest", not defined);

  // neutron capture    
  fProcessMap.Add("NeutronCaptureAtRest", kPNCapture);
  // fProcessMap.Add("LCapture", hadron capture not defined);

  // hadronic elastic incoherent scattering
  fProcessMap.Add("LElastic", kPHElastic);

  // hadronic inelastic scattering
  fProcessMap.Add("inelastic", kPHInhelastic);

  // muon nuclear interaction
  fProcessMap.Add("MuNucl", kPMuonNuclear);

  // exceeded time of flight cut
  // kPTOFlimit
  
  // nuclear photofission
  // kPPhotoFission

  // Rayleigh scattering
  fProcessMap.Add("Rayleigh Scattering", kPRayleigh);

  // no mechanism is active, usually at the entrance of a new volume
  // kPNull

  // particle has fallen below energy threshold and tracking stops
  // kPStop
  
  // Cerenkov photon absorption
  fProcessMap.Add("Absorption", kPLightAbsorption);

  // Cerenkov photon reflection/refraction
  // kPLightScattering, kPLightReflection, kPLightRefraction
  // has to be inquired from the G4OpBoundary process

  // synchrotron radiation
  fProcessMap.Add("SynchrotronRadiation", kPSynchrotron);
}  


G4int TG4PhysicsManager::GetPDGEncoding(G4ParticleDefinition* particle)
{
// Returns the PDG code of particle;
// if standard PDG code is not defined the TDatabasePDG
// is used.
// ---

  // get PDG encoding from G4 particle definition
  G4int pdgEncoding = particle->GetPDGEncoding();

  if (pdgEncoding == 0) {
    // get PDG encoding from TDatabasePDG
  
    // get particle name from the name map
    G4String g4name = particle->GetParticleName();
    G4String tname = fParticleNameMap.GetSecond(g4name);
    if (tname == "Undefined") {
      G4String text = "TG4PhysicsManager::GetPDGEncoding: \n";
      text = text + "    Particle " + g4name;
      text = text + " was not found in the name map.";
      TG4Globals::Exception(text);
    }  
  
    // get particle from TDatabasePDG
    TDatabasePDG* pdgDB = TDatabasePDG::Instance();
    TParticlePDG* particle = pdgDB->GetParticle(tname);
    if (!particle) {
      G4String text = "TG4PhysicsManager::GetPDGEncoding: \n";
      text = text + "    Particle " + tname;
      text = text + " was not found in TDatabasePDG.";
      TG4Globals::Exception(text);
    }  
    
    // get PDG encoding
    pdgEncoding = particle->PdgCode();
  }
    
  return pdgEncoding;  
}  
     
AliMCProcess TG4PhysicsManager::GetMCProcess(const G4String& g4ProcessName)
{
// Returns the AliMCProcess code of process specified by name.
// ---

  G4int processCode = fProcessMap.GetSecond(g4ProcessName);
  
  if (processCode == 0) return kPNoProcess;
  
  return (AliMCProcess)processCode; 
}

G4int TG4PhysicsManager::GetPDGEncoding(G4String particleName)
{
// Returns the PDG code of particle specified by name.
// ---

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4ParticleDefinition* particle = 0;  
  particle = particleTable->FindParticle(particleName);
  if (!particle) {    
    G4String text = "TG4PhysicsManager::GetPDGEncoding:\n";
    text = text + "   G4ParticleTable::FindParticle() " + particleName;
    text = text + " failed.";
    TG4Globals::Exception(text);
  }	

  return GetPDGEncoding(particle);
}  
  
void TG4PhysicsManager::SetCut(TG3Cut cut, G4double cutValue)
{  
// Sets kinetic energy cut (in a G3-like way).
// ---

  if (!fCutVector) {
    // create vector of kinetic energy cut values  
    fCutVector = new TG4CutVector();
  }  
  fCutVector->SetG3Cut(cut, cutValue);
  SwitchIsCutVector(cut);
}  

void TG4PhysicsManager::SetProcess(TG3Flag flag, G4int flagValue)
{
// Sets control process flag (in a G3-like way).
// ---

  if (!fFlagVector) {
    // create vector of control process flag values
    fFlagVector = new TG4FlagVector;
  }  
  fFlagVector->SetG3Flag(flag, flagValue);
} 


Float_t TG4PhysicsManager::Xsec(char* ch, Float_t p1, Int_t i1, Int_t i2)
{
// Not yet implemented -> gives exception.
// ---

  TG4Globals::Exception(
    "TG4PhysicsManager::Xsec: not yet implemented.");

  return 0.;
}    
  
Int_t TG4PhysicsManager::IdFromPDG(Int_t pdgID) const
{
// G4 does not use the integer particle identifiers
// Id <-> PDG is identity.
// ---

  return pdgID;
}  

Int_t TG4PhysicsManager::PDGFromId(Int_t mcID) const
{
// G4 does not use integer particle identifiers
// Id <-> PDG is identity.
// ---

  return mcID;
}  

void  TG4PhysicsManager::DefineParticles()
{
  // ======
  // Taken from TGeant3
  //
  // Use ENDF-6 mapping for ions = 10000*z+10*a+iso
  // and add 1 000 000
  // and numbers above 5 000 000 for special applications
  //

  const Int_t kion=10000000;
  const Int_t kspe=50000000;

  const Double_t kGeV=0.9314943228;
  const Double_t kHslash = 1.0545726663e-27;
  const Double_t kErgGeV = 1/1.6021773349e-3;
  const Double_t kHshGeV = kHslash*kErgGeV;
  const Double_t kYearsToSec = 3600*24*365.25;

  TDatabasePDG *pdgDB = TDatabasePDG::Instance();

  pdgDB->AddParticle("Deuteron","Deuteron",2*kGeV+8.071e-3,kTRUE,
		     0,1,"Ion",kion+10020);
		     
  pdgDB->AddParticle("Triton","Triton",3*kGeV+14.931e-3,kFALSE,
		     kHshGeV/(12.33*kYearsToSec),1,"Ion",kion+10030);

  pdgDB->AddParticle("Alpha","Alpha",4*kGeV+2.424e-3,kTRUE,
		     kHshGeV/(12.33*kYearsToSec),2,"Ion",kion+20040);

  pdgDB->AddParticle("HE3","HE3",3*kGeV+14.931e-3,kFALSE,
		     0,2,"Ion",kion+20030);

  pdgDB->AddParticle("Cherenkov","Cherenkov",0,kFALSE,
		     0,0,"Special",kspe+50);

  pdgDB->AddParticle("FeedbackPhoton","FeedbackPhoton",0,kFALSE,
		     0,0,"Special",kspe+51);


  // To do: define the PDG database extension
  // in a common part.
  //
  // AliMC::ExtendPDGDatabase(); 
  //
  // end of "common" implementation
  // ======

  // map G4 particle names to TDatabasePDG names
  // (the map is built only for particles that have not
  //  defined standard PDG encoding)
  
  fParticleNameMap.Add("deuteron","Deuteron");
  fParticleNameMap.Add("triton",  "Triton");
  fParticleNameMap.Add("alpha",   "Alpha");
  fParticleNameMap.Add("He3",     "HE3");
  fParticleNameMap.Add("opticalphoton","Cherenkov");
  // fParticleNameMap.Add("???","FeedbackPhoton");
  fParticleNameMap.Add("geantino", "Rootino");
  
  // map G4 particle names to TDatabasePDG encodings
  fParticlePDGMap.Add("deuteron", GetPDGEncoding("deuteron"));
  fParticlePDGMap.Add("triton", GetPDGEncoding("triton"));
  fParticlePDGMap.Add("alpha", GetPDGEncoding("alpha"));
  fParticlePDGMap.Add("He3", GetPDGEncoding("He3") );
  fParticlePDGMap.Add("opticalphoton", GetPDGEncoding("opticalphoton"));
  // fParticlePDGMap.Add("???","FeedbackPhoton");
  fParticleNameMap.Add("geantino", GetPDGEncoding("geantino"));

  // add verbose
  G4cout << "Particle maps have been filled." << G4endl;
  //fParticleNameMap.PrintAll();
  //fParticlePDGMap.PrintAll();
}    

void TG4PhysicsManager::SwitchIsCutVector(TG3Cut cut)
{
// Updates the vector of booleans (fIsCutVector) for the specified cut.
// ---

  switch (cut) {
    case kCUTGAM: 
           (*fIsCutVector)[kGamma] = true; 
           break;
    case kBCUTE:
           (*fIsCutVector)[kGamma] = true; 
           break;
    case kBCUTM:
           (*fIsCutVector)[kGamma] = true; 
           break;
    case kCUTELE:
           (*fIsCutVector)[kElectron] = true; 
           break;
    case kDCUTE:
           (*fIsCutVector)[kElectron] = true; 
           break;
    case kDCUTM:
           (*fIsCutVector)[kElectron] = true; 
           break;
    case kCUTNEU:
           (*fIsCutVector)[kNeutralHadron] = true; 
           break;
    case kCUTHAD:
           (*fIsCutVector)[kChargedHadron] = true; 
           break;
    case kCUTMUO:
           (*fIsCutVector)[kMuon] = true; 
           break;
    default:
           break;
  }
}

void TG4PhysicsManager::SwitchIsFlagVector(TG3Flag flag)
{
// Updates the vector of booleans (fIsFlagVector) for the specified flag.
// ---

  switch (flag) {
    case kPAIR: 
           // gamma
           (*fIsFlagVector)[kGamma] = true; 
           break;
    case kCOMP:
           // gamma
           (*fIsFlagVector)[kGamma] = true; 
	   break;
    case kPHOT:
           // gamma
           (*fIsFlagVector)[kGamma] = true; 
	   break;
    case kPFIS:
           // gamma
           (*fIsFlagVector)[kGamma] = true; 
	   break;
    case kDRAY:	
           // all charged particles
           (*fIsFlagVector)[kElectron] = true; 
           (*fIsFlagVector)[kEplus] = true; 
           (*fIsFlagVector)[kChargedHadron] = true; 
           (*fIsFlagVector)[kMuon] = true; 
	   break;
    case kANNI:
	   // e+ only
           (*fIsFlagVector)[kEplus] = true; 
	   break;
    case kBREM:
	   // e-/e+, muons
           (*fIsFlagVector)[kElectron] = true; 
           (*fIsFlagVector)[kEplus] = true; 
           (*fIsFlagVector)[kMuon] = true; 
	   break;
    case kHADR:
	   // hadrons
           (*fIsFlagVector)[kNeutralHadron] = true; 
           (*fIsFlagVector)[kChargedHadron] = true; 
	   break;
    case kMUNU:
           // muons
           (*fIsFlagVector)[kMuon] = true; 
	   break;
    case kDCAY:
           // any
           (*fIsFlagVector)[kAny] = true; 
	   break;
    case kLOSS:
           // all charged particles
           (*fIsFlagVector)[kElectron] = true; 
           (*fIsFlagVector)[kEplus] = true; 
           (*fIsFlagVector)[kChargedHadron] = true; 
           (*fIsFlagVector)[kMuon] = true; 
	   break;
    case kMULS:
           // all charged particles
           (*fIsFlagVector)[kElectron] = true; 
           (*fIsFlagVector)[kEplus] = true; 
           (*fIsFlagVector)[kChargedHadron] = true; 
           (*fIsFlagVector)[kMuon] = true; 
	   break;
    default:
          break;
  }
}

TG3Cut TG4PhysicsManager::GetG3Cut(G4String cutName)
{
// Retrieves corresponding TG3Cut constant from the cutName.
// ---

  if      (cutName == fG3CutNameVector[kCUTGAM]) return kCUTGAM; 
  else if (cutName == fG3CutNameVector[kBCUTE])  return kBCUTE; 
  else if (cutName == fG3CutNameVector[kBCUTM])  return kBCUTM; 
  else if (cutName == fG3CutNameVector[kCUTELE]) return kCUTELE;
  else if (cutName == fG3CutNameVector[kDCUTE])  return kDCUTE;
  else if (cutName == fG3CutNameVector[kDCUTM])  return kDCUTM;
  else if (cutName == fG3CutNameVector[kCUTNEU]) return kCUTNEU;
  else if (cutName == fG3CutNameVector[kCUTHAD]) return kCUTHAD;
  else if (cutName == fG3CutNameVector[kCUTMUO]) return kCUTMUO;
  else return kNoG3Cuts;
}

TG3Flag TG4PhysicsManager::GetG3Flag(G4String flagName)
{
// Retrieves corresponding TG3Flag constant from the flagName.
// ---

  if      (flagName == fG3FlagNameVector[kPAIR]) return kPAIR;
  else if (flagName == fG3FlagNameVector[kCOMP]) return kCOMP;
  else if (flagName == fG3FlagNameVector[kPHOT]) return kPHOT;
  else if (flagName == fG3FlagNameVector[kPFIS]) return kPFIS;
  else if (flagName == fG3FlagNameVector[kDRAY]) return kDRAY;
  else if (flagName == fG3FlagNameVector[kANNI]) return kANNI;
  else if (flagName == fG3FlagNameVector[kBREM]) return kBREM;
  else if (flagName == fG3FlagNameVector[kHADR]) return kHADR;
  else if (flagName == fG3FlagNameVector[kMUNU]) return kMUNU;
  else if (flagName == fG3FlagNameVector[kDCAY]) return kDCAY;
  else if (flagName == fG3FlagNameVector[kLOSS]) return kLOSS;
  else if (flagName == fG3FlagNameVector[kMULS]) return kMULS;
  else return kNoG3Flags;
}

// public methods

void TG4PhysicsManager::BuildPhysics()
{
// Empty function - not needed in G4.
// (Physics is built within /run/initialize.)

  TG4Globals::Warning(
    "TG4PhysicsManager::BuildPhysics: is empty function in G4 MC.");
}    

void TG4PhysicsManager::SetCut(const char* cutName, Float_t cutValue)
{
// Sets the specified cut.
// ---

  if (fLock) LockException();
  TG3Cut g3Cut = GetG3Cut(cutName);
  if (g3Cut != kNoG3Cuts)
    SetCut(g3Cut, cutValue);
  else {   
    G4String text = "TG4PhysicsManager::SetCut:\n";
    text = text + "    Parameter " + cutName;
    text = text + " is not implemented.";
    TG4Globals::Warning(text);
  }  
}  
  
void TG4PhysicsManager::SetProcess(const char* flagName, Int_t flagValue)
{
// Sets the specified process control.
// ---

  if (fLock) LockException();
  TG3Flag g3Flag = GetG3Flag(flagName);
  if (g3Flag != kNoG3Flags)
    SetProcess(g3Flag, flagValue);
  else {   
    G4String text = "TG4PhysicsManager::SetProcess:\n";
    text = text + "    Parameter " + flagName;
    text = text + " is not implemented.";
    TG4Globals::Warning(text);
  }  
}  

void TG4PhysicsManager::SetProcessActivation()
{
// (In)Activates built processes according
// to the setup in fFlagVector.
// ---

  if (fPhysicsList) {
    // temporarily excluded
    // fPhysicsList->SetProcessActivation();
  }  
  else {
    G4String text = "TG4PhysicsManager::SetProcessActivation:\n";
    text = text +   "   There is no physics list set.";
    TG4Globals::Exception(text);
  }
}       

G4int TG4PhysicsManager::GetPDGEncodingFast(G4ParticleDefinition* particle)
{
// Returns the PDG code of particle;
// if standard PDG code is not defined the preregistred
// fParticlePDGMap is used.
// ---

  // get PDG encoding from G4 particle definition
  G4int pdgEncoding = particle->GetPDGEncoding();

  if (pdgEncoding == 0) {
    // use FParticlePDGMap if standard PDG code is not defined
    G4String name = particle->GetParticleName();
    pdgEncoding = fParticlePDGMap.GetSecond(name);
  }
    
  return pdgEncoding;  
}  
     
G4bool TG4PhysicsManager::CheckCutWithCutVector(G4String name, 
                             G4double value, TG3Cut& cut)
{
// Retrieves corresponding TG3Cut from the name and 
// in case the value is different from the value in cutVector
// sets true the value of the fIsCutVector element 
// corresponding to this cut and returns true; 
// returns false otherwise.
// ---

  // convert cut name -> TG3Cut
  cut = GetG3Cut(name);

  // set switch vector element only if the value
  // is different from the value in cutVector
  if (cut !=kNoG3Cuts) {
    // get tolerance from TG4G3Defaults in GeV
    G4double tolerance = TG4G3Defaults::CutTolerance()/GeV;
    if (!(fCutVector) || (abs(value - (*fCutVector)[cut]) > tolerance)) {
      SwitchIsCutVector(cut);      
      return true;
    }  
    else return false;  
  }   	   	         
  return false;
}

G4bool TG4PhysicsManager::CheckFlagWithFlagVector(G4String name, 
                              G4double value, TG3Flag& flag)
{
// Retrieves corresponding TG3Flag from the name and 
// in case the value is different from the value in flagVector
// sets true the value of the fIsFlagVector element 
// corresponding to this flag and returns true; 
// returns false otherwise.
// ---

  // convert flag name -> TG3Flag
  flag = GetG3Flag(name);

  // set switch vector element only if the value
  // is different from the value in flagVector
  if (flag !=kNoG3Flags) {
    if (!(fFlagVector) || (abs(value - (*fFlagVector)[flag]) > 0.01)) {
      SwitchIsFlagVector(flag);      
      return true;
    }  
    else return false;  
  }   	   	         
  return false;
}

G4bool TG4PhysicsManager::CheckCutWithG3Defaults(G4String name, 
                              G4double value, TG3Cut& cut)
{
// Retrieves corresponding TG3Cut from the name and 
// in case the value is different from the G3 default value
// sets true the value of the SwitchCutVector element 
// corresponding to this cut and returns true; 
// returns false otherwise.
// ---

  // convert cut name -> TG3Cut
  cut = GetG3Cut(name);

  // set switch vector element only if the value
  // is different from G3 default
  if (cut !=kNoG3Cuts) {
    if (!TG4G3Defaults::IsDefaultCut(cut, value)) {
      SwitchIsCutVector(cut);      
      return true;
    }  
    else return false;  
  }   	   	         
  return false;
}

G4bool TG4PhysicsManager::CheckFlagWithG3Defaults(G4String name, 
                              G4double value, TG3Flag& flag)
{
// Retrieves corresponding TG3Flag from the name and 
// in case the value is different from the G3 default value
// sets true the value of the SwitchFlagVector element 
// corresponding to this flag and returns true; 
// returns false otherwise.
// ---

  // convert flag name -> TG3Flag
  flag = GetG3Flag(name);

  // set switch vector element only if the value
  // is different from G3 default
  if (flag !=kNoG3Flags) {
    if (!TG4G3Defaults::IsDefaultFlag(flag, value)) {
      SwitchIsFlagVector(flag);      
      return true;
    }  
    else return false;  
  }   	   	         
  return false;
}

void TG4PhysicsManager::SetG3DefaultCuts() 
{
// Sets G3 default values of kinetic energy cuts.
// ---

  if (fLock) LockException();
  if (!fCutVector) {
    // create vector of kinetic energy cut values  
    fCutVector = new TG4CutVector();
  }  
  fCutVector->SetG3Defaults();
}

void TG4PhysicsManager::SetG3DefaultProcesses()
{
// Sets G3 default values of control process flags.
// ---

  if (fLock) LockException();
  if (!fFlagVector) {
    // create vector of control process flag values
    fFlagVector = new TG4FlagVector;
  }  
  fFlagVector->SetG3Defaults();
}  

G4bool TG4PhysicsManager::IsSpecialCuts() const
{
// Returns true if any special cut value is set.
// ---

  for (G4int i=0; i<kNofParticlesWSP; i++)
  {  if ((*fIsCutVector)[i]) return true; }

  return false;
}

G4bool TG4PhysicsManager::IsSpecialFlags() const
{
// Returns true if any special flag value is set.
// ---

  for (G4int i=0; i<kNofParticlesWSP; i++)
  {  if ((*fIsFlagVector)[i]) return true; }

  return false;
}

TG3ParticleWSP TG4PhysicsManager::GetG3ParticleWSP(
                                  G4ParticleDefinition* particle) const 
{
// Returns TG3ParticleWSP constant for the specified particle.
// (See TG3ParticleWSP.h, too.)
// ---

  G4String name = particle->GetParticleName();     
  G4String pType = particle->GetParticleType();
    
  if (name == "gamma") {
    return kGamma;
  }  
  else if (name == "e-") {    
    return kElectron;
  }  
  else if (name == "e+") {   
    return kEplus;
  }  
  else if (( pType == "baryon" || pType == "meson" || pType == "nucleus" )) {
    if (particle->GetPDGCharge() == 0) { 
      return kNeutralHadron;
    }
    else  
      return kChargedHadron;
  }    
  else if ( name == "mu-" || name == "mu+" ) {
    return kMuon;
  }  
  else {
    return kNofParticlesWSP;
  }    
}  

void TG4PhysicsManager::GetG3ParticleWSPName(G4int particleWSP,
                                             G4String& name) const 
{
// Fills the passed name with the name/type of particle specified
// by TG3ParticleWSP constant.
// (See TG3ParticleWSP.h, too.)
// ---

  switch (particleWSP) {
    case kGamma:
      name = "Gamma";
      break;
    case kElectron:
      name = "Electron";
      break;
    case kEplus:
      name = "Eplus";
      break;
    case kNeutralHadron:
      name = "NeutralHadron";
      break;
    case kChargedHadron:
      name = "ChargedHadron";
      break;
    case kMuon:
      name = "Muon";
      break;
    case kAny:
      name = "Any";
      break;
    case kNofParticlesWSP:
      name = "NoSP";
      break;
    default:
      G4String text = "TG4PhysicsList::GetG3ParticleWSPName:\n";
      text = text + "   Wrong particleWSP."; 
      TG4Globals::Exception(text);
      break;      
  }
}  

