// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4G3PhysicsManager
// -------------------------
// See the class description in the header file.

#include "TG4G3PhysicsManager.h"
#include "TG4G3CutVector.h"
#include "TG4G3ControlVector.h"
#include "TG4G3Defaults.h"
#include "TG4G3Units.h"

#include <G4ParticleDefinition.hh>
#include <G4VProcess.hh>
#include <G4UImessenger.hh>
#include <G4ProcessTable.hh>

TG4G3PhysicsManager* TG4G3PhysicsManager::fgInstance = 0;

//_____________________________________________________________________________
TG4G3PhysicsManager::TG4G3PhysicsManager()
  : fLock(false),
    fCutVector(0),
    fControlVector(0),
    fIsCutVector(0),
    fIsControlVector(0) 
{
// 
  if (fgInstance) {
    TG4Globals::Exception(
      "TG4G3PhysicsManager: attempt to create two instances of singleton.");
  }
      
  fgInstance = this;  

  // create vector of kinetic energy cut values  
  fCutVector = new TG4G3CutVector();

  // create vector of control process control values
  fControlVector = new TG4G3ControlVector();

  // initialize fIsCutVector 
  fIsCutVector = new TG4boolVector();
  G4int i;
  for (i=0; i<kNofParticlesWSP; i++) fIsCutVector->push_back(false);

  // initialize fIsControlVector 
  fIsControlVector = new TG4boolVector;
  for (i=0; i<kNofParticlesWSP; i++) fIsControlVector->push_back(false);
}

//_____________________________________________________________________________
TG4G3PhysicsManager::TG4G3PhysicsManager(const TG4G3PhysicsManager& right) {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4G3PhysicsManager singleton.");
}

//_____________________________________________________________________________
TG4G3PhysicsManager::~TG4G3PhysicsManager() {
//
//  delete fIsCutVector;
  delete fIsControlVector;
}

// operators

//_____________________________________________________________________________
TG4G3PhysicsManager& 
TG4G3PhysicsManager::operator=(const TG4G3PhysicsManager& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4G3PhysicsManager singleton.");
    
  return *this;  
}    
          
// private methods

//_____________________________________________________________________________
void TG4G3PhysicsManager::SetCut(TG4G3Cut cut, G4double cutValue)
{  
// Sets kinetic energy cut (in a G3-like way).
// ---

  fCutVector->SetCut(cut, cutValue);
  SwitchIsCutVector(cut);
}  

//_____________________________________________________________________________
void TG4G3PhysicsManager::SetProcess(TG4G3Control control, 
                                     TG4G3ControlValue controlValue)
{
// Sets control process control (in a G3-like way).
// ---
  
  if (control == kDRAY || control == kLOSS) {
      SwitchIsCutVector(kDCUTE);
      SwitchIsCutVector(kDCUTM);
  }  

  fControlVector->SetControl(control, controlValue, *fCutVector);
} 


//_____________________________________________________________________________
void TG4G3PhysicsManager::SwitchIsCutVector(TG4G3Cut cut)
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
           (*fIsCutVector)[kEplus] = true; 
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

//_____________________________________________________________________________
void TG4G3PhysicsManager::SwitchIsControlVector(TG4G3Control control)
{
// Updates the vector of booleans (fIsControlVector) for the specified control.
// ---

  switch (control) {
    case kPAIR: 
           // gamma
           (*fIsControlVector)[kGamma] = true; 
           break;
    case kCOMP:
           // gamma
           (*fIsControlVector)[kGamma] = true; 
	   break;
    case kPHOT:
           // gamma
           (*fIsControlVector)[kGamma] = true; 
	   break;
    case kPFIS:
           // gamma
           (*fIsControlVector)[kGamma] = true; 
	   break;
    case kDRAY:	
           // all charged particles
           (*fIsControlVector)[kElectron] = true; 
           (*fIsControlVector)[kEplus] = true; 
           (*fIsControlVector)[kChargedHadron] = true; 
           (*fIsControlVector)[kMuon] = true; 
	   break;
    case kANNI:
	   // e+ only
           (*fIsControlVector)[kEplus] = true; 
	   break;
    case kBREM:
	   // e-/e+, muons
           (*fIsControlVector)[kElectron] = true; 
           (*fIsControlVector)[kEplus] = true; 
           (*fIsControlVector)[kMuon] = true; 
	   break;
    case kHADR:
	   // hadrons
           (*fIsControlVector)[kNeutralHadron] = true; 
           (*fIsControlVector)[kChargedHadron] = true; 
	   break;
    case kMUNU:
           // muons
           (*fIsControlVector)[kMuon] = true; 
	   break;
    case kDCAY:
           // any
           (*fIsControlVector)[kAny] = true; 
	   break;
    case kLOSS:
           // all charged particles
           (*fIsControlVector)[kElectron] = true; 
           (*fIsControlVector)[kEplus] = true; 
           (*fIsControlVector)[kChargedHadron] = true; 
           (*fIsControlVector)[kMuon] = true; 
	   break;
    case kMULS:
           // all charged particles
           (*fIsControlVector)[kElectron] = true; 
           (*fIsControlVector)[kEplus] = true; 
           (*fIsControlVector)[kChargedHadron] = true; 
           (*fIsControlVector)[kMuon] = true; 
	   break;
    default:
          break;
  }
}

// public methods

//_____________________________________________________________________________
void TG4G3PhysicsManager::CheckLock()
{
// Gives exception in case the physics manager is locked.
// Prevents from modifying physics setup after the physics manager is locked.
// ---

  if (fLock) {
    G4String text = "TG4PhysicsManager: \n";
    text = text + "    It is too late to change physics setup. \n";
    text = text + "    PhysicsManager has been already locked.";
    TG4Globals::Exception(text);
  }  
}

//_____________________________________________________________________________
G4VProcess* TG4G3PhysicsManager::FindProcess(G4String processName) const
{
// Finds G4VProcess with specified name.
// ---

  G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();

  G4ProcessVector* processVector 
    = processTable->FindProcesses(processName);
  G4VProcess* firstFoundProcess = 0;
  if (processVector->entries()>0) firstFoundProcess= (*processVector)[0];

  processVector->clear();
  delete processVector;
  
  return firstFoundProcess;
}

//_____________________________________________________________________________
G4bool TG4G3PhysicsManager::CheckCutWithTheVector(G4String name, 
                                 G4double value, TG4G3Cut& cut)
{
// Retrieves corresponding TG4G3Cut from the name and 
// in case the value is different from the value in cutVector
// sets true the value of the fIsCutVector element 
// corresponding to this cut and returns true; 
// returns false otherwise.
// ---

  // convert cut name -> TG4G3Cut
  cut = TG4G3CutVector::GetCut(name);

  // set switch vector element only if the value
  // is different from the value in cutVector
  if (cut !=kNoG3Cuts) {
    // get tolerance from TG4G3CutVector in G3 units
    G4double tolerance = TG4G3CutVector::Tolerance()/ TG4G3Units::Energy();
    if (abs(value - (*fCutVector)[cut]) > tolerance) {
      SwitchIsCutVector(cut);      
      return true;
    }  
    else return false;  
  }   	   	         
  return false;
}

//_____________________________________________________________________________
G4bool TG4G3PhysicsManager::CheckControlWithTheVector(G4String name, 
                                 G4double value, TG4G3Control& control,
				 TG4G3ControlValue& controlValue)
{
// Retrieves corresponding TG4G3Control from the name and 
// in case the value is different from the value in controlVector
// sets true the value of the fIsControlVector element 
// corresponding to this control and returns true; 
// returns false otherwise.
// ---

  // convert control name -> TG4G3Control
  control = TG4G3ControlVector::GetControl(name);

  // convert double value -> TG4G3ControlValue
  controlValue = TG4G3ControlVector::GetControlValue(value, control);

  // set switch vector element only if the value
  // is different from the value in controlVector
  if (control !=kNoG3Controls) {
    if (controlValue != (*fControlVector)[control]) {
      SwitchIsControlVector(control);      
      return true;
    }  
    else return false;  
  }   	   	         
  return false;
}

//_____________________________________________________________________________
G4bool TG4G3PhysicsManager::CheckCutWithG3Defaults(G4String name, 
                                 G4double value, TG4G3Cut& cut)
{
// Retrieves corresponding TG4G3Cut from the name and 
// in case the value is different from the G3 default value
// sets true the value of the SwitchCutVector element 
// corresponding to this cut and returns true; 
// returns false otherwise.
// ---

  // convert cut name -> TG4G3Cut
  cut = TG4G3CutVector::GetCut(name);

  // set switch vector element only if the value
  // is different from G3 default
  if (cut !=kNoG3Cuts) {
    if (!fG3Defaults.IsDefaultCut(cut, value)) {
      SwitchIsCutVector(cut);      
      return true;
    }  
    else return false;  
  }   	   	         
  return false;
}

//_____________________________________________________________________________
G4bool TG4G3PhysicsManager::CheckControlWithG3Defaults(G4String name, 
                                 G4double value, TG4G3Control& control,
                                 TG4G3ControlValue& controlValue)
{
// Retrieves corresponding TG4G3Control from the name and 
// in case the value is different from the G3 default value
// sets true the value of the SwitchControlVector element 
// corresponding to this control and returns true; 
// returns false otherwise.
// ---

  // convert control name -> TG4G3Control
  control = TG4G3ControlVector::GetControl(name);

  // convert double value -> TG4G3ControlValue
  controlValue = TG4G3ControlVector::GetControlValue(value, control);

  // set switch vector element only if the value
  // is different from G3 default
  if (control !=kNoG3Controls) {
    if (!fG3Defaults.IsDefaultControl(control, controlValue)) {
      SwitchIsControlVector(control);      
      return true;
    }  
    else return false;  
  }   	   	         
  return false;
}

//_____________________________________________________________________________
void TG4G3PhysicsManager::SetG3DefaultCuts() 
{
// Sets G3 default values of kinetic energy cuts.
// ---

  CheckLock();
  fCutVector->SetG3Defaults();
}

//_____________________________________________________________________________
void TG4G3PhysicsManager::SetG3DefaultControls()
{
// Sets G3 default values of control process controls.
// ---

  CheckLock();
  fControlVector->SetG3Defaults();
}  

//_____________________________________________________________________________
G4bool TG4G3PhysicsManager::IsSpecialCuts() const
{
// Returns true if any special cut value is set.
// ---


  for (G4int i=0; i<kNofParticlesWSP; i++)
  {  if ((*fIsCutVector)[i]) return true; }

  return false;
}

//_____________________________________________________________________________
G4bool TG4G3PhysicsManager::IsSpecialControls() const
{
// Returns true if any special control value is set.
// ---

  for (G4int i=0; i<kNofParticlesWSP; i++)
  {  if ((*fIsControlVector)[i]) return true; }

  return false;
}

//_____________________________________________________________________________
TG4G3ParticleWSP TG4G3PhysicsManager::GetG3ParticleWSP(
                                      G4ParticleDefinition* particle) const 
{
// Returns TG4G3ParticleWSP constant for the specified particle.
// (See TG4G3ParticleWSP.h, too.)
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

//_____________________________________________________________________________
G4String TG4G3PhysicsManager::GetG3ParticleWSPName(G4int particleWSP) const 
{
// Returns the name of the particle specified by TG4G3ParticleWSP constant.
// (See TG4G3ParticleWSP.h, too.)
// ---

  switch (particleWSP) {
    case kGamma:
      return "Gamma";
      break;
    case kElectron:
      return "Electron";
      break;
    case kEplus:
      return "Eplus";
      break;
    case kNeutralHadron:
      return "NeutralHadron";
      break;
    case kChargedHadron:
      return "ChargedHadron";
      break;
    case kMuon:
      return "Muon";
      break;
    case kAny:
      return "Any";
      break;
    case kNofParticlesWSP:
      return "NoSP";
      break;
    default:
      G4String text = "TG4G3PhysicsManager::GetG3ParticleWSPName:\n";
      text = text + "   Wrong particleWSP."; 
      TG4Globals::Exception(text);
      return "";
  }
}  

