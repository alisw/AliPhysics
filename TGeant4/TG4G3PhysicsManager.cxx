// $Id$
// Category: physics
//
// See the class description in the header file.

#include "TG4G3PhysicsManager.h"
#include "TG4G3CutVector.h"
#include "TG4G3ControlVector.h"
#include "TG4G3Defaults.h"

#include <G4ParticleDefinition.hh>
#include <G4VProcess.hh>
#include <G4UImessenger.hh>
#include <G4ProcessTable.hh>

TG4G3PhysicsManager* TG4G3PhysicsManager::fgInstance = 0;

TG4G3PhysicsManager::TG4G3PhysicsManager()
  : fLock(false),
    fCutVector(0),
    fControlVector(0) 
{ 
  if (fgInstance) {
    TG4Globals::Exception(
      "TG4G3PhysicsManager: attempt to create two instances of singleton.");
  }
      
  fgInstance = this;  

  // initialize fIsCutVector 
  fIsCutVector = new TG4boolVector;
  G4int i;
  //for (i=0; i<kNofParticlesWSP; i++) fIsCutVector->insert(false);  
  for (i=0; i<kNofParticlesWSP; i++) fIsCutVector->push_back(false);

  // initialize fIsControlVector 
  fIsControlVector = new TG4boolVector;
  //for (i=0; i<kNofParticlesWSP; i++) fIsControlVector->insert(false);
  for (i=0; i<kNofParticlesWSP; i++) fIsControlVector->push_back(false);

  // define fCutNameVector, fControlNameVector
  FillG3CutNameVector();
  FillG3ControlNameVector();
}

TG4G3PhysicsManager::TG4G3PhysicsManager(const TG4G3PhysicsManager& right) {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4G3PhysicsManager singleton.");
}

TG4G3PhysicsManager::~TG4G3PhysicsManager() {
//
  delete fIsCutVector;
  delete fIsControlVector;
}

// operators

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

void TG4G3PhysicsManager::FillG3CutNameVector()
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

void TG4G3PhysicsManager::FillG3ControlNameVector() 
{
// Defines fControlNameVector.
// ---

  fG3ControlNameVector.insert("PAIR");
  fG3ControlNameVector.insert("COMP");
  fG3ControlNameVector.insert("PHOT");
  fG3ControlNameVector.insert("PFIS");
  fG3ControlNameVector.insert("DRAY");
  fG3ControlNameVector.insert("ANNI");
  fG3ControlNameVector.insert("BREM");
  fG3ControlNameVector.insert("HADR");
  fG3ControlNameVector.insert("MUNU");
  fG3ControlNameVector.insert("DCAY");
  fG3ControlNameVector.insert("LOSS");
  fG3ControlNameVector.insert("MULS");
}

void TG4G3PhysicsManager::SetCut(TG4G3Cut cut, G4double cutValue)
{  
// Sets kinetic energy cut (in a G3-like way).
// ---

  if (!fCutVector) {
    // create vector of kinetic energy cut values  
    fCutVector = new TG4G3CutVector();
  }  
  fCutVector->SetG3Cut(cut, cutValue);
  SwitchIsCutVector(cut);
}  

void TG4G3PhysicsManager::SetProcess(TG4G3Control control, G4int controlValue)
{
// Sets control process control (in a G3-like way).
// ---

  if (!fControlVector) {
    // create vector of control process control values
    fControlVector = new TG4G3ControlVector;
  }  
  fControlVector->SetG3Control(control, controlValue);
} 


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

TG4G3Cut TG4G3PhysicsManager::GetG3Cut(G4String cutName)
{
// Retrieves corresponding TG4G3Cut constant from the cutName.
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

TG4G3Control TG4G3PhysicsManager::GetG3Control(G4String controlName)
{
// Retrieves corresponding TG4G3Control constant from the controlName.
// ---

  if      (controlName == fG3ControlNameVector[kPAIR]) return kPAIR;
  else if (controlName == fG3ControlNameVector[kCOMP]) return kCOMP;
  else if (controlName == fG3ControlNameVector[kPHOT]) return kPHOT;
  else if (controlName == fG3ControlNameVector[kPFIS]) return kPFIS;
  else if (controlName == fG3ControlNameVector[kDRAY]) return kDRAY;
  else if (controlName == fG3ControlNameVector[kANNI]) return kANNI;
  else if (controlName == fG3ControlNameVector[kBREM]) return kBREM;
  else if (controlName == fG3ControlNameVector[kHADR]) return kHADR;
  else if (controlName == fG3ControlNameVector[kMUNU]) return kMUNU;
  else if (controlName == fG3ControlNameVector[kDCAY]) return kDCAY;
  else if (controlName == fG3ControlNameVector[kLOSS]) return kLOSS;
  else if (controlName == fG3ControlNameVector[kMULS]) return kMULS;
  else return kNoG3Controls;
}

// public methods


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

G4bool TG4G3PhysicsManager::CheckControlWithTheVector(G4String name, 
                                 G4double value, TG4G3Control& control)
{
// Retrieves corresponding TG4G3Control from the name and 
// in case the value is different from the value in controlVector
// sets true the value of the fIsControlVector element 
// corresponding to this control and returns true; 
// returns false otherwise.
// ---

  // convert control name -> TG4G3Control
  control = GetG3Control(name);

  // set switch vector element only if the value
  // is different from the value in controlVector
  if (control !=kNoG3Controls) {
    if (!(fControlVector) || (abs(value - (*fControlVector)[control]) > 0.01)) {
      SwitchIsControlVector(control);      
      return true;
    }  
    else return false;  
  }   	   	         
  return false;
}

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

G4bool TG4G3PhysicsManager::CheckControlWithG3Defaults(G4String name, 
                                 G4double value, TG4G3Control& control)
{
// Retrieves corresponding TG4G3Control from the name and 
// in case the value is different from the G3 default value
// sets true the value of the SwitchControlVector element 
// corresponding to this control and returns true; 
// returns false otherwise.
// ---

  // convert control name -> TG4G3Control
  control = GetG3Control(name);

  // set switch vector element only if the value
  // is different from G3 default
  if (control !=kNoG3Controls) {
    if (!TG4G3Defaults::IsDefaultControl(control, value)) {
      SwitchIsControlVector(control);      
      return true;
    }  
    else return false;  
  }   	   	         
  return false;
}

void TG4G3PhysicsManager::SetG3DefaultCuts() 
{
// Sets G3 default values of kinetic energy cuts.
// ---

  CheckLock();
  if (!fCutVector) {
    // create vector of kinetic energy cut values  
    fCutVector = new TG4G3CutVector();
  }  
  fCutVector->SetG3Defaults();
}

void TG4G3PhysicsManager::SetG3DefaultControls()
{
// Sets G3 default values of control process controls.
// ---

  CheckLock();
  if (!fControlVector) {
    // create vector of control process control values
    fControlVector = new TG4G3ControlVector;
  }  
  fControlVector->SetG3Defaults();
}  

G4bool TG4G3PhysicsManager::IsSpecialCuts() const
{
// Returns true if any special cut value is set.
// ---

  for (G4int i=0; i<kNofParticlesWSP; i++)
  {  if ((*fIsCutVector)[i]) return true; }

  return false;
}

G4bool TG4G3PhysicsManager::IsSpecialControls() const
{
// Returns true if any special control value is set.
// ---

  for (G4int i=0; i<kNofParticlesWSP; i++)
  {  if ((*fIsControlVector)[i]) return true; }

  return false;
}

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

void TG4G3PhysicsManager::GetG3ParticleWSPName(G4int particleWSP,
                                               G4String& name) const 
{
// Fills the passed name with the name/type of particle specified
// by TG4G3ParticleWSP constant.
// (See TG4G3ParticleWSP.h, too.)
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
      G4String text = "TG4G3PhysicsManager::GetG3ParticleWSPName:\n";
      text = text + "   Wrong particleWSP."; 
      TG4Globals::Exception(text);
      break;      
  }
}  

