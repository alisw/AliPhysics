// $Id$
// Category: physics
//
// See the class description in the header file.

#include "TG4PhysicsManager.h"
#include "TG4ParticlesManager.h"
#include "TG4G3PhysicsManager.h"
#include "TG4PhysicsConstructorEM.h"
#include "TG4PhysicsConstructorOptical.h"
#include "TG4PhysicsConstructorHadron.h"
#include "TG4PhysicsConstructorSpecialCuts.h"
#include "TG4PhysicsConstructorSpecialControls.h"
#include "TG4G3Cut.h"
#include "TG4G3Control.h"

#include <G4ParticleDefinition.hh>
#include <G4VProcess.hh>
#include <G4VModularPhysicsList.hh>

#include <TDatabasePDG.h>

TG4PhysicsManager* TG4PhysicsManager::fgInstance = 0;

TG4PhysicsManager::TG4PhysicsManager(G4VModularPhysicsList* physicsList)
  : fPhysicsList(physicsList),
    fSetEMPhysics(true),
    fSetOpticalPhysics(false),
    fSetHadronPhysics(false),
    fSetSpecialCutsPhysics(false),
    fSetSpecialControlsPhysics(false)

{ 
//
  if (fgInstance) {
    TG4Globals::Exception(
      "TG4PhysicsManager: attempt to create two instances of singleton.");
  }
      
  fgInstance = this;  
  
  // create particles manager
  fParticlesManager = new TG4ParticlesManager();
      
  // create G3 physics manager
  fG3PhysicsManager = new TG4G3PhysicsManager();

  // fill process name map
  FillProcessMap();
}

TG4PhysicsManager::TG4PhysicsManager(){
//
  delete fParticlesManager;
  delete fG3PhysicsManager;
}

TG4PhysicsManager::TG4PhysicsManager(const TG4PhysicsManager& right) {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4PhysicsManager singleton.");
}

TG4PhysicsManager::~TG4PhysicsManager() {
//
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
  fProcessMap.Add("Transportation", kPNull);

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


// public methods

void TG4PhysicsManager::BuildPhysics()
{
// Empty function - not needed in G4.
// (Physics is built within /run/initialize.)
// ---

  TG4Globals::Warning(
    "TG4PhysicsManager::BuildPhysics: is empty function in G4 MC.");
}    

void TG4PhysicsManager::CreatePhysicsConstructors()
{
// Creates the selected physics constructors
// and registeres them in the modular physics list.
// ---

  // electromagnetic physics
  if (fSetEMPhysics) 
    fPhysicsList->RegisterPhysics(new TG4PhysicsConstructorEM());

  // optical physics
  if (fSetOpticalPhysics) 
    fPhysicsList->RegisterPhysics(new TG4PhysicsConstructorOptical());

  // hadron physics
  if (fSetHadronPhysics) 
    fPhysicsList->RegisterPhysics(new TG4PhysicsConstructorHadron());

  if (fSetSpecialCutsPhysics) 
    fPhysicsList->RegisterPhysics(new TG4PhysicsConstructorSpecialCuts());

  if (fSetSpecialControlsPhysics) 
    fPhysicsList->RegisterPhysics(new TG4PhysicsConstructorSpecialControls());

         // all created physics constructors are deleted
	 // in the G4VModularPhysicsList destructor
}    

void TG4PhysicsManager::SetCut(const char* cutName, Float_t cutValue)
{
// Sets the specified cut.
// ---

  fG3PhysicsManager->CheckLock();
  TG4G3Cut g3Cut = fG3PhysicsManager->GetG3Cut(cutName);
  if (g3Cut != kNoG3Cuts)
    fG3PhysicsManager->SetCut(g3Cut, cutValue);
  else {   
    G4String text = "TG4PhysicsManager::SetCut:\n";
    text = text + "    Parameter " + cutName;
    text = text + " is not implemented.";
    TG4Globals::Warning(text);
  }  
}  
  
void TG4PhysicsManager::SetProcess(const char* controlName, Int_t controlValue)
{
// Sets the specified process control.
// ---

  fG3PhysicsManager->CheckLock();
  TG4G3Control control = fG3PhysicsManager->GetG3Control(controlName);
  if (control != kNoG3Controls)
    fG3PhysicsManager->SetProcess(control, controlValue);
  else {   
    G4String text = "TG4PhysicsManager::SetProcess:\n";
    text = text + "    Parameter " + controlName;
    text = text + " is not implemented.";
    TG4Globals::Warning(text);
  }  
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

  fParticlesManager->MapParticles();
}    


void TG4PhysicsManager::SetProcessActivation()
{
// (In)Activates built processes according
// to the setup in fControlVector.
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


AliMCProcess TG4PhysicsManager::GetMCProcess(const G4VProcess* process)
{
// Returns the AliMCProcess code of the specified G4 process.
// ---
 
  if (!process) return kPNoProcess;

  G4String name = process->GetProcessName();
  G4int code = fProcessMap.GetSecond(name);
  
  if (code == 0) return kPNoProcess;
  
  return (AliMCProcess)code; 
}

