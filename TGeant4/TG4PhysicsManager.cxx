// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsManager
// -----------------------
// See the class description in the header file.

#include "TG4PhysicsManager.h"
#include "TG4ModularPhysicsList.h"
#include "TG4ParticlesManager.h"
#include "TG4G3PhysicsManager.h"
#include "TG4PhysicsConstructorGeneral.h"
#include "TG4PhysicsConstructorEM.h"
#include "TG4PhysicsConstructorMuon.h"
#include "TG4PhysicsConstructorHadron.h"
#include "TG4PhysicsConstructorIon.h"
#include "TG4PhysicsConstructorOptical.h"
#include "TG4PhysicsConstructorSpecialCuts.h"
#include "TG4PhysicsConstructorSpecialControls.h"
#include "TG4GeometryServices.h"
#include "TG4G3Cut.h"
#include "TG4G3Control.h"
#include "TG4G3Units.h"
#include "TG4Limits.h"
#include "AliDecayer.h"

#include <G4ParticleDefinition.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4VProcess.hh>
#include <G3MedTable.hh>

#include <TDatabasePDG.h>

TG4PhysicsManager* TG4PhysicsManager::fgInstance = 0;

//_____________________________________________________________________________
TG4PhysicsManager::TG4PhysicsManager(TG4ModularPhysicsList* physicsList)
  : fPhysicsList(physicsList),
    fDecayer(0),
    fSetEMPhysics(true),
    fSetMuonPhysics(true),
    fSetHadronPhysics(false),
    fSetOpticalPhysics(false),
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
  // FillProcessMap();
}

//_____________________________________________________________________________
TG4PhysicsManager::TG4PhysicsManager(){
//
  delete fDecayer;
  delete fParticlesManager;
  delete fG3PhysicsManager;
}

//_____________________________________________________________________________
TG4PhysicsManager::TG4PhysicsManager(const TG4PhysicsManager& right) {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4PhysicsManager singleton.");
}

//_____________________________________________________________________________
TG4PhysicsManager::~TG4PhysicsManager() {
//
}

// operators

//_____________________________________________________________________________
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

//_____________________________________________________________________________
void TG4PhysicsManager::FillProcessMap()
{
// Fills fProcessMCMap.
// The default G4 process names are used in the map.
// Not used - the map is filled in physics constructors
// only with processes that are really built.
// ---

  // multiple scattering
  fProcessMCMap.Add("msc",  kPMultipleScattering);
  fProcessMCMap.Add("Imsc", kPMultipleScattering);
    
  // continuous energy loss
  // !! including delta rays
  fProcessMCMap.Add("eIoni",  kPEnergyLoss);
  fProcessMCMap.Add("IeIoni", kPEnergyLoss);
  fProcessMCMap.Add("LowEnergyIoni", kPEnergyLoss);
  fProcessMCMap.Add("hIoni",  kPEnergyLoss);
  fProcessMCMap.Add("IhIoni", kPEnergyLoss);
  fProcessMCMap.Add("hLowEIoni", kPEnergyLoss);
  fProcessMCMap.Add("MuIoni", kPEnergyLoss);
  fProcessMCMap.Add("IMuIonisation", kPEnergyLoss);
  fProcessMCMap.Add("ionIoni",  kPEnergyLoss);
  fProcessMCMap.Add("ionLowEIoni",  kPEnergyLoss);
  fProcessMCMap.Add("PAIonisation",  kPEnergyLoss);
  
  // bending in mag. field
  // kPMagneticFieldL

  // particle decay
  fProcessMCMap.Add("Decay", kPDecay);
  
  // photon pair production or
  // muon direct pair production
  fProcessMCMap.Add("conv", kPPair);
  fProcessMCMap.Add("LowEnConversion", kPPair);
  fProcessMCMap.Add("MuPairProd", kPPair);
  fProcessMCMap.Add("IMuPairProduction", kPPair);

  // Compton scattering
  fProcessMCMap.Add("compt", kPCompton);
  fProcessMCMap.Add("LowEnCompton", kPCompton);
  fProcessMCMap.Add("polarCompt", kPCompton);

  // photoelectric effect
  fProcessMCMap.Add("phot", kPPhotoelectric);
  fProcessMCMap.Add("LowEnPhotoElec", kPPhotoelectric);

  // bremsstrahlung
  fProcessMCMap.Add("eBrem", kPBrem);
  fProcessMCMap.Add("IeBrem", kPBrem);
  fProcessMCMap.Add("MuBrems", kPBrem);
  fProcessMCMap.Add("IMuBremsstrahlung", kPBrem);
  fProcessMCMap.Add("LowEnBrem", kPBrem);

  // delta-ray production
  // kPDeltaRay
  // has to be distinguished from kPEnergyLoss on flight
  
  // positron annihilation
  fProcessMCMap.Add("annihil", kPAnnihilation);
  fProcessMCMap.Add("Iannihil", kPAnnihilation);

  // hadronic interaction
  // kPHadronic

  // nuclear evaporation
  // kPEvaporation
  
  // nuclear fission
  // kPNuclearFission

  // nuclear absorption
  fProcessMCMap.Add("PionMinusAbsorptionAtRest", kPNuclearAbsorption);
  fProcessMCMap.Add("PiMinusAbsorptionAtRest", kPNuclearAbsorption);
  fProcessMCMap.Add("KaonMinusAbsorption", kPNuclearAbsorption);         
  fProcessMCMap.Add("KaonMinusAbsorptionAtRest", kPNuclearAbsorption);         
  
  // antiproton annihilation
  fProcessMCMap.Add("AntiProtonAnnihilationAtRest", kPPbarAnnihilation);
  // fProcessMCMap.Add("AntiNeutronAnnihilationAtRest", not defined);

  // neutron capture    
  fProcessMCMap.Add("NeutronCaptureAtRest", kPNCapture);
  // fProcessMCMap.Add("LCapture", hadron capture not defined);

  // hadronic elastic incoherent scattering
  fProcessMCMap.Add("LElastic", kPHElastic);

  // hadronic inelastic scattering
  fProcessMCMap.Add("inelastic", kPHInhelastic);

  // muon nuclear interaction
  fProcessMCMap.Add("MuNucl", kPMuonNuclear);

  // exceeded time of flight cut
  // kPTOFlimit
  
  // nuclear photofission
  // kPPhotoFission

  // Rayleigh scattering
  fProcessMCMap.Add("Rayleigh Scattering", kPRayleigh);

  // no mechanism is active, usually at the entrance of a new volume
  fProcessMCMap.Add("Transportation", kPNull);

  // particle has fallen below energy threshold and tracking stops
  // kPStop
  
  // Cerenkov photon absorption
  fProcessMCMap.Add("Absorption", kPLightAbsorption);

  // Cerenkov photon reflection/refraction
  // kPLightScattering, kPLightReflection, kPLightRefraction
  // has to be inquired from the G4OpBoundary process

  // synchrotron radiation
  fProcessMCMap.Add("SynchrotronRadiation", kPSynchrotron);
}  

//_____________________________________________________________________________
void TG4PhysicsManager::GstparCut(G4int itmed, TG4G3Cut par, G4double parval)
{
// Sets special tracking medium parameter. 
// It is applied to all logical volumes that use the specified 
// tracking medium.
// ---

  // get medium from table
  G3MedTableEntry* medium = G3Med.get(itmed);
  if (!medium) {
    G4String text = "TG4PhysicsManager::GstparCut: \n";
    text = text + "    Medium not found."; 
    G4Exception(text);
  }  

  // get/create user limits
  TG4Limits* limits 
    = TG4GeometryServices::Instance()->GetLimits(medium->GetLimits());
    
  if (!limits) {
    limits = new TG4Limits(*fG3PhysicsManager->GetCutVector(),
                           *fG3PhysicsManager->GetControlVector());
    medium->SetLimits(limits);

    // add verbose 
    G4cout << "TG4PhysicsManager::GstparCut: new TG4Limits() for medium " 
           << itmed << " has been created." << G4endl;  
  }	   

  // add units
  if (par == kTOFMAX) parval *= TG4G3Units::Time();
  else                parval *= TG4G3Units::Energy();

  // set parameter
  limits->SetG3Cut(par, parval);
}


//_____________________________________________________________________________
void TG4PhysicsManager::GstparControl(G4int itmed, TG4G3Control par, 
                                      TG4G3ControlValue parval)
{
// Sets special tracking medium parameter. 
// It is applied to all logical volumes that use the specified 
// tracking medium.
// ---

  // get medium from table
  G3MedTableEntry* medium = G3Med.get(itmed);
  if (!medium) {
    G4String text = "TG4PhysicsManager::GstparControl: \n";
    text = text + "    Medium not found."; 
    G4Exception(text);
  }  

  // get/create user limits
  TG4Limits* limits 
    = TG4GeometryServices::Instance()->GetLimits(medium->GetLimits());

  if (!limits) {
    limits = new TG4Limits(*fG3PhysicsManager->GetCutVector(),
                           *fG3PhysicsManager->GetControlVector());
    medium->SetLimits(limits);

    // add verbose 
    G4cout << "TG4PhysicsManager::GstparControl: new TG4Limits() for medium" 
           << itmed << " has been created." << G4endl;  
  }	   
  
  // set parameter
  limits->SetG3Control(par, parval);
}

// public methods

//_____________________________________________________________________________
void  TG4PhysicsManager::Gstpar(Int_t itmed, const char *param, Float_t parval) 
{ 
// Passes the tracking medium parameter to TG4Limits.
// The tracking medium parameter is set only in case
// its value is different from the "global" physics setup.
// (If: CheckCut/ControlWithG3Defaults is used checking
//  is performed with respect to G3 default values.)
// When any cut/control parameter is set in limits
// the physics manager is locked and the physics setup
// cannot be changed.
// Applying this TG4Limits to particles and physical
// processes is still in development.
//
//  Geant3 desription:
//  ==================
//  To change the value of cut  or mechanism "CHPAR"
//  to a new value PARVAL  for tracking medium ITMED
//  The  data   structure  JTMED   contains  the   standard  tracking
//  parameters (CUTS and flags to control the physics processes)  which
//  are used  by default  for all  tracking media.   It is  possible to
//  redefine individually  with GSTPAR  any of  these parameters  for a
//  given tracking medium. 
//  ITMED     tracking medium number 
//  CHPAR     is a character string (variable name) 
//  PARVAL    must be given as a floating point.
// ---

  G4String name = TG4GeometryServices::Instance()->CutName(param); 
  TG4G3Cut cut;
  if (fG3PhysicsManager->CheckCutWithTheVector(name, parval, cut)) {
      GstparCut(itmed, cut, parval);
      fG3PhysicsManager->Lock();
  }  
  else {
    TG4G3Control control;
    TG4G3ControlValue controlValue; 
    if (fG3PhysicsManager
         ->CheckControlWithTheVector(name, parval, control, controlValue)) {
      GstparControl(itmed, control, controlValue);
      fG3PhysicsManager->Lock();
    } 
    else if (cut==kNoG3Cuts && control==kNoG3Controls) { 
      G4String text = "TG4PhysicsManager::Gstpar:";
      text = text + name;
      text = text + " parameter is not yet implemented.";
      TG4Globals::Warning(text);
    }	
  }   
} 
 
//_____________________________________________________________________________
void TG4PhysicsManager::CreatePhysicsConstructors()
{
// Creates the selected physics constructors
// and registeres them in the modular physics list.
// ---

  // general physics
  fPhysicsList->RegisterPhysics(new TG4PhysicsConstructorGeneral());

  // electromagnetic physics
  if (fSetEMPhysics) 
    fPhysicsList->RegisterPhysics(new TG4PhysicsConstructorEM());

  // muon physics
  if (fSetMuonPhysics) 
    fPhysicsList->RegisterPhysics(new TG4PhysicsConstructorMuon());

  // hadron physics
  if (fSetHadronPhysics) {
    fPhysicsList->RegisterPhysics(new TG4PhysicsConstructorIon());
    fPhysicsList->RegisterPhysics(new TG4PhysicsConstructorHadron());
  }  

  // optical physics
  if (fSetOpticalPhysics) 
    fPhysicsList->RegisterPhysics(new TG4PhysicsConstructorOptical());

  if (fSetSpecialCutsPhysics) 
    fPhysicsList->RegisterPhysics(new TG4PhysicsConstructorSpecialCuts());

  if (fSetSpecialControlsPhysics) 
    fPhysicsList->RegisterPhysics(new TG4PhysicsConstructorSpecialControls());

         // all created physics constructors are deleted
	 // in the TG4ModularPhysicsList destructor
}    

//_____________________________________________________________________________
void TG4PhysicsManager::SetCut(const char* cutName, Float_t cutValue)
{
// Sets the specified cut.
// ---

  fG3PhysicsManager->CheckLock();
  TG4G3Cut g3Cut = TG4G3CutVector::GetCut(cutName);

  if (g3Cut == kNoG3Cuts) {
    G4String text = "TG4PhysicsManager::SetCut:\n";
    text = text + "    Parameter " + cutName;
    text = text + " is not implemented.";
    TG4Globals::Warning(text);
    return;
  }  
  
  // add units
  if (g3Cut == kTOFMAX) cutValue *= TG4G3Units::Time();
  else                  cutValue *= TG4G3Units::Energy();

  fG3PhysicsManager->SetCut(g3Cut, cutValue);    
}  
  
//_____________________________________________________________________________
void TG4PhysicsManager::SetProcess(const char* controlName, Int_t value)
{
// Sets the specified process control.
// ---

  fG3PhysicsManager->CheckLock();
  TG4G3Control control = TG4G3ControlVector::GetControl(controlName);
  
  if (control != kNoG3Controls) {
    TG4G3ControlValue controlValue 
      = TG4G3ControlVector::GetControlValue(value, control);
    fG3PhysicsManager->SetProcess(control, controlValue);
  }  
  else {   
    G4String text = "TG4PhysicsManager::SetProcess:\n";
    text = text + "    Parameter " + controlName;
    text = text + " is not implemented.";
    TG4Globals::Warning(text);
  }  
}  

//_____________________________________________________________________________
Float_t TG4PhysicsManager::Xsec(char* ch, Float_t p1, Int_t i1, Int_t i2)
{
// Not yet implemented -> gives exception.
// ---

  TG4Globals::Exception(
    "TG4PhysicsManager::Xsec: not yet implemented.");

  return 0.;
}    
  
//_____________________________________________________________________________
Int_t TG4PhysicsManager::IdFromPDG(Int_t pdgID) const
{
// G4 does not use the integer particle identifiers
// Id <-> PDG is identity.
// ---

  return pdgID;
}  

//_____________________________________________________________________________
Int_t TG4PhysicsManager::PDGFromId(Int_t mcID) const
{
// G4 does not use integer particle identifiers
// Id <-> PDG is identity.
// ---

  return mcID;
}  

//_____________________________________________________________________________
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

//_____________________________________________________________________________
void TG4PhysicsManager::SetProcessActivation()
{
// (In)Activates built processes according
// to the setup in TG4G3PhysicsManager::fControlVector.
// ---

  fPhysicsList->SetProcessActivation();
}       


//_____________________________________________________________________________
AliMCProcess TG4PhysicsManager::GetMCProcess(const G4VProcess* process)
{
// Returns the AliMCProcess code of the specified G4 process.
// ---
 
  if (!process) return kPNoProcess;

  return fProcessMCMap.GetMCProcess(process);
}

//_____________________________________________________________________________
AliMCProcess TG4PhysicsManager::GetOpBoundaryStatus(const G4VProcess* process)
{
// Returns the AliMCProcess code according to the OpBoundary process
// status.
// ---
 
  if (!process) return kPNoProcess;

#ifdef TGEANT4_DEBUG
  G4OpBoundaryProcess* opBoundary
    = dynamic_cast<G4OpBoundaryProcess*>(process);
    
  if (!opBoundary) 
    TG4Globals::Exception(
      "TG4PhysicsManager::GetOpBoundaryStatus: Wrong process type.");
    return kPNoProcess;
  }
  
  return opBoundary;  
#else
  G4OpBoundaryProcess* opBoundary = (G4OpBoundaryProcess*)process;
#endif  

  switch (opBoundary->GetStatus()) {
    // reflection
    case FresnelReflection: 
    case TotalInternalReflection:
    case LambertianReflection: 
    case LobeReflection:
    case SpikeReflection: 
    case BackScattering:
       return kPLightReflection;
       ;;

    // refraction
    case FresnelRefraction: 
       return kPLightReflection;
       ;;

    // absorption
    case Absorption:
    case Detection: 
       return kPLightAbsorption;
       ;;
  }
  
  // should not happen
  return kPNoProcess;
}

