// $Id$
// Category: run
//
// See the class description in the header file.

#include "TGeant4.h"
#include "TG4Messenger.h"
#include "TG4GeometryManager.h" 
#include "TG4PhysicsManager.h" 
#include "TG4StepManager.h" 
#include "TG4VisManager.h"
#include "TG4RunManager.h"
#include "TG4Globals.h"

#include "AliDecayer.h"

TGeant4::TGeant4(const char* name, const char* title,
                 TG4VRunConfiguration* configuration, int argc, char** argv)
  : AliMC(name, title),
    fVisManager(0)
{
  // create run manager
  fRunManager = new TG4RunManager(configuration, argc, argv);
  // add verbose level
  //G4cout << "TG4RunManager has been created." << endl;

  // create geometry manager
  fGeometryManager = new TG4GeometryManager();
  // add verbose level
  //G4cout << "TG4GeometryManager has been created." << endl;
  
  // create physics manager  
  fPhysicsManager = new TG4PhysicsManager();
  // add verbose level
  //G4cout << "TG4GeometryManager has been created." << endl;
  
  // create step manager 
  fStepManager = new TG4StepManager();
  // add verbose level
  //G4cout << "TG4StepManager has been created." << endl;
  
#ifdef G4VIS_USE
  // create visualization manager
  fVisManager = new TG4VisManager();
  fVisManager->Initialize();
#endif

  // create messenger
  fMessenger = 
    new TG4Messenger(fGeometryManager, fPhysicsManager, fStepManager);
}
    
TGeant4::TGeant4(const char* name, const char* title,
                 TG4VRunConfiguration* configuration)
  : AliMC(name, title),
    fVisManager(0)
{
  // create run manager
  fRunManager = new TG4RunManager(configuration);
  // add verbose level
  //G4cout << "TG4RunManager has been created." << endl;

  // create geometry manager
  fGeometryManager = new TG4GeometryManager();
  // add verbose level
  //G4cout << "TG4GeometryManager has been created." << endl;
  
  // create physics manager  
  fPhysicsManager = new TG4PhysicsManager();
  // add verbose level
  //G4cout << "TG4GeometryManager has been created." << endl;
  
  // create step manager 
  fStepManager = new TG4StepManager();
  // add verbose level
  //G4cout << "TG4StepManager has been created." << endl;
  
#ifdef G4VIS_USE
  // create visualization manager
  fVisManager = new TG4VisManager();
  fVisManager->Initialize();
#endif

  // create messenger
  fMessenger = 
    new TG4Messenger(fGeometryManager, fPhysicsManager, fStepManager);
}
    
TGeant4::TGeant4() {
//
}

TGeant4::TGeant4(const TGeant4& right) {
// 
  TG4Globals::Exception("TGeant4 is protected from copying.");
}

TGeant4::~TGeant4() {
//
  delete fRunManager;
  delete fGeometryManager;
  delete fPhysicsManager;
  delete fStepManager;
  // fVisManager is deleted with G4RunManager destructor
  delete fMessenger;
}

// operators

TGeant4& TGeant4::operator=(const TGeant4& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception("TGeant4 is protected from assigning.");
    
  return *this;  
}    
          

// methods for building/management of geometry
// ------------------------------------------------

void TGeant4::FinishGeometry() {
//
  fGeometryManager->Ggclos();
} 

void TGeant4::Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
  		         Float_t &dens, Float_t &radl, Float_t &absl,
		         Float_t* ubuf, Int_t& nbuf) {
//
  fGeometryManager
    ->Gfmate(imat, name, a, z, dens, radl, absl, ubuf, nbuf);
} 

void TGeant4::Material(Int_t& kmat, const char* name, Float_t a, 
                     Float_t z, Float_t dens, Float_t radl, Float_t absl,
                     Float_t* buf, Int_t nwbuf) {
//
  fGeometryManager
    ->Material(kmat, name, a, z, dens, radl, absl, buf, nwbuf); 
} 

void TGeant4::Mixture(Int_t& kmat, const char *name, Float_t *a, 
                     Float_t *z, Float_t dens, Int_t nlmat, Float_t *wmat) {
//
   fGeometryManager
     ->Mixture(kmat, name, a, z, dens, nlmat, wmat); 
} 

void TGeant4::Medium(Int_t& kmed, const char *name, Int_t nmat, 
                     Int_t isvol, Int_t ifield, Float_t fieldm, Float_t tmaxfd, 
                     Float_t stemax, Float_t deemax, Float_t epsil, 
		     Float_t stmin, Float_t* ubuf, Int_t nbuf) { 
//
  fGeometryManager
    ->Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, 
        epsil, stmin, ubuf, nbuf);
} 

void TGeant4::Matrix(Int_t& krot, Float_t thetaX, Float_t phiX, 
                     Float_t thetaY, Float_t phiY, Float_t thetaZ, 
		     Float_t phiZ) {
//		     
  fGeometryManager
    ->Matrix(krot, thetaX, phiX, thetaY, phiY, thetaZ, phiZ); 
} 

void TGeant4::Gstpar(Int_t itmed, const char *param, Float_t parval) {
//
  fGeometryManager->Gstpar(itmed, param, parval); 
}    

Int_t TGeant4::Gsvolu(const char *name, const char *shape, Int_t nmed,  
                         Float_t *upar, Int_t np)  {
//
  return fGeometryManager->Gsvolu(name, shape, nmed, upar, np); 
}
 
void TGeant4::Gsdvn(const char *name, const char *mother, Int_t ndiv, 
                        Int_t iaxis) {
//
  fGeometryManager->Gsdvn(name, mother, ndiv, iaxis); 
} 

void TGeant4::Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
                         Int_t iaxis, Float_t c0i, Int_t numed) {
//
  fGeometryManager->Gsdvn2(name, mother, ndiv, iaxis, c0i, numed); 
} 

void TGeant4::Gsdvt(const char *name, const char *mother, Float_t step, 
                        Int_t iaxis, Int_t numed, Int_t ndvmx) {
//			
  fGeometryManager->Gsdvt(name, mother, step, iaxis, numed, ndvmx); 
} 

void TGeant4::Gsdvt2(const char *name, const char *mother, Float_t step, 
                         Int_t iaxis, Float_t c0, Int_t numed, Int_t ndvmx) { 
//
  fGeometryManager->Gsdvt2(name, mother, step, iaxis, c0, numed, ndvmx); 
} 

void TGeant4::Gsord(const char *name, Int_t iax) {
//
  fGeometryManager->Gsord(name, iax); 
} 

void TGeant4::Gspos(const char *name, Int_t nr, const char *mother,  
                        Float_t x, Float_t y, Float_t z, Int_t irot, 
                        const char *konly) {
//
  fGeometryManager->Gspos(name, nr, mother, x, y, z, irot, konly); 
} 

void TGeant4::Gsposp(const char *name, Int_t nr, const char *mother,  
                         Float_t x, Float_t y, Float_t z, Int_t irot,
                         const char *konly, Float_t *upar, Int_t np)  {
//
  fGeometryManager->Gsposp(name, nr, mother, x, y, z, irot, konly, upar, np); 
} 

void TGeant4::SetCerenkov(Int_t itmed, Int_t npckov, Float_t *ppckov,
                  Float_t *absco, Float_t *effic, Float_t *rindex) {
//
  fGeometryManager->SetCerenkov(itmed, npckov, ppckov, absco, effic, rindex);
}  
    
void TGeant4::WriteEuclid(const char* fileName, const char* topVol, 
                          Int_t number, Int_t nlevel) {
//
  fGeometryManager->WriteEuclid(fileName, topVol, number, nlevel); 
} 
		               
Int_t TGeant4::VolId(const Text_t* volName) const {
//
  return fGeometryManager->VolId(volName); 
} 

const char* TGeant4::VolName(Int_t id) const {
//
  return fGeometryManager->VolName(id); 
}
 
Int_t TGeant4::NofVolumes() const {
//
  return fGeometryManager->NofVolumes(); 
} 

Int_t TGeant4::VolId2Mate(Int_t id) const {
//
  return fGeometryManager->VolId2Mate(id); 
} 



// methods for physics management
// ------------------------------------------------
 
void TGeant4::BuildPhysics() {
//
  fPhysicsManager->BuildPhysics();
}  

void TGeant4::SetCut(const char* cutName, Float_t cutValue) { 
//
  fPhysicsManager->SetCut(cutName, cutValue);
}  

void TGeant4::SetProcess(const char* flagName, Int_t flagValue) {
//
  fPhysicsManager->SetProcess(flagName, flagValue);
}  
 
Float_t TGeant4::Xsec(char* reac, Float_t energy, Int_t part, Int_t mate) {
//
  return fPhysicsManager->Xsec(reac, energy, part, mate);
}  

void TGeant4::SetExternalDecayer(AliDecayer* decayer) {
//
  TG4Globals:: Warning("TGeant4::SetExternalDecayer() is not implemented."); 
}

AliDecayer* TGeant4::Decayer() const {
//
  TG4Globals:: Warning("TGeant4::Decayer() is not implemented.");
  return 0; 
}
  

Int_t TGeant4::IdFromPDG(Int_t pdgID) const { 
//
  return fPhysicsManager->IdFromPDG(pdgID);
}  

Int_t TGeant4::PDGFromId(Int_t mcID) const {
//
  return fPhysicsManager->PDGFromId(mcID);
}  

void TGeant4::DefineParticles() { 
//
  fPhysicsManager->DefineParticles();
}  

// methods for step management
// ------------------------------------------------
// inlined (in TGeant4.icc)

// methods for visualization
// ------------------------------------------------

#ifdef G4VIS_USE
void TGeant4::DrawOneSpec(const char* name) {
//
  fVisManager->DrawOneSpec(name); 
} 

void TGeant4::Gsatt(const char* name, const char* att, Int_t val) {
//
  fVisManager->Gsatt(name, att, val); 
} 

void TGeant4::Gdraw(const char* name, Float_t theta, Float_t phi,
		        Float_t psi, Float_t u0, Float_t v0,
		        Float_t ul, Float_t vl) {
//
  fVisManager->Gdraw(name, theta, phi, psi, u0, v0, ul, vl); 
} 

#else
void TGeant4::DrawOneSpec(const char* name) {
//
  TG4Globals:: Warning("TGeant4::DrawOneSpec(): no visualization available."); 
} 

void TGeant4::Gsatt(const char* name, const char* att, Int_t val) {
//
  TG4Globals:: Warning("TGeant4::Gsatt(): no visualization available."); 
} 

void TGeant4::Gdraw(const char* p1, Float_t theta, Float_t phi,
		        Float_t psi, Float_t u0, Float_t v0,
		        Float_t ul, Float_t vl) {
//
  TG4Globals:: Warning("TGeant4::Gdraw(): no visualization available."); 
} 

#endif //G4VIS_USE

// methods for run control
// ------------------------------------------------

void TGeant4::Init() { 
//
  fRunManager->Initialize();
}  
  
void TGeant4::ProcessEvent() { 
//
  fRunManager->ProcessEvent();
}  

void TGeant4::ProcessRun(Int_t nofEvents) {
//
  fRunManager->ProcessRun(nofEvents);
}  

void TGeant4::StartGeantUI() {
//
  fRunManager->StartGeantUI();
}  

void TGeant4::StartRootUI() {
//
  fRunManager->StartRootUI();
}  

void TGeant4::ProcessGeantMacro(const char* macroName) {
//
  fRunManager->ProcessGeantMacro(macroName);
}  
  
void TGeant4::ProcessGeantCommand(const char* command) {
//
  fRunManager->ProcessGeantCommand(command);
}  
  
Int_t TGeant4::CurrentEvent() const {
//
  return fRunManager->CurrentEvent(); 
} 

// Geant3 specific methods
// !!! need to be transformed to common interface
// ------------------------------------------------
    
void TGeant4::Gdopt(const char* name, const char* value) {
//
  TG4Globals:: Warning("TGeant4::Gdopt(..) is not implemented."); 
}

void TGeant4::SetClipBox(const char *name, Float_t xmin, Float_t xmax,
		     Float_t ymin, Float_t ymax, Float_t zmin, Float_t zmax) { 
//
  TG4Globals:: Warning("TGeant4::SetClipBox(..) is not implemented."); 
}

void TGeant4::DefaultRange() { 
//
  TG4Globals:: Warning("TGeant4::DefaultRange() is not implemented."); 
}

void TGeant4::Gdhead(Int_t isel, const char* name, Float_t chrsiz) { 
//
  TG4Globals:: Warning("TGeant4::Gdhead(..) is not implemented."); 
}

void TGeant4::Gdman(Float_t u, Float_t v, const char* type) { 
//
  TG4Globals:: Warning("TGeant4::Gdman(..) is not implemented."); 
}

void TGeant4::SetColors() { 
//
  TG4Globals:: Warning("TGeant4::SetColours() is not implemented."); 
}

void TGeant4::Gtreve() { 
//
  TG4Globals:: Warning("TGeant4::Gtreve() is not implemented."); 
}

void TGeant4::GtreveRoot() { 
//
  TG4Globals:: Warning("TGeant4::GtreveRoot() is not implemented."); 
}

void TGeant4::Gckmat(Int_t itmed, char* natmed) { 
//
  TG4Globals:: Warning("TGeant4::Gckmat(..) is not implemented."); 
}

void TGeant4::InitLego() { 
//
  TG4Globals:: Warning("TGeant4::InitLego() is not implemented."); 
}

void TGeant4::Gfpart(Int_t ipart, char *name, Int_t& itrtyp,  
		     Float_t& amass, Float_t& charge, Float_t& tlife) { 
//
  TG4Globals:: Warning("TGeant4::Gfpart(..) is not implemented."); 
}

void TGeant4::Gspart(Int_t ipart, const char *name, Int_t itrtyp,  
		     Float_t amass, Float_t charge, Float_t tlife) {  
//
  TG4Globals:: Warning("TGeant4::Gspart(..) is not implemented."); 
}
