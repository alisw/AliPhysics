/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
//
// Author: I. Hrivnacova
//
// Class TFluka
// ------------------------
// See the class description in the header file.

#include "TFluka.h"

#include "TG4Globals.h"
#include "TG4GeometryManager.h" 

#include "AliDecayer.h"

TFluka::TFluka(const char* name, const char* title)
  : TVirtualMC(name, title)
{
  // create geometry manager
  fGeometryManager = new TG4GeometryManager();
  // add verbose level
  //G4cout << "TG4GeometryManager has been created." << endl;
}
    
TFluka::TFluka() {
//
}

TFluka::TFluka(const TFluka& right) {
// 
  TG4Globals::Exception("TFluka is protected from copying.");
}

TFluka::~TFluka() {
//
  delete fGeometryManager;
}

// operators

TFluka& TFluka::operator=(const TFluka& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception("TFluka is protected from assigning.");
    
  return *this;  
}    
          

// methods for building/management of geometry
// ------------------------------------------------

void TFluka::FinishGeometry() {
//
  fGeometryManager->Ggclos();
} 

void TFluka::Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
  		         Float_t &dens, Float_t &radl, Float_t &absl,
		         Float_t* ubuf, Int_t& nbuf) {
//
  fGeometryManager
    ->Gfmate(imat, name, a, z, dens, radl, absl, ubuf, nbuf);
} 

void TFluka::Material(Int_t& kmat, const char* name, Float_t a, 
                     Float_t z, Float_t dens, Float_t radl, Float_t absl,
                     Float_t* buf, Int_t nwbuf) {
//
  fGeometryManager
    ->Material(kmat, name, a, z, dens, radl, absl, buf, nwbuf); 
} 

void TFluka::Mixture(Int_t& kmat, const char *name, Float_t *a, 
                     Float_t *z, Float_t dens, Int_t nlmat, Float_t *wmat) {
//
   fGeometryManager
     ->Mixture(kmat, name, a, z, dens, nlmat, wmat); 
} 

void TFluka::Medium(Int_t& kmed, const char *name, Int_t nmat, 
                     Int_t isvol, Int_t ifield, Float_t fieldm, Float_t tmaxfd, 
                     Float_t stemax, Float_t deemax, Float_t epsil, 
		     Float_t stmin, Float_t* ubuf, Int_t nbuf) { 
//
  fGeometryManager
    ->Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, 
        epsil, stmin, ubuf, nbuf);
} 

void TFluka::Matrix(Int_t& krot, Float_t thetaX, Float_t phiX, 
                     Float_t thetaY, Float_t phiY, Float_t thetaZ, 
		     Float_t phiZ) {
//		     
  fGeometryManager
    ->Matrix(krot, thetaX, phiX, thetaY, phiY, thetaZ, phiZ); 
} 

void TFluka::Gstpar(Int_t itmed, const char *param, Float_t parval) {
//
  fGeometryManager->Gstpar(itmed, param, parval); 
}    

Int_t TFluka::Gsvolu(const char *name, const char *shape, Int_t nmed,  
                         Float_t *upar, Int_t np)  {
//
  return fGeometryManager->Gsvolu(name, shape, nmed, upar, np); 
}
 
void TFluka::Gsdvn(const char *name, const char *mother, Int_t ndiv, 
                        Int_t iaxis) {
//
  fGeometryManager->Gsdvn(name, mother, ndiv, iaxis); 
} 

void TFluka::Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
                         Int_t iaxis, Float_t c0i, Int_t numed) {
//
  fGeometryManager->Gsdvn2(name, mother, ndiv, iaxis, c0i, numed); 
} 

void TFluka::Gsdvt(const char *name, const char *mother, Float_t step, 
                        Int_t iaxis, Int_t numed, Int_t ndvmx) {
//			
  fGeometryManager->Gsdvt(name, mother, step, iaxis, numed, ndvmx); 
} 

void TFluka::Gsdvt2(const char *name, const char *mother, Float_t step, 
                         Int_t iaxis, Float_t c0, Int_t numed, Int_t ndvmx) { 
//
  fGeometryManager->Gsdvt2(name, mother, step, iaxis, c0, numed, ndvmx); 
} 

void TFluka::Gsord(const char *name, Int_t iax) {
//
  fGeometryManager->Gsord(name, iax); 
} 

void TFluka::Gspos(const char *name, Int_t nr, const char *mother,  
                        Float_t x, Float_t y, Float_t z, Int_t irot, 
                        const char *konly) {
//
  fGeometryManager->Gspos(name, nr, mother, x, y, z, irot, konly); 
} 

void TFluka::Gsposp(const char *name, Int_t nr, const char *mother,  
                         Float_t x, Float_t y, Float_t z, Int_t irot,
                         const char *konly, Float_t *upar, Int_t np)  {
//
  fGeometryManager->Gsposp(name, nr, mother, x, y, z, irot, konly, upar, np); 
} 

void TFluka::SetCerenkov(Int_t itmed, Int_t npckov, Float_t *ppckov,
                  Float_t *absco, Float_t *effic, Float_t *rindex) {
//
  fGeometryManager->SetCerenkov(itmed, npckov, ppckov, absco, effic, rindex);
}  
    
void TFluka::WriteEuclid(const char* fileName, const char* topVol, 
                          Int_t number, Int_t nlevel) {
//
  fGeometryManager->WriteEuclid(fileName, topVol, number, nlevel); 
} 
		               
//======================================================================
//
// NOT IMPLEMENTED METHODS
//
//======================================================================

Int_t TFluka::VolId2Mate(Int_t id) const {
//
  TG4Globals:: Warning("TFluka::VolId2Mate(..) is not implemented."); 
  return 0; 
} 

Int_t TFluka::VolId(const Text_t* volName) const {
//
  TG4Globals:: Warning("TFluka::VolId(..) is not implemented."); 
  return 0; 
} 

const char* TFluka::VolName(Int_t id) const {
//
  TG4Globals:: Warning("TFluka::VolName(..) is not implemented."); 
  return 0; 
}
 
Int_t TFluka::NofVolumes() const {
//
  TG4Globals:: Warning("TFluka::NofVolumes(..) is not implemented."); 
  return 0; 
} 


// methods for physics management
// ------------------------------------------------
 
void TFluka::BuildPhysics() {
//
  TG4Globals:: Warning("TFluka::BuildPhysics(..) is not implemented."); 
}  

void TFluka::SetCut(const char* cutName, Float_t cutValue) { 
//
  TG4Globals:: Warning("TFluka::SetCut(..) is not implemented."); 
}  

void TFluka::SetProcess(const char* flagName, Int_t flagValue) {
//
  TG4Globals:: Warning("TFluka::SetProcess(..) is not implemented."); 
}  
 
Float_t TFluka::Xsec(char* reac, Float_t energy, Int_t part, Int_t mate) {
//
  TG4Globals:: Warning("TFluka::Xsec(..) is not implemented."); 
  return 0;
}  

void TFluka::SetExternalDecayer(AliDecayer* decayer) {
//
  TG4Globals:: Warning("TFluka::SetExternalDecayer(..) is not implemented."); 
}

AliDecayer* TFluka::Decayer() const {
//
  TG4Globals:: Warning("TFluka::Decayer(..) is not implemented."); 
  return 0; 
}
  

Int_t TFluka::IdFromPDG(Int_t pdgID) const { 
//
  TG4Globals:: Warning("TFluka::IdFromPDG(..) is not implemented."); 
  return 0;
}  

Int_t TFluka::PDGFromId(Int_t mcID) const {
//
  TG4Globals:: Warning("TFluka::PDGFromId(..) is not implemented."); 
  return 0;
}  

void TFluka::DefineParticles() { 
//
  TG4Globals:: Warning("TFluka::DefineParticles(..) is not implemented."); 
}  

// methods for step management
// ------------------------------------------------

void TFluka::StopTrack()
{ 
  TG4Globals:: Warning("TFluka::StopTrack(..) is not implemented."); 
} 

void TFluka::StopEvent()   
{ 
  TG4Globals:: Warning("TFluka::StopEvent(..) is not implemented."); 
} 

void TFluka::SetMaxStep(Float_t step)
{ 
  TG4Globals:: Warning("TFluka::SetMaxStep(..) is not implemented."); 
} 

void TFluka::SetMaxNStep(Int_t number)
{ 
  TG4Globals:: Warning("TFluka::SetMaxNStep(..) is not implemented."); 
} 

void TFluka::SetUserDecay(Int_t number)
{ 
  TG4Globals:: Warning("TFluka::SetUserDecay(..) is not implemented."); 
} 

Int_t TFluka::CurrentVolID(Int_t& copyNo) const
{ 
  TG4Globals:: Warning("TFluka::CurrentVolID(..) is not implemented."); 
  return 0;
} 

Int_t TFluka::CurrentVolOffID(Int_t off, Int_t& copyNo) const
{ 
  TG4Globals:: Warning("TFluka::CurrentVolOffID(..) is not implemented."); 
  return 0;
} 

const char* TFluka::CurrentVolName() const
{ 
  TG4Globals:: Warning("TFluka::CurrentVolName(..) is not implemented."); 
  return 0;
} 

const char* TFluka::CurrentVolOffName(Int_t off) const
{ 
  TG4Globals:: Warning("TFluka::CurrentVolOffName(..) is not implemented."); 
  return 0;
} 

Int_t TFluka::CurrentMaterial(Float_t &a, Float_t &z, 
                    Float_t &dens, Float_t &radl, Float_t &absl) const  
{ 
  TG4Globals:: Warning("TFluka::CurrentMaterial(..) is not implemented."); 
  return 0;
} 

void TFluka::Gmtod(Float_t* xm, Float_t* xd, Int_t iflag)
{ 
  TG4Globals:: Warning("TFluka::Gmtod(..) is not implemented."); 
} 

void TFluka::Gdtom(Float_t* xd, Float_t* xm, Int_t iflag)
{ 
  TG4Globals:: Warning("TFluka::Gdtom(..) is not implemented."); 
} 
    
Float_t TFluka::MaxStep() const
{ 
  TG4Globals:: Warning("TFluka::MaxStep(..) is not implemented."); 
  return 0;
} 

Int_t TFluka::GetMaxNStep() const
{ 
  TG4Globals:: Warning("TFluka::GetMaxNStep(..) is not implemented."); 
  return 0;
} 

Int_t TFluka::GetMedium() const
{ 
  TG4Globals:: Warning("TFluka::GetMedium(..) is not implemented."); 
  return 0;
} 

void TFluka::TrackPosition(TLorentzVector& position) const
{ 
  TG4Globals:: Warning("TFluka::TrackPosition(..) is not implemented."); 
} 

void TFluka::TrackMomentum(TLorentzVector& momentum) const
{ 
  TG4Globals:: Warning("TFluka::TrackMomentum(..) is not implemented."); 
} 

void TFluka::TrackVertexPosition(TLorentzVector& position) const
{ 
  TG4Globals:: Warning("TFluka::TrackVertexPosition(..) is not implemented."); 
} 

void TFluka::TrackVertexMomentum(TLorentzVector& momentum) const
{ 
  TG4Globals:: Warning("TFluka::TrackVertexMomentum(..) is not implemented."); 
} 

Float_t TFluka::TrackStep() const
{ 
  TG4Globals:: Warning("TFluka::TrackStep(..) is not implemented."); 
  return 0;
} 

Float_t TFluka::TrackLength() const 
{ 
  TG4Globals:: Warning("TFluka::TrackLength(..) is not implemented."); 
  return 0; 
} 

Float_t TFluka::TrackTime() const
{ 
  TG4Globals:: Warning("TFluka::TrackTime(..) is not implemented."); 
  return 0; 
} 

Float_t TFluka::Edep() const
{ 
  TG4Globals:: Warning("TFluka::Edep(..) is not implemented."); 
  return 0; 
} 

Int_t TFluka::TrackPid() const
{ 
  TG4Globals:: Warning("TFluka::TrackPid(..) is not implemented."); 
  return 0; 
} 

Float_t TFluka::TrackCharge() const
{ 
  TG4Globals:: Warning("TFluka::TrackCharge(..) is not implemented."); 
  return 0; 
} 

Float_t TFluka::TrackMass() const
{ 
  TG4Globals:: Warning("TFluka::TrackMass(..) is not implemented."); 
  return 0; 
} 

Float_t TFluka::Etot() const
{ 
  TG4Globals:: Warning("TFluka::Etot(..) is not implemented."); 
  return 0; 
} 

Bool_t  TFluka::IsTrackInside() const
{   
  TG4Globals:: Warning("TFluka::IsTrackInside(..) is not implemented."); 
  return false; 
} 

Bool_t  TFluka::IsTrackEntering() const
{   
  TG4Globals:: Warning("TFluka::IsTrackEntering(..) is not implemented."); 
  return false; 
} 

Bool_t  TFluka::IsTrackExiting() const
{   
  TG4Globals:: Warning("TFluka::IsTrackExiting(..) is not implemented."); 
  return false; 
} 

Bool_t  TFluka::IsTrackOut() const
{   
  TG4Globals:: Warning("TFluka::IsTrackOut(..) is not implemented."); 
  return false; 
} 

Bool_t  TFluka::IsTrackDisappeared() const
{   
  TG4Globals:: Warning("TFluka::IsTrackDisappeared(..) is not implemented."); 
  return false; 
} 

Bool_t  TFluka::IsTrackStop() const
{  
  TG4Globals:: Warning("TFluka::IsTrackStop(..) is not implemented."); 
  return false; 
} 

Bool_t  TFluka::IsTrackAlive() const
{   
  TG4Globals:: Warning("TFluka::IsTrackAlive(..) is not implemented."); 
  return false; 
} 

Bool_t TFluka::IsNewTrack() const
{
  TG4Globals:: Warning("TFluka::IsNewTrack(..) is not implemented."); 
  return false;
}  

Int_t TFluka::NSecondaries() const
{ 
  TG4Globals:: Warning("TFluka::NSecondaries(..) is not implemented."); 
  return 0; 
} 

void  TFluka::GetSecondary(Int_t isec, Int_t& particleId, 
                         TLorentzVector& position, TLorentzVector& momentum)
{ 
  TG4Globals:: Warning("TFluka::GetSecondary(..) is not implemented."); 
} 

TMCProcess TFluka::ProdProcess(Int_t isec) const 
{ 
  TG4Globals:: Warning("TFluka::ProdProcess(..) is not implemented."); 
  return kPNoProcess; 
} 

Int_t TFluka::StepProcesses(TArrayI &proc) const
{ 
  TG4Globals:: Warning("TFluka::StepProcesses(..) is not implemented."); 
  return 0; 
} 

// methods for visualization
// ------------------------------------------------

void TFluka::DrawOneSpec(const char* name) {
//
  TG4Globals:: Warning("TFluka::DrawOneSpec(): no visualization available."); 
} 

void TFluka::Gsatt(const char* name, const char* att, Int_t val) {
//
  TG4Globals:: Warning("TFluka::Gsatt(): no visualization available."); 
} 

void TFluka::Gdraw(const char* p1, Float_t theta, Float_t phi,
		        Float_t psi, Float_t u0, Float_t v0,
		        Float_t ul, Float_t vl) {
//
  TG4Globals:: Warning("TFluka::Gdraw(): no visualization available."); 
} 


// methods for run control
// ------------------------------------------------

void TFluka::Init() { 
//
  TG4Globals:: Warning("TFluka::Init(..) is not implemented."); 
}  
  
void TFluka::ProcessEvent() { 
//
  TG4Globals:: Warning("TFluka::ProcessEvent(..) is not implemented."); 
}  

void TFluka::ProcessRun(Int_t nofEvents) {
//
  TG4Globals:: Warning("TFluka::ProcessRun(..) is not implemented."); 
}  

Int_t TFluka::CurrentEvent() const {
//
  TG4Globals:: Warning("TFluka::CurrentEvent(..) is not implemented."); 
  return 0;
} 

// Geant3 specific methods
// !!! need to be transformed to common interface
// ------------------------------------------------
    
void TFluka::Gdopt(const char* name, const char* value) {
//
  TG4Globals:: Warning("TFluka::Gdopt(..) is not implemented."); 
}

void TFluka::SetClipBox(const char *name, Float_t xmin, Float_t xmax,
		     Float_t ymin, Float_t ymax, Float_t zmin, Float_t zmax) { 
//
  TG4Globals:: Warning("TFluka::SetClipBox(..) is not implemented."); 
}

void TFluka::DefaultRange() { 
//
  TG4Globals:: Warning("TFluka::DefaultRange() is not implemented."); 
}

void TFluka::Gdhead(Int_t isel, const char* name, Float_t chrsiz) { 
//
  TG4Globals:: Warning("TFluka::Gdhead(..) is not implemented."); 
}

void TFluka::Gdman(Float_t u, Float_t v, const char* type) { 
//
  TG4Globals:: Warning("TFluka::Gdman(..) is not implemented."); 
}

void TFluka::SetColors() { 
//
  TG4Globals:: Warning("TFluka::SetColours() is not implemented."); 
}

void TFluka::Gtreve() { 
//
  TG4Globals:: Warning("TFluka::Gtreve() is not implemented."); 
}

void TFluka::GtreveRoot() { 
//
  TG4Globals:: Warning("TFluka::GtreveRoot() is not implemented."); 
}

void TFluka::Gckmat(Int_t itmed, char* natmed) { 
//
  TG4Globals:: Warning("TFluka::Gckmat(..) is not implemented."); 
}

void TFluka::InitLego() { 
//
  TG4Globals:: Warning("TFluka::InitLego() is not implemented."); 
}

void TFluka::Gfpart(Int_t ipart, char *name, Int_t& itrtyp,  
		     Float_t& amass, Float_t& charge, Float_t& tlife) { 
//
  TG4Globals:: Warning("TFluka::Gfpart(..) is not implemented."); 
}

void TFluka::Gspart(Int_t ipart, const char *name, Int_t itrtyp,  
		     Float_t amass, Float_t charge, Float_t tlife) {  
//
  TG4Globals:: Warning("TFluka::Gspart(..) is not implemented."); 
}
