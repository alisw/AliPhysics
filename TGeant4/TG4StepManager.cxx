// $Id$
// Category: event
//
// See the class description in the header file.

#include "TG4StepManager.h"
#include "TG4GeometryServices.h"
#include "TG4ParticlesManager.h"
#include "TG4PhysicsManager.h"
#include "TG4VSensitiveDetector.h"
#include "TG4Limits.h"
#include "TG4Globals.h"
#include "TG4G3Units.h"

#include <G4SteppingManager.hh>
#include <G4UserLimits.hh>
#include <G4UImanager.hh>
#include <G4AffineTransform.hh>
#include <G4TransportationManager.hh>
#include <G4Navigator.hh>
#include <G4VProcess.hh>
#include <G4ProcessManager.hh>
#include <G4ProcessVector.hh>

#include <Randomize.hh>
#include <TLorentzVector.h>

TG4StepManager* TG4StepManager::fgInstance = 0;

TG4StepManager::TG4StepManager() 
  : fTrack(0),
    fStep(0),
    fStepStatus(kNormalStep),
    fSteppingManager(0)
{
// 
  if (fgInstance) {
    TG4Globals::Exception(
      "TG4StepManager: attempt to create two instances of singleton.");
  }
      
  fgInstance = this;  
}

TG4StepManager::TG4StepManager(const TG4StepManager& right) {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4StepManager singleton.");
}

TG4StepManager::~TG4StepManager() {
//
}

// operators

TG4StepManager& TG4StepManager::operator=(const TG4StepManager& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4StepManager singleton.");
    
  return *this;  
}    
          
// private methods

void TG4StepManager::CheckTrack() const
{
// Gives exception in case the track is not defined.
// ---

  if (!fTrack) 
    TG4Globals::Exception("TG4StepManager: Track is not defined.");
}     


void TG4StepManager::CheckStep(const G4String& method) const
{
// Gives exception in case the step is not defined.
// ---

  if (!fStep) {
    G4String text = "TG4StepManager::";
    text = text + method + ": Step is not defined.";
    TG4Globals::Exception(text);
  }
}     


void TG4StepManager::CheckSteppingManager() const
{
// Gives exception in case the step is not defined.
// ---

  if (!fSteppingManager) 
    TG4Globals::Exception("TG4StepManager: Stepping manager is not defined.");
}     


void TG4StepManager::SetTLorentzVector(G4ThreeVector xyz, G4double t, 
                                       TLorentzVector& lv) const				       
{
// Fills TLorentzVector with G4ThreeVector and G4double.
// ---

   lv[0] = xyz.x();  				       
   lv[1] = xyz.y();  				       
   lv[2] = xyz.z();  				       
   lv[3] = t;
}     				       

G4VPhysicalVolume* TG4StepManager::GetCurrentOffPhysicalVolume(G4int off) const 
{
// Returns the physical volume of the off-th mother's
// of the current volume.
// ---
 
  G4VPhysicalVolume* physVolume = GetCurrentPhysicalVolume(); 

  G4VPhysicalVolume* mother = physVolume; 

  Int_t level = off;
  while (level > 0) { 
    if (mother) mother = mother->GetMother();
    level--;
  }
    
  if (!mother) {
    G4String text = "TG4StepManager::CurrentVolOff: \n";
    text = text + "    Volume ";
    text = text + physVolume->GetName();
    text = text + " has not defined mother in the required level.";
    TG4Globals::Warning(text);
  }  

  return mother;  
}     

// public methods

void TG4StepManager::StopTrack()
{
// Stops the current track and skips to the next.
// ?? do we want to invoke rest processes?
// ?? do we want to stop secondaries too?
//   possible "stop" track status from G4:
//   fStopButAlive,      // Invoke active rest physics processes and
//                       // and kill the current track afterward
//   fStopAndKill,       // Kill the current track
//   fKillTrackAndSecondaries,
//                       // Kill the current track and also associated
//                       // secondaries.
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif  
  
  fTrack->SetTrackStatus(fStopAndKill);
  // fTrack->SetTrackStatus(fStopButAlive);
  // fTrack->SetTrackStatus(fKillTrackAndSecondaries);
}

void TG4StepManager::StopEvent()
{
// Aborts the current event processing.
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif  
  
  fTrack->SetTrackStatus(fKillTrackAndSecondaries);
          //StopTrack();   // cannot be used as it keeps secondaries
  G4UImanager::GetUIpointer()->ApplyCommand("/event/abort");
  G4UImanager::GetUIpointer()->ApplyCommand("/alStacking/clearStack");
}

void TG4StepManager::SetMaxStep(Float_t step)
{
// Maximum step allowed in the current logical volume.
// ---

  G4LogicalVolume* curLogVolume 
    = GetCurrentPhysicalVolume()->GetLogicalVolume();
  G4UserLimits* userLimits 
    = curLogVolume->GetUserLimits();

  if (userLimits == 0) {
    // create new limits
    userLimits = new TG4Limits();
    
    // set limits to all logical volumes
    // corresponding to the current "G3" volume 
    TG4GeometryServices* geometryServices = TG4GeometryServices::Instance();
    G4int nofLV = geometryServices->SetUserLimits(userLimits, curLogVolume);
  }

  // set max step
  userLimits->SetMaxAllowedStep(step*TG4G3Units::Length()); 
}

void TG4StepManager::SetMaxNStep(Int_t maxNofSteps)
{
// Not yet implemented.
// ---

  TG4Globals::Warning(
    "TG4StepManager::SetMaxNStep(..) is not yet implemented.");
}

void TG4StepManager::SetUserDecay(Int_t pdg)
{
// Not yet implemented.
// ---

  TG4Globals::Exception(
    "TG4StepManager::SetUserDecay(..) is not yet implemented.");
}

G4VPhysicalVolume* TG4StepManager::GetCurrentPhysicalVolume() const 
{
// Returns the current physical volume.
// According to fStepStatus the volume from track vertex,
// pre step point or post step point is returned.
// ---

  G4VPhysicalVolume* physVolume; 
  if (fStepStatus == kNormalStep) {

#ifdef TGEANT4_DEBUG
    CheckStep("GetCurrentPhysicalVolume");
#endif    

    physVolume = fStep->GetPreStepPoint()->GetPhysicalVolume();
  }  
  else if (fStepStatus == kBoundary) {

#ifdef TGEANT4_DEBUG
    CheckStep("GetCurrentPhysicalVolume");
#endif 

    physVolume = fStep->GetPostStepPoint()->GetPhysicalVolume();
  }  
  else {

#ifdef TGEANT4_DEBUG
    CheckTrack();
#endif 

    G4ThreeVector position = fTrack->GetPosition();
    G4Navigator* navigator =
      G4TransportationManager::GetTransportationManager()->
        GetNavigatorForTracking();
    physVolume
     = navigator->LocateGlobalPointAndSetup(position);  
  }   
    
  return physVolume;
}     

Int_t TG4StepManager::CurrentVolID(Int_t& copyNo) const
{
// Returns the current sensitive detector ID
// and the copy number of the current physical volume.
// ---

  G4VPhysicalVolume* physVolume = GetCurrentPhysicalVolume(); 
  copyNo = physVolume->GetCopyNo() + 1;

  // sensitive detector ID
  TG4GeometryServices* geometryServices = TG4GeometryServices::Instance();
  return geometryServices->GetVolumeID(physVolume->GetLogicalVolume());
} 

Int_t TG4StepManager::CurrentVolOffID(Int_t off, Int_t&  copyNo) const
{ 
// Returns the off-th mother's of the current volume
// the sensitive detector ID and the copy number.
// ---

  if (off == 0) return CurrentVolID(copyNo);

  G4VPhysicalVolume* mother = GetCurrentOffPhysicalVolume(off); 

  if (mother) {
    copyNo = mother->GetCopyNo() + 1;

    // sensitive detector ID
    TG4GeometryServices* geometryServices = TG4GeometryServices::Instance();
    return geometryServices->GetVolumeID(mother->GetLogicalVolume());
  }
  else {
    copyNo = 0;
    return 0;
  }  
}

const char* TG4StepManager::CurrentVolName() const
{
// Returns the current physical volume name.
// ---

  return GetCurrentPhysicalVolume()->GetName();
}

const char* TG4StepManager::CurrentVolOffName(Int_t off) const
{ 
// Returns the off-th mother's physical volume name.
// ---

  if (off == 0) return CurrentVolName();

  G4VPhysicalVolume* mother = GetCurrentOffPhysicalVolume(off); 

  if (mother) 
    return mother->GetName();
  else 
    return 0;
}

Int_t TG4StepManager::CurrentMaterial(Float_t &a, Float_t &z, Float_t &dens, 
                          Float_t &radl, Float_t &absl) const
{
// Returns the parameters of the current material during transport
// the return value is the number of elements in the mixture.
// ---

  G4VPhysicalVolume* physVolume = GetCurrentPhysicalVolume(); 
    
  G4Material* material 
    = physVolume->GetLogicalVolume()->GetMaterial();

  if (material) {
    G4int nofElements = material->GetNumberOfElements();
    TG4GeometryServices* geometryServices = TG4GeometryServices::Instance();
    a = geometryServices->GetEffA(material);
    z = geometryServices->GetEffZ(material);
      
    // density 
    dens = material->GetDensity();
    dens /= TG4G3Units::MassDensity();      
      
    // radiation length
    radl = material->GetRadlen();
    radl /= TG4G3Units::Length();
      
    absl = 0.;  // this parameter is not defined in Geant4
    return nofElements;
  }
  else {
    TG4Globals::Exception(
     "TG4StepManager::CurrentMaterial(..): material is not defined.");
    return 0;
  }
}

void TG4StepManager::Gmtod(Float_t* xm, Float_t* xd, Int_t iflag) 
{ 
// Transforms a position from the world reference frame
// to the current volume reference frame.
//
//  Geant3 desription:
//  ==================
//       Computes coordinates XD (in DRS) 
//       from known coordinates XM in MRS 
//       The local reference system can be initialized by
//         - the tracking routines and GMTOD used in GUSTEP
//         - a call to GMEDIA(XM,NUMED)
//         - a call to GLVOLU(NLEVEL,NAMES,NUMBER,IER) 
//             (inverse routine is GDTOM) 
//
//        If IFLAG=1  convert coordinates 
//           IFLAG=2  convert direction cosinus
//
// ---

  G4AffineTransform affineTransform;

  if (fStepStatus == kVertex) {
    G4Navigator* navigator =
      G4TransportationManager::GetTransportationManager()->
        GetNavigatorForTracking();
	
    affineTransform = navigator->GetGlobalToLocalTransform();
  }
  else {

#ifdef TGEANT4_DEBUG
    CheckStep("Gmtod");
#endif
 
    affineTransform
      = fStep->GetPreStepPoint()->GetTouchable()->GetHistory()
        ->GetTopTransform();
  }	

  G4ThreeVector theGlobalPoint(xm[0],xm[1],xm[2]); 
  G4ThreeVector theLocalPoint;
  if(iflag == 1) 
       theLocalPoint = affineTransform.TransformPoint(theGlobalPoint);
  else if ( iflag == 2)
       theLocalPoint = affineTransform.TransformAxis(theGlobalPoint);
  else 
    TG4Globals::Exception(
      "TG4StepManager::Gmtod(..,iflag): iflag is not in 1..2");

  xd[0] = theLocalPoint.x();
  xd[1] = theLocalPoint.y();
  xd[2] = theLocalPoint.z();
} 
 
void TG4StepManager::Gdtom(Float_t* xd, Float_t* xm, Int_t iflag) 
{ 
// Transforms a position from the current volume reference frame
// to the world reference frame.
//
//  Geant3 desription:
//  ==================
//  Computes coordinates XM (Master Reference System
//  knowing the coordinates XD (Detector Ref System)
//  The local reference system can be initialized by
//    - the tracking routines and GDTOM used in GUSTEP
//    - a call to GSCMED(NLEVEL,NAMES,NUMBER)
//        (inverse routine is GMTOD)
// 
//   If IFLAG=1  convert coordinates
//      IFLAG=2  convert direction cosinus
//
// ---

  G4AffineTransform affineTransform;

  if (fStepStatus == kVertex) {
    G4Navigator* navigator =
      G4TransportationManager::GetTransportationManager()->
        GetNavigatorForTracking();
	
    affineTransform = navigator->GetLocalToGlobalTransform();
  }
  else {

#ifdef TGEANT4_DEBUG
    CheckStep("Gdtom");
#endif

    // check this
     
    affineTransform
      = fStep->GetPreStepPoint()->GetTouchable()->GetHistory()
        ->GetTopTransform().Inverse();
  }	
  
  G4ThreeVector theLocalPoint(xd[0],xd[1],xd[2]); 
  G4ThreeVector theGlobalPoint;
  if(iflag == 1)
       theGlobalPoint = affineTransform.TransformPoint(theLocalPoint);
  else if( iflag == 2)
       theGlobalPoint = affineTransform.TransformAxis(theLocalPoint);
  else 
    TG4Globals::Warning(
      "TG4StepManager::Gdtom(...,iflag): iflag is not in 1..2");

  xm[0] = theGlobalPoint.x();
  xm[1] = theGlobalPoint.y();
  xm[2] = theGlobalPoint.z();
} 
 
Float_t TG4StepManager::MaxStep() const
{   
// Returns maximum step allowed in the current logical volume
// by User Limits.
// ---

  G4LogicalVolume* curLogVolume
    = GetCurrentPhysicalVolume()->GetLogicalVolume();

  // check this
  G4UserLimits* userLimits 
    = curLogVolume->GetUserLimits();

  G4double maxStep;
  if (userLimits == 0) { 
    G4String text = "User Limits are not defined for log volume ";
    text = text + curLogVolume->GetName();
    TG4Globals::Warning(text);
    return FLT_MAX;
  }
  else { 
    const G4Track& trackRef = *(fTrack);
    maxStep = userLimits->GetMaxAllowedStep(trackRef); 
    maxStep /= TG4G3Units::Length(); 
    return maxStep;
  }  
}

Int_t TG4StepManager::GetMaxNStep() const
{   
// Not yet implemented.
// ---

  TG4Globals::Warning(
    "Method GetMaxNStep is not yet implemented in TG4StepManager.");
  return 0; 
}

void TG4StepManager::TrackPosition(TLorentzVector& position) const
{ 
// Current particle position (in the world reference frame)
// and the local time since the current track is created
// (position of the PostStepPoint).
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif

  // get position
  // check if this is == to PostStepPoint position !!
  G4ThreeVector positionVector = fTrack->GetPosition();
  positionVector *= 1./(TG4G3Units::Length());   
     
  // local time   
  G4double time = fTrack->GetLocalTime();
  time /= TG4G3Units::Time();
    
  SetTLorentzVector(positionVector, time, position);
}

Int_t TG4StepManager::GetMedium() const
{   
// Returns the second index of the current material (corresponding to
// G3 tracking medium index).
// --- 

  // current material
  G4Material* curMaterial
    = GetCurrentPhysicalVolume()->GetLogicalVolume()->GetMaterial();

  // medium index  
  TG4GeometryServices* geometryServices = TG4GeometryServices::Instance();
  return geometryServices->GetMediumId(curMaterial);
}

void TG4StepManager::TrackMomentum(TLorentzVector& momentum) const
{  
// Current particle "momentum" (px, py, pz, Etot).
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif

  G4ThreeVector momentumVector = fTrack->GetMomentum(); 
  momentumVector *= 1./(TG4G3Units::Energy());   

  G4double energy = fTrack->GetDynamicParticle()->GetTotalEnergy();
  energy /= TG4G3Units::Energy();  

  SetTLorentzVector(momentumVector, energy, momentum);
}

void TG4StepManager::TrackVertexPosition(TLorentzVector& position) const
{ 
// The vertex particle position (in the world reference frame)
// and the local time since the current track is created.
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif

  // position
  G4ThreeVector positionVector = fTrack->GetVertexPosition();
  positionVector *= 1./(TG4G3Units::Length());   
     
  // local time 
  // to be checked  
  G4double time = fTrack->GetLocalTime();
  time /= TG4G3Units::Time();
      
  SetTLorentzVector(positionVector, time, position);
}

void TG4StepManager::TrackVertexMomentum(TLorentzVector& momentum) const
{  
// The vertex particle "momentum" (px, py, pz, Ekin)
// to do: change Ekin -> Etot 
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif

  G4ThreeVector momentumVector = fTrack->GetVertexMomentumDirection(); 
  momentumVector *= 1./(TG4G3Units::Energy());   

  G4double energy = fTrack->GetVertexKineticEnergy();
  energy /= TG4G3Units::Energy();  

  SetTLorentzVector(momentumVector, energy, momentum);
}

Float_t TG4StepManager::TrackStep() const
{   
// Returns the current step length.
// ---

  G4double length;
  if (fStepStatus == kNormalStep) {

#ifdef TGEANT4_DEBUG
    CheckStep("TrackStep");    
#endif

    length = fStep->GetStepLength();
    length /= TG4G3Units::Length();
  }  
  else 
    length = 0;

  return length;
}

Float_t TG4StepManager::TrackLength() const
{
// Returns the length of the current track from its origin.
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif

  G4double length = fTrack->GetTrackLength();
  length /= TG4G3Units::Length();
  return length;
}

Float_t TG4StepManager::TrackTime() const
{
// Returns the local time since the current track is created.
// Comment:
// in Geant4: there is also defined proper time as
// the proper time of the dynamical particle of the current track.
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif
  
  G4double time = fTrack->GetLocalTime();
  time /= TG4G3Units::Time();
  return time;
}

Float_t TG4StepManager::Edep() const
{   
// Returns total energy deposit in this step.
// ---

  G4double energyDeposit;
  if (fStepStatus == kNormalStep) {

#ifdef TGEANT4_DEBUG
    CheckStep("Edep");
#endif

    energyDeposit = fStep->GetTotalEnergyDeposit();
    energyDeposit /= TG4G3Units::Energy();
  }
  else   
    energyDeposit = 0;

  return energyDeposit;
}

Int_t TG4StepManager::TrackPid() const
{   
// Returns the current particle PDG encoding.
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif

  G4ParticleDefinition* particle
    = fTrack->GetDynamicParticle()->GetDefinition();
    
  // ask TG4ParticlesManager to get PDG encoding 
  // (in order to get PDG from extended TDatabasePDG
  // in case the standard PDG code is not defined)
  G4int pdgEncoding 
    = TG4ParticlesManager::Instance()->GetPDGEncodingFast(particle);

  return pdgEncoding;
}

Float_t TG4StepManager::TrackCharge() const
{   
// Returns the current particle charge.
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif

  G4double charge
    = fTrack->GetDynamicParticle()->GetDefinition()
      ->GetPDGCharge();
  charge /= TG4G3Units::Charge();	
  return charge;
}

Float_t TG4StepManager::TrackMass() const
{   
// Returns current particle rest mass.
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif

  G4double mass
    = fTrack->GetDynamicParticle()->GetDefinition()
      ->GetPDGMass();
  mass /= TG4G3Units::Mass();	
  return mass;
}

Float_t TG4StepManager::Etot() const
{   
// Returns total energy of the current particle.
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif

  G4double energy
    = fTrack->GetDynamicParticle()->GetTotalEnergy();
  energy /= TG4G3Units::Energy();  
  return energy;
}

Bool_t TG4StepManager::IsTrackInside() const
{   
// Returns true if particle does not cross geometrical boundary
// and is not in vertex.
// ---

  if (fStepStatus == kNormalStep  && !(IsTrackExiting()) ) {
    // track is always inside during a normal step
    return true; 
  }    

  return false;    
}

Bool_t TG4StepManager::IsTrackEntering() const
{   
// Returns true if particle cross a geometrical boundary
// or is in vertex.
// ---

  if (fStepStatus != kNormalStep) {
    // track is entering during a vertex or boundary step
    return true;  
  }
  
  return false;  
}

Bool_t TG4StepManager::IsTrackExiting() const
{   
// Returns true if particle cross a geometrical boundary.
// ---

  if (fStepStatus == kNormalStep) {

#ifdef TGEANT4_DEBUG
    CheckStep("IsTrackExiting");
#endif    

    if (fStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) 
       return true;  
  }
  
  return false;  
}

Bool_t TG4StepManager::IsTrackOut() const
{   
// Returns true if particle cross the world boundary
// at post-step point.
// ---

  if (fStepStatus == kVertex) return false;

#ifdef TGEANT4_DEBUG
  CheckStep("IsTrackCut");
#endif

  // check
  G4StepStatus status
    = fStep->GetPostStepPoint()->GetStepStatus();
  if (status == fWorldBoundary)
    return true; 
  else
    return false;
}

Bool_t TG4StepManager::IsTrackStop() const
{   
// Returns true if particle has stopped 
// or has been killed, suspended or postponed to next event.
//
// Possible track status from G4:
//   fAlive,             // Continue the tracking
//   fStopButAlive,      // Invoke active rest physics processes and
//                            // and kill the current track afterward
//   fStopAndKill,       // Kill the current track
//   fKillTrackAndSecondaries,
//                       // Kill the current track and also associated
//                       // secondaries.
//   fSuspend,           // Suspend the current track
//   fPostponeToNextEvent
//                       // Postpones the tracking of thecurrent track 
//                       // to the next event.
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif

  // check
  G4TrackStatus status
     = fTrack->GetTrackStatus();
  if ((status == fStopAndKill) ||  
      (status == fKillTrackAndSecondaries) ||
      (status == fSuspend) ||
      (status == fPostponeToNextEvent)) {
    return true; 
  }
  else
    return false; 
}

Bool_t TG4StepManager::IsTrackDisappeared() const
{ 
// Returns true if particle has disappeared 
// (due to any physical process)
// or has been killed, suspended or postponed to next event.
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif

  // check
  G4TrackStatus status
     = fTrack->GetTrackStatus();
  if ((status == fStopButAlive) ||  
      (status == fKillTrackAndSecondaries) ||
      (status == fSuspend) ||
      (status == fPostponeToNextEvent)) {
    return true; 
  }
  else
    return false;
}

Bool_t TG4StepManager::IsTrackAlive() const
{   
// Returns true if particle continues tracking.
// ---

#ifdef TGEANT4_DEBUG
  CheckTrack();
#endif

  G4TrackStatus status
     = fTrack->GetTrackStatus();
  if (status == fAlive)
    return true; 
  else
    return false; 
}

Bool_t TG4StepManager::IsNewTrack() const
{
// Returns true when track performs the first step.
// ---

  if (fStepStatus == kVertex)
    return true;
  else  
    return false;
}

Int_t TG4StepManager::NSecondaries() const
{
// Returns the number of secondary particles generated 
// in the current step.
// ---

#ifdef TGEANT4_DEBUG
  CheckSteppingManager();
#endif

  G4int nofSecondaries = 0;
  nofSecondaries += fSteppingManager->GetfN2ndariesAtRestDoIt();
  nofSecondaries += fSteppingManager->GetfN2ndariesAlongStepDoIt();
  nofSecondaries += fSteppingManager->GetfN2ndariesPostStepDoIt();

  return nofSecondaries;
}

void TG4StepManager::GetSecondary(Int_t index, Int_t& particleId, 
                          TLorentzVector& position, TLorentzVector& momentum)
{
// Fills the parameters (particle encoding, position, momentum)
// of the generated secondary particle which is specified by index.
// !! Check if indexing of secondaries is same !!
// ---

#ifdef TGEANT4_DEBUG
  CheckSteppingManager();
#endif

  G4int nofSecondaries = NSecondaries();
  G4TrackVector* secondaryTracks = fSteppingManager->GetSecondary();

  if (secondaryTracks){
    if (index < nofSecondaries) {

      // the index of the first secondary of this step
      G4int startIndex 
        = secondaryTracks->entries() - nofSecondaries;
             // (the secondaryTracks vector contains secondaries 
             // produced by the track at previous steps, too)
      G4Track* track 
        = (*secondaryTracks)[startIndex + index]; 
   
      // particle encoding
      particleId 
        = track->GetDynamicParticle()->GetDefinition()->GetPDGEncoding();
 
      // position & time
      G4ThreeVector positionVector = track->GetPosition();
      positionVector *= 1./(TG4G3Units::Length());
      G4double time = track->GetLocalTime();
      time /= TG4G3Units::Time();
      SetTLorentzVector(positionVector, time, position);

      // momentum & energy
      G4ThreeVector momentumVector = track->GetMomentum();	
      G4double energy = track->GetDynamicParticle()->GetTotalEnergy();
      energy /= TG4G3Units::Energy();
      SetTLorentzVector(momentumVector, energy, momentum);
    }
    else {
      TG4Globals::Exception(
        "TG4StepManager::GetSecondary(): wrong secondary track index.");
    }
  }
  else {
    TG4Globals::Exception(
      "TG4StepManager::GetSecondary(): secondary tracks vector is empty");
  }
}

AliMCProcess TG4StepManager::ProdProcess(Int_t isec) const
{
// The process that has produced the secondary particles specified 
// with isec index in the current step.
// ---

  G4int nofSecondaries = NSecondaries();
  if (fStepStatus == kVertex || !nofSecondaries) return kPNoProcess;

#ifdef TGEANT4_DEBUG
  CheckStep("ProdProcess");
#endif

  G4TrackVector* secondaryTracks = fSteppingManager->GetSecondary();
 
  // should never happen
  if (!secondaryTracks) {
    TG4Globals::Exception(
      "TG4StepManager::ProdProcess(): secondary tracks vector is empty.");

    return kPNoProcess;  
  }    

  if (isec < nofSecondaries) {

    // the index of the first secondary of this step
    G4int startIndex 
      = secondaryTracks->entries() - nofSecondaries;
           // the secondaryTracks vector contains secondaries 
           // produced by the track at previous steps, too

    // the secondary track with specified isec index
    G4Track* track = (*secondaryTracks)[startIndex + isec]; 
   
    const G4VProcess* kpProcess = track->GetCreatorProcess(); 
  
    AliMCProcess mcProcess 
     = TG4PhysicsManager::Instance()->GetMCProcess(kpProcess);
  
    // distinguish kPDeltaRay from kPEnergyLoss  
    if (mcProcess == kPEnergyLoss) mcProcess = kPDeltaRay;
  
    return mcProcess;
  }
  else {
    TG4Globals::Exception(
      "TG4StepManager::GetSecondary(): wrong secondary track index.");

    return kPNoProcess;  
  }
}


Int_t TG4StepManager::StepProcesses(TArrayI &proc) const
{
// Fills the array of processes that were active in the current step
// and returns the number of them.
// TBD: Distinguish between kPDeltaRay and kPEnergyLoss
// ---

 if (fStepStatus == kVertex) {
   G4cout << "kVertex" << G4endl;
   G4int nofProcesses = 1;
   proc.Set(nofProcesses);
   proc[0] = kPNull;
   return nofProcesses;
 }  
   
#ifdef TGEANT4_DEBUG
  CheckSteppingManager();
  CheckStep("StepProcesses");
#endif

  // along step processes
  G4ProcessManager* processManager
    = fStep->GetTrack()->GetDefinition()->GetProcessManager();
  G4ProcessVector* alongStepProcessVector 
    = processManager->GetAlongStepProcessVector();
  G4int nofProcesses = alongStepProcessVector->entries();
  
  // process defined step
  const G4VProcess* kpLastProcess 
    = fStep->GetPostStepPoint()->GetProcessDefinedStep();

  // fill the array of processes 
  proc.Set(nofProcesses);
  TG4PhysicsManager* physicsManager = TG4PhysicsManager::Instance();
  G4int i;  
  for (i=0; i<nofProcesses-1; i++) {
    G4VProcess* g4Process = (*alongStepProcessVector)[i];    
    // do not fill transportation along step process
    if (g4Process->GetProcessName() != "Transportation") {
      physicsManager->GetMCProcess(g4Process);   
      proc[i] = physicsManager->GetMCProcess(g4Process);
    }  
  }  
  proc[nofProcesses-1] = physicsManager->GetMCProcess(kpLastProcess);
    
  return nofProcesses;  
}
