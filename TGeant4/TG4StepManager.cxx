// $Id$
// Category: event
//
// See the class description in the header file.

#include "TG4StepManager.h"
#include "TG4GeometryManager.h"
#include "TG4PhysicsManager.h"
#include "TG4VSensitiveDetector.h"
#include "TG4Globals.h"
#include "TG3Units.h"

#include <G4Step.hh>
#include <G4SteppingManager.hh>
#include <G4UserLimits.hh>
#include <G4ParticleTable.hh>
#include <G4UImanager.hh>
#include <G4AffineTransform.hh>
//#include <G4TransportationManager.hh>
//#include <G4Navigator.hh>

#include <Randomize.hh>

#include <TLorentzVector.h>

TG4StepManager* TG4StepManager::fgInstance = 0;

TG4StepManager::TG4StepManager() 
  : fStep(0),
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

  if (fStep) {
    // check
    fStep->GetTrack()->SetTrackStatus(fStopAndKill);
    // fStep->GetTrack()->SetTrackStatus(fStopButAlive);
    // fStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
  }
}

void TG4StepManager::StopEvent()
{
// Aborts the current event processing.
// ---

  if (fStep) {
    fStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
      //StopTrack();   // cannot be used as it keeps secondaries
    G4UImanager::GetUIpointer()->ApplyCommand("/event/abort");
    G4UImanager::GetUIpointer()->ApplyCommand("/alStacking/clearStack");
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
  }
}

void TG4StepManager::Rndm(Float_t* array, const Int_t size) const
{   
// Random numbers array of the specified size.
// ---

  G4double* const kpDoubleArray = new G4double[size];
  RandFlat::shootArray(size,kpDoubleArray);
  for (G4int i=0; i<size; i++) { 
    array[i] = kpDoubleArray[i]; 
  } 
  delete [] kpDoubleArray;
}
 
void TG4StepManager::SetMaxStep(Float_t step)
{
// Maximum step allowed in the current logical volume.
// The maximum value is kept for following tracks - is it ok ??
// ---

  if (fStep) {
    // check this
    G4LogicalVolume* curLogVolume
      = fStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume();
    G4UserLimits* userLimits 
      = curLogVolume->GetUserLimits();
    if (userLimits == 0)
    { 
      userLimits = new G4UserLimits(step); 
      curLogVolume->SetUserLimits(userLimits);
    }
    else
    { userLimits->SetMaxAllowedStep(step); }  
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
  }
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

Int_t TG4StepManager::CurrentVolID(Int_t& copyNo) const
{
// Returns the current sensitive detector ID
// and the copy number of the current physical volume.
// ---

  if (fStep) {
    // check this
    G4VPhysicalVolume* physVolume 
      = fStep->GetPreStepPoint()->GetPhysicalVolume();
    copyNo = physVolume->GetCopyNo();

    // sensitive detector ID
    G4VSensitiveDetector* sd
        = physVolume->GetLogicalVolume()->GetSensitiveDetector();
    if (sd) {
      TG4VSensitiveDetector* tsd = dynamic_cast<TG4VSensitiveDetector*>(sd);
      if (tsd)
        return tsd->GetID();
      else {
        TG4Globals::Exception(
          "TG4StepManager::CurrentVol: Unknown sensitive detector type");
        return 0;
      }   	
    }  
    else {
      G4String text = "TG4StepManager::CurrentVol: \n";
      text = text + "    Volume " + physVolume->GetName();
      text = text + " has not a sensitive detector.";
      TG4Globals::Exception(text);
      return 0;
    }      	
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return 0;
  }
}

Int_t TG4StepManager::CurrentVolOffID(Int_t off, Int_t&  copyNo) const
{ 
// Returns the off-th mother's of the current volume
// the sensitive detector ID and the copy number.
// ---

  if (off == 0) return CurrentVolID(copyNo);

  if (fStep) {
    // check this !!
    // mother of physical volume may not be set ?!!
 
    G4VPhysicalVolume* physVolume 
      = fStep->GetPreStepPoint()->GetPhysicalVolume();

    G4VPhysicalVolume* mother = 0;
    Int_t level = off;
    while (level > 0) { 
      mother = physVolume->GetMother();
      level--;
    }
    
    if (mother) {
      copyNo = mother->GetCopyNo();

      // sensitive detector ID
      G4VSensitiveDetector* sd
        = mother->GetLogicalVolume()->GetSensitiveDetector();
      if (sd) {
        TG4VSensitiveDetector* tsd = dynamic_cast<TG4VSensitiveDetector*>(sd);
        if (tsd)
          return tsd->GetID();
        else {
          TG4Globals::Exception(
            "TG4StepManager::CurrentVolOff: Unknown sensitive detector type");
          return 0;
        }   	
      }  
      else {
        G4String text = "TG4StepManager::CurrentVolOff: \n";
        text = text + "    Volume " + mother->GetName();
        text = text + " has not a sensitive detector.";
        TG4Globals::Exception(text);
        return 0;
      } 
    }
    else {       	
      G4String text = "TG4StepManager::CurrentVolOff: Volume ";
      text = text + physVolume->GetName();
      text = text + " has not defined mother in required level.";
      TG4Globals::Exception(text);
      return 0;
    }  
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return 0;
  }
}

const char* TG4StepManager::CurrentVolName() const
{
// Returns the current physical volume name.
// ---

  if (fStep) {
    G4VPhysicalVolume* physVolume 
      = fStep->GetPreStepPoint()->GetPhysicalVolume();
    return physVolume->GetName();
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return 0;
  }
}

const char* TG4StepManager::CurrentVolOffName(Int_t off) const
{ 
// Returns the off-th mother's physical volume name.
// ---

  if (fStep) {
    G4VPhysicalVolume* physVolume 
      = fStep->GetPreStepPoint()->GetPhysicalVolume();

    G4VPhysicalVolume* mother = 0;
    Int_t level = off;
    while (level > 0) { 
      mother = physVolume->GetMother();
      level--;
    }
    
    if (mother) {
      return mother->GetName();
    }
    else {
      TG4Globals::Exception("TG4StepManager::CurrentVolOff: wrong usage.");
      return 0;
    }
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return 0;
  }
}

Int_t TG4StepManager::CurrentMaterial(Float_t &a, Float_t &z, Float_t &dens, 
                          Float_t &radl, Float_t &absl) const
{
// Returns the parameters of the current material during transport
// the return value is the number of elements in the mixture.
// ---

  if (fStep) {
    // check this
    G4LogicalVolume* curLogVolume 
      = fStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume();
    G4Material* material
      = curLogVolume->GetMaterial();

    // this if may be redundant - check
    if (material)
    {
      G4int nofElements = material->GetNumberOfElements();
      TG4GeometryManager* pGeometryManager = TG4GeometryManager::Instance();
      a = pGeometryManager->GetEffA(material);
      z = pGeometryManager->GetEffZ(material);
      
      // density 
      dens = material->GetDensity();
      dens /= TG3Units::MassDensity();      
      
      // radiation length
      radl = material->GetRadlen();
      radl /= TG3Units::Length();
      
      absl = 0.;  // this parameter is not defined in Geant4
      return nofElements;
    }
    else {
      TG4Globals::Exception(
       "TG4StepManager::CurrentMaterial(..): material is not defined.");
      return 0;
    }
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
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

  if (fStep) {
 
    G4ThreeVector theGlobalPoint(xm[0],xm[1],xm[2]); 

    //const G4NavigationHistory* history
    //   =  fStep->GetPreStepPoint()->GetTouchable()->GetHistory();
    G4AffineTransform affineTransform
      = fStep->GetPreStepPoint()->GetTouchable()->GetHistory()
        ->GetTopTransform();

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
     
 /* 
   // does not work ???
   G4ThreeVector direction(0,0,0);
   G4bool RelativeSearch = true;
   G4Navigator* theNavigator =
     G4TransportationManager::GetTransportationManager()->
     GetNavigatorForTracking();

   G4VPhysicalVolume* pPhysVol;
   pPhysVol
     = LocateGlobalPointAndSetup(theGlobalPoint, &direction, RelativeSearch);  
   //LocateGlobalPointWithinVolume(theGlobalPoint);

   G4AffineTransform at
     = theNavigator->GetGlobalToLocalTransform();
 */  
  } 
  else { 
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
  }  
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

  if (fStep) {
    // check this
    
   G4ThreeVector theLocalPoint(xd[0],xd[1],xd[2]); 

   G4AffineTransform affineTransform
      = fStep->GetPreStepPoint()->GetTouchable()->GetHistory()
        ->GetTopTransform().Inverse();
  
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
  else { 
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
  }  
} 
 
Float_t TG4StepManager::MaxStep() const
{   
// Returns maximum step allowed in the current logical volume
// by User Limits.
// ---

  if (fStep) {
    // check this
    G4LogicalVolume* curLogVolume
      = fStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume();
    G4UserLimits* userLimits 
      = curLogVolume->GetUserLimits();
    G4double maxStep;
    if (userLimits == 0)
    { 
      G4String text = "User Limits are not defined for log volume ";
      text = text + curLogVolume->GetName();
      TG4Globals::Warning(text);
      return FLT_MAX;
    }
    else
    { 
      const G4Track& trackRef = *(fStep->GetTrack());
      maxStep = userLimits->GetMaxAllowedStep(trackRef); 
      maxStep /= TG3Units::Length(); 
      return maxStep;
    }  
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return FLT_MAX;
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

  if (fStep) {
    // position
    G4ThreeVector positionVector
       = fStep->GetPostStepPoint()->GetPosition();
    positionVector *= 1./(TG3Units::Length());   
     
    // local time   
    G4double time
      = fStep->GetTrack()->GetLocalTime();
    time /= TG3Units::Time();
      
    SetTLorentzVector(positionVector, time, position);
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
  }
}

Int_t TG4StepManager::GetMedium() const
{   
// Returns the second index of the current material (corresponding to
// G3 tracking medium index).
// --- 

  if (fStep) {
    // current material
    G4Material* curMaterial
      = fStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial();

    // medium index  
    TG4GeometryManager* pGeometryManager = TG4GeometryManager::Instance();
    return pGeometryManager->GetMediumId(curMaterial);
  }  
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return 0;
  }
}

void TG4StepManager::TrackMomentum(TLorentzVector& momentum) const
{  
// Current particle "momentum" (px, py, pz, Etot).
// ---

  if (fStep) {
    G4ThreeVector momentumVector
      = fStep->GetTrack()->GetMomentum(); 

    G4double energy
      = fStep->GetTrack()->GetDynamicParticle()->GetTotalEnergy();
    energy /= TG3Units::Energy();  

    SetTLorentzVector(momentumVector, energy, momentum);
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
  }
}

void TG4StepManager::TrackVertexPosition(TLorentzVector& position) const
{ 
// The vertex particle position (in the world reference frame)
// and the local time since the current track is created.
// ---

  if (fStep) {
    // position
    G4ThreeVector positionVector
       = fStep->GetTrack()->GetVertexPosition();
    positionVector *= 1./(TG3Units::Length());   
     
    // local time 
    // to be checked  
    G4double time
      = fStep->GetTrack()->GetLocalTime();
    time /= TG3Units::Time();
      
    SetTLorentzVector(positionVector, time, position);
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
  }
}

void TG4StepManager::TrackVertexMomentum(TLorentzVector& momentum) const
{  
// The vertex particle "momentum" (px, py, pz, Ekin)
// to do: change Ekin -> Etot 
// ---

  if (fStep) {
    G4ThreeVector momentumVector
      = fStep->GetTrack()->GetVertexMomentumDirection(); 

    G4double energy
      = fStep->GetTrack()->GetVertexKineticEnergy();
    energy /= TG3Units::Energy();  

    SetTLorentzVector(momentumVector, energy, momentum);
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
  }
}

Float_t TG4StepManager::TrackStep() const
{   
// Returns the current step length.
// ---

  if (fStep) {
    // check
    G4double length = fStep->GetStepLength();
    length /= TG3Units::Length();
    return length;
  }
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return 0.;
  }
}

Float_t TG4StepManager::TrackLength() const
{
// Returns the length of the current track from its origin.
// ---

  if (fStep) {
    G4double length
      = fStep->GetTrack()->GetTrackLength();
    length /= TG3Units::Length();
    return length;
  }
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return 0;
  }
}

Float_t TG4StepManager::TrackTime() const
{
// Returns the local time since the current track is created.
// Comment:
// in Geant4: there is also defined proper time as
// the proper time of the dynamical particle of the current track.
// ---

  if (fStep) {
    G4double time
      = fStep->GetTrack()->GetLocalTime();
    time /= TG3Units::Time();
    return time;
  }
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return 0;
  }
}

Float_t TG4StepManager::Edep() const
{   
// Returns total energy deposit in this step.
// ---

  if (fStep) {
    // G4double deltaEnergy = fStep->GetDeltaEnergy();
    G4double deltaEnergy = fStep->GetTotalEnergyDeposit();
    deltaEnergy /= TG3Units::Energy();
    return deltaEnergy;
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return 0.;
  }
}

Int_t TG4StepManager::TrackPid() const
{   
// Returns the current particle PDG encoding.
// ---

  if (fStep) {
    G4ParticleDefinition* particle
      = fStep->GetTrack()->GetDynamicParticle()->GetDefinition();
      
    // ask TG4PhysicsManager to get PDG encoding 
    // (in order to get PDG from extended TDatabasePDG
    // in case the standard PDG code is not defined)
    TG4PhysicsManager* pPhysicsManager = TG4PhysicsManager::Instance();
    G4int pdgEncoding = pPhysicsManager->GetPDGEncodingFast(particle);

    return pdgEncoding;
  }
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return 0;
  }
}

Float_t TG4StepManager::TrackCharge() const
{   
// Returns the current particle charge.
// ---

  if (fStep) {
    G4double charge
      = fStep->GetTrack()->GetDynamicParticle()->GetDefinition()
        ->GetPDGCharge();
    charge /= TG3Units::Charge();	
    return charge;
  }
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return 0;
  }
}

Float_t TG4StepManager::TrackMass() const
{   
// Returns current particle rest mass.
// ---

  if (fStep) {
    G4double mass
      = fStep->GetTrack()->GetDynamicParticle()->GetDefinition()
        ->GetPDGMass();
    mass /= TG3Units::Mass();	
    return mass;
  }
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return 0;
  }
}

Float_t TG4StepManager::Etot() const
{   
// Returns total energy of the current particle.
// ---

  if (fStep) {
    G4double energy
      = fStep->GetTrack()->GetDynamicParticle()->GetTotalEnergy();
    energy /= TG3Units::Energy();  
    return energy;
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return 0.;
  }
}

Bool_t TG4StepManager::IsTrackInside() const
{   
// Returns true if particle does not cross geometrical boundary
// at both pre-step and post-step points.
// ---

  if ( !(IsTrackEntering()) && !(IsTrackExiting()) )
  { return true; }
  else
  { return false; }
}

Bool_t TG4StepManager::IsTrackEntering() const
{   
// Returns true if particle cross a geometrical boundary
// at pre-step point.
// ---

  if (fStep) {
    G4StepStatus status
       = fStep->GetPreStepPoint()->GetStepStatus();
    if (status == fGeomBoundary)
    { return true; }
    else
    { return false; }
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return false;
  }
}

Bool_t TG4StepManager::IsTrackExiting() const
{   
// Returns true if particle cross a geometrical boundary
// at post-step point.
// ---

  if (fStep) {
    // check
    G4StepStatus status
       = fStep->GetPostStepPoint()->GetStepStatus();
    if (status == fGeomBoundary)
    { return true; }
    else
    { return false; }
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return false;
  }
}

Bool_t TG4StepManager::IsTrackOut() const
{   
// Returns true if particle cross the Hall frame 
// at post-step point.
// ---

  if (fStep) {
    // check
    G4String volName
       = fStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
    if (volName == "pHallFrame")
    { return true; }
    else
    { return false; }
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return false;
  }
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

  if (fStep) {
    // check
    G4TrackStatus status
       = fStep->GetTrack()->GetTrackStatus();
    if ((status == fStopAndKill) ||  
        (status == fKillTrackAndSecondaries) ||
        (status == fSuspend) ||
        (status == fPostponeToNextEvent))
    { return true; }
    else
    { return false; }
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return false;
  }
}

Bool_t TG4StepManager::IsTrackDisappeared() const
{ 
// Returns true if particle has disappeared 
// (due to any physical process)
// or has been killed, suspended or postponed to next event.
// ---

  if (fStep) {
    // check
    G4TrackStatus status
       = fStep->GetTrack()->GetTrackStatus();
    if ((status == fStopButAlive) ||  
        (status == fKillTrackAndSecondaries) ||
        (status == fSuspend) ||
        (status == fPostponeToNextEvent))
    { return true; }
    else
    { return false; }
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return false;
  }
}

Bool_t TG4StepManager::IsTrackAlive() const
{   
// Returns true if particle continues tracking.
// ---

  if (fStep) {
    G4TrackStatus status
       = fStep->GetTrack()->GetTrackStatus();
    if (status == fAlive)
    { return true; }
    else
    { return false; }
  } 
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return false;
  }
}

Bool_t TG4StepManager::IsNewTrack() const
{
// Returns true when track performs the first step.
// ---

  if (fStep->GetTrack()->GetCurrentStepNumber() == 1)
    return true;
  else  
    return false;
}

Int_t TG4StepManager::NSecondaries() const
{
// Returns the number of secondary particles generated 
// in the current step.
// ---

  if (fSteppingManager) 
  {
    G4int nofSecondaries = 0;
    nofSecondaries += fSteppingManager->GetfN2ndariesAtRestDoIt();
    nofSecondaries += fSteppingManager->GetfN2ndariesAlongStepDoIt();
    nofSecondaries += fSteppingManager->GetfN2ndariesPostStepDoIt();

    return nofSecondaries;
  }
  else {   
    TG4Globals::Exception("TG4StepManager: SteppingManager is not defined.");
    return 0;
  }
}

void TG4StepManager::GetSecondary(Int_t index, Int_t& particleId, 
                          TLorentzVector& position, TLorentzVector& momentum)
{
// Fills the parameters (particle encoding, position, momentum)
// of the generated secondary particle which is specified by index.
// !! Check if indexing of secondaries is same !!
// ---

  if (fSteppingManager) 
  {
    G4int nofSecondaries = NSecondaries();
    G4TrackVector* secondaryTracks = fSteppingManager->GetSecondary();

    if (secondaryTracks)
    {
      if (index < nofSecondaries)
      {
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
	positionVector *= 1./(TG3Units::Length());
        G4double time = track->GetLocalTime();
	time /= TG3Units::Time();
        SetTLorentzVector(positionVector, time, position);

        // momentum & energy
        G4ThreeVector momentumVector = track->GetMomentum();	
        G4double energy = track->GetDynamicParticle()->GetTotalEnergy();
	energy /= TG3Units::Energy();
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
  else {   
    TG4Globals::Exception("TG4StepManager: SteppingManager is not defined.");
  }
}

const char* TG4StepManager::ProdProcess() const
{
// Returns the name of the process that defined current step
// (and may produce the secondary particles).
// ---

  if (fStep) 
  {
    const G4VProcess* curProcess 
      = fStep->GetPostStepPoint()->GetProcessDefinedStep(); 

    G4String g4Name = curProcess->GetProcessName(); 
    return g4Name;
  }
  else {   
    TG4Globals::Exception("TG4StepManager: Step is not defined.");
    return "";
  }
}
