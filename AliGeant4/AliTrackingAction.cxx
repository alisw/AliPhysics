// $Id$
// Category: event
//
// See the class description in the header file.

#include "AliTrackingAction.h"
#include "AliTrackingActionMessenger.h"
#include "AliTrackInformation.h"
#include "AliRun.h"
#include "AliGlobals.h"  
#include "TG4StepManager.h"
#include "TG4PhysicsManager.h"

#include <G4Track.hh>
#include <G4TrackVector.hh>
#include <G4VUserTrackInformation.hh>
#include <G4TrackingManager.hh>
#include <G4SteppingManager.hh>

/*
#include <TTree.h>
#include <TParticle.h>
#include <TClonesArray.h>
*/

// static data members
AliTrackingAction* AliTrackingAction::fgInstance = 0;

AliTrackingAction::AliTrackingAction()
  : fPrimaryTrackID(0),
    fVerboseLevel(2),
    fSavePrimaries(true),
    fTrackCounter(0)
{
//
  if (fgInstance) { 
    AliGlobals::Exception("AliTrackingAction constructed twice."); 
  }

  fMessenger = new AliTrackingActionMessenger(this);
  fgInstance = this;
}

AliTrackingAction::AliTrackingAction(const AliTrackingAction& right) {
//
  AliGlobals::Exception("AliTrackingAction is protected from copying.");
}

AliTrackingAction::~AliTrackingAction() {
//
  delete fMessenger;
}

// operators

AliTrackingAction& 
AliTrackingAction::operator=(const AliTrackingAction &right)
{
  // check assignement to self
  if (this == &right) return *this;
  
  AliGlobals::Exception("AliTrackingAction is protected from assigning.");

  return *this;
}

// private methods

AliTrackInformation* AliTrackingAction::GetTrackInformation(
                                           const G4Track* track,
                                           const G4String& method) const
{
// Returns user track information.
// ---
 
  G4VUserTrackInformation* trackInfo = track->GetUserInformation();
  if (!trackInfo) return 0;  

  AliTrackInformation* aliTrackInfo
    = dynamic_cast<AliTrackInformation*>(trackInfo);
  if (!aliTrackInfo) { 
     G4String text = "AliTrackingAction::" + method + ":\n";
     text = text + "   Unknown track information type";
     AliGlobals::Exception(text);
  }
  
  return aliTrackInfo;
}    
  
// public methods

void AliTrackingAction::PrepareNewEvent()
{
// Called by G4 kernel at the beginning of event.
// ---

  fTrackCounter = 0;

  // set g4 stepping manager pointer
  TG4StepManager* stepManager = TG4StepManager::Instance();
  stepManager->SetSteppingManager(fpTrackingManager->GetSteppingManager());
}

void AliTrackingAction::PreTrackingAction(const G4Track* aTrack)
{
// Called by G4 kernel before starting tracking.
// ---

  // track index in the particles array
  G4int trackID = aTrack->GetTrackID();
  G4int parentID = aTrack->GetParentID();
  Int_t trackIndex;
  if (parentID==0) { 
    // in AliRoot (from V3.0) track numbering starts from 0
    trackIndex = trackID-1; 
  } 
  else { 
    trackIndex = gAlice->GetNtrack();
  }
  
  // set track index to track information
  AliTrackInformation* trackInfo
    = GetTrackInformation(aTrack, "PreTrackingAction");
  if (!trackInfo) {
    // create track information and set it to G4Track
    // if it does not yet exist
    trackInfo = new AliTrackInformation(trackIndex);
    fpTrackingManager->SetUserTrackInformation(trackInfo);
        // the track information is deleted together with its
        // G4Track object  
  }
  else
    trackInfo->SetTrackParticleID(trackIndex);	

  // set current track number
  gAlice->SetCurrentTrack(trackIndex);

  if (parentID == 0) {  
    // finish previous primary track
    FinishPrimaryTrack();
    fPrimaryTrackID = aTrack->GetTrackID();
  }
  else { 
    // save secondary particles info 
    SaveTrack(aTrack);
  }
  
  // aliroot pre track actions
  gAlice->PreTrack();
}

void AliTrackingAction::PostTrackingAction(const G4Track* aTrack)
{
// Called by G4 kernel after finishing tracking.
// ---

  fTrackCounter++;
  
  // set parent track particle index to all secondary tracks 
  G4TrackVector* secondaryTracks 
    = fpTrackingManager->GetSteppingManager()->GetSecondary();
  if (secondaryTracks){
    G4int i;
    for (i=0; i<secondaryTracks->entries(); i++) {
      G4Track* track = (*secondaryTracks)[i]; 

      if (track->GetUserInformation()) {
        // this should never happen
	G4String text = "AliTrackingAction::PostTrackingAction:\n";
	text = text + "    Inconsistent track information."; 
        AliGlobals::Exception(text);
      }	
      
      // get parent track index
      AliTrackInformation* aliParentInfo
        = GetTrackInformation(aTrack, "PostTrackingAction");
      G4int parentParticleID 
        = aliParentInfo->GetTrackParticleID();

      // create track information and set it to the G4Track
      AliTrackInformation* trackInfo 
        = new AliTrackInformation(-1, parentParticleID);
      track->SetUserInformation(trackInfo);
        // the track information is deleted together with its
        // G4Track object  
    } 	
  }
      
  // aliroot post track actions
  gAlice->PostTrack();
}

void AliTrackingAction::FinishPrimaryTrack()
{
// Calls AliRun::PurifyKine and fills trees of hits
// after finishing tracking of each primary track.
// !! This method has to be also called from AlEventAction::EndOfEventAction() 
// for storing the last primary track of the current event.
// --- 

  if (fPrimaryTrackID>0) {

    // verbose
    if (fVerboseLevel == 3) { 
      G4cout << "$$$ Primary track " << fPrimaryTrackID << G4endl;
    } 
    else if ( fVerboseLevel == 2 &&  fPrimaryTrackID % 10 == 0 ) {
      G4cout << "$$$ Primary track " << fPrimaryTrackID  << G4endl;
    } 
    else if ( fVerboseLevel == 1 &&  fPrimaryTrackID % 100 == 0 ) {
      G4cout << "$$$ Primary track " << fPrimaryTrackID  << G4endl;
    } 

    // aliroot finish primary track         
    gAlice->FinishPrimary();
  }
  fPrimaryTrackID = 0;
}  

void AliTrackingAction::SaveTrack(const G4Track* track)
{
// Get all needed parameters from G4track and pass them
// to AliRun::SetTrack() that creates corresponding TParticle
// in the AliRun::fParticles array.
// ----

  // parent particle index 
  G4int parentID = track->GetParentID();
  G4int motherIndex;
  if (parentID == 0) { 
    motherIndex = -1; 
  }
  else {
    motherIndex 
      = GetTrackInformation(track,"SaveTrack")->GetParentParticleID();
  }
     
  //G4cout << "SaveTrack: TrackID = " << track->GetTrackID()
  //       << "  Parent ID = " << track->GetParentID() 
  //       << "  Index = " << gAlice->CurrentTrack() 
  //       << "  Parent Index = " << motherIndex
  //       << G4endl;

  // PDG code
  G4int pdg = track->GetDefinition()->GetPDGEncoding();

  // track kinematics  
  G4ThreeVector momentum = track->GetMomentum(); 
  
  G4double px = momentum.x()/GeV;
  G4double py = momentum.y()/GeV;
  G4double pz = momentum.z()/GeV;
  G4double e = track->GetTotalEnergy()/GeV;

  G4ThreeVector position = track->GetPosition(); 
  G4double vx = position.x()/cm;
  G4double vy = position.y()/cm;
  G4double vz = position.z()/cm;
  // time of production - check if ekvivalent with G4
  G4double t = track->GetGlobalTime();

  G4ThreeVector polarization = track->GetPolarization(); 
  G4double polX = polarization.x();
  G4double polY = polarization.y();
  G4double polZ = polarization.z();

  // production process
  AliMCProcess mcProcess;  
  const G4VProcess* kpProcess = track->GetCreatorProcess();
  if (!kpProcess) {
    mcProcess = kPPrimary;
  }
  else {  
    TG4PhysicsManager* pPhysicsManager = TG4PhysicsManager::Instance();
    mcProcess = pPhysicsManager->GetMCProcess(kpProcess);  
    // distinguish kPDeltaRay from kPEnergyLoss  
    if (mcProcess == kPEnergyLoss) mcProcess = kPDeltaRay;
  }  
  
  G4int ntr;
  // create particle 
  gAlice->SetTrack(1, motherIndex, pdg, px, py, pz, e, vx, vy, vz, t,
                   polX, polY, polZ, mcProcess, ntr);  
                   
}

