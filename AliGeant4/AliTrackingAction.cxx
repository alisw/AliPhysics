// $Id$
// Category: event
//
// See the class description in the header file.

#include "AliTrackingAction.h"
#include "AliTrackingActionMessenger.h"
#include "AliSensitiveDetector.h"
#include "AliRun.h"
#include "AliGlobals.h"  
#include "TG4StepManager.h"

#include <G4TrackingManager.hh>
#include <G4Track.hh>
#include <G4Event.hh>
#include <G4SDManager.hh>
#include <G4VSensitiveDetector.hh>
#include <G4VHitsCollection.hh>

#include <TParticle.h>

// static data members
AliTrackingAction* AliTrackingAction::fgInstance = 0;

// static methods
AliTrackingAction* AliTrackingAction::Instance() {
// 
  return fgInstance; 
}

AliTrackingAction::AliTrackingAction()
  : fParticles(0),
    fPrimaryTrackID(0),
    fVerboseLevel(2),
    fSavePrimaries(true),
    fPrimariesCounter(0),
    fParticlesCounter(0),
    fLastParticleIndex(-1)
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

// public methods

void AliTrackingAction::PrepareNewEvent()
{
// Called by G4 kernel at the beginning of event.
// ---

  // aliroot
  if (!fParticles) fParticles = gAlice->Particles();

  // set the particles and primaries counter to the
  // number of tracks already stored in AliRun::fParticles
  // (fParticles can already contain primary particles
  //  saved by AliGenerator)
  G4int nofTracks = gAlice->GetNtrack();
  fPrimariesCounter = nofTracks;
  fParticlesCounter = nofTracks;
  fLastParticleIndex = nofTracks-1;

  // set g4 stepping manager pointer
  G4SteppingManager* pG4StepManager 
    = fpTrackingManager->GetSteppingManager();
  TG4StepManager* pStepManager = TG4StepManager::Instance();
  pStepManager->SetSteppingManager(pG4StepManager);
}

void AliTrackingAction::PreTrackingAction(const G4Track* aTrack)
{
// Called by G4 kernel before starting tracking.
// ---

  // aliroot
  // track index in the particles array
  G4int trackID = aTrack->GetTrackID();
  G4int parentID = aTrack->GetParentID();
  Int_t trackIndex;
  if (parentID==0) { 
    trackIndex = trackID; 
  } 
  else { 
    trackIndex = GetNofSavedTracks(); 
  }

  // in AliRoot (from V3.0) track numbering starts from 0
  gAlice->SetCurrentTrack(--trackIndex);

  if (parentID == 0) {  
    // save primary track
    SaveAndDestroyTrack();
    fPrimaryTrackID = aTrack->GetTrackID();
  }
  else { 
    // save secondary particles info 
    // improve this later with retrieving the generation process
    // (primary particles are stored 
    //  by AlStackingAction in ClassifyNewTrack() method)
    G4String origin = "secondary"; 
    SaveParticle(aTrack, origin);
  }
}

void AliTrackingAction::PostTrackingAction(const G4Track* aTrack)
{
// Called by G4 kernel after finishing tracking.
// ---

  G4String particleName 
    = aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName();
  if (particleName == "opticalphoton") {
    G4cout << "$$$ Track " <<  aTrack->GetTrackID()
           << " is optical photon." << G4endl;
  }	                
}

void AliTrackingAction::SaveAndDestroyTrack()
{
// Calls AliRun::PurifyKine and fills trees of hits
// after finishing tracking of each primary track.
// !! This method has to be also called from AlEventAction::EndOfEventAction() 
// for storing the last primary track of the current event.
// --- 

  if (fPrimaryTrackID>0)
  {
     if (fVerboseLevel == 3) { 
       G4cout << "$$$ Primary track " << fPrimaryTrackID << G4endl;
     } 
     else if ( fVerboseLevel == 2 &&  fPrimaryTrackID % 10 == 0 ) {
         G4cout << "$$$ Primary track " << fPrimaryTrackID  << G4endl;
     } 
     else if ( fVerboseLevel == 1 &&  fPrimaryTrackID % 100 == 0 ) {
         G4cout << "$$$ Primary track " << fPrimaryTrackID  << G4endl;
     } 
     
     // aliroot
     G4int lastSavedTrack 
       = gAlice->PurifyKine(fLastParticleIndex, fParticlesCounter);
     G4int nofPuredSavedTracks 
       = gAlice->GetNtrack();
     fLastParticleIndex = lastSavedTrack;
     fParticlesCounter = nofPuredSavedTracks;

     if(gAlice->TreeH()) gAlice->TreeH()->Fill();
     gAlice->ResetHits();
   }
   else { 
     fLastParticleIndex = fPrimariesCounter-1; 
   }
   fPrimaryTrackID = 0;
}  

void AliTrackingAction::SaveParticle(const G4Track* track, 
                                     G4String processName)
{
// Converts G4track to TParticle and saves it in AliRun::fParticles
// array.
// ----

  fParticlesCounter++;

  // track history
  G4int firstDaughter = -1; 
  G4int lastDaughter = -1;      
  G4int parentID = track->GetParentID();
  G4int motherIndex1;
  G4int motherIndex2 = -1;
  if (parentID == 0) { 
    motherIndex1 = -1; 
    fPrimariesCounter++;    
  }
  else {
    // set first/last child for already defined particles
    motherIndex1 = GetParticleIndex(parentID);
    // aliroot
    TParticle* parent 
      = dynamic_cast<TParticle*>(fParticles->UncheckedAt(motherIndex1));
    if (parent) {  
      if (parent->GetFirstDaughter()<0) 
        parent->SetFirstDaughter(fParticlesCounter);
      parent->SetLastDaughter(fParticlesCounter);
    }
    else {
      AliGlobals::Exception(
        "AliTrackingAction::SaveParticle: Unknown particle type");
    }   	
  };
     
  // ?? status of particle 
  // temporarily used for storing trackID from G4
  G4int ks = track->GetTrackID();
   
  // particle type
  // is this G3 ekvivalent ?? - check
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

  // aliroot
  TClonesArray& theCollectionRef = *fParticles;
  G4int nofParticles = theCollectionRef.GetEntriesFast();
  TParticle* particle 
   = new(theCollectionRef[nofParticles]) 
       TParticle(pdg, ks, motherIndex1, motherIndex1, 
         firstDaughter, lastDaughter, px, py, pz, e, vx, vy, vz, t);
  particle->SetPolarisation(polX, polY, polZ);
  particle->SetBit(kKeepBit, false); 
}

G4int AliTrackingAction::GetParticleIndex(G4int trackID)
{
// Converts the trackID into the index of the particle
// in AliRun::fParticles array.
// ---

  if (trackID <= fPrimariesCounter) { 
    return trackID-1; 
  }
  else
    for (G4int i=fLastParticleIndex+1; i<fParticlesCounter; i++) {
      // aliroot
      TParticle* particle
        = (TParticle*)fParticles->UncheckedAt(i);
      if (trackID == particle->GetStatusCode()) return i;
    }

  AliGlobals::Exception("AliTrackingAction::GetParticleIndex() has failed.");
  return 0;   
}
