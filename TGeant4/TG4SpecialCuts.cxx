// $Id$
// Category: physics
//
// See the class description in the header file.

#include "TG4SpecialCuts.h"
#include "TG4G3CutVector.h"
#include "TG4Limits.h"

#include <G4UserLimits.hh>

#include <G4EnergyLossTables.hh>

TG4SpecialCuts::TG4SpecialCuts(TG4G3ParticleWSP particle, 
                               TG4G3CutVector* cutVector,
                               const G4String& processName)
  : G4UserSpecialCuts(processName),
    fCutVector(cutVector)
{
//
  switch (particle) {
    case kGamma:
      fPtrMinEkineInCutVector = &TG4G3CutVector::GetMinEkineForGamma;
      fPtrMinEkineInLimits = &TG4Limits::GetMinEkineForGamma;
      break;
    case kElectron: case kEplus:  
      fPtrMinEkineInCutVector = &TG4G3CutVector::GetMinEkineForElectron;
      fPtrMinEkineInLimits = &TG4Limits::GetMinEkineForElectron;
      break;
    case kChargedHadron:  
      fPtrMinEkineInCutVector = &TG4G3CutVector::GetMinEkineForHadron;
      fPtrMinEkineInLimits = &TG4Limits::GetMinEkineForHadron;
      break;
    case kNeutralHadron:  
      fPtrMinEkineInCutVector = &TG4G3CutVector::GetMinEkineForNeutralHadron;
      fPtrMinEkineInLimits = &TG4Limits::GetMinEkineForNeutralHadron;
      break;
    case kMuon:  
      fPtrMinEkineInCutVector = &TG4G3CutVector::GetMinEkineForMuon;
      fPtrMinEkineInLimits = &TG4Limits::GetMinEkineForMuon;
      break;
    case kAny:
      fPtrMinEkineInCutVector = &TG4G3CutVector::GetMinEkineForOther;
      fPtrMinEkineInLimits = &TG4Limits::GetMinEkineForOther;
      break;
    case kNofParticlesWSP:
      TG4Globals::Exception("TG4SpecialCuts: Wrong particle specification.");
      break;  
  }
}

TG4SpecialCuts::TG4SpecialCuts() {
//
}

TG4SpecialCuts::TG4SpecialCuts(const TG4SpecialCuts& right) {
// 
  TG4Globals::Exception(
    "TG4SpecialCuts is protected from copying.");
}

TG4SpecialCuts::~TG4SpecialCuts() {
//
}

// operators

TG4SpecialCuts& TG4SpecialCuts::operator=(const TG4SpecialCuts& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "TG4SpecialCuts is protected from assigning.");
    
  return *this;  
}    
          
// public methods

G4double TG4SpecialCuts::PostStepGetPhysicalInteractionLength(
                           const G4Track& track, G4double previousStepSize,
			   G4ForceCondition* condition)
{
// Returns the Step-size (actual length) which is allowed 
// by this process.
// ---

  // set condition
  *condition = NotForced;

  G4double proposedStep = DBL_MAX;
  G4UserLimits* limits 
    = track.GetVolume()->GetLogicalVolume()->GetUserLimits();
  if (limits) { 
    // max track length
    proposedStep 
      = (limits->GetUserMaxTrackLength(track) - track.GetTrackLength());
    if (proposedStep < 0.) return 0.;

    // max time limit
    G4double beta 
      = (track.GetDynamicParticle()->GetTotalMomentum()) /
        (track.GetTotalEnergy());
    G4double time
      = (limits->GetUserMaxTime(track) - track.GetGlobalTime());
    G4double temp = beta*c_light*time;
    if (temp < 0.) return 0.;
    if (proposedStep > temp) proposedStep = temp;

    // min remaining range
    G4ParticleDefinition* particle = track.GetDefinition();
    G4double kinEnergy = track.GetKineticEnergy();
    G4Material* material = track.GetMaterial();
    G4double rangeNow 
      = G4EnergyLossTables::GetRange(particle, kinEnergy, material);
    temp = (rangeNow - limits->GetUserMinRange(track));
    if (temp < 0.) return 0.;
    if (proposedStep > temp) proposedStep = temp;

    // min kinetic energy (from limits)
    // the kin energy cut can be applied only in case
    // G4EnergyLossTables are defined for the particle
    if (G4EnergyLossTables::GetDEDXTable(particle)) {
      TG4Limits* tg4Limits = dynamic_cast<TG4Limits*>(limits);
      if (!tg4Limits) {
        G4String text = "TG4SpecialCuts::PostStepGetPhysicalInteractionLength:\n";
        text = text + "    Unknown limits type.";
        TG4Globals::Exception(text);
      }  
      G4double minEkine 
        = (tg4Limits->*fPtrMinEkineInLimits)(track);
      G4double minR 
        = G4EnergyLossTables::GetRange(particle, minEkine, material);
      temp = rangeNow - minR;
      if (temp < 0.) return 0.;
      if (proposedStep > temp) proposedStep = temp;  
    }  

  }
  else if (fCutVector) {
    // min kinetic energy (from cut vector)
    // the kin energy cut can be applied only in case
    // G4EnergyLossTables are defined for the particle
    G4ParticleDefinition* particle = track.GetDefinition();
    if (G4EnergyLossTables::GetDEDXTable(particle)) {
      G4double kinEnergy = track.GetKineticEnergy();
      G4Material* material = track.GetMaterial();
      G4double rangeNow 
        = G4EnergyLossTables::GetRange(particle, kinEnergy, material);
      G4double minEkine 
        = (fCutVector->*fPtrMinEkineInCutVector)(track);
      G4double minR 
        = G4EnergyLossTables::GetRange(particle, minEkine, material);
      G4double temp = rangeNow - minR;
      if (temp < 0.) return 0.;
      if (proposedStep > temp) proposedStep = temp;  
    }  
  }   
  return proposedStep;
}

G4VParticleChange* TG4SpecialCuts::PostStepDoIt(const G4Track& track, 
                                                const G4Step& step)
{
// Kills the current particle, if requested by G4UserLimits.
// ---
 
  aParticleChange.Initialize(track);
  aParticleChange.SetEnergyChange(0.) ;
  aParticleChange.SetLocalEnergyDeposit(track.GetKineticEnergy()) ;
  aParticleChange.SetStatusChange(fStopAndKill);
  return &aParticleChange;
}
