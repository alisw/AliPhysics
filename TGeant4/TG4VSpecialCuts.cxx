// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4VSpecialCuts
// ---------------------
// See the class description in the header file.

#include "TG4VSpecialCuts.h"
#include "TG4G3CutVector.h"
#include "TG4GeometryServices.h"
#include "TG4Limits.h"

#include <G4UserLimits.hh>
#include <G4EnergyLossTables.hh>

//_____________________________________________________________________________
TG4VSpecialCuts::TG4VSpecialCuts(const G4String& processName)
  : G4VProcess(processName) {
//
}

//_____________________________________________________________________________
TG4VSpecialCuts::TG4VSpecialCuts() {
//
}

//_____________________________________________________________________________
TG4VSpecialCuts::TG4VSpecialCuts(const TG4VSpecialCuts& right) {
// 
  TG4Globals::Exception(
    "TG4VSpecialCuts is protected from copying.");
}

//_____________________________________________________________________________
TG4VSpecialCuts::~TG4VSpecialCuts() {
//
}

// operators

//_____________________________________________________________________________
TG4VSpecialCuts& TG4VSpecialCuts::operator=(const TG4VSpecialCuts& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "TG4VSpecialCuts is protected from assigning.");
    
  return *this;  
}    
          
// public methods

//_____________________________________________________________________________
G4double TG4VSpecialCuts::PostStepGetPhysicalInteractionLength(
                           const G4Track& track, G4double previousStepSize,
			   G4ForceCondition* condition)
{
// Returns the Step-size (actual length) which is allowed 
// by this process.
// ---

  // set condition
  *condition = NotForced;

  G4double proposedStep = DBL_MAX;
  G4double minStep = (1.0e-9)*m;
  
  // get limits
#ifdef TGEANT4_DEBUG
  TG4Limits* limits 
     = TG4GeometryServices::Instance()
         ->GetLimits(track.GetVolume()->GetLogicalVolume()->GetUserLimits());

  if (!limits) {
    G4String text = "TG4VSpecialCuts::PostStepGetPhysicalInteractionLength:\n";
    text = text + "    " + track.GetVolume()->GetLogicalVolume()->GetName();
    text = text + " has not limits.";
    TG4Globals::Exception(text);
  }  	  
#else  
  TG4Limits* limits 
    = (TG4Limits*) track.GetVolume()->GetLogicalVolume()->GetUserLimits();
#endif    

  // max track length
  proposedStep 
    = (limits->GetUserMaxTrackLength(track) - track.GetTrackLength());
  if (proposedStep < 0.) return minStep;

  // max time limit
  G4double beta 
    = (track.GetDynamicParticle()->GetTotalMomentum()) /
      (track.GetTotalEnergy());
  G4double time
    = (limits->GetUserMaxTime(track) - track.GetGlobalTime());
  G4double temp = beta*c_light*time;
  if (temp < 0.) return minStep;
  if (proposedStep > temp) proposedStep = temp;

  G4ParticleDefinition* particle = track.GetDefinition();
  if (particle->GetPDGCharge() != 0.) {

    // min remaining range
    G4double kinEnergy = track.GetKineticEnergy();
    G4Material* material = track.GetMaterial();
    G4double rangeNow 
      = G4EnergyLossTables::GetRange(particle, kinEnergy, material);

    temp = (rangeNow - limits->GetUserMinRange(track));
    if (temp < 0.) return minStep;
    if (proposedStep > temp) proposedStep = temp;

    // min kinetic energy (from limits)
    // the kin energy cut can be applied only in case
    // G4EnergyLossTables are defined for the particle    
    if (G4EnergyLossTables::GetDEDXTable(particle)) {
      G4double minEkine = GetMinEkine(*limits, track);
      G4double minR 
        = G4EnergyLossTables::GetRange(particle, minEkine, material);
      temp = rangeNow - minR;
      if (temp < 0.) return minStep;
      if (proposedStep > temp) proposedStep = temp;  
    }
  }
  else {  
    // min kinetic energy (from limits)
    // for neutral particles
    G4double minEkine = GetMinEkine(*limits, track);
    if (track.GetKineticEnergy() <= minEkine) return minStep;
  }
    
  return proposedStep;
}

//_____________________________________________________________________________
G4VParticleChange* TG4VSpecialCuts::PostStepDoIt(const G4Track& track, 
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
