// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorMuon
// -------------------------------
// Constructor of electromagnetic physics.
// According to ExN04MuonPhysics.hh, GEANT4 tag Name: geant4-03-02

#ifndef TG4_PHYSICS_CONSTRUCTOR_MUON_H
#define TG4_PHYSICS_CONSTRUCTOR_MUON_H

#include <G4VPhysicsConstructor.hh>
#include "G4MultipleScattering.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuIonisation.hh"
#include "G4MuonMinusCaptureAtRest.hh"
#include "G4hIonisation.hh"
#include <globals.hh>

class TG4PhysicsConstructorMuon: public G4VPhysicsConstructor
{
  public:
    TG4PhysicsConstructorMuon(const G4String& name = "EM");
    virtual ~TG4PhysicsConstructorMuon();

  protected:
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    // data members
          // mu+ physics
   G4MultipleScattering   fMuPlusMultipleScattering;
   G4MuBremsstrahlung     fMuPlusBremsstrahlung ;
   G4MuPairProduction     fMuPlusPairProduction;
   G4MuIonisation         fMuPlusIonisation;

          // mu- physics
   G4MultipleScattering   fMuMinusMultipleScattering;
   G4MuBremsstrahlung     fMuMinusBremsstrahlung ;
   G4MuPairProduction     fMuMinusPairProduction;
   G4MuIonisation         fMuMinusIonisation;
   G4MuonMinusCaptureAtRest fMuMinusCaptureAtRest;

          // tau+ physics
   G4MultipleScattering   fTauPlusMultipleScattering;
   G4hIonisation          fTauPlusIonisation;

          // tau+ physics
   G4MultipleScattering   fTauMinusMultipleScattering;
   G4hIonisation          fTauMinusIonisation;
    
  private:
    // methods
    void ConstructProcessForMuonPlus();    
    void ConstructProcessForMuonMinus();    
    void ConstructProcessForTauPlus();    
    void ConstructProcessForTauMinus();    
};

#endif //TG4_PHYSICS_CONSTRUCTOR_MUON_H

