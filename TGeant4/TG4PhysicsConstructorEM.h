// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorEM
// -----------------------------
// Constructor of electromagnetic physics.
// According to the corresponding part of:
// ExN04PhysicsList.hh, GEANT4 tag Name: geant4-03-02

#ifndef TG4_PHYSICS_CONSTRUCTOR_EM_H
#define TG4_PHYSICS_CONSTRUCTOR_EM_H

#include <G4VPhysicsConstructor.hh>
#include <G4PhotoElectricEffect.hh>
#include <G4ComptonScattering.hh>
#include <G4GammaConversion.hh>
#include <G4MultipleScattering.hh>
#include <G4eIonisation.hh>
#include <G4eBremsstrahlung.hh>
#include <G4eplusAnnihilation.hh>
#include <globals.hh>

class TG4PhysicsConstructorEM: public G4VPhysicsConstructor
{
  public:
    TG4PhysicsConstructorEM(const G4String& name = "EM");
    virtual ~TG4PhysicsConstructorEM();

  protected:
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    // data members
            // Gamma physics
    G4PhotoElectricEffect fPhotoEffect;
    G4ComptonScattering   fComptonEffect;
    G4GammaConversion     fPairProduction;
  
            // Electron physics
    G4MultipleScattering  fElectronMultipleScattering;
    G4eIonisation         fElectronIonisation;
    G4eBremsstrahlung     fElectronBremsStrahlung;
  
            //Positron physics
    G4MultipleScattering  fPositronMultipleScattering;
    G4eIonisation         fPositronIonisation; 
    G4eBremsstrahlung     fPositronBremsStrahlung;  
    G4eplusAnnihilation   fAnnihilation; 
    
  private:
    // methods
    void ConstructProcessForGamma();    
    void ConstructProcessForElectron();    
    void ConstructProcessForPositron();    
};

#endif //TG4_PHYSICS_CONSTRUCTOR_H

