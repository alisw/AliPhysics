// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorIon
// -----------------------------
// Constructor of physics for ions.
// According to ExN04IonPhysics.hh, GEANT4 tag Name: geant4-03-02

#ifndef TG4_PHYSICS_CONSTRUCTOR_ION_H
#define TG4_PHYSICS_CONSTRUCTOR_ION_H

#include <G4VPhysicsConstructor.hh>
#include <G4HadronElasticProcess.hh>
#include <G4LElastic.hh>
#include <G4DeuteronInelasticProcess.hh>
#include <G4LEDeuteronInelastic.hh>
#include <G4TritonInelasticProcess.hh>
#include <G4LETritonInelastic.hh>
#include <G4AlphaInelasticProcess.hh>
#include <G4LEAlphaInelastic.hh>
#include <G4hIonisation.hh>
#include <G4MultipleScattering.hh>
#include <globals.hh>

class TG4PhysicsConstructorIon: public G4VPhysicsConstructor
{
  public:
    TG4PhysicsConstructorIon(const G4String& name = "Ion");
    // --> protected
    // TG4PhysicsConstructorIon(const TG4PhysicsConstructorIon& right);
    virtual ~TG4PhysicsConstructorIon();

  protected:
    TG4PhysicsConstructorIon(const TG4PhysicsConstructorIon& right);
    
    // operators
    TG4PhysicsConstructorIon& operator=(const TG4PhysicsConstructorIon& right);
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    // data members
          // Elastic Process
   G4HadronElasticProcess fElasticProcess;
   G4LElastic*            fElasticModel;

          // Generic Ion physics
   G4MultipleScattering   fIonMultipleScattering;
   G4hIonisation          fIonIonisation;

          // Deuteron physics
   G4MultipleScattering        fDeuteronMultipleScattering;
   G4hIonisation               fDeuteronIonisation;
   G4DeuteronInelasticProcess  fDeuteronProcess;
   G4LEDeuteronInelastic*      fDeuteronModel;

          // Triton physics
   G4MultipleScattering        fTritonMultipleScattering;
   G4hIonisation               fTritonIonisation;
   G4TritonInelasticProcess    fTritonProcess;
   G4LETritonInelastic*        fTritonModel;
  
         // Alpha physics
   G4MultipleScattering        fAlphaMultipleScattering;
   G4hIonisation               fAlphaIonisation;
   G4AlphaInelasticProcess     fAlphaProcess;
   G4LEAlphaInelastic*         fAlphaModel;

        // He3 physics
   G4MultipleScattering        fHe3MultipleScattering;
   G4hIonisation               fHe3Ionisation;
    
  private:
    // methods
    void ConstructProcessForGenericIon();    
    void ConstructProcessForDeuteron();    
    void ConstructProcessForTriton();    
    void ConstructProcessForAlpha();    
    void ConstructProcessForHe3();    
};

#endif //TG4_PHYSICS_CONSTRUCTOR_ION_H

