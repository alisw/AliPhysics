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

#include "TG4VPhysicsConstructor.h"

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

class TG4PhysicsConstructorIon: public TG4VPhysicsConstructor
{
  public:
    TG4PhysicsConstructorIon(const G4String& name = "Ion");
    TG4PhysicsConstructorIon(G4int verboseLevel, 
                             G4bool setEM, G4bool setHadron,
                             const G4String& name = "Ion");
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
   G4HadronElasticProcess fElasticProcess;         //elastic process
   G4LElastic*            fElasticModel;           //elastic model

          // Generic Ion physics
   G4MultipleScattering   fIonMultipleScattering;  //ion multiple scattering
   G4hIonisation          fIonIonisation;          //ion ionisation

          // Deuteron physics
   G4MultipleScattering        fDeuteronMultipleScattering; //D msc
   G4hIonisation               fDeuteronIonisation;//D ionisation
   G4DeuteronInelasticProcess  fDeuteronProcess;   //D inelastic process
   G4LEDeuteronInelastic*      fDeuteronModel;     //D LE inelastic model

          // Triton physics
   G4MultipleScattering        fTritonMultipleScattering; //T msc
   G4hIonisation               fTritonIonisation;  //T ionisation
   G4TritonInelasticProcess    fTritonProcess;     //T inelastic process
   G4LETritonInelastic*        fTritonModel;       //T LE inelastic model
  
         // Alpha physics
   G4MultipleScattering        fAlphaMultipleScattering; //alpha msc
   G4hIonisation               fAlphaIonisation;   //alpha ionisation
   G4AlphaInelasticProcess     fAlphaProcess;      //alpha inelastic process
   G4LEAlphaInelastic*         fAlphaModel;        //alpha LE inelastic model

        // He3 physics
   G4MultipleScattering        fHe3MultipleScattering; //He3 msc
   G4hIonisation               fHe3Ionisation;     //He3 ionisation
    
  private:
    // methods

    // Hadron processes
    void ConstructHadProcessForGenericIon();    
    void ConstructHadProcessForDeuteron();    
    void ConstructHadProcessForTriton();    
    void ConstructHadProcessForAlpha();    
    void ConstructHadProcessForHe3();    

    // EM processes
    void ConstructEMProcessForGenericIon();    
    void ConstructEMProcessForDeuteron();    
    void ConstructEMProcessForTriton();    
    void ConstructEMProcessForAlpha();    
    void ConstructEMProcessForHe3();    

    // data members
    G4bool  fSetEM;    //if true - EM processes are constructed
    G4bool  fSetHadron;//if true - hadron processes are constructed
};

#endif //TG4_PHYSICS_CONSTRUCTOR_ION_H

