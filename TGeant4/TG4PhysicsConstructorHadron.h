// $Id$
// Category: physics
//
// Author: I. Hrivnacova
//
// Class TG4PhysicsConstructorHadron
// ---------------------------------
// Constructor of hadron physics.
// According to ExN04HadronPhysics.hh, GEANT4 tag Name: geant4-03-02

#ifndef TG4_PHYSICS_CONSTRUCTOR_HADRON_H
#define TG4_PHYSICS_CONSTRUCTOR_HADRON_H

#include <G4VPhysicsConstructor.hh>
#include <g4std/vector>
#include <globals.hh>

#include <G4MultipleScattering.hh>
#include <G4hIonisation.hh>

#include <G4HadronElasticProcess.hh>
#include <G4HadronFissionProcess.hh>
#include <G4HadronCaptureProcess.hh>

#include <G4PionPlusInelasticProcess.hh>
#include <G4PionMinusInelasticProcess.hh>
#include <G4KaonPlusInelasticProcess.hh>
#include <G4KaonZeroSInelasticProcess.hh>
#include <G4KaonZeroLInelasticProcess.hh>
#include <G4KaonMinusInelasticProcess.hh>
#include <G4ProtonInelasticProcess.hh>
#include <G4AntiProtonInelasticProcess.hh>
#include <G4NeutronInelasticProcess.hh>
#include <G4AntiNeutronInelasticProcess.hh>
#include <G4LambdaInelasticProcess.hh>
#include <G4AntiLambdaInelasticProcess.hh>
#include <G4SigmaPlusInelasticProcess.hh>
#include <G4SigmaMinusInelasticProcess.hh>
#include <G4AntiSigmaPlusInelasticProcess.hh>
#include <G4AntiSigmaMinusInelasticProcess.hh>
#include <G4XiZeroInelasticProcess.hh>
#include <G4XiMinusInelasticProcess.hh>
#include <G4AntiXiZeroInelasticProcess.hh>
#include <G4AntiXiMinusInelasticProcess.hh>
#include <G4DeuteronInelasticProcess.hh>
#include <G4TritonInelasticProcess.hh>
#include <G4AlphaInelasticProcess.hh>
#include <G4OmegaMinusInelasticProcess.hh>
#include <G4AntiOmegaMinusInelasticProcess.hh>

// Low-energy Models
#include <G4LElastic.hh>   
#include <G4LFission.hh>
#include <G4LCapture.hh>

#include <G4LEPionPlusInelastic.hh>
#include <G4LEPionMinusInelastic.hh>
#include <G4LEKaonPlusInelastic.hh>
#include <G4LEKaonZeroSInelastic.hh>
#include <G4LEKaonZeroLInelastic.hh>
#include <G4LEKaonMinusInelastic.hh>
#include <G4LEProtonInelastic.hh>
#include <G4LEAntiProtonInelastic.hh>
#include <G4LENeutronInelastic.hh>
#include <G4LEAntiNeutronInelastic.hh>
#include <G4LELambdaInelastic.hh>
#include <G4LEAntiLambdaInelastic.hh>
#include <G4LESigmaPlusInelastic.hh>
#include <G4LESigmaMinusInelastic.hh>
#include <G4LEAntiSigmaPlusInelastic.hh>
#include <G4LEAntiSigmaMinusInelastic.hh>
#include <G4LEXiZeroInelastic.hh>
#include <G4LEXiMinusInelastic.hh>
#include <G4LEAntiXiZeroInelastic.hh>
#include <G4LEAntiXiMinusInelastic.hh>
#include <G4LEDeuteronInelastic.hh>
#include <G4LETritonInelastic.hh>
#include <G4LEAlphaInelastic.hh>
#include <G4LEOmegaMinusInelastic.hh>
#include <G4LEAntiOmegaMinusInelastic.hh>

// High-energy Models

#include <G4HEPionPlusInelastic.hh>
#include <G4HEPionMinusInelastic.hh>
#include <G4HEKaonPlusInelastic.hh>
#include <G4HEKaonZeroInelastic.hh>
#include <G4HEKaonZeroInelastic.hh>
#include <G4HEKaonMinusInelastic.hh>
#include <G4HEProtonInelastic.hh>
#include <G4HEAntiProtonInelastic.hh>
#include <G4HENeutronInelastic.hh>
#include <G4HEAntiNeutronInelastic.hh>
#include <G4HELambdaInelastic.hh>
#include <G4HEAntiLambdaInelastic.hh>
#include <G4HESigmaPlusInelastic.hh>
#include <G4HESigmaMinusInelastic.hh>
#include <G4HEAntiSigmaPlusInelastic.hh>
#include <G4HEAntiSigmaMinusInelastic.hh>
#include <G4HEXiZeroInelastic.hh>
#include <G4HEXiMinusInelastic.hh>
#include <G4HEAntiXiZeroInelastic.hh>
#include <G4HEAntiXiMinusInelastic.hh>
#include <G4HEOmegaMinusInelastic.hh>
#include <G4HEAntiOmegaMinusInelastic.hh>

// Stopping processes
#include <G4AntiProtonAnnihilationAtRest.hh>
#include <G4AntiNeutronAnnihilationAtRest.hh>

#ifdef TRIUMF_STOP_PIMINUS
#include <G4PionMinusAbsorptionAtRest.hh>
#else
#include <G4PiMinusAbsorptionAtRest.hh>
#endif
#ifdef TRIUMF_STOP_KMINUS
#include <G4KaonMinusAbsorption.hh>
#else
#include <G4KaonMinusAbsorptionAtRest.hh>
#endif

class TG4PhysicsConstructorHadron: public G4VPhysicsConstructor
{
  typedef G4std::vector<G4VProcess*>  ProcessVector;

  public:
    TG4PhysicsConstructorHadron(const G4String& name = "Hadron");
    // --> protected
    // TG4PhysicsConstructorHadron(const TG4PhysicsConstructorHadron& right);
    virtual ~TG4PhysicsConstructorHadron();

  protected:
    TG4PhysicsConstructorHadron(const TG4PhysicsConstructorHadron& right);
    
    // operators
    TG4PhysicsConstructorHadron& operator=(
                                const TG4PhysicsConstructorHadron& right);
  
    // methods
          // construct particle and physics
    virtual void ConstructParticle();
    virtual void ConstructProcess();

    // data members
         // Elastic Process
    G4HadronElasticProcess  fElasticProcess;
    G4LElastic*             fElasticModel;
  
         // Pi + 
    G4PionPlusInelasticProcess fPionPlusInelastic;
    G4LEPionPlusInelastic*     fLEPionPlusModel;
    G4HEPionPlusInelastic*     fHEPionPlusModel;
    G4MultipleScattering       fPionPlusMult;
    G4hIonisation              fPionPlusIonisation;

         // Pi -
    G4PionMinusInelasticProcess  fPionMinusInelastic;
    G4LEPionMinusInelastic*      fLEPionMinusModel;
    G4HEPionMinusInelastic*      fHEPionMinusModel;
    G4MultipleScattering         fPionMinusMult;
    G4hIonisation                fPionMinusIonisation;
#ifdef TRIUMF_STOP_PIMINUS
    G4PionMinusAbsorptionAtRest  fPionMinusAbsorption;
#else
    G4PiMinusAbsorptionAtRest    fPionMinusAbsorption;
#endif

         // K + 
    G4KaonPlusInelasticProcess  fKaonPlusInelastic;
    G4LEKaonPlusInelastic*      fLEKaonPlusModel;
    G4HEKaonPlusInelastic*      fHEKaonPlusModel;
    G4MultipleScattering        fKaonPlusMult;
    G4hIonisation               fKaonPlusIonisation;

         // K -
    G4KaonMinusInelasticProcess  fKaonMinusInelastic;
    G4LEKaonMinusInelastic*      fLEKaonMinusModel;
    G4HEKaonMinusInelastic*      fHEKaonMinusModel;
    G4MultipleScattering         fKaonMinusMult;
    G4hIonisation                fKaonMinusIonisation;
#ifdef TRIUMF_STOP_KMINUS
    G4KaonMinusAbsorption        fKaonMinusAbsorption;
#else
    G4PiMinusAbsorptionAtRest    fKaonMinusAbsorption;
#endif

        // K0L
    G4KaonZeroLInelasticProcess  fKaonZeroLInelastic;
    G4LEKaonZeroLInelastic*      fLEKaonZeroLModel;
    G4HEKaonZeroInelastic*       fHEKaonZeroLModel;

        // K0S
    G4KaonZeroSInelasticProcess  fKaonZeroSInelastic;
    G4LEKaonZeroSInelastic*      fLEKaonZeroSModel;
    G4HEKaonZeroInelastic*       fHEKaonZeroSModel;

        // Proton
    G4ProtonInelasticProcess  fProtonInelastic;
    G4LEProtonInelastic*      fLEProtonModel;
    G4HEProtonInelastic*      fHEProtonModel;
    G4MultipleScattering      fProtonMult;
    G4hIonisation             fProtonIonisation;
 
        // anti-proton
    G4AntiProtonInelasticProcess    fAntiProtonInelastic;
    G4LEAntiProtonInelastic*        fLEAntiProtonModel;
    G4HEAntiProtonInelastic*        fHEAntiProtonModel;
    G4MultipleScattering            fAntiProtonMult;
    G4hIonisation                   fAntiProtonIonisation;
    G4AntiProtonAnnihilationAtRest  fAntiProtonAnnihilation;
    
       // neutron
    G4NeutronInelasticProcess  fNeutronInelastic;
    G4LENeutronInelastic*      fLENeutronModel;
    G4HENeutronInelastic*      fHENeutronModel;
    G4HadronFissionProcess     fNeutronFission;
    G4LFission*                fNeutronFissionModel;
    G4HadronCaptureProcess     fNeutronCapture;
    G4LCapture*                fNeutronCaptureModel;

       // anti-neutron
    G4AntiNeutronInelasticProcess    fAntiNeutronInelastic;
    G4LEAntiNeutronInelastic*        fLEAntiNeutronModel;
    G4HEAntiNeutronInelastic*        fHEAntiNeutronModel;
    G4AntiNeutronAnnihilationAtRest  fAntiNeutronAnnihilation;
     
       // Lambda
    G4LambdaInelasticProcess  fLambdaInelastic;
    G4LELambdaInelastic*      fLELambdaModel;
    G4HELambdaInelastic*      fHELambdaModel;
  
       // AntiLambda
    G4AntiLambdaInelasticProcess  fAntiLambdaInelastic;
    G4LEAntiLambdaInelastic*      fLEAntiLambdaModel;
    G4HEAntiLambdaInelastic*      fHEAntiLambdaModel;
  
       // SigmaMinus
    G4SigmaMinusInelasticProcess  fSigmaMinusInelastic;
    G4LESigmaMinusInelastic*      fLESigmaMinusModel;
    G4HESigmaMinusInelastic*      fHESigmaMinusModel;
    G4MultipleScattering          fSigmaMinusMult;
    G4hIonisation                 fSigmaMinusIonisation;
  
       // AntiSigmaMinus
    G4AntiSigmaMinusInelasticProcess  fAntiSigmaMinusInelastic;
    G4LEAntiSigmaMinusInelastic*      fLEAntiSigmaMinusModel;
    G4HEAntiSigmaMinusInelastic*      fHEAntiSigmaMinusModel;
    G4MultipleScattering              fAntiSigmaMinusMult;
    G4hIonisation                     fAntiSigmaMinusIonisation;
   
       // SigmaPlus
    G4SigmaPlusInelasticProcess  fSigmaPlusInelastic;
    G4LESigmaPlusInelastic*      fLESigmaPlusModel;
    G4HESigmaPlusInelastic*      fHESigmaPlusModel;
    G4MultipleScattering         fSigmaPlusMult;
    G4hIonisation                fSigmaPlusIonisation;
  
       // AntiSigmaPlus
    G4AntiSigmaPlusInelasticProcess  fAntiSigmaPlusInelastic;
    G4LEAntiSigmaPlusInelastic*      fLEAntiSigmaPlusModel;
    G4HEAntiSigmaPlusInelastic*      fHEAntiSigmaPlusModel;
    G4MultipleScattering             fAntiSigmaPlusMult;
    G4hIonisation                    fAntiSigmaPlusIonisation;
  
      // XiZero
    G4XiZeroInelasticProcess  fXiZeroInelastic;
    G4LEXiZeroInelastic*      fLEXiZeroModel;
    G4HEXiZeroInelastic*      fHEXiZeroModel;
  
      // AntiXiZero
    G4AntiXiZeroInelasticProcess  fAntiXiZeroInelastic;
    G4LEAntiXiZeroInelastic*      fLEAntiXiZeroModel;
    G4HEAntiXiZeroInelastic*      fHEAntiXiZeroModel;
  
      // XiMinus
    G4XiMinusInelasticProcess  fXiMinusInelastic;
    G4LEXiMinusInelastic*      fLEXiMinusModel;
    G4HEXiMinusInelastic*      fHEXiMinusModel;
    G4MultipleScattering       fXiMinusMult;
    G4hIonisation              fXiMinusIonisation;

      // AntiXiMinus
    G4AntiXiMinusInelasticProcess  fAntiXiMinusInelastic;
    G4LEAntiXiMinusInelastic*      fLEAntiXiMinusModel;
    G4HEAntiXiMinusInelastic*      fHEAntiXiMinusModel;
    G4MultipleScattering           fAntiXiMinusMult;
    G4hIonisation                  fAntiXiMinusIonisation;
  
      // OmegaMinus
    G4OmegaMinusInelasticProcess  fOmegaMinusInelastic;
    G4LEOmegaMinusInelastic*      fLEOmegaMinusModel;
    G4HEOmegaMinusInelastic*      fHEOmegaMinusModel;
    G4MultipleScattering          fOmegaMinusMult;
    G4hIonisation                 fOmegaMinusIonisation;
   
      // AntiOmegaMinus
    G4AntiOmegaMinusInelasticProcess  fAntiOmegaMinusInelastic;
    G4LEAntiOmegaMinusInelastic*      fLEAntiOmegaMinusModel;
    G4HEAntiOmegaMinusInelastic*      fHEAntiOmegaMinusModel;
    G4MultipleScattering              fAntiOmegaMinusMult;
    G4hIonisation                     fAntiOmegaMinusIonisation;
    
      // Other
    ProcessVector  fOtherProcesses;  
    

  private:
    // methods
    void ConstructProcessForPionPlus();
    void ConstructProcessForPionMinus();
    void ConstructProcessForKaonPlus();
    void ConstructProcessForKaonMinus();
    void ConstructProcessForKaonZeroLong();
    void ConstructProcessForKaonZeroShort();
    void ConstructProcessForProton();
    void ConstructProcessForAntiProton();
    void ConstructProcessForNeutron();
    void ConstructProcessForAntiNeutron();
    void ConstructProcessForLambda();
    void ConstructProcessForAntiLambda();
    void ConstructProcessForSigmaMinus();
    void ConstructProcessForAntiSigmaMinus();
    void ConstructProcessForSigmaPlus();
    void ConstructProcessForAntiSigmaPlus();
    void ConstructProcessForXiMinus();
    void ConstructProcessForAntiXiMinus();
    void ConstructProcessForXiZero();
    void ConstructProcessForAntiXiZero();
    void ConstructProcessForOmegaMinus();
    void ConstructProcessForAntiOmegaMinus();
    void ConstructProcessForOther();
};

#endif //TG4_PHYSICS_CONSTRUCTOR_HADRON_H

