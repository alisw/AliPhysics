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
    G4HadronElasticProcess  fElasticProcess;            //hadron elastic process
    G4LElastic*             fElasticModel;              //elastic model
  
         // Pi + 
    G4PionPlusInelasticProcess fPionPlusInelastic;      //pi+ inel process
    G4LEPionPlusInelastic*     fLEPionPlusModel;        //pi+ LE inel model
    G4HEPionPlusInelastic*     fHEPionPlusModel;        //pi+ HE inel model
    G4MultipleScattering       fPionPlusMult;           //pi+ msc
    G4hIonisation              fPionPlusIonisation;     //pi+ ionisation

         // Pi -
    G4PionMinusInelasticProcess  fPionMinusInelastic;   //pi- inel process
    G4LEPionMinusInelastic*      fLEPionMinusModel;     //pi- LE inel model
    G4HEPionMinusInelastic*      fHEPionMinusModel;     //pi- HE inel model
    G4MultipleScattering         fPionMinusMult;        //pi- msc
    G4hIonisation                fPionMinusIonisation;  //pi- ionisation
#ifdef TRIUMF_STOP_PIMINUS
    G4PionMinusAbsorptionAtRest  fPionMinusAbsorption;  //pi- absorption
#else
    G4PiMinusAbsorptionAtRest    fPionMinusAbsorption;  //pi- absorption
#endif

         // K + 
    G4KaonPlusInelasticProcess  fKaonPlusInelastic;     //kaon+ inel process
    G4LEKaonPlusInelastic*      fLEKaonPlusModel;       //kaon+ LE inel model
    G4HEKaonPlusInelastic*      fHEKaonPlusModel;       //kaon+ HE inel model
    G4MultipleScattering        fKaonPlusMult;          //kaon+ msc
    G4hIonisation               fKaonPlusIonisation;    //kaon+ ionisation

         // K -
    G4KaonMinusInelasticProcess  fKaonMinusInelastic;   //kaon- inel process
    G4LEKaonMinusInelastic*      fLEKaonMinusModel;     //kaon- LE inel model
    G4HEKaonMinusInelastic*      fHEKaonMinusModel;     //kaon- HE inel model
    G4MultipleScattering         fKaonMinusMult;        //kaon- msc
    G4hIonisation                fKaonMinusIonisation;  //kaon- ionisation
#ifdef TRIUMF_STOP_KMINUS
    G4KaonMinusAbsorption        fKaonMinusAbsorption;  //kaon- absorption
#else
    G4PiMinusAbsorptionAtRest    fKaonMinusAbsorption;  //kaon- absorption
#endif

        // K0L
    G4KaonZeroLInelasticProcess  fKaonZeroLInelastic;   //kaon0 inel process
    G4LEKaonZeroLInelastic*      fLEKaonZeroLModel;     //kaon0 LE inel model
    G4HEKaonZeroInelastic*       fHEKaonZeroLModel;     //kaon0 HE inel model

        // K0S
    G4KaonZeroSInelasticProcess  fKaonZeroSInelastic;   //kaon0S inel process
    G4LEKaonZeroSInelastic*      fLEKaonZeroSModel;     //kaon0S LE inel model
    G4HEKaonZeroInelastic*       fHEKaonZeroSModel;     //kaon0S HE inel mode

        // Proton
    G4ProtonInelasticProcess  fProtonInelastic;         //p inel process
    G4LEProtonInelastic*      fLEProtonModel;           //p LE inel model
    G4HEProtonInelastic*      fHEProtonModel;           //p HE inel model
    G4MultipleScattering      fProtonMult;              //p msc
    G4hIonisation             fProtonIonisation;        //p ionisation
 
        // anti-proton
    G4AntiProtonInelasticProcess    fAntiProtonInelastic; //p_bar inel process
    G4LEAntiProtonInelastic*        fLEAntiProtonModel;   //p_bar LE inel model
    G4HEAntiProtonInelastic*        fHEAntiProtonModel;   //p_bar HE inel model
    G4MultipleScattering            fAntiProtonMult;      //p_bar msc
    G4hIonisation                   fAntiProtonIonisation;//p_bar ionisation
    G4AntiProtonAnnihilationAtRest  fAntiProtonAnnihilation;//p_bar annihilation
    
       // neutron
    G4NeutronInelasticProcess  fNeutronInelastic;       //n inel process
    G4LENeutronInelastic*      fLENeutronModel;         //n LE inel model
    G4HENeutronInelastic*      fHENeutronModel;         //n HE inel model
    G4HadronFissionProcess     fNeutronFission;         //n fission
    G4LFission*                fNeutronFissionModel;    //n fission model
    G4HadronCaptureProcess     fNeutronCapture;         //n capture
    G4LCapture*                fNeutronCaptureModel;    //n capture model

       // anti-neutron
    G4AntiNeutronInelasticProcess    fAntiNeutronInelastic;//n_bar inel process
    G4LEAntiNeutronInelastic*        fLEAntiNeutronModel;  //n_bar LE inel model
    G4HEAntiNeutronInelastic*        fHEAntiNeutronModel;  //n_bar HE inel model
    G4AntiNeutronAnnihilationAtRest  fAntiNeutronAnnihilation;//n_bar ionisation
     
       // Lambda
    G4LambdaInelasticProcess  fLambdaInelastic;         //lambda inel process
    G4LELambdaInelastic*      fLELambdaModel;           //lambda LE inel model
    G4HELambdaInelastic*      fHELambdaModel;           //lambda HE inel model
  
       // AntiLambda
    G4AntiLambdaInelasticProcess  fAntiLambdaInelastic; //lambda_bar inel process
    G4LEAntiLambdaInelastic*      fLEAntiLambdaModel;   //lambda_bar LE inel model
    G4HEAntiLambdaInelastic*      fHEAntiLambdaModel;   //lambda_bar HE inel model
  
       // SigmaMinus
    G4SigmaMinusInelasticProcess  fSigmaMinusInelastic; //sigma- inel process
    G4LESigmaMinusInelastic*      fLESigmaMinusModel;   //sigma- LE inel model
    G4HESigmaMinusInelastic*      fHESigmaMinusModel;   //sigma- HE inel model
    G4MultipleScattering          fSigmaMinusMult;      //sigma- msc
    G4hIonisation                 fSigmaMinusIonisation;//sigma- ionisation
  
       // AntiSigmaMinus
    G4AntiSigmaMinusInelasticProcess  fAntiSigmaMinusInelastic; //sigma-_bar inel process
    G4LEAntiSigmaMinusInelastic*      fLEAntiSigmaMinusModel;   //sigma-_bar LE inel model
    G4HEAntiSigmaMinusInelastic*      fHEAntiSigmaMinusModel;   //sigma-_bar HE inel model
    G4MultipleScattering              fAntiSigmaMinusMult;      //sigma-_bar msc
    G4hIonisation                     fAntiSigmaMinusIonisation;//sigma-_bar ionisation
   
       // SigmaPlus
    G4SigmaPlusInelasticProcess  fSigmaPlusInelastic;   //sigma+ inel process
    G4LESigmaPlusInelastic*      fLESigmaPlusModel;     //sigma+ LE inel model
    G4HESigmaPlusInelastic*      fHESigmaPlusModel;     //sigma+ HE inel model
    G4MultipleScattering         fSigmaPlusMult;        //sigma+ msc
    G4hIonisation                fSigmaPlusIonisation;  //sigma+ ionisation
  
       // AntiSigmaPlus
    G4AntiSigmaPlusInelasticProcess  fAntiSigmaPlusInelastic;  //sigma+_bar inel process
    G4LEAntiSigmaPlusInelastic*      fLEAntiSigmaPlusModel;    //sigma+_bar LE inel model
    G4HEAntiSigmaPlusInelastic*      fHEAntiSigmaPlusModel;    //sigma+_bar HE inel model
    G4MultipleScattering             fAntiSigmaPlusMult;       //sigma+_bar msc
    G4hIonisation                    fAntiSigmaPlusIonisation; //sigma+_bar ionisation
  
      // XiZero
    G4XiZeroInelasticProcess  fXiZeroInelastic;        //xi0 inel process
    G4LEXiZeroInelastic*      fLEXiZeroModel;          //xi0 LE inel model
    G4HEXiZeroInelastic*      fHEXiZeroModel;          //xi0 HE inel model
  
      // AntiXiZero
    G4AntiXiZeroInelasticProcess  fAntiXiZeroInelastic;//xi0_bar inel process
    G4LEAntiXiZeroInelastic*      fLEAntiXiZeroModel;  //xi0_bar LE inel model
    G4HEAntiXiZeroInelastic*      fHEAntiXiZeroModel;  //xi0_bar HE inel model
  
      // XiMinus
    G4XiMinusInelasticProcess  fXiMinusInelastic;      //xi- inel process
    G4LEXiMinusInelastic*      fLEXiMinusModel;        //xi- LE inel model
    G4HEXiMinusInelastic*      fHEXiMinusModel;        //xi- HE inel model
    G4MultipleScattering       fXiMinusMult;           //xi- msc
    G4hIonisation              fXiMinusIonisation;     //xi- ionisation

      // AntiXiMinus
    G4AntiXiMinusInelasticProcess  fAntiXiMinusInelastic; //xi-_bar inel process
    G4LEAntiXiMinusInelastic*      fLEAntiXiMinusModel;   //xi-_bar LE inel model
    G4HEAntiXiMinusInelastic*      fHEAntiXiMinusModel;   //xi-_bar HE inel model
    G4MultipleScattering           fAntiXiMinusMult;      //xi-_bar msc
    G4hIonisation                  fAntiXiMinusIonisation;//xi-_bar ionisation
  
      // OmegaMinus
    G4OmegaMinusInelasticProcess  fOmegaMinusInelastic;   //omega- inel process
    G4LEOmegaMinusInelastic*      fLEOmegaMinusModel;     //omega- LE inel model
    G4HEOmegaMinusInelastic*      fHEOmegaMinusModel;     //omega- HE inel model
    G4MultipleScattering          fOmegaMinusMult;        //omega- msc
    G4hIonisation                 fOmegaMinusIonisation;  //omega- ionisation
   
      // AntiOmegaMinus
    G4AntiOmegaMinusInelasticProcess  fAntiOmegaMinusInelastic; //omega-_bar inel process
    G4LEAntiOmegaMinusInelastic*      fLEAntiOmegaMinusModel;   //omega-_bar LE inel model
    G4HEAntiOmegaMinusInelastic*      fHEAntiOmegaMinusModel;   //omega-_bar HE inel model
    G4MultipleScattering              fAntiOmegaMinusMult;      //omega-_bar msc
    G4hIonisation                     fAntiOmegaMinusIonisation;//omega-_bar ionisation
    
      // Other
    ProcessVector  fOtherProcesses; //other process
    

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

