// $Id$
// Category: physics
//
// According to corresponding part of:
// ExN04PhysicsList.cc,v 1.7 1999/12/15 14:49:26 gunter
// GEANT4 tag Name: geant4-01-01

#include "TG4PhysicsConstructorHadron.h"

#include <G4ParticleDefinition.hh>
#include <G4ProcessManager.hh>

// Hadron Processes
#include <G4HadronElasticProcess.hh>
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

TG4PhysicsConstructorHadron::TG4PhysicsConstructorHadron(const G4String& name)
  : G4VPhysicsConstructor(name)
{
//
  SetVerboseLevel(1);
}

TG4PhysicsConstructorHadron::~TG4PhysicsConstructorHadron() {
//
}

// protected methods

void TG4PhysicsConstructorHadron::ConstructParticle()
{
// The particles are constructed in the 
// TG4ModularPhysicsList.
// ---
}

void TG4PhysicsConstructorHadron::ConstructProcess()
{
// Makes discrete physics processes for the hadrons, at present limited
// to those particles with GHEISHA interactions (INTRC > 0).
// The processes are: Elastic scattering and Inelastic scattering.
//
// F.W.Jones  09-JUL-1998
// ---

   G4HadronElasticProcess* theElasticProcess = 
                                    new G4HadronElasticProcess;
   G4LElastic* theElasticModel = new G4LElastic;
   theElasticProcess->RegisterMe(theElasticModel);

   theParticleIterator->reset();
   while ((*theParticleIterator)()) {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "pi+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4PionPlusInelasticProcess* theInelasticProcess = 
                                new G4PionPlusInelasticProcess("inelastic");
         G4LEPionPlusInelastic* theLEInelasticModel = 
                                new G4LEPionPlusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEPionPlusInelastic* theHEInelasticModel = 
                                new G4HEPionPlusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "pi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4PionMinusInelasticProcess* theInelasticProcess = 
                                new G4PionMinusInelasticProcess("inelastic");
         G4LEPionMinusInelastic* theLEInelasticModel = 
                                new G4LEPionMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEPionMinusInelastic* theHEInelasticModel = 
                                new G4HEPionMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
#ifdef TRIUMF_STOP_PIMINUS
         pmanager->AddRestProcess(new G4PionMinusAbsorptionAtRest, ordDefault);
#else
         G4String prcNam;
         pmanager->AddRestProcess(
           new G4PiMinusAbsorptionAtRest(
                prcNam="PiMinusAbsorptionAtRest"), ordDefault);
#endif
      }
      else if (particleName == "kaon+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonPlusInelasticProcess* theInelasticProcess = 
                                  new G4KaonPlusInelasticProcess("inelastic");
         G4LEKaonPlusInelastic* theLEInelasticModel = 
                                  new G4LEKaonPlusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEKaonPlusInelastic* theHEInelasticModel = 
                                  new G4HEKaonPlusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0S") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonZeroSInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroSInelasticProcess("inelastic");
         G4LEKaonZeroSInelastic* theLEInelasticModel = 
                             new G4LEKaonZeroSInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEKaonZeroInelastic* theHEInelasticModel = 
                             new G4HEKaonZeroInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0L") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonZeroLInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroLInelasticProcess("inelastic");
         G4LEKaonZeroLInelastic* theLEInelasticModel = 
                             new G4LEKaonZeroLInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEKaonZeroInelastic* theHEInelasticModel = 
                             new G4HEKaonZeroInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonMinusInelasticProcess* theInelasticProcess = 
                                 new G4KaonMinusInelasticProcess("inelastic");
         G4LEKaonMinusInelastic* theLEInelasticModel = 
                                 new G4LEKaonMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEKaonMinusInelastic* theHEInelasticModel = 
                                 new G4HEKaonMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
#ifdef TRIUMF_STOP_KMINUS
         pmanager->AddRestProcess(new G4KaonMinusAbsorption, ordDefault);
#else
         pmanager->AddRestProcess(new G4KaonMinusAbsorptionAtRest, ordDefault);
#endif
      }
      else if (particleName == "proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4ProtonInelasticProcess* theInelasticProcess = 
                                    new G4ProtonInelasticProcess("inelastic");
         G4LEProtonInelastic* theLEInelasticModel = new G4LEProtonInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEProtonInelastic* theHEInelasticModel = new G4HEProtonInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiProtonInelasticProcess* theInelasticProcess = 
                                new G4AntiProtonInelasticProcess("inelastic");
         G4LEAntiProtonInelastic* theLEInelasticModel = 
                                new G4LEAntiProtonInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiProtonInelastic* theHEInelasticModel = 
                                new G4HEAntiProtonInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "neutron") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4NeutronInelasticProcess* theInelasticProcess = 
                                    new G4NeutronInelasticProcess("inelastic");
         G4LENeutronInelastic* theLEInelasticModel = 
                                    new G4LENeutronInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HENeutronInelastic* theHEInelasticModel = 
                                    new G4HENeutronInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }  
      else if (particleName == "anti_neutron") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiNeutronInelasticProcess* theInelasticProcess = 
                               new G4AntiNeutronInelasticProcess("inelastic");
         G4LEAntiNeutronInelastic* theLEInelasticModel = 
                               new G4LEAntiNeutronInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiNeutronInelastic* theHEInelasticModel = 
                               new G4HEAntiNeutronInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "lambda") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4LambdaInelasticProcess* theInelasticProcess = 
                                    new G4LambdaInelasticProcess("inelastic");
         G4LELambdaInelastic* theLEInelasticModel = new G4LELambdaInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HELambdaInelastic* theHEInelasticModel = new G4HELambdaInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_lambda") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiLambdaInelasticProcess* theInelasticProcess = 
                                new G4AntiLambdaInelasticProcess("inelastic");
         G4LEAntiLambdaInelastic* theLEInelasticModel = 
                                new G4LEAntiLambdaInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiLambdaInelastic* theHEInelasticModel = 
                                new G4HEAntiLambdaInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "sigma+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4SigmaPlusInelasticProcess* theInelasticProcess = 
                                 new G4SigmaPlusInelasticProcess("inelastic");
         G4LESigmaPlusInelastic* theLEInelasticModel = 
                                 new G4LESigmaPlusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HESigmaPlusInelastic* theHEInelasticModel = 
                                 new G4HESigmaPlusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "sigma-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4SigmaMinusInelasticProcess* theInelasticProcess = 
                                 new G4SigmaMinusInelasticProcess("inelastic");
         G4LESigmaMinusInelastic* theLEInelasticModel = 
                                 new G4LESigmaMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HESigmaMinusInelastic* theHEInelasticModel = 
                                 new G4HESigmaMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_sigma+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiSigmaPlusInelasticProcess* theInelasticProcess = 
                             new G4AntiSigmaPlusInelasticProcess("inelastic");
         G4LEAntiSigmaPlusInelastic* theLEInelasticModel = 
                                 new G4LEAntiSigmaPlusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiSigmaPlusInelastic* theHEInelasticModel = 
                                 new G4HEAntiSigmaPlusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_sigma-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiSigmaMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiSigmaMinusInelasticProcess("inelastic");
         G4LEAntiSigmaMinusInelastic* theLEInelasticModel = 
                                 new G4LEAntiSigmaMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiSigmaMinusInelastic* theHEInelasticModel = 
                                 new G4HEAntiSigmaMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "xi0") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4XiZeroInelasticProcess* theInelasticProcess = 
                            new G4XiZeroInelasticProcess("inelastic");
         G4LEXiZeroInelastic* theLEInelasticModel = 
                                 new G4LEXiZeroInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEXiZeroInelastic* theHEInelasticModel = 
                                 new G4HEXiZeroInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "xi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4XiMinusInelasticProcess* theInelasticProcess = 
                            new G4XiMinusInelasticProcess("inelastic");
         G4LEXiMinusInelastic* theLEInelasticModel = 
                                 new G4LEXiMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEXiMinusInelastic* theHEInelasticModel = 
                                 new G4HEXiMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_xi0") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiXiZeroInelasticProcess* theInelasticProcess = 
                            new G4AntiXiZeroInelasticProcess("inelastic");
         G4LEAntiXiZeroInelastic* theLEInelasticModel = 
                                 new G4LEAntiXiZeroInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiXiZeroInelastic* theHEInelasticModel = 
                                 new G4HEAntiXiZeroInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_xi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiXiMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiXiMinusInelasticProcess("inelastic");
         G4LEAntiXiMinusInelastic* theLEInelasticModel = 
                                 new G4LEAntiXiMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiXiMinusInelastic* theHEInelasticModel = 
                                 new G4HEAntiXiMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "deuteron") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4DeuteronInelasticProcess* theInelasticProcess = 
                            new G4DeuteronInelasticProcess("inelastic");
         G4LEDeuteronInelastic* theLEInelasticModel = 
                                 new G4LEDeuteronInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "triton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4TritonInelasticProcess* theInelasticProcess = 
                            new G4TritonInelasticProcess("inelastic");
         G4LETritonInelastic* theLEInelasticModel = 
                                 new G4LETritonInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "alpha") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AlphaInelasticProcess* theInelasticProcess = 
                            new G4AlphaInelasticProcess("inelastic");
         G4LEAlphaInelastic* theLEInelasticModel = 
                                 new G4LEAlphaInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "omega-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4OmegaMinusInelasticProcess* theInelasticProcess = 
                            new G4OmegaMinusInelasticProcess("inelastic");
         G4LEOmegaMinusInelastic* theLEInelasticModel = 
                                 new G4LEOmegaMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEOmegaMinusInelastic* theHEInelasticModel = 
                                 new G4HEOmegaMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_omega-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiOmegaMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiOmegaMinusInelasticProcess("inelastic");
         G4LEAntiOmegaMinusInelastic* theLEInelasticModel = 
                                 new G4LEAntiOmegaMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiOmegaMinusInelastic* theHEInelasticModel = 
                                 new G4HEAntiOmegaMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
   }

   if (verboseLevel>0)
     G4cout << "### Hadron physics constructed." << G4endl;
}

