#ifndef ALIANALYSISTASKSE_H

#include <Riostream.h>
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TKey.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TH1I.h>
#include <TArrayF.h>
#include <TArrayD.h>
#include <TVector2.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TList.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TMath.h>
#include <TSystem.h>
#include <TInterpreter.h>
#include "AliAnalysisTask.h"
#include "AliCentrality.h"
#include "AliStack.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliParticleContainer.h"
#include "AliInputEventHandler.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#endif

#include <string>
#include <time.h>
#include <TRandom3.h>
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"
#include "AliLog.h"
#include <AliEmcalJet.h>
#include <AliPicoTrack.h>
#include "AliVEventHandler.h"
#include "AliVParticle.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisUtils.h"
#include "AliRhoParameter.h"
#include "TVector3.h"
#include "AliVVertex.h"
#include "AliExternalTrackParam.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "AliGenDPMjetEventHeader.h"
#include "AliJetContainer.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisTaskV0Mnorm.h"
#include "AliHeader.h"
#include "AliRunLoader.h"
#include "AliVVZERO.h"
#include "AliAODZDC.h"
#include "AliVZDC.h"
#include "AliAnalysisEmcalJetHelperEA.h"
#include "AliAnalysisTaskEmcalEmbeddingHelper.h"

//#include "AliEmcalDownscaleFactorsOCDB.h"
//#include "AliEmcalAnalysisFactory.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskV0Mnorm)

using namespace PWGJE::EMCALJetTasks;
using namespace std;

// ANALYSIS OF V0M IN PP 13 TeV
// Author Filip Krizek   (12 Mar 2024)

//________________________________________________________________________________________

AliAnalysisTaskV0Mnorm::AliAnalysisTaskV0Mnorm():
AliAnalysisTaskEmcalJet("AliAnalysisTaskV0Mnorm", kTRUE),
fUseDefaultVertexCut(1),
fUsePileUpCut(1),
fMyParticleContainerName(""),
fParticleContainerPartLevel(0x0),
fHelperEA(0x0),
fHelperClass(0), 
fInitializedLocal(0),  
hEA_correlations(nullptr)
{
   //default constructor

   fHelperEA = new PWGJE::EMCALJetTasks::AliAnalysisEmcalJetHelperEA();
}

//________________________________________________________________________
AliAnalysisTaskV0Mnorm::AliAnalysisTaskV0Mnorm(const char *name):
AliAnalysisTaskEmcalJet(name,kTRUE),
fUseDefaultVertexCut(1),
fUsePileUpCut(1),
fMyParticleContainerName(""),
fParticleContainerPartLevel(0x0),
fHelperEA(0x0),
fHelperClass(0), //FFF
fInitializedLocal(0),  //FFF
hEA_correlations(nullptr)
{
   //Constructor


   fHelperEA = new PWGJE::EMCALJetTasks::AliAnalysisEmcalJetHelperEA();

   DefineOutput(1, TList::Class());
}

/**
 * This function adds the task to the analysis manager. Often, this function is called
 * by an AddTask C macro. However, by compiling the code, it ensures that we do not
 * have to deal with difficulties caused by CINT.
 */

//_____________________________________________________________________________________
AliAnalysisTaskV0Mnorm*  AliAnalysisTaskV0Mnorm::AddTaskV0Mnorm(
  const char* mcpariclearraynamePartMC,
  Bool_t      useVertexCut,
  Bool_t      usePileUpCut,
  const char* suffix
){
   // Get the pointer to the existing analysis manager via the static access method.
   //==============================================================================


   // #### DEFINE MANAGER AND DATA CONTAINER NAMES
   AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
   if(!manager){
      ::Error("AliAnalysisTaskV0Mnorm.cxx", "No analysis manager to connect to.");
      return NULL;
   }

   //__________________________________________________________________________________
   // #### DEFINE MY ANALYSIS TASK

   TString myContName("V0normCorrelation");
   myContName.Append(suffix);

   AliAnalysisTaskV0Mnorm *task = new AliAnalysisTaskV0Mnorm(myContName.Data());
   //for PYTHIA
   task->SetIsPythia(kTRUE);  //NECESSARY IN ORDER TO FILL XSEC AND TRIALS
   task->SetMakeGeneralHistograms(kTRUE); //NECESSARY IN ORDER TO FILL XSEC AND TRIALS

   //inspired by AliAnalysisTaskEmcalQGTagging
   //_____________________________________________
   //TRACK/PARTICLE CONTAINTERS
   AliParticleContainer *trackContTrue    = 0x0; //mc particle container on  particle level for jets

   trackContTrue = task->AddMCParticleContainer(mcpariclearraynamePartMC); //particle level MC particles
   trackContTrue->SetClassName("AliAODMCParticle");
   trackContTrue->SetMinPt(1e-3); // KA: old value was 0.15 GeV/c
   trackContTrue->SetEtaLimits(-5.1,5.1); //V0 eta range


   // #### Task configuration
   task->SetUsePileUpCut(usePileUpCut);
   task->SetUseDefaultVertexCut(useVertexCut);
   task->SelectCollisionCandidates(AliVEvent::kAny);

   task->SetMCParticleContainerName(mcpariclearraynamePartMC);

   task->SetDebugLevel(0); //No debug messages 0

   // output container
   AliAnalysisDataContainer *contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChJetSpectra%s", AliAnalysisManager::GetCommonFileName(), myContName.Data()));


   // #### ADD ANALYSIS TASK
   manager->AddTask(task);
   manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
   manager->ConnectOutput(task, 1, contHistos);

   return task;
}
//_____________________________________________________________________________________

Bool_t AliAnalysisTaskV0Mnorm::IsEventInAcceptance(AliVEvent* event){
   //EVENT SELECTION RECONSTRUCTED DATA

   if(!event) return kFALSE;

   //incomplete DAQ events rejection Run2 data 2015
   // https://twiki.cern.ch/twiki/bin/view/ALICE/PWGPPEvSelRun2pp
   Bool_t bIncompleteDAQ = event->IsIncompleteDAQ();
   if(bIncompleteDAQ){
      return kFALSE;
   }
   //___________________________________________________
   //TEST PILE UP
   if(fUsePileUpCut){
      if(!fHelperClass || fHelperClass->IsPileUpEvent(event)){
         return kFALSE;
      }

      if(!fHelperClass || fHelperClass->IsSPDClusterVsTrackletBG(event)){
         return kFALSE;
      }

      if(event->IsPileupFromSPDInMultBins()){
         return kFALSE;
      }
   }


   //___________________________________________________
   //VERTEX CUT

   if(fUseDefaultVertexCut){
      if(!fHelperClass || !fHelperClass->IsVertexSelected2013pA(event)){  
         return kFALSE;
      }
   }

   if(TMath::Abs(event->GetPrimaryVertex()->GetZ()) > 10.){  //VERTEX CUT
      return kFALSE;
   }


  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskV0Mnorm::ExecOnceLocal(){
   // Initialization of jet containers done in  AliAnalysisTaskEmcalJet::ExecOnce()
   //Read arrays of jets and tracks
   fInitializedLocal = kTRUE;

   // Initialize helper class (for vertex selection & pile up correction)
   fHelperClass = new AliAnalysisUtils();
   fHelperClass->SetCutOnZVertexSPD(kFALSE); // kFALSE: no cut; kTRUE: |zvtx-SPD - zvtx-TPC|<0.5cm

   return;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskV0Mnorm::FillHistograms(){
   // executed in each event
   //called in AliAnalysisTaskEmcal::UserExec(Option_t *)
   //   Analyze the event and Fill histograms

   if(!InputEvent()){
      AliError("??? Event pointer == 0 ???");
      return kFALSE;
   }

   //Execute only once:  Get tracks, jets from arrays if not already given
   if(!fInitializedLocal) ExecOnceLocal();


   AliGenEventHeader* mcHeader = NULL;
   AliAODMCHeader* aodMCH = NULL;

   //+++++++++++++++++++++++++++++ check MC z vertex position ++++++++++++++++++++++++++
   if(MCEvent()){

      mcHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
      if(!mcHeader){
         // Check if AOD
          aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

          if(aodMCH){
             for(UInt_t i = 0; i<aodMCH->GetNCocktailHeaders(); i++){
               mcHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
               if(mcHeader) break;
            }
         }
      }
   }

   if(mcHeader){
      TArrayF pyVtx;
      mcHeader->PrimaryVertex(pyVtx);
      if(TMath::Abs(pyVtx[2]) > 10.) return kTRUE; //skip events with particle level vertex out of +-10 cm
   }

   //_________________________________________________________
   //READ  TRACK AND JET CONTAINERS
   //Container operations   http://alidoc.cern.ch/AliPhysics/master/READMEcontainers.html#emcalContainerIterateTechniques

   fParticleContainerPartLevel = GetParticleContainer(0); //pythia particle level particles

   //________________________________________________________________
   //DATA ANALYSIS DETECTOR LEVEL

   //Check Reconstructed event vertex and pileup
   if(!IsEventInAcceptance(InputEvent())) return kFALSE; //post data is in UserExec

   //_________________________________________________________________

   AliVParticle *mcParticle = NULL; //mc particle
   //Detector level
   double EA_V0A_det = .0, EA_V0C_det = .0, EA_V0M_det = .0, meanV0M_det = .0, EA_V0Mnorm_det = .0;

   AliVVZERO *vzeroAOD = InputEvent()->GetVZEROData();
   if(vzeroAOD){
      EA_V0A_det     = vzeroAOD->GetMTotV0A();
      EA_V0C_det     = vzeroAOD->GetMTotV0C();
      EA_V0M_det     = EA_V0A_det + EA_V0C_det;
      meanV0M_det    = fHelperEA->GetV0MDetLevel();
      EA_V0Mnorm_det = EA_V0M_det/meanV0M_det;
   }

   //Particle level
   double EA_V0A_part = .0, EA_V0C_part = .0, EA_V0M_part = .0, meanV0M_part = .0, EA_V0Mnorm_part = .0;
   if(fParticleContainerPartLevel){
      for(auto mcPartIterator : fParticleContainerPartLevel->accepted_momentum()){
         mcParticle = mcPartIterator.second;
         if(!mcParticle) continue;

         if(mcParticle->Charge()){
            if((static_cast<AliAODMCParticle*>(mcParticle))->IsPhysicalPrimary()){
               //get particle level charged particles multiplicities in V0A and V0C
               if(-3.7 < mcParticle->Eta() && mcParticle->Eta() < -1.7) EA_V0C_part++;
               if( 2.8 < mcParticle->Eta() && mcParticle->Eta() < 5.1)  EA_V0A_part++;
            }
         }
      }
   } 

   // Coincidence condition (neglegible impact at high multiplicity)
   if(EA_V0C_part > .0 && EA_V0A_part > .0){ 
      if(EA_V0C_det > .0 && EA_V0A_det > .0){

         //combined V0 multiplicities particle level
         EA_V0M_part     = EA_V0A_part + EA_V0C_part;
         meanV0M_part    = fHelperEA->GetV0MPartLevel();
         EA_V0Mnorm_part = EA_V0M_part/meanV0M_part; 
         
         hEA_correlations->Fill(EA_V0Mnorm_part, EA_V0Mnorm_det);
      }
   }

   return kTRUE; 
}
//________________________________________________________________________
void AliAnalysisTaskV0Mnorm::Terminate(Option_t *){
   //Treminate
   PostData(1, fOutput);

   // Mandatory
   fOutput = dynamic_cast<AliEmcalList*> (GetOutputData(1)); // '1' refers to the output slot
   if(!fOutput) {
      printf("ERROR: Output list not available\n");
      return;
   }
}

//________________________________________________________________________
AliAnalysisTaskV0Mnorm::~AliAnalysisTaskV0Mnorm(){
   // Destructor. Clean-up the output list, but not the histograms that are put inside
   // (the list is owner and will clean-up these histograms). Protect in PROOF case.
   if(fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
   }
   delete fHelperClass;
   delete fHelperEA;
}
//________________________________________________________________________
void AliAnalysisTaskV0Mnorm::UserCreateOutputObjects(){
  // called once to create user defined output objects like histograms, plots etc.
  // and to put it on the output list.
  // Note: Saving to file with e.g. OpenFile(0) is must be before creating other objects.
  //fOutput TList defined in the mother class

   AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   //__________________________________________________________

   // EA correlations   
   hEA_correlations = new TH2D("EA_correlations_MC", "EA_correlations_MC", 300, 0.0, 15.0, 300, 0.0, 15.0);
   fOutput->Add((TH2D*) hEA_correlations);


   // =========== Switch on Sumw2 for all histos ===========
   for(Int_t i=0; i<fOutput->GetEntries(); i++){
      TH1 *h1 = dynamic_cast<TH1*>(fOutput->At(i));
      if(h1){
         h1->Sumw2();
         continue;
      }
      THnSparse *hn = dynamic_cast<THnSparse*>(fOutput->At(i));
      if(hn){
         hn->Sumw2();
      }
   }
   TH1::AddDirectory(oldStatus);


   PostData(1, fOutput);
}
//________________________________________________________________________
Bool_t AliAnalysisTaskV0Mnorm::Run(){
   // Run analysis code here, if needed. It will be executed before FillHistograms().

   return kTRUE;
}

