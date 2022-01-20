#ifndef ALIANALYSISTASKSE_H

#include <Riostream.h>
#include <TROOT.h>
#include <TFile.h>
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
#include <TH1I.h>
#include <TArrayF.h>
#include <TArrayD.h>
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
#endif

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
#include "TVector2.h"
#include "TLorentzVector.h"
#include "AliVVertex.h"
#include "AliExternalTrackParam.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "AliGenDPMjetEventHeader.h"
#include "AliJetContainer.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisTaskHJetSpectra.h"
#include "AliHeader.h" 
#include "AliRunLoader.h"  
#include "AliVVZERO.h"
#include "AliAODZDC.h" 
#include "AliVZDC.h"
#include "AliMultSelection.h"
#include "AliCentralitySelectionTask.h"

using namespace std;

// ANALYSIS OF HIGH PT HADRON TRIGGER ASSOCIATED SPECTRUM OF RECOIL JETS IN P+PB
// Author Filip Krizek   (7.Oct. 2015)

ClassImp(AliAnalysisTaskHJetSpectra)
//________________________________________________________________________________________

AliAnalysisTaskHJetSpectra::AliAnalysisTaskHJetSpectra(): 
AliAnalysisTaskEmcalJet("AliAnalysisTaskHJetSpectra", kTRUE),  
 fCollisionSystem(0),  fTypeOfAnal(0),
  fUseDefaultVertexCut(1), fUsePileUpCut(1),  
 fSignalJetRadius(0.4), fSignalJetRadiusSquared(fSignalJetRadius*fSignalJetRadius),
fSignalJetEtaWindow(0.9 - fSignalJetRadius), fTrackEtaWindow(0.9), fMinTrackPt(0.150), fMinJetArea(0.0),  
 fRandom(0), fHelperClass(0), fInitializedLocal(0),
 fDphiCut(TMath::Pi()-0.6), 
fHistEvtSelection(0x0),
fhTrackPhi(0x0), fhTrackEta(0x0), fhJetPhiDet(0x0), fhJetEtaDet(0x0),
fhJetEtaPart(0x0),fhJetPtPart(0x0),
fhVertexZ(0x0), fhVertexXAccept(0x0), fhVertexYAccept(0x0), fhVertexZAccept(0x0), 
fhVertexZMC(0x0), fhVertexZAcceptMC(0x0),
fhCentralityV0A(0x0),fhCentralityZNA(0x0), 
fhVzeroATotMult(0x0), fhZNAEnergy(0x0), 
fhJetPtPartVsJetPtDet(0x0),
fhPhysPrimaryPtDet(0x0),
fhPhysPrimaryPtPart(0x0),
fhPtTrkSecOrFakeDet(0x0),
fhJetAreaVsPt(0x0),                             //dch
fZVertexCut(10.0),fCutPhi(0.6),
fMultSelection(0x0),
fMyTrackContainerName(""),
fMyParticleContainerName(""),
fMyJetContainerName(""),
fMyJetParticleContainerName(""),
fMyKTJetContainerName(""),
fMyKTJetParticleContainerName("")
{
   //default constructor
  
   fTTlow[kRef] = 6.0;
   fTThigh[kRef]= 7.0;
   fTTlow[kSig] = 12.0;
   fTThigh[kSig]= 50.0;

   for(Int_t it=0; it<kTT;it++){ 
      fhTTDet_V0A[it]=NULL;
      fhTTDet_ZNA[it]=NULL;
      fhTTMultDet_V0A[it] = NULL; 
      fhTTMultDet_ZNA[it] = NULL;
      fhTTPart[it]=NULL;
      fhTTMultPart[it]=NULL;



      fhCentralityTT_V0A[it] = NULL;
      fhCentralityTT_ZNA[it] = NULL;
      fTTH_Det[it].resize(0);
      fTTH_Part[it].resize(0);

      fHJetSpecDet_V0A[it]=NULL;
      fHJetSpecDet_ZNA[it]=NULL;
      fHJetSpecPart[it] = NULL;

      fRhoDet[it] = 0;
      fhRhoTT[it] = NULL;
      fhDeltaPt_V0A[it]=NULL;
      fhDeltaPt_ZNA[it]=NULL;

      fhCentralityTT_V0A[it]=NULL;
      fhCentralityTT_ZNA[it]=NULL;

      fTTH_Det[it].resize(0);  
      fTTH_Part[it].resize(0); 
   }


   for(Int_t i=0; i<999; i++){
      frhovec[i] = 0.;
   }

}

//________________________________________________________________________
AliAnalysisTaskHJetSpectra::AliAnalysisTaskHJetSpectra(const char *name) : 
AliAnalysisTaskEmcalJet(name,kTRUE),  
fCollisionSystem(0), fTypeOfAnal(0),
  fUseDefaultVertexCut(1), fUsePileUpCut(1),
fSignalJetRadius(0.4), fSignalJetRadiusSquared(fSignalJetRadius*fSignalJetRadius), 
fSignalJetEtaWindow(0.9 - fSignalJetRadius), fTrackEtaWindow(0.9), fMinTrackPt(0.150), fMinJetArea(0.0),   
fRandom(0), fHelperClass(0), fInitializedLocal(0), 
fDphiCut(TMath::Pi()-0.6), 
fHistEvtSelection(0x0), 
fhTrackPhi(0x0), fhTrackEta(0x0), fhJetPhiDet(0x0), fhJetEtaDet(0x0),
fhJetEtaPart(0x0),fhJetPtPart(0x0),
fhVertexZ(0x0), fhVertexXAccept(0x0), fhVertexYAccept(0x0), fhVertexZAccept(0x0), 
fhVertexZMC(0x0), fhVertexZAcceptMC(0x0),
fhCentralityV0A(0x0), fhCentralityZNA(0x0),
fhVzeroATotMult(0x0), fhZNAEnergy(0x0),
fhJetPtPartVsJetPtDet(0x0),
fhPhysPrimaryPtDet(0x0),
fhPhysPrimaryPtPart(0x0),
fhPtTrkSecOrFakeDet(0x0),
fhJetAreaVsPt(0x0),                             //dch
fZVertexCut(10.0),fCutPhi(0.6),
fMultSelection(0x0),
fMyTrackContainerName(""),
fMyParticleContainerName(""),
fMyJetContainerName(""),
fMyJetParticleContainerName(""),
fMyKTJetContainerName(""),
fMyKTJetParticleContainerName("")
{
//Constructor

   fTTlow[kRef] = 6.0;
   fTThigh[kRef]= 7.0;
   fTTlow[kSig] = 12.0;
   fTThigh[kSig]= 50.0;

   for(Int_t it=0; it<kTT;it++){ 
      fhTTDet_V0A[it]=NULL;
      fhTTDet_ZNA[it]=NULL;
      fhTTMultDet_V0A[it] = NULL; 
      fhTTMultDet_ZNA[it] = NULL;
      fhTTPart[it]=NULL;
      fhTTMultPart[it]=NULL;



      fhCentralityTT_V0A[it] = NULL;
      fhCentralityTT_ZNA[it] = NULL;
      fTTH_Det[it].resize(0);
      fTTH_Part[it].resize(0);

      fHJetSpecDet_V0A[it]=NULL;
      fHJetSpecDet_ZNA[it]=NULL;
      fHJetSpecPart[it] = NULL;

      fRhoDet[it] = 0;
      fhRhoTT[it] = NULL;
      fhDeltaPt_V0A[it]=NULL;
      fhDeltaPt_ZNA[it]=NULL;

      fhCentralityTT_V0A[it]=NULL;
      fhCentralityTT_ZNA[it]=NULL;

      fTTH_Det[it].resize(0);  
      fTTH_Part[it].resize(0); 
   }


   for(Int_t i=0; i<999; i++){
      frhovec[i] = 0.;
   }


   DefineOutput(1, TList::Class());
}
//________________________________________________________________________
void  AliAnalysisTaskHJetSpectra::SetTT(Double_t tlr, Double_t thr,Double_t tls, Double_t ths){
   //set trigger track pT bins
   fTTlow[kRef]  = tlr; 
   fTThigh[kRef] = thr;  
   fTTlow[kSig]  = tls; 
   fTThigh[kSig] = ths;  
}
//________________________________________________________________________
 AliAnalysisTaskHJetSpectra* AliAnalysisTaskHJetSpectra::AddTaskHJetSpectra(
      Int_t collisionSystem, 
      Int_t typeOfAnal,
      const char* jetarrayname, 
      const char* jetarraynamePartMC, 
      const char* trackarrayname, 
      const char* mcpariclearraynamePartMC,
      const char* ktjetarrayname,
      const char* ktjetarraynamePartMC, 
      Double_t    jetRadius,
      UInt_t      trigger,
      Double_t    trackEtaWindow,
      Bool_t      useVertexCut,
      Bool_t      usePileUpCut,
      Double_t    acut
   ){

   Double_t jetEtaRange   = TMath::Abs(trackEtaWindow - jetRadius);
   Double_t jetRadiuskt   = 0.4;  //for all kt jets use fixed jet radius
   Double_t jetEtaRangekt = TMath::Abs(trackEtaWindow - jetRadiuskt);

   // #### DEFINE MANAGER AND DATA CONTAINER NAMES
   AliAnalysisManager *manager = AliAnalysisManager::GetAnalysisManager();
   if(!manager){
      ::Error("AliAnalysisTaskEA.cxx", "No analysis manager to connect to.");
      return NULL;
   }


   TString myContName("");
   myContName = Form("JetAnalysisR%02d_Acut%02d", TMath::Nint(jetRadius*10), TMath::Nint(acut*10));

   AliAnalysisTaskHJetSpectra *task = new AliAnalysisTaskHJetSpectra(myContName.Data());
    if(typeOfAnal == AliAnalysisTaskHJetSpectra::kEff){  
      task->SetIsPythia(kTRUE);  //for PYTHIA NECESSARY IN ORDER TO FILL XSEC AND TRIALS
      task->SetMakeGeneralHistograms(kTRUE); //NECESSARY IN ORDER TO FILL XSEC AND TRIALS
   }

   AliTrackContainer    *trackCont        = 0x0; // detector level track container (or tracks in  combined events when embedding )
   AliParticleContainer *trackContTrue    = 0x0; //mc particle container on  particle level for jets

   trackCont = task->AddTrackContainer(trackarrayname);  //detector level tracks (or combined tracks if embedding)
   trackCont->SetMinPt(0.15);
   trackCont->SetEtaLimits(-trackEtaWindow, trackEtaWindow);

   if(typeOfAnal == AliAnalysisTaskHJetSpectra::kEff){  
      trackContTrue = task->AddMCParticleContainer(mcpariclearraynamePartMC); //particle level MC particles
      trackContTrue->SetClassName("AliAODMCParticle");
      trackContTrue->SetMinPt(0.15);
      trackContTrue->SetEtaLimits(-trackEtaWindow, trackEtaWindow);
   } 


   //_____________________________________________
   //JET CONTAINERS
   AliJetContainer *jetContRec    = 0x0; //AKT jet container with detector level tracks 
   AliJetContainer *jetContTrue   = 0x0; //AKT jet container with mc particle level jets pythia
   AliJetContainer *jetContRecKT  = 0x0; //KT jet container with detector level tracks 
   AliJetContainer *jetContTrueKT = 0x0; //KT jet container with mc particle level jets pythia

   //AKT DETECTOR LEVEL JET    (or combined event jet container when embedding)
   jetContRec   = task->AddJetContainer(jetarrayname,"TPC",jetRadius);

   if(jetContRec) {
      jetContRec->ConnectParticleContainer(trackCont);
      jetContRec->SetPercAreaCut(acut);
      jetContRec->SetMinPt(0.150);
      jetContRec->SetMaxTrackPt(100.);
      jetContRec->SetJetAcceptanceType(AliEmcalJet::kUser);
      jetContRec->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
   }


   //KT DETECTOR LEVEL JET    (or combined event jet container when embedding)
   jetContRecKT   = task->AddJetContainer(ktjetarrayname,"TPC",jetRadiuskt);

   if(jetContRecKT) {
      jetContRecKT->ConnectParticleContainer(trackCont);
      jetContRecKT->SetPercAreaCut(acut);
      jetContRecKT->SetMinPt(0.);
      jetContRecKT->SetMaxTrackPt(100.);
      jetContRecKT->SetJetAcceptanceType(AliEmcalJet::kUser);
      jetContRecKT->SetJetEtaLimits(-jetEtaRangekt,jetEtaRangekt);
   }

   if(typeOfAnal == AliAnalysisTaskHJetSpectra::kEff){ 
      jetContTrue = task->AddJetContainer(jetarraynamePartMC,"TPC",jetRadius);

      if(jetContTrue){
         jetContTrue->ConnectParticleContainer(trackContTrue);
         jetContTrue->SetPercAreaCut(acut);
         jetContTrue->SetMinPt(0.15);
         jetContTrue->SetJetAcceptanceType(AliEmcalJet::kUser);
         jetContTrue->SetJetEtaLimits(-jetEtaRange,jetEtaRange);
      }

      //KT JETS PARTICLE LEVEL
      jetContTrueKT = task->AddJetContainer(ktjetarraynamePartMC,"TPC",jetRadiuskt);

      if(jetContTrueKT){
         jetContTrueKT->ConnectParticleContainer(trackContTrue);
         jetContTrueKT->SetMinPt(0.);
         jetContTrueKT->SetJetAcceptanceType(AliEmcalJet::kUser);
         jetContTrueKT->SetJetEtaLimits(-jetEtaRangekt,jetEtaRangekt);
      }
   }


   // #### Task configuration
   task->SetUsePileUpCut(usePileUpCut);
   task->SetUseDefaultVertexCut(useVertexCut);
   task->SetAcceptanceWindows(trackEtaWindow, jetRadius);
   task->SelectCollisionCandidates(trigger);
   task->SetAnalysisType(collisionSystem, typeOfAnal); 

   task->SetTrackContainerName(trackarrayname);
   task->SetMCParticleContainerName(mcpariclearraynamePartMC);

   task->SetJetContainerName(jetarrayname);
   task->SetMCPartJetContainerName(jetarraynamePartMC);
   task->SetKTJetContainerName(ktjetarrayname);
   task->SetKTMCPartJetContainerName(ktjetarraynamePartMC);

   task->SetSignalJetMinArea(acut);

   task->SetDebugLevel(0); //No debug messages 0

   // output container
   AliAnalysisDataContainer *contHistos = manager->CreateContainer(myContName.Data(), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChJetSpectra%s", AliAnalysisManager::GetCommonFileName(), myContName.Data()));


   // #### ADD ANALYSIS TASK
   manager->AddTask(task);
   manager->ConnectInput(task, 0, manager->GetCommonInputContainer());
   manager->ConnectOutput(task, 1, contHistos);

   return task;


}
//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetConePt(Double_t eta, Double_t phi, Double_t radius, AliParticleContainer *recTrkCont){
   //sum up pt inside a cone
   Double_t tmpConePt = 0.0;

   if(!recTrkCont) return 0.0;

   AliVParticle* tmpTrack=NULL; 
   for(auto trackIterator : recTrkCont->accepted_momentum()){
      // trackIterator is a std::map of AliTLorentzVector and AliVTrack
      tmpTrack = trackIterator.second;  // Get the full track
      if(!tmpTrack) continue; 
      if(IsTrackInAcceptance(tmpTrack, kFALSE)){  //check if reconstructed track is in acceptance 
         if(GetDeltaR(tmpTrack->Phi(), phi, tmpTrack->Eta(), eta) < radius){
            tmpConePt = tmpConePt + tmpTrack->Pt();
         }
      }
   }
   return tmpConePt;
}

//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetSimPrimaryVertex(){
   //get generator level primary vertex

   AliGenEventHeader* mcHeader = NULL; 
   AliAODMCHeader* aodMCH = NULL;

   if(MCEvent()){
      if(fTypeOfAnal == kEff){ 
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
   }


   if(mcHeader){
      TArrayF pyVtx;   //primary vertex x,y,z
      mcHeader->PrimaryVertex(pyVtx);
      return (Double_t) (pyVtx[2]);
   }
   AliWarning(Form("In task %s: Pythia Vertex failed!", GetName()));
   return 9999.0;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHJetSpectra::IsMCEventInAcceptance(AliVEvent* event){
   //EVENT SELECTION PURE MC

   if(!event) return kFALSE;
   if(!MCEvent()) return kFALSE; 

   Double_t vtxMC = GetSimPrimaryVertex();
   fhVertexZMC->Fill(vtxMC); //Fill BEFORE vertex cut

   if(TMath::Abs(vtxMC) > fZVertexCut){
      return kFALSE;
   }
   fhVertexZAcceptMC->Fill(vtxMC);//Fill AFTER vertex cut

   return kTRUE;
   
}
//________________________________________________________________________
Bool_t AliAnalysisTaskHJetSpectra::IsEventInAcceptance(AliVEvent* event){
   //EVENT SELECTION RECONSTRUCTED DATA

   if(!event) return kFALSE;

   //___________________________________________________
   //TEST PILE UP
   if(fUsePileUpCut){
      if(!fHelperClass || fHelperClass->IsPileUpEvent(event)){ 
         fHistEvtSelection->Fill(1.5); //count events rejected by pileup
         return kFALSE;
      }
   }

   //___________________________________________________
   //BEFORE VERTEX CUT
   fhVertexZ->Fill(event->GetPrimaryVertex()->GetZ()); 

   if(fUseDefaultVertexCut){
      if(!fHelperClass || !fHelperClass->IsVertexSelected2013pA(event)){
         fHistEvtSelection->Fill(2.5); //count events rejected by vertex cut 
         return kFALSE;
      }
   }else{
      if(TMath::Abs(event->GetPrimaryVertex()->GetZ()) > fZVertexCut){
         fHistEvtSelection->Fill(2.5); //count events rejected by vertex cut 
         return kFALSE;
      }
      if(event->GetPrimaryVertex()->GetNContributors()<1){
         fHistEvtSelection->Fill(2.5); //count events rejected by vertex cut 
         return kFALSE;
      }
   }
   //___________________________________________________
   //AFTER VERTEX CUT
   fhVertexXAccept->Fill(event->GetPrimaryVertex()->GetX());
   fhVertexYAccept->Fill(event->GetPrimaryVertex()->GetY());
   fhVertexZAccept->Fill(event->GetPrimaryVertex()->GetZ());

   //___________________________________________________
   if(fTypeOfAnal == kRec){   
      fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
      if(!fMultSelection ||
         !fMultSelection->GetThisEventIsNotPileup() ||
         !fMultSelection->GetThisEventIsNotPileupInMultBins() ||
         !fMultSelection->GetThisEventHasNoInconsistentVertices() ||
         !fMultSelection->GetThisEventPassesTrackletVsCluster()){

         fHistEvtSelection->Fill(3.5); //count events rejected by multiplicity selection
         return kFALSE;
      }
   }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHJetSpectra::IsTrackInAcceptance(AliVParticle* track, Bool_t isGen){
   // Check if the track pt and eta range 
   if(!track) return kFALSE;

   if(isGen){ //pure MC select charged primary tracks 
      //Apply only for kine level or MC containers   
      if(!track->Charge()) return kFALSE;
      if(fTypeOfAnal == kEff){
         if(!(static_cast<AliAODMCParticle*>(track))->IsPhysicalPrimary()) return kFALSE;
      }    
   }
   if(TMath::Abs(track->Eta()) <= fTrackEtaWindow){ //APPLY TRACK ETA CUT
      if(track->Pt() >= fMinTrackPt){   //APPLY TRACK CUT
         return kTRUE;
      }
   }
   return kFALSE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskHJetSpectra::IsSignalJetInAcceptance(AliEmcalJet *jet, Bool_t suppressGhost){   
   //select jets in acceptance 
   if(!jet) return kFALSE;
   if(TMath::Abs(jet->Eta()) <= fSignalJetEtaWindow){
      if(jet->Area() >= fMinJetArea){
         if(suppressGhost){
            if(jet->Pt() >= fMinTrackPt) return kTRUE;
         }else{
            return kTRUE;
         }
      }
   }  
   return kFALSE;
}


//________________________________________________________________________
void AliAnalysisTaskHJetSpectra::ExecOnceLocal(){
   // Initialization of jet containers done in  AliAnalysisTaskEmcalJet::ExecOnce()
   //Read arrays of jets and tracks
   fInitializedLocal = kTRUE; 

   // Initialize helper class (for vertex selection & pile up correction)
   fHelperClass = new AliAnalysisUtils();
   fHelperClass->SetCutOnZVertexSPD(kFALSE); // kFALSE: no cut; kTRUE: |zvtx-SPD - zvtx-TPC|<0.5cm

   return;
}


//________________________________________________________________________
Double_t  AliAnalysisTaskHJetSpectra::GetDeltaPt(Double_t rho, Double_t ttPhi, Double_t ttEta, AliParticleContainer *recTrkCont){ 

   //delta pt = pT in random cone - rho * Area of random cone
   // processes real reconstructed data. Exclude region around jet with TT   


   // Define random cone Eta+Phi
   Double_t tmpRandConeEta = -fSignalJetEtaWindow  + fRandom->Rndm()*2*fSignalJetEtaWindow;
   Double_t tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();

   while(GetDeltaR( tmpRandConePhi, ttPhi, tmpRandConeEta, ttEta)<2*fSignalJetRadius){
      tmpRandConeEta = -fSignalJetEtaWindow  + fRandom->Rndm()*2*fSignalJetEtaWindow;
      tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();
   }
 
   Double_t conePt = GetConePt(tmpRandConeEta, tmpRandConePhi, fSignalJetRadius, recTrkCont);
    
   return  conePt - (rho*fSignalJetRadiusSquared*TMath::Pi());
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHJetSpectra::FillHistograms(){  
   // executed in each event 
   //called in AliAnalysisTaskEmcal::UserExec(Option_t *)
   //   Analyze the event and Fill histograms

   if(!InputEvent()){
      AliError("??? Event pointer == 0 ???");
      return kFALSE;
   }

   
   if(!fInitializedLocal) ExecOnceLocal(); //Executed only once 

   //_________________________________________________________________
   //  FILL EVENT STATISTICS
   fHistEvtSelection->Fill(0.5); //Count input event

   if(fTypeOfAnal == kEff){   //Check MC event vertex
      if(!IsMCEventInAcceptance(InputEvent())) return kFALSE; //post data is in UserExec
   }

   //Check Reconstructed event vertex
   if(!IsEventInAcceptance(InputEvent())) return kFALSE; //post data is in UserExec

   if(fTypeOfAnal == kRec){   //Select minimum bias trigger in real data
      UInt_t triggerMask = fInputHandler->IsEventSelected();
      if(!(triggerMask & AliVEvent::kINT7)){
         fHistEvtSelection->Fill(4.5);
	 return kFALSE;
      }
   }

   fHistEvtSelection->Fill(5.5); //Count Accepted input event
   // END EVENT SELECTION

   //___________________
   //VZERO AND ZNA CENTRALITY
   Double_t multVzero   = -1.0;
   Double_t energyZdcNA = -1.0;
   if(InputEvent()->GetVZEROData()){
      multVzero = InputEvent()->GetVZEROData()->GetMTotV0A();
   }

   //ZDC ENERGY
   AliVZDC *aodZDC = dynamic_cast<AliVZDC*> (InputEvent()->GetZDCData());
   if(aodZDC){
      const Double_t *ZNAtower = aodZDC->GetZNATowerEnergy();
      energyZdcNA = ZNAtower[0]*4.*82./208./12.96; //ZNA energy in TeV 
   }
   // Get centrality
   Double_t centralityPercentileV0A = 1.0;  //initialize default centrality for MC and pp
   Double_t centralityPercentileZNA = 1.0;

   if(fTypeOfAnal == kRec  && fCollisionSystem != kpp){
      fMultSelection = (AliMultSelection*) InputEvent()->FindListObject("MultSelection");
      if(fMultSelection){
         centralityPercentileV0A = fMultSelection->GetMultiplicityPercentile("V0A");
         centralityPercentileZNA = fMultSelection->GetMultiplicityPercentile("ZNA");
         
         fhCentralityV0A->Fill(centralityPercentileV0A);
         fhCentralityZNA->Fill(centralityPercentileZNA);
         fhVzeroATotMult->Fill(centralityPercentileV0A, multVzero);
         fhZNAEnergy->Fill(centralityPercentileZNA, energyZdcNA);
      }
   }

   //_________________________________________________________________
   // JET+TRACK CONTAINERS
   AliJetContainer *jetContDet    = NULL; //AKTjet container detector level 
   AliJetContainer *jetContPart   = NULL; //AKT jet container particle level 
   AliJetContainer *jetContDetKT  = NULL; //KT jet container detector level 
   AliJetContainer *jetContPartKT = NULL; //KT jet containe particle level 

   AliEmcalJet  *jetDet  = NULL;  //jet pointer detector level jet
   AliEmcalJet  *jetPart = NULL;  //jet pointer particle level jet

   AliParticleContainer *trkContDet = NULL; //track container at detector level 
   AliParticleContainer *parContPart = NULL; //particle container at particle level

   AliVParticle *pointerTrackDet = NULL; // detector level track
   AliVParticle *pointerParticle = NULL; // particle level particle
   //_________________________________________________________
   //READ JET TRACK CONTAINERS
   jetContDet   = GetJetContainer(fMyJetContainerName.Data()); //AKT detector level jet container
   trkContDet   = GetTrackContainer(fMyTrackContainerName.Data()); //container with detector level tracks
   jetContDetKT = GetJetContainer(fMyKTJetContainerName.Data());//KT detector level jet container

   if(fTypeOfAnal == kEff){ 
      //MC + EMB DATA
      jetContPart   = GetJetContainer(fMyJetParticleContainerName.Data());      //AKT particle level jet container 
      parContPart   = GetParticleContainer(fMyParticleContainerName.Data()); //particle level paticle container
      jetContPartKT = GetJetContainer(fMyKTJetParticleContainerName.Data());   //KT particle level jet container
   }

   Double_t criteria = fRandom->Rndm(); //dch: get random [0;1], if c < 0.1, fill reference, otherwise signal histo for tts
   Bool_t bFillSig = (criteria>0.1) ? kTRUE : kFALSE;    //dch if kTRUE => fill signal,  if kFALSE => fill ref

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //   INCLUSIVE JETS RECONSTRUCTED DATA 
   if(jetContDet){
      for(auto jetIteratorDet : jetContDet->accepted_momentum()){
         jetDet = jetIteratorDet.second;  // Get the full track
         if(!jetDet) continue;
         if(!IsSignalJetInAcceptance(jetDet,kTRUE)) continue; //cuts on eta, pT ,area
                  
         fhJetPhiDet->Fill(jetDet->Pt(), jetDet->Phi());
         fhJetEtaDet->Fill(jetDet->Pt(), jetDet->Eta());
	 fhJetAreaVsPt->Fill(jetDet->Pt(), jetDet->Area());
      }
   }


   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //LOOP OVER TRACKS  SEARCH FOR TRIGGER CANDIDATES AT PARTICLE LEVEL 
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   TLorentzVector myTT;  //trigger track

   for(Int_t it = 0; it < kTT; it++){
      fTTH_Det[it].resize(0);  //vector of trigger track candidates at detector level
      fTTH_Part[it].resize(0); //vector of trigger track candidates at particle level
   }


   Int_t indexSingleRndTrigPart[kTT] = {-1,-1}; //index of single random trigger


   if(fTypeOfAnal == kEff){

      if(parContPart){ 
      
         for(auto mcPartIterator : parContPart->accepted_momentum()){
            pointerParticle = mcPartIterator.second;  // Get the pointer to mc particle object
            if(!pointerParticle) continue;
       
            if(IsTrackInAcceptance(pointerParticle, kTRUE)){ //select charged physics primary particles 
  
               for(Int_t it=0; it<kTT; it++){ 
                  if((fTTlow[it] <= pointerParticle->Pt()) && (pointerParticle->Pt() < fTThigh[it])){
		     myTT.SetPtEtaPhiM(pointerParticle->Pt(), pointerParticle->Eta(), pointerParticle->Phi(), 0.);
                     fTTH_Part[it].push_back(myTT);  //fill vector of TT candiates at particle level
                  }
               }
            }
         }//end of loop over MC particles

         for(Int_t it=0; it<kTT; it++){ 
            //MULTIPLICITY OF TRIGGER TRACK CANDIDATES AT PARTICLE LEVEL
            fhTTMultPart[it]->Fill(fTTH_Part[it].size()); 
         
            //SELECT TRIGGER PARTICLE AT PARTICLE LEVEL
            if(fTTH_Part[it].size()>0){
               indexSingleRndTrigPart[it] = fRandom->Integer(fTTH_Part[it].size()); //Integer 0 ... ntriggers-1
            }
         }
      }
   }
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //LOOP OVER TRACKS  SEARCH FOR TRIGGER CANDIDATES IN DETECTOR LEVEL TRACKS
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   Int_t indexSingleRndTrigDet[kTT] = {-1,-1}; //index of single random trigger
 
   if(trkContDet){

      for(auto trackIterator : trkContDet->accepted_momentum()){
         // trackIterator is a std::map of AliTLorentzVector and AliVTrack
         pointerTrackDet = trackIterator.second;  // Get the full track

         if(IsTrackInAcceptance(pointerTrackDet, kFALSE)){  //check if reconstructed track is in acceptance 
            fhTrackPhi->Fill(pointerTrackDet->Pt(), pointerTrackDet->Phi()); 
            fhTrackEta->Fill(pointerTrackDet->Pt(), pointerTrackDet->Eta());

            for(Int_t it=0; it<kTT; it++){ //select trigger candidates
               if(fTTlow[it] <= pointerTrackDet->Pt() && pointerTrackDet->Pt() < fTThigh[it]){
		  myTT.SetPtEtaPhiM(pointerTrackDet->Pt(), pointerTrackDet->Eta(), pointerTrackDet->Phi(), 0.);
                  fTTH_Det[it].push_back(myTT);
               }
            }
         }
      }
   }
   
   for(Int_t it=0; it<kTT; it++){
      fhTTMultDet_V0A[it]->Fill(centralityPercentileV0A, fTTH_Det[it].size()); 
      fhTTMultDet_ZNA[it]->Fill(centralityPercentileZNA, fTTH_Det[it].size()); 
   
      
      //  SELECT SINGLE INCLUSIVE REC LEVEL TRIGGER 
      if(fTTH_Det[it].size()>0){ 
         indexSingleRndTrigDet[it] =  fRandom->Integer(fTTH_Det[it].size()); //index of the selected candidate
      }
   }

   //___________________________________________________________
   //_________________________________________________________
   //EVALUATE SINGLE PARTICLE EFFICIENCY + FILL RESPONSE MATRIX  for PYTHIA PP EVENTS
   Bool_t bRecPrim = kFALSE; //tags the reconstructed detector level physical primary particles

   if(fTypeOfAnal == kEff){

      //1) FILL HISTOS FOR SINGLE PARTICLE EFFICIENCY
      if(parContPart){
         for(auto mcPartIterator : parContPart->accepted_momentum()){
            pointerParticle = mcPartIterator.second;  // Get the pointer to mc particle object
            if(!pointerParticle) continue;
 
            if(IsTrackInAcceptance(pointerParticle, kTRUE)){
               //pT spectrum of particle level physical primary particles
               fhPhysPrimaryPtPart->Fill(pointerParticle->Pt());
            }  
         }

         //single particle efficiency and contamination
         if(trkContDet && parContPart){ 
            for(auto trackIterator : trkContDet->accepted_momentum()){
               pointerTrackDet = trackIterator.second;  // Get the full track
               if(!pointerTrackDet) continue;
               if(!IsTrackInAcceptance(pointerTrackDet, kFALSE)) continue; //reconstructed level tracks
               bRecPrim = kFALSE; //not yet matched to generator level physical primary

               for(auto mcPartIterator : parContPart->accepted_momentum()){
                  pointerParticle = mcPartIterator.second;  // Get the pointer to mc particle object
                  if(!pointerParticle) continue;
                  if(!IsTrackInAcceptance(pointerParticle, kTRUE)) continue; //gen level physical primary
                  if(TMath::Abs(pointerTrackDet->GetLabel()) == TMath::Abs(pointerParticle->GetLabel())){ 
                     //has the same label as reconstr track
 
                     bRecPrim = kTRUE;
                     fhPhysPrimaryPtDet->Fill(pointerParticle->Pt()); //this is well recontr phys primary
                   
                     break;
                  }//same label with rec particle
               }//loop over gen tracks
               if(!bRecPrim){
                  fhPtTrkSecOrFakeDet->Fill(pointerTrackDet->Pt()); //matchnig to phys primary not found, this is fake or second.
               }
            }//loop over detector level tracks
         }//detector level track array exists
      }// particle level array exists
   
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //2) FILL JET RESPONSE MATRIX
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Double_t ptRecCorr; //REC jet pt corrected for rho

      //Response matrix normalization - spectrum of all generator level jets in acceptance
      if(jetContPart){
         for(auto jetIteratorPart : jetContPart->accepted_momentum()){
            jetPart = jetIteratorPart.second;  // Get the full track
            if(!jetPart) continue;
            if(!IsSignalJetInAcceptance(jetPart,kTRUE)) continue; //cuts on eta, pT ,area

            fhJetPtPart->Fill(jetPart->Pt());
            fhJetEtaPart->Fill(jetPart->Pt(), jetPart->Eta());
         }
      }
 
      //Find closest gen level+rec level  jets
      if(jetContDet){
         for(auto jetIteratorDet : jetContDet->accepted_momentum()){
            jetDet = jetIteratorDet.second;  // Get the full track
            if(!jetDet) continue;
            if(!IsSignalJetInAcceptance(jetDet,kTRUE)) continue; //cuts on eta, pT ,area

            //Find closest gen level+rec level  jets
            jetPart = jetDet->ClosestJet();

            if(!jetPart){ //did not find matching generator level jet
               continue;
            }
            if(jetPart->Pt()<1e-3){
                //if(fDebug>20)  Printf("SKIP MATCH WITH GHOST JET");
                continue; 
            }
            //corresponding generator level jet found
            if(!IsSignalJetInAcceptance(jetPart,kTRUE)) continue; //cuts on eta, pT ,area

            fhJetPtPartVsJetPtDet->Fill(jetDet->Pt(), jetPart->Pt()); //response matrix

         } 
      }//rec jet container exists
   }//analyze efficiency mode (response matrix + single particle efficiency) 

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Double_t areaJet,  pTJet; 
   Double_t dphi, deta;

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //H-JET CORRELATIONS IN MC TRUTH
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(fTypeOfAnal == kEff){
      for(Int_t it=0; it< kTT; it++){
	 if(bFillSig  && it != kSig ) continue; //dch bFillSig is kTRUE => we want to fill Signal TT   but the it counter is at reference TT => skip 
	 if(!bFillSig && it == kSig ) continue; //dch bFillSig is kFALSE => we want to fill Ref TT  but the it counter is at signal TT => skip 


         if(parContPart && fTTH_Part[it].size() >0){
	    myTT = (fTTH_Part[it])[indexSingleRndTrigPart[it]]; //TT particle level
         
            fhTTPart[it]->Fill((Float_t) myTT.Pt()); //trigger pT gen 
         
            //JET LOOP
	    if(jetContPart){
               for(auto jetIteratorPart : jetContPart->accepted_momentum()){
                  jetPart = jetIteratorPart.second;  // Get the full track
                  if(!jetPart) continue;
                  if(!IsSignalJetInAcceptance(jetPart,kTRUE)) continue; //cuts on eta, pT ,area

                  if(TMath::Abs(TVector2::Phi_mpi_pi(jetPart->Phi() - myTT.Phi())) < fDphiCut) continue;  //Dphi cut between trigger and assoc

		  fHJetSpecPart[it]->Fill(jetPart->Pt());
               }//JET LOOP
            }//container exists
         }//TT exists
      }//TT loop
   }
   
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // RECONSTRUCTED DATA ANALYSIS

   // CALCULATE RHO

   for(Int_t it=0;it<kTT;it++){
      if(fTypeOfAnal == kEff || fCollisionSystem == kpp){ 
         fRhoDet[it] = 0.0; 
      }else{
         if(fTTH_Det[it].size()>0){
            fRhoDet[it] = EstimateBgKT( jetContDetKT, trkContDet, (fTTH_Det[it])[indexSingleRndTrigDet[it]]);
         }
      }
   } 

   if(fCollisionSystem != kpp){
      for(Int_t it=0; it<kTT; it++){
         if(fTTH_Det[it].size()>0){
         
            fhRhoTT[it]->Fill(centralityPercentileV0A, (Float_t) fRhoDet[it]);
      
            fhCentralityTT_V0A[it]->Fill((Float_t) centralityPercentileV0A);
            fhCentralityTT_ZNA[it]->Fill((Float_t) centralityPercentileZNA);
         }
      }
   }


   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // CALCULATE DELTA PT IN RECONSTRUCTED DATA WITH RANDOM CONE METHOD
   // Exclude region around TT from delta pt calculation
   Double_t deltapt, phiTT = 0., etaTT = -1000.;       
   Bool_t  bRecJetCloseToTT = kFALSE;
   if(fCollisionSystem != kpp){
      for(Int_t it=0; it<kTT;it++){
         if(fTTH_Det[it].size()>0){ //get phi and eta of the TT or the reconstructed jet that contains TT
            myTT = (fTTH_Det[it])[indexSingleRndTrigDet[it]];  
            phiTT = myTT.Phi(); // TT is proxy for jet
            etaTT = myTT.Eta();
     
            if(jetContDet){
               for(auto jetIteratorDet : jetContDet->accepted_momentum()){
                  jetDet = jetIteratorDet.second;  // Get the full track
                  if(!jetDet) continue;
                  if(!IsSignalJetInAcceptance(jetDet,kTRUE)) continue; //cuts on eta, pT ,area
      
                  //find jet which containes TT as a constituent
                  bRecJetCloseToTT = kFALSE;
                 
                  if(jetDet->Pt() > myTT.Pt()*0.5){ //jet containing TT has pT larger than pT of TT
                     for(Int_t iq=0; iq < jetDet->GetNumberOfTracks(); iq++) {
                        pointerTrackDet = (AliVParticle*) (jetDet->TrackAt(iq,trkContDet->GetArray())); //matched rec and emb tracks
                        if(!pointerTrackDet) continue;
      
           	     if(TMath::Abs(pointerTrackDet->Eta() - myTT.Eta()) < 1e-3){
                           if(TMath::Abs(TVector2::Phi_mpi_pi(pointerTrackDet->Phi() - myTT.Phi())) < 1e-3){
                              phiTT = jetDet->Phi(); // jet which contains TT
                              etaTT = jetDet->Eta(); // jet which contains TT
      
                              bRecJetCloseToTT = kTRUE;
                              break;
                           }
                        }
           	  }
                  }
                  if(bRecJetCloseToTT) break; 
               }
      
              
      	    deltapt = GetDeltaPt((Double_t) (fRhoDet[it]),  phiTT, etaTT, trkContDet);
               
               fhDeltaPt_V0A[it]->Fill(centralityPercentileV0A, deltapt); 
               fhDeltaPt_ZNA[it]->Fill(centralityPercentileZNA, deltapt);
            }       
         }//number TT gt 0
      }//trigger loop 
   }
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //  H+JET IN RECONSTRUCTED DATA  
 
   for(Int_t it=0; it<kTT;it++){
      if(bFillSig  && it != kSig ) continue; //dch bFillSig is kTRUE => we want to fill Signal TT   but the it counter is at reference TT => skip 
      if(!bFillSig && it == kSig ) continue; //dch bFillSig is kFALSE => we want to fill Ref TT  but the it counter is at signal TT => skip 

      if(fTTH_Det[it].size()>0){
         myTT = (fTTH_Det[it])[indexSingleRndTrigDet[it]]; 
      
         fhTTDet_V0A[it]->Fill(centralityPercentileV0A, (Float_t) myTT.Pt()); //trigger p 
         fhTTDet_ZNA[it]->Fill(centralityPercentileZNA, (Float_t) myTT.Pt()); //trigger p 
 
      	 //   INCLUSIVE JETS RECONSTRUCTED DATA 
         if(jetContDet){
            for(auto jetIteratorDet : jetContDet->accepted_momentum()){
               jetDet = jetIteratorDet.second;  // Get the full track
               if(!jetDet) continue;
               if(!IsSignalJetInAcceptance(jetDet,kTRUE)) continue; //cuts on eta, pT ,area
            
                   if(TMath::Abs(TVector2::Phi_mpi_pi(jetDet->Phi() - myTT.Phi())) < fDphiCut) continue;  //Dphi cut between trigger and assoc
                   if(fTypeOfAnal == kEff ||  fCollisionSystem == kpp){
                      pTJet = jetDet->Pt();
		   }else{
                      pTJet = jetDet->Pt() - jetDet->Area()*fRhoDet[it];
		   }
      
    
                   fHJetSpecDet_V0A[it]->Fill(centralityPercentileV0A, pTJet);
                   fHJetSpecDet_ZNA[it]->Fill(centralityPercentileZNA, pTJet);
      
            }//jet loop
         }//jet container exists
      } //TT exists
   } //TT loop

   return kTRUE;
}

//________________________________________________________________________
/*Double_t AliAnalysisTaskHJetSpectra::GetNcoll(Double_t centr){
   //Get Ncoll for given centrality
   if(fCollisionSystem == kpPb){
      //convert centrality percentle to Ncoll (2014-Sep-26-analysis_note-AnalysisNoteDijetspPb.pdf table 1)
      //  PLATI JEN PRO pA !!!!!!!!!!   CO PP?
      if(centr < 0.0 || centr > 100.) return -1;

      if(centr < 5.0)     return 14.7;  //0-5%
      else if(centr < 10) return 13.0; //5-10%
      else if(centr < 20) return 11.7; //10-20%
      else if(centr < 40) return 9.38; //20-40%
      else if(centr < 60) return 6.49; //40-60%
      else if(centr < 80) return 3.96; //60-80%
      else                return 1.52; //80-100%
   }
     

   return -1.; //pp, pbpb

}*/
//________________________________________________________________________
void AliAnalysisTaskHJetSpectra::Terminate(Option_t *){
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
AliAnalysisTaskHJetSpectra::~AliAnalysisTaskHJetSpectra(){
   // Destructor. Clean-up the output list, but not the histograms that are put inside
   // (the list is owner and will clean-up these histograms). Protect in PROOF case.
   if(fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
   }
   delete fRandom;
   delete fHelperClass;
 
} 

//________________________________________________________________________
void AliAnalysisTaskHJetSpectra::UserCreateOutputObjects(){
  // called once to create user defined output objects like histograms, plots etc. 
  // and to put it on the output list.
  // Note: Saving to file with e.g. OpenFile(0) is must be before creating other objects.
  //fOutput TList defined in the mother class
   AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

   fRandom = new TRandom3(0);

   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   TString name;

   //__________________________________________________________
   // Event statistics
   fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 6, 0, 6);
   fHistEvtSelection->GetXaxis()->SetBinLabel(1,"All events");
   fHistEvtSelection->GetXaxis()->SetBinLabel(2,"Rejected by pile up");
   fHistEvtSelection->GetXaxis()->SetBinLabel(3,"Rejected by vertex");
   fHistEvtSelection->GetXaxis()->SetBinLabel(4,"Without mult object");
   fHistEvtSelection->GetXaxis()->SetBinLabel(5,"Rejected by trigger");
   fHistEvtSelection->GetXaxis()->SetBinLabel(6,"Accepted events");

   fOutput->Add(fHistEvtSelection);
   //___________________________________________________________
   // Hard trigger counter
   for(Int_t it =0; it<kTT; it++){
      name = "fhTTDet_V0A"; 
      name.Append(Form("_TT%d%d", TMath::Nint(fTTlow[it]), TMath::Nint(fTThigh[it])));      
      fhTTDet_V0A[it] = new TH2D(name.Data(),"# of triggers",11,0,110, 50,0.0,50.0);
      fOutput->Add((TH2D*) fhTTDet_V0A[it]);
 
      name = "fhTTDet_ZNA"; 
      name.Append(Form("_TT%d%d", TMath::Nint(fTTlow[it]), TMath::Nint(fTThigh[it])));      
      fhTTDet_ZNA[it] = new TH2D(name.Data(),"# of triggers",11,0,110, 50,0.0,50.0);
      fOutput->Add((TH2D*) fhTTDet_ZNA[it]);
      
      name = "fhTTMultDet_V0A";
      name.Append(Form("_TT%d%d", TMath::Nint(fTTlow[it]), TMath::Nint(fTThigh[it])));      
      fhTTMultDet_V0A[it] = new TH2F(name.Data(),"multiplicity of triggers in event V0A",11,0,110,10,0.0,10.0);
      fOutput->Add((TH2F*) fhTTMultDet_V0A[it]);

      name = "fhTTMultDet_ZNA";
      name.Append(Form("_TT%d%d", TMath::Nint(fTTlow[it]), TMath::Nint(fTThigh[it])));      
      fhTTMultDet_ZNA[it] = new TH2F(name.Data(),"multiplicity of triggers in event ZNA",11,0,110,10,0.0,10.0);
      fOutput->Add((TH2F*) fhTTMultDet_ZNA[it]);
   }

   if(fTypeOfAnal == kEff){
      for(Int_t it=0; it<kTT;it++){ 
         name = "fhTTPart";
         name.Append(Form("TT%d%d",TMath::Nint(fTTlow[it]), TMath::Nint(fTThigh[it]))); 
         fhTTPart[it] = new TH1F(name.Data(), name.Data(),100, 0, 100);
         fOutput->Add((TH1D*) fhTTPart[it]);
        
        
         name =Form("fhTTMultPartTT%d%d",TMath::Nint(fTTlow[it]), TMath::Nint(fTThigh[it])); 
         fhTTMultPart[it] = new TH1F(name.Data(), name.Data(),10, 0, 10);
         fOutput->Add((TH1D*) fhTTMultPart[it]);
      }
   }
   //___________________________________________________________
   // trigger associated jet spectra (jet pT not corrected for UE)

   for(Int_t it =0; it<kTT; it++){
      name = Form("fHJetSpecTT%d%d_V0A", TMath::Nint(fTTlow[it]), TMath::Nint(fTThigh[it]));      
      fHJetSpecDet_V0A[it] = new TH2D(name.Data(),name.Data(),110,0,110,170,-20,150);
      fOutput->Add((TH2D*) fHJetSpecDet_V0A[it]);

      name = Form("fHJetSpecTT%d%d_ZNA", TMath::Nint(fTTlow[it]), TMath::Nint(fTThigh[it]));      
      fHJetSpecDet_ZNA[it] = new TH2D(name.Data(),name.Data(),110,0,110,170,-20,150);
      fOutput->Add((TH2D*) fHJetSpecDet_ZNA[it]);
   } 

   if(fTypeOfAnal == kEff){
      for(Int_t it=0; it<kTT;it++){ 
     	 name = Form("fHJetSpecPartTT%d%d", TMath::Nint(fTTlow[it]), TMath::Nint(fTThigh[it]));
         fHJetSpecPart[it] = new TH1D(name.Data(), name.Data(), 250, 0, 250);
         fOutput->Add((TH1D*) fHJetSpecPart[it]);
      }
   } 
   //____________________________________________________________________
   //UE from cell median  
   if(fCollisionSystem != kpp){
      for(Int_t it=0; it<kTT; it++){
         //rho in events with TT 
         name = Form("fhRhoTT%d%d",TMath::Nint(fTTlow[it]),TMath::Nint(fTThigh[it])); 
         fhRhoTT[it] = new TH2F(name.Data(),name.Data(),110,0,100, 80, 0.0, 40.0);
         fOutput->Add((TH2F*) fhRhoTT[it]);
      }
   }
   //_______________________________________________________________________
   // Delta pt distributions  
   if(fCollisionSystem != kpp){ 
      for(Int_t it =0; it< kTT; it++){
         name = Form("fhDeltaPtTT%d%d_V0A",TMath::Nint(fTTlow[it]),TMath::Nint(fTThigh[it])); 
         fhDeltaPt_V0A[it] = new TH2D(name.Data(),name.Data(), 110, 0, 110, 150, -50, 100);
         fOutput->Add((TH2D*) fhDeltaPt_V0A[it]);
      
         name = Form("fhDeltaPtTT%d%d_ZNA",TMath::Nint(fTTlow[it]),TMath::Nint(fTThigh[it])); 
         fhDeltaPt_ZNA[it] = new TH2D(name.Data(),name.Data(), 110, 0, 110, 150, -50, 100);
         fOutput->Add((TH2D*) fhDeltaPt_ZNA[it]);
      }
   }
   //_______________________________________________________________________
   //inclusive azimuthal and pseudorapidity histograms

   name = "fhTrackPhi"; 
   fhTrackPhi = new TH2F(name.Data(),"azim dist trig had vs pT,trk", 50, 0, 50, 50, 0,TMath::TwoPi());
   fOutput->Add((TH2F*) fhTrackPhi);

   name = "fhTrackEta";
   fhTrackEta = new TH2F(name.Data(),"Eta dist trig had vs pT,trk", 100, 0, 100, 40,-0.9,0.9);
   fOutput->Add((TH2F*) fhTrackEta);

   fhJetPhiDet = new TH2F("fhJetPhiDet","Azim dist jets vs pTjet", 50, 0, 100, 50,0,TMath::TwoPi());
   fOutput->Add((TH2F*) fhJetPhiDet);

   fhJetEtaDet = new TH2F("fhJetEtaDet","Eta dist inclusive jets vs pTjet", 100,0, 100, 40,-0.9,0.9);
   fOutput->Add((TH2F*) fhJetEtaDet);

   if(fTypeOfAnal== kEff ){
      name = "fhJetEtaPart";
      fhJetEtaPart = (TH2F*) fhJetEtaDet->Clone(name.Data());
      fOutput->Add((TH2F*) fhJetEtaPart);

      name = "fhJetPtPart";
      fhJetPtPart = new TH1D(name.Data(), name.Data(), 250, 0, 250);
      fOutput->Add((TH1D*) fhJetPtPart);
   }

   fhJetAreaVsPt = new TH2D("fhJetAreaVsPt","fhJetAreaVsPt",100, 0, 100, 100, 0, 1); //dch 
   fOutput->Add(fhJetAreaVsPt);

   //-------------------------
   fhVertexZ = new TH1F("fhVertexZ","z vertex",40,-20,20);
   fOutput->Add(fhVertexZ);

   fhVertexXAccept = new TH1F("fhVertexXAccept","vertex after cut",600,-3,3);
   fOutput->Add(fhVertexXAccept);
 
   fhVertexYAccept = (TH1F*) fhVertexXAccept->Clone("fhVertexYAccept");
   fOutput->Add(fhVertexYAccept);
 
   fhVertexZAccept = new TH1F("fhVertexZAccept","z vertex after cut",40,-20,20);
   fOutput->Add(fhVertexZAccept);

   if(fTypeOfAnal== kEff ){
      fhVertexZMC = new TH1F("fhVertexZMC","z vertex",40,-20,20);
      fOutput->Add(fhVertexZMC);
 
      fhVertexZAcceptMC = new TH1F("fhVertexZAcceptMC","z vertex after cut",40,-20,20);
      fOutput->Add(fhVertexZAcceptMC);
   }
   //-------------------------
   if(fCollisionSystem != kpp){

      fhCentralityV0A = new TH1F("hCentralityV0A","hCentralityV0A",110,0,110);
      fOutput->Add(fhCentralityV0A); 
      
      fhCentralityZNA = new TH1F("hCentralityZNA","hCentralityZNA",110,0,110);
      fOutput->Add(fhCentralityZNA);
      
      for(Int_t it =0; it<kTT;it++){
         name = Form("fhCentralityTT%d%d_V0A", TMath::Nint(fTTlow[it]), TMath::Nint(fTThigh[it])); 
         fhCentralityTT_V0A[it] = (TH1F*) fhCentralityV0A->Clone(name.Data());
         fOutput->Add(fhCentralityTT_V0A[it]);
      }
      
      for(Int_t it =0; it<kTT;it++){
         name = Form("fhCentralityTT%d%d_ZNA", TMath::Nint(fTTlow[it]), TMath::Nint(fTThigh[it])); 
         fhCentralityTT_ZNA[it] = (TH1F*) fhCentralityZNA->Clone(name.Data());
         fOutput->Add(fhCentralityTT_ZNA[it]);
      }
   }

   if(fTypeOfAnal == kRec  && fCollisionSystem != kpp){
      //-----------------------------------------------------
      //  vzero multiplicity
      name = "fhVzeroATotMult_V0A"; 
      fhVzeroATotMult = new TH2F(name.Data(),"hVzeroATotMult",11,0,110,200,0,1000);
      fOutput->Add((TH2F*) fhVzeroATotMult);
      
      //-----------------------------------------------------
      //  ZDC ZNA energy
      name = "fhZNAEnergy_ZNA";
      fhZNAEnergy = new TH2F(name.Data(),"fhZNAEnergy",11,0,110,200,0,3000);
      fOutput->Add((TH2F*)fhZNAEnergy);
   }
   

   //-----------------------------------------------------
   if(fTypeOfAnal== kEff ){
      name = "fhJetPtPartVsJetPtDet";
      fhJetPtPartVsJetPtDet = new TH2D(name.Data(), name.Data(), 250, 0, 250, 250, 0, 250); 
      fOutput->Add((TH2D*) fhJetPtPartVsJetPtDet);
      
      name = "fhPhysPrimaryPtDet"; 
      fhPhysPrimaryPtDet = new TH1D(name.Data(),"",100,0,100);
      fOutput->Add((TH1D*) fhPhysPrimaryPtDet);
 

      name = "fhPhysPrimaryPtPart";
      fhPhysPrimaryPtPart = new TH1D(name.Data(),"",100,0,100);
      fOutput->Add((TH1D*) fhPhysPrimaryPtPart);

      name = "fhPtTrkSecOrFakeDetMB";
      fhPtTrkSecOrFakeDet = new TH1D(name.Data(),"",100,0,100);
      fOutput->Add((TH1D*) fhPtTrkSecOrFakeDet);

   }


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
Bool_t AliAnalysisTaskHJetSpectra::RetrieveEventObjects() {
   //
   // retrieve event objects
   //
    if(!AliAnalysisTaskEmcalJet::RetrieveEventObjects())  return kFALSE;
 
   return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskHJetSpectra::Run()
{
   // Run analysis code here, if needed. It will be executed before FillHistograms().
   
   return kTRUE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::EstimateBgKT(AliJetContainer *jetCont, AliParticleContainer *trkCont, TLorentzVector myTT){
   //Estimate rho from KT jet median. Ignore jet that contains TT
   Double_t rhoKT = 0.0;
 
   if(!jetCont) return rhoKT;   

   AliEmcalJet*  jet        = NULL;
   AliVParticle* constTrack = NULL;
   Bool_t bKTJetCloseToTT = kFALSE;
   Int_t nJetAcc = 0;
   Double_t jetpt;
   Double_t sumEmbPt;

   jetCont->ResetCurrentID();
   while((jet = jetCont->GetNextAcceptJet())){ //loop over KT jets
      if(!jet) continue;
      if(!IsSignalJetInAcceptance(jet,kFALSE)) continue;

      bKTJetCloseToTT = kFALSE;

      if(jet->Pt() > myTT.Pt()*0.5){ //jet containing TT has pT larger than pT of TT
         for(Int_t iq=0; iq < jet->GetNumberOfTracks(); iq++) {
            constTrack = (AliVParticle*) (jet->TrackAt(iq,trkCont->GetArray())); //matched rec and emb tracks
            if(!constTrack) continue;
            if(TMath::Abs(constTrack->Eta() - myTT.Eta()) < 1e-3){
               if(TMath::Abs(TVector2::Phi_mpi_pi(constTrack->Phi() - myTT.Phi())) < 1e-3){
                  bKTJetCloseToTT = kTRUE; 
                  break;
               } 
            }
         }
      }
      if(bKTJetCloseToTT) continue; //skip the jet that contains TT 

      jetpt = jet->Pt();
      if(jetpt <0.005) jetpt = 0.; //set pt of ghost jets identical to zero
      if(jet->Area() < 0.2) continue; //skip jets which have too small area
      frhovec[nJetAcc] = jetpt/jet->Area();
      nJetAcc++;
   }

   if(nJetAcc>0){
      rhoKT = TMath::Median(nJetAcc, frhovec);
   }
 
  return rhoKT; 
}
//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetDeltaR(Double_t phi1, Double_t phi2, Double_t eta1, Double_t eta2){
   //angular distance between two jets
   Double_t dphi = TVector2::Phi_mpi_pi(phi1-phi2);
   Double_t deta = eta1 - eta2;
   return sqrt(dphi*dphi + deta*deta); 
}
