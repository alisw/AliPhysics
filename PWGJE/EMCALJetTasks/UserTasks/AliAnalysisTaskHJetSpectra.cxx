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
#include "AliVVertex.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "AliJetContainer.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisTaskHJetSpectra.h"
#include "AliHeader.h" //KINE
#include "AliRunLoader.h"   //KINE
using namespace std;

// ANALYSIS OF HIGH PT HADRON TRIGGER ASSOCIATED SPECTRUM OF RECOIL JETS IN P+PB
// Author Filip Krizek   (17.May. 2014)

//TODO: Not accessing the particles when using MC
//TODO: FillHistogram can be done better with virtual TH1(?)
ClassImp(AliAnalysisTaskHJetSpectra)
//________________________________________________________________________________________

AliAnalysisTaskHJetSpectra::AliAnalysisTaskHJetSpectra(): 
AliAnalysisTaskEmcalJet("AliAnalysisTaskHJetSpectra", kTRUE),  
 fCollisionSystem(0), fTypeOfData(0), fTypeOfAnal(0),
  fUseDefaultVertexCut(1), fUsePileUpCut(1),  
  fRhoTaskName(), fRhoTaskNameMC(),
fPerpConeRadius(0.4),fPerpConeRadiusSquared(fPerpConeRadius*fPerpConeRadius),
 fSignalJetRadius(0.4), fSignalJetRadiusSquared(fSignalJetRadius*fSignalJetRadius),
fSignalJetEtaWindow(0.9 - fSignalJetRadius), fTrackEtaWindow(0.9), fMinTrackPt(0.150), fMinJetArea(0.0),  
fCentralityType("V0A"), fCentPercMin(0.), fCentPercMax(100.), fMinFractionShared(0.5), 
fCrossSection(0.0), fTrials(0.0), fImpParam(-1.0), fRandom(0), fHelperClass(0), fInitializedLocal(0),
fTTlow(8.0), fTThigh(9.0), fTTtype(0), fDphiCut(TMath::Pi()-0.6), fUseDoubleBinPrecision(0),
fHistEvtSelection(0x0), fh1Ntriggers(0x0), fh1NtriggersGen(0x0), 
fhJetPhi(0x0),  fhJetPhiGen(0x0), fhTrackPhi(0x0), fhJetEta(0x0), fhJetEtaGen(0x0), fhTrackEta(0x0), fhTrackPt(0x0), fhTrackPtGen(0x0), fhVertexZ(0x0), fhVertexZAccept(0x0), fhVertexZMC(0x0), fhVertexZAcceptMC(0x0),
 fhDphiTriggerJetAccept(0x0),
fhCentrality(0x0), fhCentralityV0M(0x0), fhCentralityV0A(0x0), fhCentralityV0C(0x0), fhCentralityZNA(0x0),
fh1Xsec(0x0), fh1Trials(0x0), fh1PtHard(0x0), fhImpactParameter(0x0), fhImpactParameterTT(0x0),
fhPtTrkTruePrimRec(0x0), fhPtTrkTruePrimGen(0x0), fhPtTrkSecOrFakeRec(0x0),
fRhoRec(kRho),fRhoMC(kRho),
fNofRandomCones(1),
fZVertexCut(10.0),
fTrigTracksGen(),
fTrigTracks()
{
    //default constructor
   for(Int_t ir=0; ir<kRho; ir++){
      fHJetSpec[ir]=NULL;
      fhDphiTriggerJet[ir]=NULL;
      fhDphiTriggerJetGen[ir]=NULL;
      fhJetPtGen[ir]=NULL;
      fhJetPtGenVsJetPtRec[ir]=NULL;
      fhJetPtResolutionVsPtGen[ir]=NULL;

   }

   for(Int_t ir=0; ir<kRho-1; ir++){
      fhRhoTT[ir]=NULL;
      fhRhoIncl[ir]=NULL;
      fARhoTT[ir]=NULL;
      fhDeltaPt[ir]=NULL; 
      fhDeltaPtEmb[ir]=NULL; 
      fhDeltaPtEmb2D[ir]=NULL;
      fhDeltaPtIncl[ir]=NULL;
   }

}

//________________________________________________________________________
AliAnalysisTaskHJetSpectra::AliAnalysisTaskHJetSpectra(const char *name) : 
AliAnalysisTaskEmcalJet(name,kTRUE),  
fCollisionSystem(0), fTypeOfData(0), fTypeOfAnal(0),
  fUseDefaultVertexCut(1), fUsePileUpCut(1),
   fRhoTaskName(), fRhoTaskNameMC(),
fPerpConeRadius(0.4),  fPerpConeRadiusSquared(fPerpConeRadius*fPerpConeRadius), 
fSignalJetRadius(0.4), fSignalJetRadiusSquared(fSignalJetRadius*fSignalJetRadius), 
fSignalJetEtaWindow(0.9 - fSignalJetRadius), fTrackEtaWindow(0.9), fMinTrackPt(0.150), fMinJetArea(0.0),   
fCentralityType("V0A"), fCentPercMin(0.), fCentPercMax(100.),  fMinFractionShared(0.5), 
fCrossSection(0.0), fTrials(0.0), fImpParam(-1.0), fRandom(0), fHelperClass(0), fInitializedLocal(0), 
fTTlow(8.0), fTThigh(9.0), fTTtype(0), fDphiCut(TMath::Pi()-0.6), fUseDoubleBinPrecision(0),
fHistEvtSelection(0x0), fh1Ntriggers(0x0), fh1NtriggersGen(0x0), 
fhJetPhi(0x0), fhJetPhiGen(0x0), fhTrackPhi(0x0), fhJetEta(0x0), fhJetEtaGen(0x0), fhTrackEta(0x0), fhTrackPt(0x0), fhTrackPtGen(0x0), fhVertexZ(0x0), fhVertexZAccept(0x0), fhVertexZMC(0x0), fhVertexZAcceptMC(0x0),
fhDphiTriggerJetAccept(0x0),
fhCentrality(0x0), fhCentralityV0M(0x0), fhCentralityV0A(0x0), fhCentralityV0C(0x0), fhCentralityZNA(0x0),
fh1Xsec(0x0), fh1Trials(0x0), fh1PtHard(0x0), fhImpactParameter(0x0), fhImpactParameterTT(0x0),
 fhPtTrkTruePrimRec(0x0), fhPtTrkTruePrimGen(0x0), fhPtTrkSecOrFakeRec(0x0),
fRhoRec(kRho),fRhoMC(kRho),
fNofRandomCones(1),
fZVertexCut(10.0),
fTrigTracksGen(),
fTrigTracks()
{
//Constructor
   for(Int_t ir=0; ir<kRho; ir++){
      fHJetSpec[ir]=NULL;
      fhDphiTriggerJet[ir]=NULL;
      fhDphiTriggerJetGen[ir]=NULL;
      fhJetPtGen[ir]=NULL;
      fhJetPtGenVsJetPtRec[ir]=NULL;
      fhJetPtResolutionVsPtGen[ir]=NULL;

   }

   for(Int_t ir=0; ir<kRho-1; ir++){
      fhRhoTT[ir]=NULL;
      fhRhoIncl[ir]=NULL;
      fARhoTT[ir]=NULL;
      fhDeltaPt[ir]=NULL;
      fhDeltaPtEmb[ir]=NULL; 
      fhDeltaPtEmb2D[ir]=NULL;
      fhDeltaPtIncl[ir]=NULL;
   }

   DefineOutput(1, TList::Class());
}
  
//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetConePt(Double_t eta, Double_t phi, Double_t radius, Int_t icont){
   //sum up pt inside a cone
   Double_t tmpConePt = 0.0;

   TClonesArray *trkArray = GetParticleArray(icont); 
   if(!trkArray) return 0.0;

   for(Int_t i = 0; i < trkArray->GetEntries(); i++){
      AliVParticle* tmpTrack = static_cast<AliVParticle*>(trkArray->At(i));
      if(!tmpTrack) continue; 
      if(IsTrackInAcceptance(tmpTrack, icont)){ 
         if(GetDeltaR(tmpTrack->Phi(), phi, tmpTrack->Eta(), eta) < radius){
            tmpConePt = tmpConePt + tmpTrack->Pt();
         }
      }
   }
   return tmpConePt;
}

//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetPtHard(){

   //Get pt hard from pythia header
   //Fill Xsection, trials histograms

   AliGenPythiaEventHeader* pythiaHeader = NULL; 

   if(fTypeOfAnal == kKine){ //KINE 
      AliRunLoader *rl = AliRunLoader::Instance();
      if(rl)  pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(rl->GetHeader()->GenEventHeader());
      if(pythiaHeader){
         fh1Xsec->Fill("<#sigma>", pythiaHeader->GetXsection());
         fh1Trials->Fill("#sum{ntrials}", pythiaHeader->Trials());

         return pythiaHeader->GetPtHard();
      }

   } else {
      
      if(MCEvent()){ 
         pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader()); 
         if(!pythiaHeader){
            // Check if AOD
            AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));
      
            if(aodMCH){
               for(UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++){
                  pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
                  if(pythiaHeader) break;
               }
            }
         }
      }
      
      if(pythiaHeader){
      
         TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
      
         if(tree){
            TFile *curfile = tree->GetCurrentFile();
            if(!curfile) {
               Error("Notify","No current file");
               return 0.0;
            }
            Float_t xsection = 0.0;
            Float_t trials  = 1.0;
      
            AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,trials);
      
            fCrossSection = (Double_t) xsection;//pythiaHeader->GetXsection();
      

            if(fCrossSection > 0.){ //save cross-section and the number of trials
               fTrials = (Double_t) trials; //pythiaHeader->Trials();
               fh1Xsec->Fill("<#sigma>", fCrossSection);
               fh1Trials->Fill("#sum{ntrials}", fTrials);
            }
         }

         //if(fDebug>20)  Printf("XSECTION 1=%f   2=%f   TRIALS 1=%f  2=%f", 
         //                       fCrossSection, pythiaHeader->GetXsection(), fTrials,  pythiaHeader->Trials());

         return pythiaHeader->GetPtHard();
      }
   }
   AliWarning(Form("In task %s: GetPtHard() failed!", GetName()));
   return -1.0;
}
//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetImpactParameter(){
   //get impact parameter from hijing
   AliGenHijingEventHeader* hijingHeader = dynamic_cast<AliGenHijingEventHeader*>(MCEvent()->GenEventHeader());
   if(MCEvent()){ 
      if(!hijingHeader){
         // Check if AOD
         AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

         if(aodMCH){
            for(UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++){
               hijingHeader = dynamic_cast<AliGenHijingEventHeader*>(aodMCH->GetCocktailHeader(i));
               if(hijingHeader) break;
            }
         }
      }
   }
   if(hijingHeader){
      return (Double_t) (hijingHeader->ImpactParameter());
   }
   AliWarning(Form("In task %s: GetImpactParameter() failed!", GetName()));
   return -1.0;
}
//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetSimPrimaryVertex(){
   //get generator level primary vertex

   AliGenEventHeader* mcHeader = NULL; 
   AliAODMCHeader* aodMCH = NULL;

   if(fTypeOfAnal == kKine){ //KINE
      AliRunLoader *rl = AliRunLoader::Instance();
      if(rl)  mcHeader = dynamic_cast<AliGenPythiaEventHeader*>(rl->GetHeader()->GenEventHeader());

   } else {

      if(MCEvent()){
         if(fTypeOfData == kPythia){ 
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
         
         if(fTypeOfData == kHijing){ 
            mcHeader = dynamic_cast<AliGenHijingEventHeader*>(MCEvent()->GenEventHeader());
            if(!mcHeader){
               // Check if AOD
                aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));
         
                if(aodMCH){
                   for(UInt_t i = 0; i<aodMCH->GetNCocktailHeaders(); i++){
                     mcHeader = dynamic_cast<AliGenHijingEventHeader*>(aodMCH->GetCocktailHeader(i));
                     if(mcHeader) break;
                  }
               }
            }
         }
      }
   }


   if(mcHeader){
      
      TArrayF pyVtx(3);
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

   //___________________________________________________

   if(fTypeOfAnal!= kKine  &&   !MCEvent()) return kFALSE; 

   
   Double_t vtxMC = GetSimPrimaryVertex();
   fhVertexZMC->Fill(vtxMC); //Fill BEFORE vertex cut

   if(TMath::Abs(vtxMC) > fZVertexCut){
      fHistEvtSelection->Fill(3); //count events rejected by vertex cut 
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
         fHistEvtSelection->Fill(2); //count events rejected by pileup
         return kFALSE;
      }
   }
   //___________________________________________________
   //BEFORE VERTEX CUT
   fhVertexZ->Fill(event->GetPrimaryVertex()->GetZ()); 

   if(fUseDefaultVertexCut){
      if(!fHelperClass || !fHelperClass->IsVertexSelected2013pA(event)){
         fHistEvtSelection->Fill(3); //count events rejected by vertex cut 
         return kFALSE;
      }
   }else{
      if(TMath::Abs(event->GetPrimaryVertex()->GetZ()) > fZVertexCut){
         fHistEvtSelection->Fill(3); //count events rejected by vertex cut 
         return kFALSE;
      }
   }
   //___________________________________________________
   //AFTER VERTEX CUT
   fhVertexZAccept->Fill(event->GetPrimaryVertex()->GetZ());

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHJetSpectra::IsTrackInAcceptance(AliVParticle* track, Int_t icont){
   // Check if the track pt and eta range 
   if(!track) return kFALSE;

   if(fTypeOfAnal == kKine || icont == kContainerTwo){ //pure MC select charged primary tracks 
      //Apply only for kine level or MC containers   
      if((!track->Charge()) || (!(static_cast<AliAODMCParticle*>(track))->IsPhysicalPrimary()))
         return kFALSE;
   }
   if(TMath::Abs(track->Eta()) <= fTrackEtaWindow){ //APPLY TRACK ETA CUT
      if(track->Pt() >= fMinTrackPt){   //APPLY TRACK CUT
         return kTRUE;
      }
   }
   return kFALSE;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskHJetSpectra::IsSignalJetInAcceptance(AliEmcalJet *jet){   
   //select jets in acceptance 
   if(!jet) return kFALSE;
   if(TMath::Abs(jet->Eta()) <= fSignalJetEtaWindow){
      if(jet->Pt() >= fMinTrackPt){
         if(jet->Area() >= fMinJetArea){
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
void  AliAnalysisTaskHJetSpectra::GetDeltaPt(Int_t nrho, TArrayD &rho, Double_t *dpt,
                                           Double_t ttPhi, Double_t ttEta,
                                           Double_t leadingJetExclusionProbability){

   //delta pt = random cone - rho
   // processes real reconstructed data

   for(Int_t ir=0;ir<nrho;ir++){
      dpt[ir] = -10000.0;   // Set an invalid delta pt
   }

   // Define random cone Eta+Phi
   Bool_t coneValid = kTRUE;
   Double_t tmpRandConeEta = -fSignalJetEtaWindow  + fRandom->Rndm()*2*fSignalJetEtaWindow;
   Double_t tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();

   if(ttEta > -2.0){  //TT exists => make sure that one generates RC far away from the trigger track
      while(GetDeltaR( tmpRandConePhi, ttPhi, tmpRandConeEta, ttEta)<fSignalJetRadius){
         tmpRandConeEta = -fSignalJetEtaWindow  + fRandom->Rndm()*2*fSignalJetEtaWindow;
         tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();
      }
   }

   //if(fDebug>20)  Printf("RND CONE RCphi=%f   RCeta=%f   TTphi=%f   TTeta=%f", 
    //                     tmpRandConePhi,tmpRandConeEta,ttPhi,ttEta);

   //Based on Ncoll probability that   Reject a jet based on the probability  check for overlap if demanded
   /*

   //WHAT TO DO HERE ???  SHOULD RANDOM CONES BE CORRELAtED WITH JETs and rejected based on 1./Ncoll ?

   if(leadingJetExclusionProbability>0){ //skips pp
      
      AliEmcalJet* tmpLeading = NULL; 
      Double_t lpt = -1.0; 
      // Get leading jet (regardless of pT)

      AliEmcalJet* tmpJet = NULL;
      AliJetContainer *jetContRec = GetJetContainer(kContainerOne);
      if(jetContRec){
         jetContRec->ResetCurrentID();
         while((tmpJet = jetContRec->GetNextAcceptJet())) {
            if(!tmpJet) continue;
            if((TMath::Abs(tmpJet->Eta()) <= fSignalJetEtaWindow) && (tmpJet->Area() >= fMinJetArea)){
               if(tmpJet->Pt() > lpt){
                  tmpLeading = tmpJet;
                  lpt =  tmpJet->Pt();
               }
            }
         }
      
         if(tmpLeading){
            Double_t tmpDeltaPhi    = RelativePhi(tmpRandConePhi, tmpLeading->Phi());
            Double_t tmpDeltaEta = tmpLeading->Eta()-tmpRandConeEta;
         
            // Check, if cone has overlap with jet
            if(sqrt(tmpDeltaPhi*tmpDeltaPhi + tmpDeltaEta*tmpDeltaEta) <= fSignalJetRadius){
               // probability to exclude the RC
               // Only exclude cone with a given probability
               if(fRandom->Rndm()<=leadingJetExclusionProbability)  coneValid = kFALSE;
            }
         }
      }
      if(fRandom->Rndm()<=leadingJetExclusionProbability)  coneValid = kFALSE;
   }*/
 
   // Get the cones' pt and calculate delta pt
   if(coneValid){
      Double_t conePt = GetConePt(tmpRandConeEta,tmpRandConePhi, fSignalJetRadius, kContainerOne);
      for(Int_t ir=0; ir < nrho; ir++){
         dpt[ir] =  conePt - (rho[ir]*fSignalJetRadiusSquared*TMath::Pi());
      }
   }
 
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

   //Execute only once:  Get tracks, jets from arrays if not already given 
   if(!fInitializedLocal) ExecOnceLocal(); 

  //_________________________________________________________________
   // fill MC information

   if(fTypeOfData == kPythia){   //FK//How does this work with Embedding ?
      fh1PtHard->Fill(GetPtHard());  //Fills cross section
   }

   if(fTypeOfData == kHijing){
      fImpParam = GetImpactParameter(); 
      fhImpactParameter->Fill(fImpParam);
   }
   //_________________________________________________________________
   //  FILL EVENT STATISTICS
   fHistEvtSelection->Fill(1); //Count input event

   if(fTypeOfData != kReal){   //Check MC event vertex
      if(!IsMCEventInAcceptance(InputEvent())) return kFALSE; //post data is in UserExec
   }

   //Check Reconstructed event vertex
   if(fTypeOfAnal != kKine){   //Check MC event vertex
      if(!IsEventInAcceptance(InputEvent())) return kFALSE; //post data is in UserExec
   }

   //___________________
   // Get centrality
   Double_t centralityPercentile    = -1.0; //KINE

   if(fCollisionSystem != kpp){   //KINE Check MC event vertex
      AliCentrality* tmpCentrality = InputEvent()->GetCentrality();
      if(!tmpCentrality){
         fHistEvtSelection->Fill(4);
         return kFALSE; //post data is in UserExec
      }
      Double_t centralityPercentileV0A = -1.0;
      Double_t centralityPercentileV0C = -1.0;
      Double_t centralityPercentileV0M = -1.0;
      Double_t centralityPercentileZNA = -1.0;
      if(tmpCentrality != NULL){
         centralityPercentile    = tmpCentrality->GetCentralityPercentile(fCentralityType.Data());
         centralityPercentileV0A = tmpCentrality->GetCentralityPercentile("V0A");
         centralityPercentileV0C = tmpCentrality->GetCentralityPercentile("V0C");
         centralityPercentileV0M = tmpCentrality->GetCentralityPercentile("V0M");
         centralityPercentileZNA = tmpCentrality->GetCentralityPercentile("ZNA");
      }
      
      fhCentrality->Fill(centralityPercentile);
      
      if((centralityPercentile < fCentPercMin) || (centralityPercentile > fCentPercMax)){ //cut on centrality
         AliWarning(Form("Centrality value not valid (c=%E)",centralityPercentile)); 
         fHistEvtSelection->Fill(4);
         return kFALSE;
      }
      fhCentralityV0M->Fill(centralityPercentileV0M); 
      fhCentralityV0A->Fill(centralityPercentileV0A);
      fhCentralityV0C->Fill(centralityPercentileV0C); 
      fhCentralityZNA->Fill(centralityPercentileZNA);
   } 
   fHistEvtSelection->Fill(0); //Count Accepted input event

   // END EVENT SELECTION
   //_________________________________________________________________
   // JET+TRACK CONTAINERS
   AliJetContainer *jetContRec = NULL; //jet container from reconstruced tracks
   AliJetContainer *jetContGen = NULL; // jet container from  MC particles
   AliEmcalJet* jetGen = NULL;
   AliEmcalJet* jetRec = NULL;
   TClonesArray *trkArrayRec = NULL; //track array of real reconstructed tracks 
   TClonesArray *trkArrayGen = NULL; //track array of MC particles

   //_________________________________________________________
   //READ JET TRACK CONTAINERS
   if(fTypeOfAnal == kRec){ 
      jetContRec  = GetJetContainer(kContainerOne);
      trkArrayRec = GetParticleArray(kContainerOne); //reconstructed particle container
   }

   if(fTypeOfAnal == kEff || fTypeOfAnal == kEmb || fTypeOfAnal == kEmbSingl){
      jetContRec  = GetJetContainer(kContainerOne);
      trkArrayRec = GetParticleArray(kContainerOne); //reconstructed particle container
      jetContGen  = GetJetContainer(kContainerTwo); //GENERATOR LEVEL JETS
      trkArrayGen = GetParticleArray(kContainerTwo); //true MC particle container
   }

   if(fTypeOfAnal == kKine){   // Kine written to the 0th container !!!!!!!!!!!!!!!! 
      jetContGen  = GetJetContainer(kContainerOne);
      trkArrayGen = GetParticleArray(kContainerOne); //process kine particle container as reco???
   }
 
   //if(fDebug>20)  Printf("POINTER TO CONTAINERS   JETrec=%p  TRKrec=%p   JETgen=%p   TRKgen=%p", 
   //                       jetContRec,trkArrayRec,jetContGen, trkArrayGen);


   //___________________________________________________________
   // CALCULATE RHO
   fRhoRec.Reset(0.);  
   fRhoMC.Reset(0.);

   if(fTypeOfAnal == kRec || fTypeOfAnal == kEff || fTypeOfAnal == kEmb || fTypeOfAnal == kEmbSingl){
      fRhoRec[kConeRho] = EstimateBgCone(kContainerOne); //container ID=0 reconstructed tracks
      fRhoRec[kCMSRho]  = GetRhoVal(kContainerOne);   //funkce AliAnalysisTaskEmcalJet.h 
      fRhoRec[kZeroRho] = 0.0; 
   }

   if(fTypeOfAnal == kEff || fTypeOfAnal == kEmb){ //rho in MC events 
      fRhoMC[kConeRho] = EstimateBgCone(kContainerTwo); //container ID=1 mc particles
      fRhoMC[kCMSRho]  = GetRhoVal(kContainerTwo); //funkce AliAnalysisTaskEmcalJet.h
      fRhoMC[kZeroRho] = 0.0;
   }

   if(fTypeOfAnal == kEmbSingl){ //embedding single track
      fRhoMC[kConeRho] = 0.0; 
      fRhoMC[kCMSRho]  = 0.0; 
      fRhoMC[kZeroRho] = 0.0;
   }

   if(fTypeOfAnal == kKine){ //rho in KINE MC events 
      fRhoMC[kConeRho] = EstimateBgCone(kContainerOne); //container ID=1 mc particles
      fRhoMC[kCMSRho]  = GetRhoVal(kContainerOne); //funkce AliAnalysisTaskEmcalJet.h
      fRhoMC[kZeroRho] = 0.0;
   }


         
   //_________________________________________________________
   //Evaluate Single particle Efficiency + Fill Response Matrix
   if(fTypeOfAnal == kEff){

      //1) FILL HISTOS FOR SINGLE PARTICLE EFFICIENCY
      if(trkArrayGen){

         for(Int_t i = 0; i < trkArrayGen->GetEntries(); i++){
            AliVParticle* tmpTrackGen = static_cast<AliVParticle*>(trkArrayGen->At(i));
            if(!tmpTrackGen) continue;
            if(IsTrackInAcceptance(tmpTrackGen, kContainerTwo)){
               //pT spectrum of generator level particles
               fhPtTrkTruePrimGen->Fill(tmpTrackGen->Pt(),tmpTrackGen->Eta()); 
            }  
         }

         //single particle efficiency and contamination
         Bool_t bRecPrim = kFALSE; //tags the reconstructed primary particles
         if(trkArrayRec && trkArrayGen){ 

            for(Int_t j = 0; j < trkArrayRec->GetEntries(); j++){ // loop over reconstructed tracks 
               AliVParticle* tmpTrackRec = static_cast<AliVParticle*>(trkArrayRec->At(j));
               if(!tmpTrackRec) continue;
               if(!IsTrackInAcceptance(tmpTrackRec, kContainerOne)) continue; //reconstructed level tracks
               bRecPrim = kFALSE; //not yet matched to generator level physical primary

               for(Int_t i = 0; i < trkArrayGen->GetEntries(); i++){
                  AliVParticle* tmpTrackGen = static_cast<AliVParticle*>(trkArrayGen->At(i));
                  if(!tmpTrackGen) continue;
                  if(!IsTrackInAcceptance(tmpTrackGen, kContainerTwo)) continue; //gen level physical primary
                  if(TMath::Abs(tmpTrackRec->GetLabel()) == TMath::Abs(tmpTrackGen->GetLabel())){ 
                     //has the same label as reconstr track
 
                     bRecPrim = kTRUE;
                     fhPtTrkTruePrimRec->Fill(tmpTrackGen->Pt(),tmpTrackGen->Eta()); //this is well recontr phys primary
                     break;
                  }//same label with rec particle
               }//loop over gen tracks
               if(!bRecPrim) fhPtTrkSecOrFakeRec->Fill(tmpTrackRec->Pt(),tmpTrackRec->Eta()); //matchnig to phys primary not found, this is fake or second.
            }//loop over rec tracks
         }//rec track array exists
      }//gen particle array exists
   
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      //2) FILL JET RESPONSE MATRIX
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++
      Double_t ptGenCorr; //GEN jet pt corrected for rho
      Double_t ptRecCorr; //REC jet pt corrected for rho

      //Response matrix normalization - spectrum of all generator level jets in acceptance
      if(jetContGen){
         jetContGen->ResetCurrentID();
         while((jetGen = jetContGen->GetNextAcceptJet())){ 
            if(!jetGen) continue;
            if(!IsSignalJetInAcceptance(jetGen)) continue; //cuts on eta, pT ,area

            //if(fDebug>20) Printf("GEN JET phi=%f  eta=%f  pt=%f", jetGen->Phi(), jetGen->Eta(), jetGen->Pt());

            for(Int_t ir=0; ir< kRho; ir++){
               ptGenCorr = jetGen->Pt() - jetGen->Area()*fRhoMC[ir]; // correct for rho
               fhJetPtGen[ir]->Fill(ptGenCorr);
            }
         }
      }
 
      //Find closest gen level+rec level  jets ala AliAnalysisHelperJetTasks::GetClosestJets 
      if(jetContRec){
         jetContRec->ResetCurrentID();
         while((jetRec = jetContRec->GetNextAcceptJet())) {
            if(!jetRec) continue;
            if(!IsSignalJetInAcceptance(jetRec)) continue; //cuts on eta, pT ,area

            //if(fDebug>20) Printf("REC JET phi=%f  eta=%f  pt=%f",jetRec->Phi(), jetRec->Eta(), jetRec->Pt());

            jetGen = 0x0;
            jetGen = jetRec->ClosestJet();

            if(!jetGen){ //did not find matching generator level jet
               //if(fDebug>20)  Printf("NO MATCH (NO SUCH GEN JET)");

               continue;
            }
            if(jetGen->Pt()<1e-3){
                //if(fDebug>20)  Printf("SKIP MATCH WITH GHOST JET");
                continue; 
            }
            //corresponding generator level jet found
            //if(fDebug>20)  Printf("MATCHED WITH  phi=%f  eta=%f  pt=%f",jetGen->Phi(), jetGen->Eta(), jetGen->Pt());
            
            if(!IsSignalJetInAcceptance(jetGen)) continue; //cuts on eta, pT ,area

            //check fraction of tracks from generator level jet in rec level jet
            //if(fDebug>20)  Printf("FRACTIONH SHARED = %f ", GetFractionSharedPt(jetRec,jetContRec,jetGen, jetContGen));

            if(GetFractionSharedPt(jetRec,jetContRec,jetGen, jetContGen) < fMinFractionShared) continue;

            //if(fDebug>20)  Printf("PASSED MIN FRACTION CRITERION ");

            for(Int_t ir=0; ir< kRho; ir++){
               ptGenCorr = jetGen->Pt() - jetGen->Area()*fRhoMC[ir]; //perp cone bg correction to pt 
               ptRecCorr = jetRec->Pt() - jetRec->Area()*fRhoRec[ir]; //perp cone bg correction to pt 

               fhJetPtGenVsJetPtRec[ir]->Fill(ptRecCorr, ptGenCorr); //response matrix

               if(ptGenCorr >0){
                  fhJetPtResolutionVsPtGen[ir]->Fill(ptGenCorr,(ptRecCorr-ptGenCorr)/ptGenCorr); //jet pT resolution
               }
            }
         } 
      }//rec jet comtainer exists
   }//analyze efficiency mode (response matrix + single particle efficiency) 

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Double_t areaJet,  pTJet; 
   Double_t tmpArray[3];              
   Double_t dphi, dfi;
   Bool_t bFirstCycle = kTRUE; 
   Int_t  ic;  //container

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //H-JET CORRELATIONS IN MC TRUTH
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(fTypeOfAnal == kEff || fTypeOfAnal == kKine){

      //std::vector<Int_t> fTrigTracksGen is a list of trigger particle indices in true MC
      fTrigTracksGen.clear(); //list of trigger particle indices true MC

      if(trkArrayGen){ 
         for(Int_t i = 0; i < trkArrayGen->GetEntries(); i++){ // loop over reconstructed tracks 
            AliVParticle* track = static_cast<AliVParticle*>(trkArrayGen->At(i));

            if(!track) continue;
            ic = (fTypeOfAnal == kEff ) ? kContainerTwo : kContainerOne; 
       
            if(IsTrackInAcceptance(track, ic)){  
               fhTrackPtGen->Fill(track->Pt());  //inclusive pT spectrum of tracks
                
               if((fTTlow <= track->Pt()) && (track->Pt() < fTThigh)){
                  fTrigTracksGen.push_back(i);  //trigger candidates

                  //if(fDebug>20)  Printf("GEN TT candidate  index = %d  phi=%f  eta=%f pT=%f ", 
                  //            i, track->Phi(), track->Eta(), track->Pt());

               }
            }
         }
      }
 
      Int_t ntriggersGen = (Int_t) fTrigTracksGen.size();
      Int_t indexSingleRndTrigGen = -1; //index of single random trigger
   
      if(ntriggersGen>0){
         if(fTTtype==0){ //select single inclusive trigger
            indexSingleRndTrigGen = fRandom->Integer(ntriggersGen); //Integer 0 ... ntriggers-1
         }
      }

      if(trkArrayGen){
         bFirstCycle = kTRUE; 
         for(Int_t it=0; it<ntriggersGen; it++){ //loop over trigger configurations
        
            if(fTTtype==0){
               if(it != indexSingleRndTrigGen) continue;
            }
         
            AliVParticle* triggerHadronGen = static_cast<AliVParticle*>(trkArrayGen->At(fTrigTracksGen[it]));
            if(!triggerHadronGen) continue;
          
            fh1NtriggersGen->Fill((Float_t) triggerHadronGen->Pt()); //trigger pT gen 
         
            if( fTypeOfData == kHijing){ //impact parameter for triggered events
               fhImpactParameterTT->Fill(fImpParam);
            }
         
         
            //JET LOOP
            if(jetContGen){
               jetContGen->ResetCurrentID();
               while((jetGen = jetContGen->GetNextAcceptJet())) {
                  if(!jetGen){
                     AliError(Form("%s: Could not receive gen jet", GetName()));
                     continue;
                  }
                  if(!IsSignalJetInAcceptance(jetGen)) continue;
                 
                  areaJet = jetGen->Area();
                  pTJet   = jetGen->Pt();
                 
                  if(bFirstCycle){
                     fhJetPhiGen->Fill( pTJet, RelativePhi(jetGen->Phi(),0.0));
                     fhJetEtaGen->Fill( pTJet, jetGen->Eta());
                  }
                 
                  dphi = RelativePhi(triggerHadronGen->Phi(), jetGen->Phi());
                 
                  dfi = dphi; //-0.5*pi to 1.5*Pi
                  if(dfi < -0.5*TMath::Pi()) dfi += TMath::TwoPi();
                  if(dfi >  1.5*TMath::Pi()) dfi -= TMath::TwoPi();
                 
                  for(Int_t ir=0; ir< kRho;ir++){ 
                     fhDphiTriggerJetGen[ir]->Fill((Float_t) (pTJet - areaJet*fRhoMC[ir]), (Float_t) dfi); //Rongrong's analysis 
                  }
                  //-------------------------
                 
                  if(TMath::Abs(dphi) < fDphiCut) continue;  //Dphi cut between trigger and assoc
                 
                 //Centrality, A, pTjet
                  tmpArray[0] =  areaJet;
                  for(Int_t ir=0; ir< kRho;ir++){ 
                     tmpArray[1] =  pTJet - areaJet*fRhoMC[ir];
                     fHJetSpecGen[ir]->Fill(tmpArray);
                  }
               }//JET LOOP
               bFirstCycle = kFALSE;
            }//container exists
         }
      }
   }
   
   //++++++++++++++++++++++++++++++++++++++++++++++++++++

   if(fTypeOfAnal == kKine) return kTRUE; 

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //LOOP OVER TRACKS  SEARCH FOR TRIGGER CANDIDATES IN REC TRACKS
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   //std::vector<Int_t> fTrigTracks; //list pf trigger particle indices
   fTrigTracks.clear(); //list pf trigger particle indices

   if(trkArrayRec){ 
      for(Int_t i = 0; i < trkArrayRec->GetEntries(); i++){ // loop over reconstructed tracks 
         AliVParticle* track = static_cast<AliVParticle*>(trkArrayRec->At(i));

         if(!track) continue;

         if(fTypeOfAnal == kEmb || fTypeOfAnal == kEmbSingl){
            //if(fDebug>99)  Printf("TRACK LABEL %d", track->GetLabel());
            //in embed evts search for TT only among real tracks; do not consider embedded tracks as trigger
            if(TMath::Abs(track->GetLabel()) == 99999) continue;//9999 set in AddTaskHJetSpectra.C
         }


         if(IsTrackInAcceptance(track, kContainerOne)){  //rec   (Same analysis for gen level?)
            //Fill some inclusive spectra
            fhTrackPhi->Fill(track->Pt(), RelativePhi(track->Phi(),0.0)); // phi = -pi,pi
            fhTrackEta->Fill(track->Pt(), track->Eta());
            fhTrackPt->Fill(track->Pt());

            if(fTTlow <= track->Pt() && track->Pt() < fTThigh){
               fTrigTracks.push_back(i);  //trigger candidates

               //if(fDebug>20)  Printf("REC TT candidate  index = %d  phi=%f  eta=%f pT=%f ", 
               //               i, track->Phi(), track->Eta(), track->Pt());
            }
         }
      }
   }
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //   SELECT SINGLE INCLUSIVE TRIGGER 
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Int_t ntriggers = (Int_t) fTrigTracks.size();
   Int_t indexSingleRndTrig = -1; //index of single random trigger

   if(ntriggers>0){
      if(fTTtype == 0){ //select single inclusive trigger
         indexSingleRndTrig = fRandom->Integer(ntriggers); //Integer 0 ... ntriggers-1
      }

      //if(fDebug>20)  Printf("TT index = %d   size =%d", indexSingleRndTrig, (int)fTrigTracks.size());
   }

   //___________________________________________________________
   if((fTypeOfAnal==kEmb || fTypeOfAnal == kEmbSingl) && indexSingleRndTrig > -1){ 
       //delta pT using embedded pythia events
       //delta pT analyzed only in events with REAL EVENT TT present  !!!!!!!!!!!  (condition above)
       // PWGJE/EMCALJetTasks/UserTasks/AliAnalysisTaskDeltaPtJEmb.cxx
       //AliJetResponseMaker
      //get deltaPt from embedding
      //++++++++++++++++++++++++++++++++++++++++
      //ANALYZE DELTA PT FOR ALL EMB MC JETS
      //++++++++++++++++++++++++++++++++++++++++
      AliVParticle* triggerHadron = static_cast<AliVParticle*>(trkArrayRec->At(fTrigTracks[indexSingleRndTrig]));
 
      AliEmcalJet  *jetEmb     = NULL;
      AliVParticle *constTrack = NULL;
      Bool_t bEmbJetCloseToTT  = kFALSE;

      if(jetContRec && triggerHadron){
         jetContRec->ResetCurrentID();
         while((jetRec = jetContRec->GetNextAcceptJet())) { //loop over reconstructed jets
            if(!jetRec) continue;
            if(!IsSignalJetInAcceptance(jetRec)) continue; //apply cuts on eta, pT ,area

            //skip the jet that contains TT 
            bEmbJetCloseToTT = kFALSE;
            for(Int_t iq=0; iq < jetRec->GetNumberOfTracks(); iq++) {
               constTrack = static_cast<AliVParticle*> (jetRec->TrackAt(iq,trkArrayRec)); //matched rec and emb tracks
               if(!constTrack) continue;
               if(constTrack != triggerHadron) continue;
               /*if(TMath::Abs(constTrack->Pt()  - triggerHadron->Pt())>0.01) continue;
               if(TMath::Abs(constTrack->Eta() - triggerHadron->Eta())>0.01) continue;
               if(TMath::Abs(RelativePhi(constTrack->Phi(), triggerHadron->Phi()))>0.01) continue;*/
               //if(fDebug>21)  Printf("EMB FIND TT COMPARE TRACK PT %f %f", constTrack->Pt(), triggerHadron->Pt());
               bEmbJetCloseToTT = kTRUE;
               break;   
            }        

            if(bEmbJetCloseToTT) continue;//skip the jet that contains TT track 

            jetEmb = NULL;
            jetEmb = jetRec->ClosestJet();
            if(!jetEmb){
               //if(fDebug>21) Printf("embedded jet does not exists, returning");
               continue;
            }
            if(jetEmb->Pt()<1e-3){
                //if(fDebug>21)  Printf("SKIP MATCH WITH EMBEDDED GHOST JET");
                continue; 
            }
            if((fTypeOfAnal==kEmb) &&  (jetEmb->Pt()<6.0))  continue; // some hard cut on pT of emb jet ???
    
            if(!IsSignalJetInAcceptance(jetEmb)) continue; //apply cuts on eta, pT ,area on the embedded jet
          
           //Check fraction of tracks from generator level jet in rec level jet
            //At least 50% of embedded jet has to be in the reconstructed jet
            if(GetFractionSharedPt(jetRec,jetContRec,jetEmb, jetContGen) < fMinFractionShared) continue;
 
            for(Int_t ir=0; ir < kRho-1; ir++){
               //1 Dim  distribution
               fhDeltaPtEmb[ir]->Fill(
                       jetRec->Pt() - jetRec->Area() * fRhoRec[ir] - jetEmb->Pt());

               //2 Dim distribution
               fhDeltaPtEmb2D[ir]->Fill(jetEmb->Pt(),                
                       jetRec->Pt() - jetRec->Area() * fRhoRec[ir] - jetEmb->Pt());
             }
          }
       }
   }//end of embedding

   if(fTypeOfAnal == kEmb || fTypeOfAnal == kEmbSingl) return kTRUE; 

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //   delta pT in reconstructed data
   for(Int_t ir=0; ir < kRho-1; ir++){
      fhRhoIncl[ir]->Fill((Float_t) fRhoRec[ir]); 
   }


   Double_t deltapt[kRho-1], phiTT = -1000., etaTT = -1000.;       
   Double_t ncoll = -1.0; 
   if(centralityPercentile>=0.)  ncoll = GetNcoll(centralityPercentile);

   if(indexSingleRndTrig > -1){ //get phi and eta of the TT
      AliVParticle* hadronTT = static_cast<AliVParticle*>(trkArrayRec->At(fTrigTracks[indexSingleRndTrig]));
      if(hadronTT){ 
         phiTT = hadronTT->Phi();
         etaTT = hadronTT->Eta();
      }
   }

   for(Int_t irc=0; irc<fNofRandomCones; irc++){ 

      //generate certain number of random cones per event
      GetDeltaPt(kRho, fRhoRec, &deltapt[0], phiTT, etaTT, 1.0/ncoll); //FK//????? prob exlude RC ??? 1/Ncoll 
       

      for(Int_t ir=0; ir < kRho-1; ir++){
         //fill delta pt histograms in inclusive events 
         fhDeltaPtIncl[ir]->Fill(deltapt[ir]); 
  
         if(ntriggers>0){
            //fill delta pt histograms in events with TT (trigger track) 
            fhDeltaPt[ir]->Fill( deltapt[ir]); 
         }
      }
   }

   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //  h+jet in reconstructed data  

   if(ntriggers>0){

      for(Int_t ir=0; ir < kRho-1; ir++){
         //Estimate UE density in events with TT
         //Fill once per event
         fhRhoTT[ir]->Fill( (Float_t) fRhoRec[ir]); 
      }

      //TRIGGER PARTICLE LOOP
      bFirstCycle=kTRUE; 
      if(trkArrayRec){

         for(Int_t it=0; it<ntriggers; it++){ //loop over trigger configurations

            if(fTTtype==0){
               if(it != indexSingleRndTrig) continue;
            }

            AliVParticle* triggerHadron = static_cast<AliVParticle*>(trkArrayRec->At(fTrigTracks[it]));
            if(!triggerHadron) continue;
          
            fh1Ntriggers->Fill((Float_t) triggerHadron->Pt()); //trigger p 
         
         
            //JET LOOP
           if(jetContRec){
              jetContRec->ResetCurrentID();
              while((jetRec = jetContRec->GetNextAcceptJet())) {
                 if(!jetRec){
                     AliError(Form("%s: Could not receive jet", GetName()));
                     continue;
                  }
                  if(!IsSignalJetInAcceptance(jetRec)) continue;
                 
                  areaJet = jetRec->Area();
                  pTJet   = jetRec->Pt();
                 
                  if(bFirstCycle){
                     fhJetPhi->Fill( pTJet, RelativePhi(jetRec->Phi(),0.0));
                     fhJetEta->Fill( pTJet, jetRec->Eta());
                  }
                 
                  dphi = RelativePhi(triggerHadron->Phi(), jetRec->Phi());
                 
                  dfi = dphi; //-0.5*pi to 1.5*Pi
                  if(dfi < -0.5*TMath::Pi()) dfi += TMath::TwoPi();
                  if(dfi >  1.5*TMath::Pi()) dfi -= TMath::TwoPi();
                 
                  for(Int_t ir=0; ir< kRho;ir++){ 
                     fhDphiTriggerJet[ir]->Fill((Float_t) (pTJet - areaJet*fRhoRec[ir]), (Float_t) dfi); //Rongrong's analysis 
                  }
                  //-------------------------
                 
                  if(TMath::Abs(dphi) < fDphiCut) continue;  //Dphi cut between trigger and assoc
                  fhDphiTriggerJetAccept->Fill(dfi); //Accepted
                 
                 //Centrality, A, pTjet
                  tmpArray[0] =  areaJet;
                  for(Int_t ir=0; ir< kRho;ir++){ 
                     tmpArray[1] =  pTJet - areaJet*fRhoRec[ir];
                     fHJetSpec[ir]->Fill(tmpArray);
                 
                     if(ir<kRho-1){
                        fARhoTT[ir]->Fill((Float_t) (areaJet*fRhoRec[ir]));
                     }
                  }
               }//JET LOOP

               bFirstCycle = kFALSE;
            }//container exists
         }
      }
   }

   return kTRUE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetNcoll(Double_t centr){
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

}
//________________________________________________________________________
void AliAnalysisTaskHJetSpectra::Terminate(Option_t *){
   //Treminate 
   PostData(1, fOutput);

   // Mandatory
   fOutput = dynamic_cast<TList*> (GetOutputData(1)); // '1' refers to the output slot
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
   AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

   //Label of background
   TString bgtype[]={"Cone","CMS","Zero"};

   fRandom = new TRandom3(0);

   fOutput = new TList();
   fOutput->SetOwner(); // otherwise it produces leaks in merging
   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   TString name;

   Bool_t bHistRec =  (fTypeOfAnal != kEmb  && fTypeOfAnal != kEmbSingl && fTypeOfAnal != kKine);
   Bool_t bNotKine =  (fTypeOfAnal != kKine);

   //__________________________________________________________
   // Event statistics
   fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 6, -0.5, 5.5);
   fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
   fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
   fHistEvtSelection->GetXaxis()->SetBinLabel(3,"pile up (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");

   fOutput->Add(fHistEvtSelection);
   //___________________________________________________________
   // Hard trigger counter
   fh1Ntriggers = new TH1D("fh1Ntriggers","# of triggers",50,0.0,50.0);
   if(bHistRec)   fOutput->Add(fh1Ntriggers);

   //___________________________________________________________
   // trigger associated jet spectra (jet pT not corrected for UE)
   Int_t bw = (fUseDoubleBinPrecision==0) ? 1 : 2; //make larger bin width

   //jet associated to given TT 
   //A, pTjet  
   const Int_t    dimSpec   = 2;
   const Int_t    nBinsSpec[dimSpec]  = { 50, bw*160};
   const Double_t lowBinSpec[dimSpec] = { 0.0,  -20.0};
   const Double_t hiBinSpec[dimSpec]  = { 2.0,  300.0};
   for(Int_t ir=0; ir< kRho; ir++){
      fHJetSpec[ir] = new THnSparseF(
                   Form("fHJetSpec%s",bgtype[ir].Data()),
                   Form("Recoil jet spectrum [A,pTjet-A*rho%s]",bgtype[ir].Data()),
                   dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);
      if(bHistRec)  fOutput->Add(fHJetSpec[ir]);
   }
   
   //____________________________________________________________________
   //UE from cell median  [Centrality, rho, pTUe ]

   for(Int_t ir=0; ir< kRho-1; ir++){ //Skip Zero bg
      //rho in events with TT 
      fhRhoTT[ir] = new TH1F(Form("fhRho%s",bgtype[ir].Data()),
                           Form("Rho%s",bgtype[ir].Data()),40, 0.0, 20.0);
      if(bHistRec) fOutput->Add(fhRhoTT[ir]);

      //rho in inclusive events
      fhRhoIncl[ir] = (TH1F*) fhRhoTT[ir]->Clone(Form("fhRhoIncl%s",bgtype[ir].Data()));
      if(bHistRec) fOutput->Add(fhRhoIncl[ir]);
 
      // rho times area in events with TT 
      fARhoTT[ir] = new TH1F(Form("fARho%s",bgtype[ir].Data()),
                            Form("Area times rho %s",bgtype[ir].Data()),40, 0.0, 20.0);
      if(bHistRec) fOutput->Add(fARhoTT[ir]);
   }
   //_______________________________________________________________________
   // Delta pt distributions   

   for(Int_t ir=0; ir< kRho-1; ir++){
      //events with TT tracks
      fhDeltaPt[ir] = new TH1D(Form("fhDeltaPt%s",bgtype[ir].Data()),
                               Form("DeltaPt%s",bgtype[ir].Data()),  150, -50, 100);
      if(bHistRec) fOutput->Add(fhDeltaPt[ir]);

      //inclusive events
      fhDeltaPtIncl[ir] = (TH1D*) fhDeltaPt[ir]->Clone(Form("fhDeltaPtIncl%s",bgtype[ir].Data()));
      if(bHistRec) fOutput->Add(fhDeltaPtIncl[ir]);

      if(fTypeOfAnal == kEmb || fTypeOfAnal == kEmbSingl){
         //Embedded PYTHIA jets
         fhDeltaPtEmb[ir] = (TH1D*) fhDeltaPt[ir]->Clone(Form("fhDeltaPtEmb%s",bgtype[ir].Data()));
         fOutput->Add(fhDeltaPtEmb[ir]);

         fhDeltaPtEmb2D[ir] = new TH2D(Form("fhDeltaPtEmb2D%s",bgtype[ir].Data()),
                                       Form("fhDeltaPtEmb2D%s",bgtype[ir].Data()), 125,0, 250, 150, -50, 100);
         fOutput->Add(fhDeltaPtEmb2D[ir]);
      }
   }


   //_______________________________________________________________________
   //inclusive azimuthal and pseudorapidity histograms
   fhJetPhi   = new TH2F("fhJetPhi","Azim dist jets vs pTjet", 50, 0, 100, 50,-TMath::Pi(),TMath::Pi());
   if(bHistRec)  fOutput->Add(fhJetPhi);
   //-------------------------
   fhTrackPhi = new TH2F("fhTrackPhi","azim dist trig had vs pT,trk", 50, 0, 50, 50,-TMath::Pi(),TMath::Pi());
   if(bNotKine) fOutput->Add(fhTrackPhi);
   //-------------------------
   fhJetEta   = new TH2F("fhJetEta","Eta dist jets vs pTjet", 50,0, 100, 40,-0.9,0.9);
   if(bHistRec) fOutput->Add(fhJetEta);
   //-------------------------
   fhTrackEta = new TH2F("fhTrackEta","Eta dist trig had vs pT,trk", 50, 0, 50, 40,-0.9,0.9);
   if(bNotKine) fOutput->Add(fhTrackEta);
   //-------------------------
   fhTrackPt = new TH1F("fhTrackPt","pT,trk ", 50, 0, 50);
   if(bNotKine) fOutput->Add(fhTrackPt);
   //-------------------------
   fhVertexZ = new TH1F("fhVertexZ","z vertex",40,-20,20);
   if(bNotKine) fOutput->Add(fhVertexZ);
   //-------------------------
   fhVertexZAccept = new TH1F("fhVertexZAccept","z vertex after cut",40,-20,20);
   if(bNotKine) fOutput->Add(fhVertexZAccept);
   //-------------------------

   if(fTypeOfData!=kReal){
      fhVertexZMC = new TH1F("fhVertexZMC","z vertex",40,-20,20);
      fOutput->Add(fhVertexZMC);
      //-------------------------
      fhVertexZAcceptMC = new TH1F("fhVertexZAcceptMC","z vertex after cut",40,-20,20);
      fOutput->Add(fhVertexZAcceptMC);
   }
   //-------------------------
 
   //fhContribVtx = new TH1F("fhContribVtx","contrib to vtx",200,0,200);
   //fOutput->Add(fhContribVtx);
   //-------------------------
   //fhContribVtxAccept = new TH1F("fhContribVtxAccept","contrib to vtx after cut",200,0,200);
   //fOutput->Add(fhContribVtxAccept);
   //------------------------
   for(Int_t ir=0; ir < kRho; ir++){ //Rongrong's analysis
 
      fhDphiTriggerJet[ir] = new TH2F(Form("fhDphiTriggerJet%s",bgtype[ir].Data()),"Deltaphi trig-jet",75,-50,100, 100, -0.5*TMath::Pi(),1.5*TMath::Pi());
      if(bHistRec) fOutput->Add(fhDphiTriggerJet[ir]);
   }
   //-------------------------

   fhDphiTriggerJetAccept = new TH1F("fhDphiTriggerJetAccept","Deltaphi trig-jet after cut",50, -0.5*TMath::Pi(),1.5*TMath::Pi());
   if(bHistRec)  fOutput->Add(fhDphiTriggerJetAccept);
   //-------------------------
   fhCentrality = new TH1F("fhCentrality","Centrality",100,0,100);
   if(bNotKine) fOutput->Add(fhCentrality);
   //-------------------------
   fhCentralityV0M = new TH1F("hCentralityV0M","hCentralityV0M",100,0,100);
   if(bNotKine)fOutput->Add(fhCentralityV0M); 
   //-------------------------
   fhCentralityV0A = new TH1F("hCentralityV0A","hCentralityV0A",100,0,100);
   if(bNotKine) fOutput->Add(fhCentralityV0A); 
   //-------------------------
   fhCentralityV0C = new TH1F("hCentralityV0C","hCentralityV0C",100,0,100);
   if(bNotKine) fOutput->Add(fhCentralityV0C);
   //-------------------------
   fhCentralityZNA = new TH1F("hCentralityZNA","hCentralityZNA",100,0,100);
   if(bNotKine) fOutput->Add(fhCentralityZNA);

   if(fTypeOfAnal > kRec){ 
      fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
      fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
      fOutput->Add(fh1Xsec);
      
      fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
      fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
      fOutput->Add(fh1Trials);
      
      fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",500,0,500);
      fOutput->Add(fh1PtHard);
     
      if(fTypeOfData == kHijing){ 
         fhImpactParameter = new TH1D("fhImpactParameter","impact parameter distribution from HIJING",50,0,10);
         fOutput->Add(fhImpactParameter);
         
         fhImpactParameterTT = new TH1D("fhImpactParameterTT","b versus TT",50,0,10);
         fOutput->Add(fhImpactParameterTT);
      }
   }

   if(fTypeOfAnal== kEff ){
      for(Int_t ir=0; ir < kRho; ir++){
         fhJetPtGen[ir] = new TH1D(Form("fhJetPtGen%s",bgtype[ir].Data()),
                                   Form("Jet pT Gen %s",bgtype[ir].Data()),bw*100,-20,200);
         fOutput->Add(fhJetPtGen[ir]);
         
         fhJetPtGenVsJetPtRec[ir] = new TH2D(Form("fhJetPtGenVsJetPtRec%s",bgtype[ir].Data()),
                                             "", bw*160,-20,300, bw*160,-20,300); 
         fOutput->Add(fhJetPtGenVsJetPtRec[ir]);
         
         fhJetPtResolutionVsPtGen[ir] = new TH2D(Form("fhJetPtResolutionVsPtGen%s",bgtype[ir].Data()), 
                                                 "Resolution", 20,0,100, 35,-1.,0.4);
         fOutput->Add(fhJetPtResolutionVsPtGen[ir]);
      }
      
      Double_t bins [] = {0, 0.2,0.4,0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 5., 6., 8., 10., 20., 50.};
      Int_t nbins = sizeof(bins)/sizeof(Double_t)-1;
      
      fhPtTrkTruePrimRec = new TH2D("fhPtTrkTruePrimRec","",nbins, bins, 18,-0.9,0.9);
      fOutput->Add(fhPtTrkTruePrimRec);
      
      fhPtTrkTruePrimGen = (TH2D*) fhPtTrkTruePrimRec->Clone("fhPtTrkTruePrimGen");
      fOutput->Add(fhPtTrkTruePrimGen);
      
      fhPtTrkSecOrFakeRec = (TH2D*) fhPtTrkTruePrimRec->Clone("fhPtTrkSecOrFakeRec");
      fOutput->Add(fhPtTrkSecOrFakeRec);
   }

   if(fTypeOfAnal == kEff || fTypeOfAnal == kKine){
      fh1NtriggersGen = (TH1D*) fh1Ntriggers->Clone("fh1NtriggersGen");
      fOutput->Add(fh1NtriggersGen);
     //-------------------------
      for(Int_t ir=0; ir< kRho; ir++){
         fHJetSpecGen[ir] = (THnSparseF*)fHJetSpec[ir]->Clone(Form("fHJetSpecGen%s",bgtype[ir].Data()));
         fOutput->Add(fHJetSpecGen[ir]);
      }

     //-------------------------
      fhJetPhiGen = (TH2F*)  fhJetPhi->Clone("fhJetPhiGen");
      fOutput->Add(fhJetPhiGen);
      //-------------------------
      fhJetEtaGen   = (TH2F*) fhJetEta->Clone("fhJetEtaGen");
      fOutput->Add(fhJetEtaGen);
      //-------------------------
      fhTrackPtGen = (TH1F*) fhTrackPt->Clone("fhTrackPtGen");
      fOutput->Add(fhTrackPtGen);
      //------------------------- Rongrong's analysis
      for(Int_t ir=0; ir< kRho; ir++){
         name = Form("%sGen",fhDphiTriggerJet[ir]->GetName());
         fhDphiTriggerJetGen[ir] = (TH2F*) fhDphiTriggerJet[ir]->Clone(name.Data());
         fOutput->Add(fhDphiTriggerJetGen[ir]);
      } 
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

Double_t AliAnalysisTaskHJetSpectra::RelativePhi(Double_t mphi,Double_t vphi){
   //Get relative azimuthal angle of two particles -pi to pi
   if      (vphi < -TMath::Pi()) vphi += TMath::TwoPi();
   else if (vphi > TMath::Pi())  vphi -= TMath::TwoPi();

   if      (mphi < -TMath::Pi()) mphi += TMath::TwoPi();
   else if (mphi > TMath::Pi())  mphi -= TMath::TwoPi();

   Double_t dphi = mphi - vphi;
   if      (dphi < -TMath::Pi()) dphi += TMath::TwoPi();
   else if (dphi > TMath::Pi())  dphi -= TMath::TwoPi();

   return dphi;//dphi in [-Pi, Pi]
}


//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::EstimateBgCone(Int_t icont){
   //Estimate background rho by means of integrating track pT outside identified jet cones
   Double_t rhoPerpCone = 0.0;
   
   Double_t pTleading  = -1.0;
   Double_t phiLeading = 1000.;
   Double_t etaLeading = 1000.;

   AliEmcalJet* jet = NULL;
   AliJetContainer *jetCont = GetJetContainer(icont);

   if(!jetCont) return 0.0;
   jetCont->ResetCurrentID(); //loop over jets
   while((jet = jetCont->GetNextAcceptJet())) {
      if(!jet) continue;

      if(!IsSignalJetInAcceptance(jet)) continue;

      if(pTleading < jet->Pt()){
         pTleading  = jet->Pt();
         phiLeading = jet->Phi();
         etaLeading = jet->Eta();
      }
   } 
   if(pTleading < 0.0) return 0.0;
   Double_t etawin = fTrackEtaWindow - fPerpConeRadius;

   if(TMath::Abs(etaLeading) > etawin){
      //make sure that the perp cone will be in acceptance
      etaLeading = (etaLeading>0) ? etawin : -etawin;
   }

   Double_t phileftcone  = phiLeading + TMath::Pi()/2;
   Double_t phirightcone = phiLeading - TMath::Pi()/2;

   rhoPerpCone +=  GetConePt(etaLeading, phileftcone,  fPerpConeRadius, icont);
   rhoPerpCone +=  GetConePt(etaLeading, phirightcone, fPerpConeRadius, icont);

   //normalize total pT by two times cone are 
   rhoPerpCone = rhoPerpCone/(2*TMath::Pi()*fPerpConeRadiusSquared);


   return rhoPerpCone;
}
//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetDeltaR(Double_t phi1, Double_t phi2, Double_t eta1, Double_t eta2){
   //angular distance between two jets
   Double_t dphi = RelativePhi(phi1,phi2);
   Double_t deta = eta1 - eta2;
   return sqrt(dphi*dphi + deta*deta); 

}
              
//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetFractionSharedPt(AliEmcalJet *jRec, AliJetContainer *jconRec, AliEmcalJet *jGen, AliJetContainer *jconGen){

   //get fraction of pT shared by reconstructed and generated level jet
   if(!jRec)    return -1.0;
   if(!jconRec)  return -1.0;
   if(!jGen)    return -1.0;
   if(!jconGen)  return -1.0;

   Double_t fraction = 0., sumPt = 0.;
   Double_t jetPt2 = jGen->Pt();
   //Int_t idxGen, idxRec;
   AliVParticle *pgen, *prec;
   if(jetPt2>0){

      for(Int_t ig=0; ig< jGen->GetNumberOfTracks(); ig++) {
         //idxGen = (Int_t) jGen->TrackAt(ig);
         pgen = static_cast<AliVParticle*>(jGen->TrackAt(ig, jconGen->GetParticleContainer()->GetArray()));

      
         for(Int_t ir=0; ir< jRec->GetNumberOfTracks(); ir++){
            //idxRec = (Int_t) jRec->TrackAt(ir);
            prec = static_cast<AliVParticle*>(jRec->TrackAt(ir, jconRec->GetParticleContainer()->GetArray()));

            if(TMath::Abs(prec->GetLabel()) == TMath::Abs(pgen->GetLabel())){

               if(fTypeOfAnal == kEmb ||  fTypeOfAnal == kEmbSingl){
                  //All embedded tracks have the same label check also spatial coordinates
                 if(TMath::Abs(prec->Eta() - pgen->Eta()) > 1e-4) continue;
                  if(TMath::Abs(RelativePhi(prec->Phi(), pgen->Phi())) > 1e-4) continue;
                  if(TMath::Abs(prec->Pt() - pgen->Pt()) > 1e-4) continue;
                  //if(fDebug>20){
                     //Printf("fraction TRACK REC eta = %f; phi = %f; pt = %f", prec->Eta(), prec->Phi(), prec->Pt());
                     //Printf("fraction TRACK GEN eta = %f; phi = %f; pt = %f", pgen->Eta(), pgen->Phi(), pgen->Pt());
                  //}
 
               }

               sumPt +=  pgen->Pt();
               break;
            }
         }
      }


      fraction = sumPt/jetPt2;
   } else{
     fraction = -1;
   }

   //if(fDebug>20) Printf("fraction return = %f ",fraction);

   return fraction;

} 
