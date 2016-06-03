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

#include "AliJetContainer.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisTaskHJetSpectra.h"
#include "AliHeader.h" 
#include "AliRunLoader.h"  
#include "AliVVZERO.h"
#include "AliAODZDC.h" 
#include "AliVZDC.h"

using namespace std;

// ANALYSIS OF HIGH PT HADRON TRIGGER ASSOCIATED SPECTRUM OF RECOIL JETS IN P+PB
// Author Filip Krizek   (7.Oct. 2015)

ClassImp(AliAnalysisTaskHJetSpectra)
//________________________________________________________________________________________

AliAnalysisTaskHJetSpectra::AliAnalysisTaskHJetSpectra(): 
AliAnalysisTaskEmcalJet("AliAnalysisTaskHJetSpectra", kTRUE),  
 fCollisionSystem(0), fTypeOfData(0), fTypeOfAnal(0),
  fUseDefaultVertexCut(1), fUsePileUpCut(1),  
  fRhoTaskName(), fRhoTaskNameMC(),
 fSignalJetRadius(0.4), fSignalJetRadiusSquared(fSignalJetRadius*fSignalJetRadius),
fSignalJetEtaWindow(0.9 - fSignalJetRadius), fTrackEtaWindow(0.9), fMinTrackPt(0.150), fMinJetArea(0.0),  
fCentralityType("V0A"), fMinFractionShared(0.5), 
fCrossSection(0.0), fTrials(0.0), fImpParam(-1.0), fRandom(0), fHelperClass(0), fInitializedLocal(0),
fTTlow(8.0), fTThigh(9.0), fTTtype(0), fDphiCut(TMath::Pi()-0.6), fUseDoubleBinPrecision(0),
fTTLowToSkip(-1.0), fTTHighToSkip(-1.0), fHistEvtSelection(0x0),fhKTAreaPt(0x0),
 fhVertexZ(0x0), fhVertexZAccept(0x0), fhVertexZMC(0x0), fhVertexZAcceptMC(0x0),
 fhDphiTriggerJetAccept(0x0),
fhCentralityV0M(0x0), fhCentralityV0A(0x0), fhCentralityV0C(0x0), fhCentralityZNA(0x0),
fhCentralityV0MTT(0x0), fhCentralityV0ATT(0x0), fhCentralityV0CTT(0x0), fhCentralityZNATT(0x0),
fRhoRec(kRho),fRhoMC(kRho), fCentralityBins(kCAll),
fNofRandomCones(1),
fZVertexCut(10.0),fCutPhi(0.6),
fpyVtx(3)
{
   //default constructor
   for(Int_t ic =0; ic<kCAll; ic++){
      fh1Ntriggers[ic]=NULL;
      fh1TriggerMult[ic]=NULL;
      fh1NtriggersGen[ic]=NULL;
      fh1TriggerMultGen[ic]=NULL;
      fhJetPhi[ic]=NULL;
      fhJetPhiGen[ic]=NULL;
      fhTrackPhi[ic]=NULL;
      fhJetEta[ic]=NULL;
      fhJetEtaGen[ic]=NULL;
      fhTrackEta[ic]=NULL;
      fhTrackPt[ic]=NULL;
      fhTrackPtGen[ic]=NULL;
      fhCentrality[ic]=NULL;
      fhVzeroATotMult[ic]=NULL;
      fhVzeroATotMultTT[ic]=NULL;
      fhZNAEnergy[ic]=NULL;
      fhZNAEnergyTT[ic]=NULL;
      fhTrackMultiplicity[ic]=NULL;
      fhTrackMultiplicityTT[ic]=NULL;
      fhZNAVzeroATrack[ic]=NULL;
      fhZNAVzeroATrackTT[ic]=NULL;
      fhImpactParameter[ic]=NULL;
      fhImpactParameterTT[ic]=NULL;

      for(Int_t ir=0; ir<kRho; ir++){
        fHJetSpec[ic][ir]=NULL;
        fhDphiTriggerJet[ic][ir]=NULL;
        fhDphiTriggerJetGen[ic][ir]=NULL;
        fhJetPtGen[ic][ir]=NULL;
        fhJetPtGenVsJetPtRec[ic][ir]=NULL;
        fhJetPtResolutionVsPtGen[ic][ir]=NULL;
      
      }
      
      for(Int_t ir=0; ir<kRho-1; ir++){
         fhRhoTT[ic][ir]=NULL;
         fhRhoIncl[ic][ir]=NULL;
         fARhoTT[ic][ir]=NULL;
         fhDeltaPt[ic][ir]=NULL; 
         fhDeltaPtEmb[ic][ir]=NULL; 
         fhDeltaPtEmb2D[ic][ir]=NULL;
         fhDeltaPtEmbPerp[ic][ir]=NULL; 
         fhDeltaPtEmbPerp2D[ic][ir]=NULL;
         fhDeltaPtEmbBc2Bc[ic][ir]=NULL;
         fhDeltaPtEmbBc2Bc2D[ic][ir]=NULL;
         fhDeltaPtIncl[ic][ir]=NULL;
      }

      fhPtTrkTruePrimRec[ic]=NULL;
      fhPtTrkTruePrimGen[ic]=NULL;
      fhPtTrkSecOrFakeRec[ic]=NULL;
   }
   //Centrality bin borders
   fCentralityBins[0]=0.;
   fCentralityBins[1]=20.;
   fCentralityBins[2]=50.;
   fCentralityBins[3]=100.;
   fCentralityBins[4]=1e6; //centrality overflow bin
 
   ficb[0]=-1;
   ficb[1]=-1;
   ftmpArray[0]=0.;              
   ftmpArray[1]=0.;              
   ftmpArrayX[0]=0.;              
   ftmpArrayX[1]=0.;              
   ftmpArrayX[2]=0.;              

   for(Int_t i=0; i<999; i++){
     frhovec[i] = 0.;
     fTrigTracksGen[i]=0x0;
     fTrigTracks[i]=0x0;
   }
}

//________________________________________________________________________
AliAnalysisTaskHJetSpectra::AliAnalysisTaskHJetSpectra(const char *name) : 
AliAnalysisTaskEmcalJet(name,kTRUE),  
fCollisionSystem(0), fTypeOfData(0), fTypeOfAnal(0),
  fUseDefaultVertexCut(1), fUsePileUpCut(1),
   fRhoTaskName(), fRhoTaskNameMC(),
fSignalJetRadius(0.4), fSignalJetRadiusSquared(fSignalJetRadius*fSignalJetRadius), 
fSignalJetEtaWindow(0.9 - fSignalJetRadius), fTrackEtaWindow(0.9), fMinTrackPt(0.150), fMinJetArea(0.0),   
fCentralityType("V0A"),   fMinFractionShared(0.5), 
fCrossSection(0.0), fTrials(0.0), fImpParam(-1.0), fRandom(0), fHelperClass(0), fInitializedLocal(0), 
fTTlow(8.0), fTThigh(9.0), fTTtype(0), fDphiCut(TMath::Pi()-0.6), fUseDoubleBinPrecision(0),
fTTLowToSkip(-1.0), fTTHighToSkip(-1.0),
fHistEvtSelection(0x0), fhKTAreaPt(0x0), 
 fhVertexZ(0x0), fhVertexZAccept(0x0), fhVertexZMC(0x0), fhVertexZAcceptMC(0x0),
fhDphiTriggerJetAccept(0x0),
fhCentralityV0M(0x0), fhCentralityV0A(0x0), fhCentralityV0C(0x0), fhCentralityZNA(0x0),
fhCentralityV0MTT(0x0), fhCentralityV0ATT(0x0), fhCentralityV0CTT(0x0), fhCentralityZNATT(0x0),
/*fh1Xsec(0x0), fh1Trials(0x0), fh1PtHard(0x0),*/
fRhoRec(kRho), fRhoMC(kRho), fCentralityBins(kCAll),
fNofRandomCones(1),
fZVertexCut(10.0),fCutPhi(0.6),
fpyVtx(3)
{
//Constructor

   for(Int_t ic =0; ic<kCAll; ic++){
      fh1Ntriggers[ic]=NULL;
      fh1TriggerMult[ic]=NULL;
      fh1NtriggersGen[ic]=NULL;
      fh1TriggerMultGen[ic]=NULL;
      fhJetPhi[ic]=NULL;
      fhJetPhiGen[ic]=NULL;
      fhTrackPhi[ic]=NULL;
      fhJetEta[ic]=NULL;
      fhJetEtaGen[ic]=NULL;
      fhTrackEta[ic]=NULL;
      fhTrackPt[ic]=NULL;
      fhTrackPtGen[ic]=NULL;
      fhCentrality[ic]=NULL;
      fhVzeroATotMult[ic]=NULL;
      fhVzeroATotMultTT[ic]=NULL;
      fhZNAEnergy[ic]=NULL;
      fhZNAEnergyTT[ic]=NULL;
      fhTrackMultiplicity[ic]=NULL;
      fhTrackMultiplicityTT[ic]=NULL;
      fhZNAVzeroATrack[ic]=NULL;
      fhZNAVzeroATrackTT[ic]=NULL; 
      fhImpactParameter[ic]=NULL;
      fhImpactParameterTT[ic]=NULL;

      for(Int_t ir=0; ir<kRho; ir++){
         fHJetSpec[ic][ir]=NULL;
         fhDphiTriggerJet[ic][ir]=NULL;
         fhDphiTriggerJetGen[ic][ir]=NULL;
         fhJetPtGen[ic][ir]=NULL;
         fhJetPtGenVsJetPtRec[ic][ir]=NULL;
         fhJetPtResolutionVsPtGen[ic][ir]=NULL;
      
      }
      
      for(Int_t ir=0; ir<kRho-1; ir++){
         fhRhoTT[ic][ir]=NULL;
         fhRhoIncl[ic][ir]=NULL;
         fARhoTT[ic][ir]=NULL;
         fhDeltaPt[ic][ir]=NULL;
         fhDeltaPtEmb[ic][ir]=NULL; 
         fhDeltaPtEmb2D[ic][ir]=NULL;
         fhDeltaPtEmbPerp[ic][ir]=NULL; 
         fhDeltaPtEmbPerp2D[ic][ir]=NULL;
         fhDeltaPtEmbBc2Bc[ic][ir]=NULL;
         fhDeltaPtEmbBc2Bc2D[ic][ir]=NULL;
         fhDeltaPtIncl[ic][ir]=NULL;
      }

      fhPtTrkTruePrimRec[ic]=NULL;
      fhPtTrkTruePrimGen[ic]=NULL;
      fhPtTrkSecOrFakeRec[ic]=NULL;
   }

   //Centrality bin borders
   fCentralityBins[0]=0.;
   fCentralityBins[1]=20.;
   fCentralityBins[2]=50.;
   fCentralityBins[3]=100.;
   fCentralityBins[4]=1e6; //centrality overflow bin
 
   ficb[0]=-1;
   ficb[1]=-1;
   ftmpArray[0]=0.;              
   ftmpArray[1]=0.;              
   ftmpArrayX[0]=0.;              
   ftmpArrayX[1]=0.;              
   ftmpArrayX[2]=0.;              

   for(Int_t i=0; i<999; i++){
      frhovec[i] = 0.;
      fTrigTracksGen[i]=0x0;
      fTrigTracks[i]=0x0;
   }
   DefineOutput(1, TList::Class());
}
  
//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetConePt(Double_t eta, Double_t phi, Double_t radius, AliParticleContainer *recTrkCont, Bool_t isGen){
   //sum up pt inside a cone
   Double_t tmpConePt = 0.0;

   if(!recTrkCont) return 0.0;

   AliVParticle* tmpTrack=NULL; 

   if(recTrkCont){
      recTrkCont->ResetCurrentID();
      while((tmpTrack = recTrkCont->GetNextAcceptParticle())){
         if(!tmpTrack) continue; 
         if(IsTrackInAcceptance(tmpTrack, isGen)){ 
            if(GetDeltaR(tmpTrack->Phi(), phi, tmpTrack->Eta(), eta) < radius){
               tmpConePt = tmpConePt + tmpTrack->Pt();
            }
         }
      }
   }
   return tmpConePt;
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
      
      mcHeader->PrimaryVertex(fpyVtx);
      return (Double_t) (fpyVtx[2]);
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
void  AliAnalysisTaskHJetSpectra::GetDeltaPt(Int_t nrho, TArrayD &rho, Double_t *dpt,
                                           Double_t ttPhi, Double_t ttEta, AliParticleContainer *recTrkCont, Bool_t isGen){ 
                                           //Double_t leadingJetExclusionProbability)

   //delta pt = pT in random cone - rho * Area of random cone
   // processes real reconstructed data. Exclude region around jet with TT   

   for(Int_t ir=0;ir<nrho;ir++){
      dpt[ir] = -10000.0;   // Set an invalid delta pt
   }

   // Define random cone Eta+Phi
   Bool_t coneValid = kTRUE;
   Double_t tmpRandConeEta = -fSignalJetEtaWindow  + fRandom->Rndm()*2*fSignalJetEtaWindow;
   Double_t tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();

   if(ttEta > -2.0){  //make sure that one generates RC far away from the trigger track jet
      while(GetDeltaR( tmpRandConePhi, ttPhi, tmpRandConeEta, ttEta)<fSignalJetRadius+0.1){
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
      Double_t conePt = GetConePt(tmpRandConeEta,tmpRandConePhi, fSignalJetRadius, recTrkCont, isGen);
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

   //if(fTypeOfAnal == kKine && fTypeOfData == kPythia){   
   //   fh1PtHard->Fill(GetPtHard());  //Fills cross section
   //}

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
   Double_t centralityPercentileV0A = -1.0;
   Double_t centralityPercentileV0C = -1.0;
   Double_t centralityPercentileV0M = -1.0;
   Double_t centralityPercentileZNA = -1.0;
   Double_t multVzero   = -1.0;
   Double_t energyZdcNA = -1.0;
   ficb[0] = 0;  //0=min bias 
   ficb[1] = -1; //-1 assigned centrality bin

   if(InputEvent()->GetVZEROData()){
      multVzero = InputEvent()->GetVZEROData()->GetMTotV0A();
   }

   //if(InputEvent()->GetZDCData()){
      // AliESDZDC.h :   C --> 1     A--->2
      //Double_t GetZNAEnergy() const {return (Double_t) fZDCN2Energy;}
      //Double_t GetZNCEnergy() const {return (Double_t) fZDCN1Energy;}
   //   energyZdcNA =  InputEvent()->GetZDCData()->GetZNAEnergy();
   //energyZdcNA =  InputEvent()->GetZDCN2Energy(); //should be ZNA
   //}
   AliVZDC *aodZDC = dynamic_cast<AliVZDC*> (InputEvent()->GetZDCData());
   if(aodZDC){
      const Double_t *ZNAtower = aodZDC->GetZNATowerEnergy();
      energyZdcNA = ZNAtower[0]*4.*82./208./12.96; //ZNA energy in TeV 
   }

   if(fCollisionSystem != kpp){   //KINE Check MC event vertex
      AliCentrality* tmpCentrality = InputEvent()->GetCentrality();
      if(!tmpCentrality){
         fHistEvtSelection->Fill(4);
         return kFALSE; //post data is in UserExec
      }else{
         centralityPercentile    = tmpCentrality->GetCentralityPercentile(fCentralityType.Data());
         centralityPercentileV0A = tmpCentrality->GetCentralityPercentile("V0A");
         centralityPercentileV0C = tmpCentrality->GetCentralityPercentile("V0C");
         centralityPercentileV0M = tmpCentrality->GetCentralityPercentile("V0M");
         centralityPercentileZNA = tmpCentrality->GetCentralityPercentile("ZNA");
      }

      if((centralityPercentile < fCentralityBins[0]) || (centralityPercentile >  fCentralityBins[4])){ //cut on centrality
         AliWarning(Form("Centrality value not valid (c=%E)",centralityPercentile)); 
         fHistEvtSelection->Fill(4);
         return kFALSE;
      }
     
      if(!(fTypeOfData == kPythia || fTypeOfAnal == kKine)){ 
         for(Int_t ic=0; ic<fCentralityBins.GetSize()-1;ic++){ 
            if(fCentralityBins[ic] <= centralityPercentile &&
               centralityPercentile < fCentralityBins[ic+1]){
               ficb[1] = ic+1; //MB is ficb[0]=0;  other centralities ficb[1]>0
            }
         }
      }
 
      for(Int_t ic=0; ic<2; ic++){
         if(ficb[ic]==-1) continue;  //MB +  CENT BIASED 
         fhCentrality[ficb[ic]]->Fill(centralityPercentile);
      } 

      //MB 
      fhCentralityV0M->Fill(centralityPercentileV0M); 
      fhCentralityV0A->Fill(centralityPercentileV0A);
      fhCentralityV0C->Fill(centralityPercentileV0C); 
      fhCentralityZNA->Fill(centralityPercentileZNA);
   }

   for(Int_t ic=0; ic<2; ic++){
      if(ficb[ic]==-1) continue;   //MB +CENT bias 
      fhVzeroATotMult[ficb[ic]]->Fill(multVzero);
      fhZNAEnergy[ficb[ic]]->Fill(energyZdcNA);

      if(fTypeOfData == kHijing){
         fImpParam = GetImpactParameter(); 
         fhImpactParameter[ficb[ic]]->Fill(fImpParam);
      }
   }
   //cout<<"energyZdcNA = "<<energyZdcNA<<endl;

   fHistEvtSelection->Fill(0); //Count Accepted input event

   // END EVENT SELECTION

   //_________________________________________________________________
   // JET+TRACK CONTAINERS
   AliJetContainer *jetContRec   = NULL; //AKTjet container from reconstruced tracks
   AliJetContainer *jetContGen   = NULL; //AKT jet container from  MC particles
   AliJetContainer *jetContRecKT = NULL; //KT jet container from reconstruced tracks
   AliJetContainer *jetContGenKT = NULL; //KT jet container from  MC particles

   AliEmcalJet  *jetGen = NULL;  //jet pointer MC jet
   AliEmcalJet  *jetRec = NULL;  //jet pointer real jet

   AliParticleContainer *trkContRec = NULL; //track contaier real reconstructed tracks 
   AliParticleContainer *parContGen = NULL; //particle container of MC particles

   AliVParticle *constTrackRec = NULL; //jet constituent
   AliVParticle *constTrackGen = NULL; //jet constituent
   //_________________________________________________________
   //READ JET TRACK CONTAINERS
   if(fTypeOfAnal == kRec){ 
      //REAL DATA
      jetContRec   = GetJetContainer(kContainerOne); //AKT jet
      trkContRec   = GetParticleContainer(kContainerOne); //reconstructed particle container 
      jetContRecKT = GetJetContainer(kContainerTwo);//KT jet
   }

   if(fTypeOfAnal == kEff || fTypeOfAnal == kEmb || fTypeOfAnal == kEmbSingl){ 
      //MC + EMB DATA
      jetContRec   = GetJetContainer(kContainerOne);  //AKT
      trkContRec   = GetParticleContainer(kContainerOne); //reconstructed particle container
      jetContGen   = GetJetContainer(kContainerTwo); //AKT generator level jets
      parContGen   = GetParticleContainer(kContainerTwo); //true MC particle container

      jetContRecKT = GetJetContainer(kContainerThree);  //KT  reconstructed jets
      jetContGenKT = GetJetContainer(kContainerFour);   //KT generator level jets
   }

   if(fTypeOfAnal == kKine){   // Kine written to the 0th container !!!!!!!!!!!!!!!!
      //KINE 
      jetContGen   = GetJetContainer(kContainerOne);  //AKT jets
      parContGen   = GetParticleContainer(kContainerOne); //kine particle container 

      jetContGenKT =  GetJetContainer(kContainerTwo); //KT jets
   }
 
   //if(fDebug>20)  Printf("POINTER TO CONTAINERS   JETrec=%p  TRKrec=%p   JETgen=%p   TRKgen=%p", 
   //                       jetContRec,trkContRec,jetContGen, parContGen);


   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //LOOP OVER TRACKS  SEARCH FOR TRIGGER CANDIDATES IN GEN TRACKS
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   Int_t indexSingleRndTrigGen = -1; //index of single random trigger
   AliVParticle* trackTTGen    = NULL; //pointer to the trigger hadron
   Int_t ntriggersGen = 0;

   if(fTypeOfAnal == kEff || fTypeOfAnal == kKine){


      if(parContGen){ 
         parContGen->ResetCurrentID();
         while((constTrackGen = (AliVParticle*)  parContGen->GetNextAcceptParticle())){ 
            if(!constTrackGen) continue;
       
            if(IsTrackInAcceptance(constTrackGen, kTRUE)){ 
               for(Int_t ic=0; ic<2; ic++){
                  if(ficb[ic]==-1) continue; 
                  fhTrackPtGen[ficb[ic]]->Fill(constTrackGen->Pt());  //inclusive pT spectrum of tracks
               }
               if((fTTlow <= constTrackGen->Pt()) && (constTrackGen->Pt() < fTThigh)){
                  if(ntriggersGen<999){
                     fTrigTracksGen[ntriggersGen] = constTrackGen;  //GEN trigger candidates
                     ntriggersGen++;
                  }
               }
            }
         }
      }
 

      for(Int_t ic=0; ic<2; ic++){
         if(ficb[ic]==-1) continue; 
         fh1TriggerMultGen[ficb[ic]]->Fill(ntriggersGen); 
      } 
 
      if(ntriggersGen>0){
         indexSingleRndTrigGen = fRandom->Integer(ntriggersGen); //Integer 0 ... ntriggers-1
         trackTTGen = (AliVParticle*) (fTrigTracksGen[indexSingleRndTrigGen]);
      }
   }
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //LOOP OVER TRACKS  SEARCH FOR TRIGGER CANDIDATES IN REC TRACKS
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   Int_t indexSingleRndTrig = -1; //index of single random trigger
   AliVParticle* trackTT = NULL;
   Int_t ntriggers = 0;
 
   if(fTypeOfAnal != kKine){
      
      if(trkContRec){
         trkContRec->ResetCurrentID();
         while((constTrackRec = (AliVParticle*) (trkContRec->GetNextAcceptParticle()))){
            if(!constTrackRec) continue;
      
            if(fTypeOfAnal == kEmb || fTypeOfAnal == kEmbSingl){
               //if(fDebug>99)  Printf("TRACK LABEL %d", track->GetLabel());
               //in embed evts search for TT only among real tracks; do not consider embedded tracks as trigger
               if(TMath::Abs(constTrackRec->GetLabel()) == 99999) continue;//9999 set in AddTaskHJetSpectra.C
            }
      
      
            if(IsTrackInAcceptance(constTrackRec, kFALSE)){  //rec  
               //Fill some inclusive spectra (Same analysis for gen level?)
               for(Int_t ic=0; ic<2; ic++){
                  if(ficb[ic]==-1) continue;
                  fhTrackPhi[ficb[ic]]->Fill(constTrackRec->Pt(), RelativePhi(constTrackRec->Phi(),0.0)); // phi = -pi,pi
                  fhTrackEta[ficb[ic]]->Fill(constTrackRec->Pt(), constTrackRec->Eta());
                  fhTrackPt[ficb[ic]]->Fill(constTrackRec->Pt());
               }

               if(fTTlow <= constTrackRec->Pt() && constTrackRec->Pt() < fTThigh){
                  if(ntriggers<999){
                     fTrigTracks[ntriggers] = constTrackRec;  //trigger candidates
                     ntriggers++; 
                  }
               }
            }
         }
      }
      
      //   SELECT SINGLE INCLUSIVE REC LEVEL TRIGGER 
      for(Int_t ic=0; ic<2; ic++){
         if(ficb[ic]==-1) continue;
         fh1TriggerMult[ficb[ic]]->Fill(ntriggers); 
      } 
      
      if(ntriggers>0){
         indexSingleRndTrig = fRandom->Integer(ntriggers); //Integer 0 ... ntriggers-1
         trackTT = (AliVParticle*) (fTrigTracks[indexSingleRndTrig]);
         //if(fDebug>20)  Printf("TT index = %d   size =%d", indexSingleRndTrig, (int)fTrigTracks.size());
      }
   }
 

   //___________________________________________________________
   // CALCULATE RHO
   fRhoRec.Reset(0.);  
   fRhoMC.Reset(0.);

   if(fTypeOfAnal == kRec || fTypeOfAnal == kEff || fTypeOfAnal == kEmb || fTypeOfAnal == kEmbSingl){
      fRhoRec[kConeRho] = EstimateBgCone( jetContRec,   trkContRec, trackTT, kFALSE); //container ID=0 reconstructed tracks
      fRhoRec[kCMSRho]  = EstimateBgKTcms(jetContRecKT, trkContRec, trackTT); 
      fRhoRec[kKtRho]   = EstimateBgKT(   jetContRecKT, trkContRec, trackTT);
      fRhoRec[kZeroRho] = 0.0; 
   }

   if(fTypeOfAnal == kEff || fTypeOfAnal == kEmb){ //rho in MC events 
      fRhoMC[kConeRho] = EstimateBgCone( jetContGen,   parContGen, trackTTGen, kTRUE); //container ID=1 mc particles
      fRhoMC[kCMSRho]  = EstimateBgKTcms(jetContGenKT, parContGen, trackTTGen);
      fRhoMC[kKtRho]   = EstimateBgKT(   jetContGenKT, parContGen, trackTTGen);
      fRhoMC[kZeroRho] = 0.0;
   }

   if(fTypeOfAnal == kEmbSingl){ //embedding single track
      fRhoMC[kConeRho] = 0.0; 
      fRhoMC[kCMSRho]  = 0.0; 
      fRhoMC[kKtRho]   = 0.0; 
      fRhoMC[kZeroRho] = 0.0;
   }

   if(fTypeOfAnal == kKine){ //rho in KINE MC events 
      fRhoMC[kConeRho] = EstimateBgCone( jetContGen,   parContGen, trackTTGen, kTRUE); //container ID=1 mc particles
      fRhoMC[kCMSRho]  = EstimateBgKTcms(jetContGenKT, parContGen, trackTTGen);
      fRhoMC[kKtRho]   = EstimateBgKT(   jetContGenKT, parContGen, trackTTGen);
      fRhoMC[kZeroRho] = 0.0;
   }

   //________________________________________________________
   //WHEN THIS IS TT REF SKIP EVENT THAT CONTAINS TT SIGNAL TRACK 
   /*
   Bool_t bStop = kFALSE;
   if(fTTLowToSkip>0 &&  fTTHighToSkip>0){
      if(fTypeOfAnal == kRec || fTypeOfAnal == kEff){
         if(trkContRec){
            trkContRec->ResetCurrentID();
            while((constTrackRec = (AliVParticle*) (trkContRec->GetNextAcceptParticle()))){
               if(!constTrackRec) continue;
               if(!IsTrackInAcceptance(constTrackRec, kFALSE)) continue; //reconstructed level tracks
               if(fTTLowToSkip <= constTrackRec->Pt() && constTrackRec->Pt() < fTTHighToSkip){
                  bStop = kTRUE;
                  break;
               }
            }
         }
      }
      
      if(fTypeOfAnal == kKine || fTypeOfAnal == kEff){
         if(parContGen){
            parContGen->ResetCurrentID();
            while((constTrackGen = (AliVParticle*) (parContGen->GetNextAcceptParticle()))){
               if(!constTrackGen) continue;
               if(!IsTrackInAcceptance(constTrackGen, kTRUE)) continue; //gen level physical primary
               if(fTTLowToSkip <= constTrackGen->Pt() && constTrackGen->Pt() < fTTHighToSkip){
                  bStop = kTRUE;
                  break;
               }
            }
         }
      }   
   }
   */
   //_________________________________________________________
   //EVALUATE SINGLE PARTICLE EFFICIENCY + FILL RESPONSE MATRIX

   if(fTypeOfAnal == kEff){

      //1) FILL HISTOS FOR SINGLE PARTICLE EFFICIENCY
      if(parContGen){
         parContGen->ResetCurrentID();
         while((constTrackGen = (AliVParticle*)(parContGen->GetNextAcceptParticle()))){
            if(!constTrackGen) continue;
            if(IsTrackInAcceptance(constTrackGen, kTRUE)){
               //pT spectrum of generator level particles
               for(Int_t ic=0; ic<2; ic++){
                  if(ficb[ic]==-1) continue;
                  fhPtTrkTruePrimGen[ficb[ic]]->Fill(constTrackGen->Pt(),constTrackGen->Eta());
               } 
            }  
         }

         //single particle efficiency and contamination
         Bool_t bRecPrim = kFALSE; //tags the reconstructed primary particles
         if(trkContRec && parContGen){ 
            trkContRec->ResetCurrentID();
            while((constTrackRec =(AliVParticle*) (trkContRec->GetNextAcceptParticle()))){
               if(!constTrackRec) continue;
               if(!IsTrackInAcceptance(constTrackRec, kFALSE)) continue; //reconstructed level tracks
               bRecPrim = kFALSE; //not yet matched to generator level physical primary

               parContGen->ResetCurrentID();
               while((constTrackGen = (AliVParticle*)(parContGen->GetNextAcceptParticle()))){
                  if(!constTrackGen) continue;
                  if(!IsTrackInAcceptance(constTrackGen, kTRUE)) continue; //gen level physical primary
                  if(TMath::Abs(constTrackRec->GetLabel()) == TMath::Abs(constTrackGen->GetLabel())){ 
                     //has the same label as reconstr track
 
                     bRecPrim = kTRUE;
                     for(Int_t ic=0; ic<2; ic++){
                        if(ficb[ic]==-1) continue;
                        fhPtTrkTruePrimRec[ficb[ic]]->Fill(constTrackGen->Pt(),constTrackGen->Eta()); //this is well recontr phys primary
                     }
                     break;
                  }//same label with rec particle
               }//loop over gen tracks
               if(!bRecPrim){
                  for(Int_t ic=0; ic<2; ic++){
                     if(ficb[ic]==-1) continue;
                     fhPtTrkSecOrFakeRec[ficb[ic]]->Fill(constTrackRec->Pt(),constTrackRec->Eta()); //matchnig to phys primary not found, this is fake or second.
                  }
               }
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
            if(!IsSignalJetInAcceptance(jetGen,kTRUE)) continue; //cuts on eta, pT ,area

            //if(fDebug>20) Printf("GEN JET phi=%f  eta=%f  pt=%f", jetGen->Phi(), jetGen->Eta(), jetGen->Pt());

            for(Int_t ir=0; ir< kRho; ir++){
               ptGenCorr = jetGen->Pt() - jetGen->Area()*fRhoMC[ir]; // correct for rho

               for(Int_t ic=0; ic<2; ic++){
                  if(ficb[ic]==-1) continue;
                  fhJetPtGen[ficb[ic]][ir]->Fill(ptGenCorr);
               }
            }
         }
      }
 
      //Find closest gen level+rec level  jets
      if(jetContRec){
         jetContRec->ResetCurrentID();
         while((jetRec = jetContRec->GetNextAcceptJet())) {
            if(!jetRec) continue;
            if(!IsSignalJetInAcceptance(jetRec,kTRUE)) continue; //cuts on eta, pT ,area

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
            
            if(!IsSignalJetInAcceptance(jetGen,kTRUE)) continue; //cuts on eta, pT ,area

            //check fraction of tracks from generator level jet in rec level jet
            //if(fDebug>20)  Printf("FRACTIONH SHARED = %f ", GetFractionSharedPt(jetRec,jetContRec,jetGen, jetContGen));

            //FK// if(GetFractionSharedPt(jetRec, jetContRec, jetGen, jetContGen) < fMinFractionShared) continue;

            //if(fDebug>20)  Printf("PASSED MIN FRACTION CRITERION ");

            for(Int_t ir=0; ir< kRho; ir++){
               ptGenCorr = jetGen->Pt() - jetGen->Area()*fRhoMC[ir]; //perp cone bg correction to pt 
               ptRecCorr = jetRec->Pt() - jetRec->Area()*fRhoRec[ir]; //perp cone bg correction to pt 

               for(Int_t ic=0; ic<2; ic++){
                  if(ficb[ic]==-1) continue;

                  fhJetPtGenVsJetPtRec[ficb[ic]][ir]->Fill(ptRecCorr, ptGenCorr); //response matrix

                  if(ptGenCorr>0){
                     fhJetPtResolutionVsPtGen[ficb[ic]][ir]->Fill(ptGenCorr,(ptRecCorr-ptGenCorr)/ptGenCorr); //jet pT resolution
                  }
               }
            }
         } 
      }//rec jet container exists
   }//analyze efficiency mode (response matrix + single particle efficiency) 


   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //if(bStop){
   //  fHistEvtSelection->Fill(5); //Rejected events
  //   return kTRUE; // SKIP EVENTS WHERE REF TT CLASS EVENT CONTAINS SIG TT TRACK 
  // }
   fHistEvtSelection->Fill(6); //Accepted events
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Double_t areaJet,  pTJet; 
   Double_t dphi, dfi;
   Bool_t bFirstCycle = kTRUE; 


   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //H-JET CORRELATIONS IN MC TRUTH
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(fTypeOfAnal == kEff || fTypeOfAnal == kKine){
      if(parContGen && ntriggersGen >0 && trackTTGen){
         bFirstCycle = kTRUE; 
         //for(Int_t it=0; it<ntriggersGen; it++){ //loop over trigger configurations
        
         //   if(fTTtype==0){
         //      if(it != indexSingleRndTrigGen) continue;
         //   }
         
            AliVParticle* triggerHadronGen = (AliVParticle*) trackTTGen; //(fTrigTracksGen[it]);
            //if(!triggerHadronGen) continue;
         
            for(Int_t ic=0; ic<2; ic++){
               if(ficb[ic]==-1) continue;
               fh1NtriggersGen[ficb[ic]]->Fill((Float_t) triggerHadronGen->Pt()); //trigger pT gen 
         
               if( fTypeOfData == kHijing){ //impact parameter for triggered events
                  fhImpactParameterTT[ficb[ic]]->Fill(fImpParam);
               }
            }
         
            //JET LOOP
            if(jetContGen){
               jetContGen->ResetCurrentID();
               while((jetGen = jetContGen->GetNextAcceptJet())) {
                  if(!jetGen){
                     AliError(Form("%s: Could not receive gen jet", GetName()));
                     continue;
                  }
                  if(!IsSignalJetInAcceptance(jetGen,kTRUE)) continue;
                 
                  areaJet = jetGen->Area();
                  pTJet   = jetGen->Pt();
                 
                  if(bFirstCycle){
                     for(Int_t ic=0; ic<2; ic++){
                        if(ficb[ic]==-1) continue;

                        fhJetPhiGen[ficb[ic]]->Fill( pTJet, RelativePhi(jetGen->Phi(),0.0));
                        fhJetEtaGen[ficb[ic]]->Fill( pTJet, jetGen->Eta());
                     }
                  }
                 
                  dphi = RelativePhi(triggerHadronGen->Phi(), jetGen->Phi());
                 
                  dfi = dphi; //-0.5*pi to 1.5*Pi
                  if(dfi < -0.5*TMath::Pi()) dfi += TMath::TwoPi();
                  if(dfi >  1.5*TMath::Pi()) dfi -= TMath::TwoPi();
                 
                  for(Int_t ir=0; ir< kRho;ir++){
                     for(Int_t ic=0; ic<2; ic++){
                        if(ficb[ic]==-1) continue;
                        fhDphiTriggerJetGen[ficb[ic]][ir]->Fill((Float_t) (pTJet - areaJet*fRhoMC[ir]), (Float_t) dfi); //Rongrong's analysis
                     } 
                  }
                  //-------------------------
                 
                  if(TMath::Abs(dphi) < fDphiCut) continue;  //Dphi cut between trigger and assoc
                 
                 //Centrality, A, pTjet
                  ftmpArray[0] =  areaJet;
                  for(Int_t ir=0; ir< kRho;ir++){ 
                     ftmpArray[1] =  pTJet - areaJet*fRhoMC[ir];
                     for(Int_t ic=0; ic<2; ic++){
                        if(ficb[ic]==-1) continue;
                        fHJetSpecGen[ficb[ic]][ir]->Fill(ftmpArray[0],ftmpArray[1]);
                     }
                  }
               }//JET LOOP
               bFirstCycle = kFALSE;
            }//container exists
        // }
      }
   }
   
   //++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(fTypeOfAnal == kKine) return kTRUE; 
   //++++++++++++++++++++++++++++++++++++++++++++++++++++

 

   if((fTypeOfAnal==kEmb || fTypeOfAnal == kEmbSingl) && trackTT){ 
       //delta pT using embedded pythia events
       //delta pT analyzed only in events with REAL EVENT TT present  !!!!!!!!!!!  (condition above)
       // PWGJE/EMCALJetTasks/UserTasks/AliAnalysisTaskDeltaPtJEmb.cxx
       //AliJetResponseMaker
      //get deltaPt from embedding
      //++++++++++++++++++++++++++++++++++++++++
      //ANALYZE DELTA PT FOR ALL EMB MC JETS
      //++++++++++++++++++++++++++++++++++++++++
 
      AliEmcalJet  *jetEmb     = NULL;
      Bool_t bEmbJetCloseToTT  = kFALSE;

      if(jetContRec){
         jetContRec->ResetCurrentID();
         while((jetRec = jetContRec->GetNextAcceptJet())) { //loop over reconstructed jets
            if(!jetRec) continue;
            if(!IsSignalJetInAcceptance(jetRec,kTRUE)) continue; //apply cuts on eta, pT ,area

            //skip the jet that contains TT 
            bEmbJetCloseToTT = kFALSE;

            if(jetRec->Pt() > trackTT->Pt()*0.5){  // if jet contains TT it has to have at least pT of TT

               for(Int_t iq=0; iq < jetRec->GetNumberOfTracks(); iq++) {
                  constTrackRec = static_cast<AliVParticle*> (jetRec->TrackAt(iq,trkContRec->GetArray())); //matched rec and emb tracks
                  if(!constTrackRec) continue;
                  if(constTrackRec != trackTT) continue;
                  //if(fDebug>21)  Printf("EMB FIND TT COMPARE TRACK PT %f %f", constTrack->Pt(), triggerHadron->Pt());
                  bEmbJetCloseToTT = kTRUE;
                  break;   
               }       
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
            //if((fTypeOfAnal==kEmb) &&  (jetEmb->Pt()<6.0))  continue; // some hard cut on pT of emb jet ???
    
            if(!IsSignalJetInAcceptance(jetEmb,kTRUE)) continue; //apply cuts on eta, pT ,area on the embedded jet
          
           //Check fraction of tracks from generator level jet in rec level jet
            //At least 50% of embedded jet has to be in the reconstructed jet
            if(GetFractionSharedPt(jetRec, jetContRec, jetEmb, jetContGen) < fMinFractionShared) continue;
 
            for(Int_t ir=0; ir < kRho-1; ir++){
               for(Int_t ic=0; ic<2; ic++){
                  if(ficb[ic]==-1) continue;
                  //1 Dim  distribution
                  fhDeltaPtEmb[ficb[ic]][ir]->Fill(
                       jetRec->Pt() - jetRec->Area() * fRhoRec[ir] - jetEmb->Pt());

                  //2 Dim distribution
                  fhDeltaPtEmb2D[ficb[ic]][ir]->Fill(jetEmb->Pt(),                
                       jetRec->Pt() - jetRec->Area() * fRhoRec[ir] - jetEmb->Pt());


                  if(trackTT){ //embedded track/jet is perp to TT
                     dphi = TMath::Abs(RelativePhi(trackTT->Phi(), jetEmb->Phi())); 
                     if(TMath::Pi()/4 <  dphi && dphi < 3*TMath::Pi()/4){
                        //1 Dim  distribution
                        fhDeltaPtEmbPerp[ficb[ic]][ir]->Fill(
                             jetRec->Pt() - jetRec->Area() * fRhoRec[ir] - jetEmb->Pt());
                        
                        //2 Dim distribution
                        fhDeltaPtEmbPerp2D[ficb[ic]][ir]->Fill(jetEmb->Pt(),                
                             jetRec->Pt() - jetRec->Area() * fRhoRec[ir] - jetEmb->Pt());

                     }else if(dphi >= 3*TMath::Pi()/4 ){
                        //embedded track back-to-back w.r.t the TT
                        //1 Dim  distribution
                        fhDeltaPtEmbBc2Bc[ficb[ic]][ir]->Fill(
                             jetRec->Pt() - jetRec->Area() * fRhoRec[ir] - jetEmb->Pt());

                        //2 Dim distribution
                        fhDeltaPtEmbBc2Bc2D[ficb[ic]][ir]->Fill(jetEmb->Pt(), 
                             jetRec->Pt() - jetRec->Area() * fRhoRec[ir] - jetEmb->Pt());
                     }
                  }
               }
            }
         }
      }
   }//end of embedding


   // RECONSTRUCTED DATA ANALYSIS

   for(Int_t ir=0; ir < kRho-1; ir++){
      for(Int_t ic=0; ic<2; ic++){
         if(ficb[ic]==-1) continue;

         fhRhoIncl[ficb[ic]][ir]->Fill((Float_t) fRhoRec[ir]); 
      }
   }
   if(ntriggers>0){

      for(Int_t ir=0; ir < kRho-1; ir++){
         //Estimate UE density in events with TT
         //Fill once per event
         for(Int_t ic=0; ic<2; ic++){
            if(ficb[ic]==-1) continue;

            fhRhoTT[ficb[ic]][ir]->Fill( (Float_t) fRhoRec[ir]); 
         }
      }
   }

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   if(fTypeOfAnal == kEmb || fTypeOfAnal == kEmbSingl) return kTRUE; 
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 

   // CALCULATE  DELTA PT IN RECONSTRUCTED DATA WITH RANDOM CONES
   Double_t deltapt[kRho-1], phiTT = 0., etaTT = -1000.;       
   //Double_t ncoll = -1.0; 
   //if(centralityPercentile>=0.)  ncoll = GetNcoll(centralityPercentile);
   Bool_t  bRecJetCloseToTT = kFALSE;
   //Exclude region around TT from delta pt calculation

   if(trackTT){ //get phi and eta of the TT or the reconstructed jet that contains TT
      phiTT = trackTT->Phi(); // TT is proxy for jet
      etaTT = trackTT->Eta();

      if(jetContRec){  //find jet that contains TT if it exists
         jetContRec->ResetCurrentID();
         while((jetRec = jetContRec->GetNextAcceptJet())) { //loop over reconstructed jets
            if(!jetRec) continue;
            if(!IsSignalJetInAcceptance(jetRec,kTRUE)) continue;

            //skip the jet that contains TT 
            bRecJetCloseToTT = kFALSE;

            if(jetRec->Pt() > trackTT->Pt()*0.5){  // if jet contains TT it has to have at leat pT of TT
               for(Int_t iq=0; iq < jetRec->GetNumberOfTracks(); iq++){
                  constTrackRec = (AliVParticle*) (jetRec->TrackAt(iq,trkContRec->GetArray())); //matched rec and emb tracks
                  if(!constTrackRec) continue;
                  if(constTrackRec != trackTT) continue;
                  phiTT = jetRec->Phi();
                  etaTT = jetRec->Eta();
                  bRecJetCloseToTT = kTRUE;
                  break;
               }
            }
            if(bRecJetCloseToTT) break;
         }
      }
   }

   for(Int_t irc=0; irc<fNofRandomCones; irc++){ 

      //generate certain number of random cones per event
      GetDeltaPt(kRho-1, fRhoRec, &deltapt[0], phiTT, etaTT, trkContRec, kFALSE);// 1.0/ncoll);//FK//????? prob exlude RC ??? 1/Ncoll 
       

      for(Int_t ir=0; ir < kRho-1; ir++){
         for(Int_t ic=0; ic<2; ic++){
            if(ficb[ic]==-1) continue;
            //fill delta pt histograms in inclusive events 
            fhDeltaPtIncl[ficb[ic]][ir]->Fill(deltapt[ir]); 
  
            if(ntriggers>0){
               //fill delta pt histograms in events with TT (trigger track) 
               fhDeltaPt[ficb[ic]][ir]->Fill( deltapt[ir]); 
            }
         }
      }
   }
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // Track Multiplicity
   
   if(trkContRec){
      
      Int_t mult   = 0;
      trkContRec->ResetCurrentID();
      while((constTrackRec = (AliVParticle*) (trkContRec->GetNextAcceptParticle()))){
         if(!constTrackRec) continue;
         if(!IsTrackInAcceptance(constTrackRec, kFALSE)) continue; //reconstructed level tracks
         mult++;
      }

      for(Int_t ic=0; ic<2; ic++){
         if(ficb[ic]==-1) continue;
         //MB  or CENT biases
         fhTrackMultiplicity[ficb[ic]]->Fill(mult);

         ftmpArrayX[0] = energyZdcNA;
         ftmpArrayX[1] = multVzero;
         ftmpArrayX[2] = mult;
         fhZNAVzeroATrack[ficb[ic]]->Fill(ftmpArrayX);
         
         if(trackTT){//TT biased + TT&&CENT biased distributions
            fhTrackMultiplicityTT[ficb[ic]]->Fill(mult);
            fhVzeroATotMultTT[ficb[ic]]->Fill(multVzero);
            fhZNAEnergyTT[ficb[ic]]->Fill(energyZdcNA);
            fhZNAVzeroATrackTT[ficb[ic]]->Fill(ftmpArrayX);
         }
      }
   
        

      if(trackTT){ //TT biased distributions
         fhCentralityV0MTT->Fill((Float_t) centralityPercentileV0M);
         fhCentralityV0ATT->Fill((Float_t) centralityPercentileV0A);
         fhCentralityV0CTT->Fill((Float_t) centralityPercentileV0C);
         fhCentralityZNATT->Fill((Float_t) centralityPercentileZNA);
      }
   }
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   //  H+JET IN RECONSTRUCTED DATA  

   if(ntriggers>0 && trackTT){

      bFirstCycle=kTRUE; 
      if(trkContRec){

       //  for(Int_t it=0; it<ntriggers; it++){ //loop over trigger configurations

       //     if(fTTtype==0){
       //        if(it != indexSingleRndTrig) continue;
      //      }

            AliVParticle* triggerHadron = (AliVParticle*)  trackTT;// (fTrigTracks[it]);
            //if(!triggerHadron) continue;
          
            for(Int_t ic=0; ic<2; ic++){
               if(ficb[ic]==-1) continue;

               fh1Ntriggers[ficb[ic]]->Fill((Float_t) triggerHadron->Pt()); //trigger p 
            }
         
            //JET LOOP
           if(jetContRec){
              jetContRec->ResetCurrentID();
              while((jetRec = jetContRec->GetNextAcceptJet())) {
                 if(!jetRec){
                     AliError(Form("%s: Could not receive jet", GetName()));
                     continue;
                  }
                  if(!IsSignalJetInAcceptance(jetRec,kTRUE)) continue;
                 
                  areaJet = jetRec->Area();
                  pTJet   = jetRec->Pt();
                 
                  if(bFirstCycle){
                     for(Int_t ic=0; ic<2; ic++){
                        if(ficb[ic]==-1) continue;
                        fhJetPhi[ficb[ic]]->Fill( pTJet, RelativePhi(jetRec->Phi(),0.0));
                        fhJetEta[ficb[ic]]->Fill( pTJet, jetRec->Eta());
                     }
                  }
                 
                  dphi = RelativePhi(triggerHadron->Phi(), jetRec->Phi());
                 
                  dfi = dphi; //-0.5*pi to 1.5*Pi
                  if(dfi < -0.5*TMath::Pi()) dfi += TMath::TwoPi();
                  if(dfi >  1.5*TMath::Pi()) dfi -= TMath::TwoPi();

                  for(Int_t ic=0; ic<2; ic++){
                     if(ficb[ic]==-1) continue;
 
                     for(Int_t ir=0; ir< kRho;ir++){ 
                        fhDphiTriggerJet[ficb[ic]][ir]->Fill((Float_t) (pTJet - areaJet*fRhoRec[ir]), (Float_t) dfi); //Rongrong's analysis
                     } 
                  }
                  //-------------------------
                 
                  if(TMath::Abs(dphi) < fDphiCut) continue;  //Dphi cut between trigger and assoc
                  fhDphiTriggerJetAccept->Fill(dfi); //Accepted
                 
                 //Centrality, A, pTjet
                  ftmpArray[0] =  areaJet;
                  for(Int_t ir=0; ir< kRho;ir++){ 
                     ftmpArray[1] =  pTJet - areaJet*fRhoRec[ir];

                     for(Int_t ic=0; ic<2; ic++){
                        if(ficb[ic]==-1) continue;
 
                        fHJetSpec[ficb[ic]][ir]->Fill(ftmpArray[0],ftmpArray[1]);
                 
                        if(ir<kRho-1){
                           fARhoTT[ficb[ic]][ir]->Fill((Float_t) (areaJet*fRhoRec[ir]));
                        }
                     }
                  } 
               }//JET LOOP

               bFirstCycle = kFALSE;
            }//container exists
         //}
      }
   }

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

   //Label of background
   TString bgtype[]={"Perp","CMS","KT","Zero"};

   fRandom = new TRandom3(0);

   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   TString name;

   Bool_t bHistRec =  (fTypeOfAnal != kEmb  && fTypeOfAnal != kEmbSingl && fTypeOfAnal != kKine);
   Bool_t bNotKine =  (fTypeOfAnal != kKine);

   Int_t icmax = (fTypeOfData == kPythia ||fTypeOfAnal == kKine) ? 1 : kCAll;
   //__________________________________________________________
   // Event statistics
   fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 7, -0.5, 6.5);
   fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
   fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
   fHistEvtSelection->GetXaxis()->SetBinLabel(3,"pile up (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(6,"sig TT (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(7,"sig TT (accepted)");

   fOutput->Add(fHistEvtSelection);
   //___________________________________________________________
   // Hard trigger counter
   for(Int_t ic =0; ic<icmax; ic++){
      name = (ic==0) ? Form("fh1NtriggersMB") : 
       Form("fh1Ntriggers%d%d",TMath::Nint(fCentralityBins[ic-1]),TMath::Nint(fCentralityBins[ic])); 
      fh1Ntriggers[ic] = new TH1D(name.Data(),"# of triggers",50,0.0,50.0);
      if(bHistRec)   fOutput->Add((TH1D*)fh1Ntriggers[ic]);
 
      name = (ic==0) ? Form("fh1TriggerMultMB") : 
       Form("fh1TriggerMult%d%d",TMath::Nint(fCentralityBins[ic-1]),TMath::Nint(fCentralityBins[ic])); 
      fh1TriggerMult[ic] = new TH1D(name.Data(),"# of triggers",50,0.0,50.0);
      if(bHistRec)   fOutput->Add((TH1D*)fh1TriggerMult[ic]);
   }
   //___________________________________________________________
   // trigger associated jet spectra (jet pT not corrected for UE)
   Int_t bw = (fUseDoubleBinPrecision==0) ? 1 : 2; //make larger bin width

   //jet associated to given TT 
   //A, pTjet  
   const Int_t    dimSpec   = 2;
   const Int_t    nBinsSpec[dimSpec]  = { 50, bw*160};
   const Double_t lowBinSpec[dimSpec] = { 0.0,  -20.0};
   const Double_t hiBinSpec[dimSpec]  = { 2.0,  300.0};

   for(Int_t ic =0; ic<icmax; ic++){
      for(Int_t ir=0; ir< kRho; ir++){
         name = (ic==0) ? Form("fHJetSpecMB%s",bgtype[ir].Data()) : 
                          Form("fHJetSpec%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 

         /*  
         fHJetSpec[ic][ir] = new THnSparseF(
                      name.Data(),
                      Form("Recoil jet spectrum [A,pTjet-A*rho%s]",bgtype[ir].Data()),
                      dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);
         if(bHistRec)  fOutput->Add((THnSparseF*) fHJetSpec[ic][ir]);
         */
         fHJetSpec[ic][ir] = new TH2D(
                      name.Data(),
                      Form("Recoil jet spectrum [A,pTjet-A*rho%s]",bgtype[ir].Data()),
                      nBinsSpec[0], lowBinSpec[0], hiBinSpec[0], nBinsSpec[1], lowBinSpec[1],hiBinSpec[1]);
         if(bHistRec)  fOutput->Add((TH2D*) fHJetSpec[ic][ir]);
      }
   } 
   //____________________________________________________________________
   //UE from cell median  [Centrality, rho, pTUe ]

   for(Int_t ic =0; ic<icmax; ic++){
      for(Int_t ir=0; ir< kRho-1; ir++){ //Skip Zero bg
         name = (ic==0) ? Form("fhRhoMB%s",bgtype[ir].Data()) : 
                          Form("fhRho%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 

         //rho in events with TT 
         fhRhoTT[ic][ir] = new TH1F(name.Data(),
                              Form("Rho%s",bgtype[ir].Data()),80, 0.0, 40.0);
         if(bNotKine) fOutput->Add((TH1F*) fhRhoTT[ic][ir]);
      
         //rho in inclusive events
         name = (ic==0) ? Form("fhRhoInclMB%s",bgtype[ir].Data()) : 
                          Form("fhRhoIncl%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 

         fhRhoIncl[ic][ir] = (TH1F*) fhRhoTT[ic][ir]->Clone(name.Data());
         if(bNotKine) fOutput->Add((TH1F*) fhRhoIncl[ic][ir]);
      
         // rho times area in events with TT 
         name = (ic==0) ? Form("fARhoMB%s",bgtype[ir].Data()) : 
                          Form("fARho%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 

         fARhoTT[ic][ir] = new TH1F(name.Data(),
                               Form("Area times rho %s",bgtype[ir].Data()),80, 0.0, 40.0);
         if(bHistRec) fOutput->Add((TH1F*) fARhoTT[ic][ir]);
      }
   }
   //_______________________________________________________________________
   // Delta pt distributions   

   for(Int_t ic =0; ic<icmax; ic++){
      for(Int_t ir=0; ir< kRho-1; ir++){
         //events with TT tracks
         name = (ic==0) ? Form("fhDeltaPtMB%s",bgtype[ir].Data()) : 
                          Form("fhDeltaPt%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 

         fhDeltaPt[ic][ir] = new TH1D(name.Data(),
                                  Form("DeltaPt%s",bgtype[ir].Data()),  150, -50, 100);
         if(bHistRec) fOutput->Add((TH1D*) fhDeltaPt[ic][ir]);

         //inclusive events
         name = (ic==0) ? Form("fhDeltaPtInclMB%s",bgtype[ir].Data()) : 
                          Form("fhDeltaPtIncl%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 

         fhDeltaPtIncl[ic][ir] = (TH1D*) fhDeltaPt[ic][ir]->Clone(name.Data());
         if(bHistRec) fOutput->Add((TH1D*) fhDeltaPtIncl[ic][ir]);

         if(fTypeOfAnal == kEmb || fTypeOfAnal == kEmbSingl){
            //Embedded PYTHIA jets
            //1D
            name = (ic==0) ? Form("fhDeltaPtEmbMB%s",bgtype[ir].Data()) : 
                             Form("fhDeltaPtEmb%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 

            fhDeltaPtEmb[ic][ir] = (TH1D*) fhDeltaPt[ic][ir]->Clone(name.Data());
            fOutput->Add((TH1D*) fhDeltaPtEmb[ic][ir]);

            //1D perp to TT
            name = (ic==0) ? Form("fhDeltaPtEmbPerpMB%s",bgtype[ir].Data()) : 
                             Form("fhDeltaPtEmbPerp%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 

            fhDeltaPtEmbPerp[ic][ir] = (TH1D*) fhDeltaPt[ic][ir]->Clone(name.Data());
            fOutput->Add((TH1D*) fhDeltaPtEmbPerp[ic][ir]);

            //1D back-to back to TT
            name = (ic==0) ? Form("fhDeltaPtEmbBc2BcMB%s",bgtype[ir].Data()) : 
                             Form("fhDeltaPtEmbBc2Bc%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 

            fhDeltaPtEmbBc2Bc[ic][ir]= (TH1D*) fhDeltaPt[ic][ir]->Clone(name.Data());
            fOutput->Add((TH1D*) fhDeltaPtEmbBc2Bc[ic][ir]);

           //2D
            name = (ic==0) ? Form("fhDeltaPtEmb2DMB%s",bgtype[ir].Data()) : 
                             Form("fhDeltaPtEmb2D%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 

            fhDeltaPtEmb2D[ic][ir] = new TH2D(name.Data(),
                Form("fhDeltaPtEmb2D%s",bgtype[ir].Data()), 125,0, 250, 150, -50, 100);
            fOutput->Add((TH2D*)fhDeltaPtEmb2D[ic][ir]);


           //2D perp to TT
            name = (ic==0) ? Form("fhDeltaPtEmbPerp2DMB%s",bgtype[ir].Data()) : 
                             Form("fhDeltaPtEmbPerp2D%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 

            fhDeltaPtEmbPerp2D[ic][ir] = (TH2D*) fhDeltaPtEmb2D[ic][ir]->Clone(name.Data());
            fOutput->Add((TH2D*)fhDeltaPtEmbPerp2D[ic][ir]);


            //2D back-to-back p to TT
            name = (ic==0) ? Form("fhDeltaPtEmbBc2Bc2DMB%s",bgtype[ir].Data()) : 
                             Form("fhDeltaPtEmbBc2Bc2D%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 

            fhDeltaPtEmbBc2Bc2D[ic][ir]=(TH2D*) fhDeltaPtEmb2D[ic][ir]->Clone(name.Data());
            fOutput->Add((TH2D*)fhDeltaPtEmbBc2Bc2D[ic][ir]);

         }
      }
   }

   //_______________________________________________________________________
   //KT jets area versus PT
   Double_t binsPt [] = {0, 0.05,0.15, 0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6.,6.5, 7., 8., 9., 10.,15., 20., 50., 100.};
   Int_t nbinsPt = sizeof(binsPt)/sizeof(Double_t)-1;
 
   fhKTAreaPt  = new TH2F("fhKTAreaPt","KT jet Area vs Pt",nbinsPt, binsPt, 50,0,2);
   if(fTypeOfAnal == kRec)  fOutput->Add(fhKTAreaPt);

   //_______________________________________________________________________
   //inclusive azimuthal and pseudorapidity histograms
   for(Int_t ic =0; ic<icmax; ic++){
       name = (ic==0) ? Form("fhJetPhiMB") : 
                        Form("fhJetPhi%d%d",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic])); 

      fhJetPhi[ic] = new TH2F(name.Data(),"Azim dist jets vs pTjet", 50, 0, 100, 50,-TMath::Pi(),TMath::Pi());
      if(bHistRec)  fOutput->Add((TH2F*)fhJetPhi[ic]);

      //-------------------------
      name = (ic==0) ? Form("fhTrackPhiMB") : 
                       Form("fhTrackPhi%d%d",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic])); 

      fhTrackPhi[ic] = new TH2F(name.Data(),"azim dist trig had vs pT,trk", 50, 0, 50, 50,-TMath::Pi(),TMath::Pi());
      if(bNotKine) fOutput->Add((TH2F*)fhTrackPhi[ic]);
      //-------------------------
      name = (ic==0) ? Form("fhJetEtaMB") : 
                       Form("fhJetEta%d%d",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic])); 

      fhJetEta[ic] = new TH2F(name.Data(),"Eta dist jets vs pTjet", 50,0, 100, 40,-0.9,0.9);
      if(bHistRec) fOutput->Add((TH2F*)fhJetEta[ic]);
      //-------------------------
      name = (ic==0) ? Form("fhTrackEtaMB") : 
                       Form("fhTrackEta%d%d",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic])); 

      fhTrackEta[ic] = new TH2F(name.Data(),"Eta dist trig had vs pT,trk", 50, 0, 50, 40,-0.9,0.9);
      if(bNotKine) fOutput->Add((TH2F*) fhTrackEta[ic]);
      //-------------------------
      name = (ic==0) ? Form("fhTrackPtMB") : 
                       Form("fhTrackPt%d%d",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic])); 

      fhTrackPt[ic] = new TH1F(name.Data(),"pT,trk ", 50, 0, 50);
      if(bNotKine) fOutput->Add((TH1F*) fhTrackPt[ic]);
   }
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
   for(Int_t ic =0; ic<icmax; ic++){
      for(Int_t ir=0; ir < kRho; ir++){ //Rongrong's analysis
         name = (ic==0) ? Form("fhDphiTriggerJetMB%s",bgtype[ir].Data()) : 
                          Form("fhDphiTriggerJet%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                  TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 

         fhDphiTriggerJet[ic][ir] = new TH2F(name.Data(),"Deltaphi trig-jet",75,-50,100, 100, -0.5*TMath::Pi(),1.5*TMath::Pi());
         if(bHistRec) fOutput->Add((TH2F*) fhDphiTriggerJet[ic][ir]);
      }
   }
   //-------------------------

   fhDphiTriggerJetAccept = new TH1F("fhDphiTriggerJetAccept","Deltaphi trig-jet after cut",50, -0.5*TMath::Pi(),1.5*TMath::Pi());
   if(bHistRec)  fOutput->Add(fhDphiTriggerJetAccept);
   //-------------------------
   for(Int_t ic =0; ic<icmax; ic++){
      name = (ic==0) ? Form("fhCentralityMB") : 
                       Form("fhCentrality%d%d",TMath::Nint(fCentralityBins[ic-1]),
                               TMath::Nint(fCentralityBins[ic])); 

      fhCentrality[ic] = new TH1F(name.Data(),"Centrality",100,0,100);
      if(bNotKine) fOutput->Add((TH1F*) fhCentrality[ic]);
   }
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
   //-------------------------
   fhCentralityV0MTT = (TH1F*) fhCentralityV0M->Clone("fhCentralityV0MTT");
   if(bHistRec) fOutput->Add(fhCentralityV0MTT);
   //-------------------------
   fhCentralityV0ATT = (TH1F*) fhCentralityV0A->Clone("fhCentralityV0ATT");
   if(bHistRec) fOutput->Add(fhCentralityV0ATT); 
   //-------------------------
   fhCentralityV0CTT = (TH1F*) fhCentralityV0C->Clone("fhCentralityV0CTT");
   if(bHistRec) fOutput->Add(fhCentralityV0CTT); 
   //-------------------------
   fhCentralityZNATT = (TH1F*) fhCentralityZNA->Clone("fhCentralityZNATT");
   if(bHistRec) fOutput->Add(fhCentralityZNATT);
   //-----------------------------------------------------
   for(Int_t ic =0; ic<icmax; ic++){
      name = (ic==0) ? Form("fhVzeroATotMultMB") : 
                       Form("fhVzeroATotMult%d%d",TMath::Nint(fCentralityBins[ic-1]),
                               TMath::Nint(fCentralityBins[ic])); 

      //  vzero multiplicity
      fhVzeroATotMult[ic] = new TH1F(name.Data(),"hVzeroATotMult",1000,0,1000);
      if(bHistRec) fOutput->Add((TH1F*) fhVzeroATotMult[ic]);

      name = (ic==0) ? Form("fhVzeroATotMultTTMB") : 
                       Form("fhVzeroATotMultTT%d%d",TMath::Nint(fCentralityBins[ic-1]),
                               TMath::Nint(fCentralityBins[ic])); 

      fhVzeroATotMultTT[ic] = (TH1F*)  fhVzeroATotMult[ic]->Clone(name.Data());  
      if(bHistRec) fOutput->Add((TH1F*) fhVzeroATotMultTT[ic]);
   }

   
   //-----------------------------------------------------
   //  ZDC ZNA energy
   
   for(Int_t ic =0; ic<icmax; ic++){
      name = (ic==0) ? Form("fhZNAEnergyMB") : 
                       Form("fhZNAEnergy%d%d",TMath::Nint(fCentralityBins[ic-1]),
                               TMath::Nint(fCentralityBins[ic])); 


      fhZNAEnergy[ic] = new TH1F(name.Data(),"fhZNAEnergy",1000,0,200);
      if(bHistRec) fOutput->Add((TH1F*)fhZNAEnergy[ic]);
 
      name = (ic==0) ? Form("fhZNAEnergyTTMB") : 
                       Form("fhZNAEnergyTT%d%d",TMath::Nint(fCentralityBins[ic-1]),
                               TMath::Nint(fCentralityBins[ic])); 

      fhZNAEnergyTT[ic] = (TH1F*)  fhZNAEnergy[ic]->Clone(name.Data());  
      if(bHistRec) fOutput->Add((TH1F*) fhZNAEnergyTT[ic]);
   }
   
   //-----------------------------------------------------
   //   track multiplicity
   
   for(Int_t ic =0; ic<icmax; ic++){
      name = (ic==0) ? Form("fhTrackMultiplicityMB") : 
                       Form("fhTrackMultiplicity%d%d",TMath::Nint(fCentralityBins[ic-1]),
                               TMath::Nint(fCentralityBins[ic])); 

      fhTrackMultiplicity[ic] = new TH1D(name.Data(),"fhTrackMultiplicity",1000,0,1000);
      if(bNotKine) fOutput->Add((TH1D*)fhTrackMultiplicity[ic]);

      name = (ic==0) ? Form("fhTrackMultiplicityTTMB") : 
                       Form("fhTrackMultiplicityTT%d%d",TMath::Nint(fCentralityBins[ic-1]),
                               TMath::Nint(fCentralityBins[ic])); 

      fhTrackMultiplicityTT[ic] = (TH1D*)fhTrackMultiplicity[ic]->Clone(name.Data());
      if(bNotKine) fOutput->Add(fhTrackMultiplicityTT[ic]);
   }

   //-----------------------------------------------------
   // ZNA energy versus Vzero mult. versus track mult. in all events 
   
   const Int_t    dimZ   = 3;
   const Int_t    nBinsZ[dimZ]  = { 100,     100,   25 };
   const Double_t lowBinZ[dimZ] = { 0.0,     0.0,   0.0};
   const Double_t hiBinZ[dimZ]  = { 200.0, 1000.0, 250};

   for(Int_t ic =0; ic<icmax; ic++){
      name = (ic==0) ? Form("fhZNAVzeroATrackMB") : 
                       Form("fhZNAVzeroATrack%d%d",TMath::Nint(fCentralityBins[ic-1]),
                                                   TMath::Nint(fCentralityBins[ic])); 

      fhZNAVzeroATrack[ic] = new THnSparseF(
                      name.Data(),
                      Form("ZNA, V0A mult, track mult"),
                      dimZ, nBinsZ,lowBinZ,hiBinZ);
      if(fTypeOfAnal ==  kRec)  fOutput->Add((THnSparseF*) fhZNAVzeroATrack[ic]);


      //-------------
      // ZNA energy versus Vzero mult. versus track mult. in events with TT
      name = (ic==0) ? Form("fhZNAVzeroATrackTTMB") : 
                       Form("fhZNAVzeroATrackTT%d%d",TMath::Nint(fCentralityBins[ic-1]),
                                                   TMath::Nint(fCentralityBins[ic])); 

      fhZNAVzeroATrackTT[ic] = (THnSparseF*) fhZNAVzeroATrack[ic]->Clone(name.Data()); 
      if(fTypeOfAnal ==  kRec)  fOutput->Add((THnSparseF*) fhZNAVzeroATrackTT[ic]);
   }
    
   //-----------------------------------------------------

  /*
   if(fTypeOfAnal == kKine){ 
      fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
      fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
      fOutput->Add(fh1Xsec);
      
      fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
      fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
      fOutput->Add(fh1Trials);
      
      //fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",1000,0,1000);
      //fOutput->Add(fh1PtHard);
   } 
*/
   if(fTypeOfData == kHijing){ 
      for(Int_t ic =0; ic<kCAll; ic++){
         name = (ic==0) ? Form("fhImpactParameterMB") : 
                          Form("fhImpactParameter%d%d",TMath::Nint(fCentralityBins[ic-1]),
                               TMath::Nint(fCentralityBins[ic])); 

         fhImpactParameter[ic] = new TH1D(name.Data(),"impact parameter distribution from HIJING",50,0,10);
         fOutput->Add((TH1D*) fhImpactParameter[ic]);
 

         name = (ic==0) ? Form("fhImpactParameterTTMB") : 
                          Form("fhImpactParameterTT%d%d",TMath::Nint(fCentralityBins[ic-1]),
                               TMath::Nint(fCentralityBins[ic])); 
           
         fhImpactParameterTT[ic] = new TH1D(name.Data(),"b versus TT",50,0,10);
         fOutput->Add((TH1D*) fhImpactParameterTT[ic]);
      }
   }

   if(fTypeOfAnal== kEff ){
      for(Int_t ic =0; ic<icmax; ic++){
         for(Int_t ir=0; ir < kRho; ir++){
            name = (ic==0) ? Form("fhJetPtGenMB%s",bgtype[ir].Data()) : 
                             Form("fhJetPtGen%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                               TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 
 
            fhJetPtGen[ic][ir] = new TH1D(name.Data(),
                                      Form("Jet pT Gen %s",bgtype[ir].Data()),bw*160,-20,300);
            fOutput->Add((TH1D*) fhJetPtGen[ic][ir]);


            name = (ic==0) ? Form("fhJetPtGenVsJetPtRecMB%s",bgtype[ir].Data()) : 
                             Form("fhJetPtGenVsJetPtRec%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                               TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 
            
            fhJetPtGenVsJetPtRec[ic][ir] = new TH2D(name.Data(),
                                                "", bw*160,-20,300, bw*160,-20,300); 
            fOutput->Add((TH2D*) fhJetPtGenVsJetPtRec[ic][ir]);


            name = (ic==0) ? Form("fhJetPtResolutionVsPtGenMB%s",bgtype[ir].Data()) : 
                             Form("fhJetPtResolutionVsPtGen%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                               TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 
 
            fhJetPtResolutionVsPtGen[ic][ir] = new TH2D(name.Data(), 
                                                    "Resolution", 20,0,100, 35,-1.,0.4);
            fOutput->Add((TH2D*) fhJetPtResolutionVsPtGen[ic][ir]);
         }
      }
      
      Double_t bins [] = {0, 0.2,0.4,0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.5, 3., 3.5, 4., 5., 6., 8., 10., 20., 50.};
      Int_t nbins = sizeof(bins)/sizeof(Double_t)-1;
      
      for(Int_t ic =0; ic<icmax; ic++){
         name = (ic==0) ? Form("fhPtTrkTruePrimRecMB") : 
                          Form("fhPtTrkTruePrimRec%d%d",TMath::Nint(fCentralityBins[ic-1]),
                            TMath::Nint(fCentralityBins[ic])); 
 
         fhPtTrkTruePrimRec[ic] = new TH2D(name.Data(),"",nbins, bins, 18,-0.9,0.9);
         fOutput->Add((TH2D*) fhPtTrkTruePrimRec[ic]);
 

         name = (ic==0) ? Form("fhPtTrkTruePrimGenMB") : 
                          Form("fhPtTrkTruePrimGen%d%d",TMath::Nint(fCentralityBins[ic-1]),
                            TMath::Nint(fCentralityBins[ic])); 
         
         fhPtTrkTruePrimGen[ic] = (TH2D*) fhPtTrkTruePrimRec[ic]->Clone(name.Data());
         fOutput->Add((TH2D*) fhPtTrkTruePrimGen[ic]);

         name = (ic==0) ? Form("fhPtTrkSecOrFakeRecMB") : 
                          Form("fhPtTrkSecOrFakeRec%d%d",TMath::Nint(fCentralityBins[ic-1]),
                           TMath::Nint(fCentralityBins[ic])); 
         
         fhPtTrkSecOrFakeRec[ic] = (TH2D*) fhPtTrkTruePrimRec[ic]->Clone(name.Data());
         fOutput->Add((TH2D*) fhPtTrkSecOrFakeRec[ic]);
      }
   }

   if(fTypeOfAnal == kEff || fTypeOfAnal == kKine){
      for(Int_t ic =0; ic<icmax; ic++){
         name = (ic==0) ? Form("fh1NtriggersGenMB") : 
                          Form("fh1NtriggersGen%d%d",TMath::Nint(fCentralityBins[ic-1]),
                            TMath::Nint(fCentralityBins[ic])); 
 
         fh1NtriggersGen[ic] = (TH1D*) fh1Ntriggers[ic]->Clone(name.Data());
         fOutput->Add((TH1D*) fh1NtriggersGen[ic]);
  

          name = (ic==0) ? Form("fh1TriggerMultGenMB") : 
                           Form("fh1TriggerMultGen%d%d",TMath::Nint(fCentralityBins[ic-1]),
                            TMath::Nint(fCentralityBins[ic])); 
      
         fh1TriggerMultGen[ic] = (TH1D*) fh1TriggerMult[ic]->Clone();
         fOutput->Add((TH1D*) fh1TriggerMultGen[ic]);
        
        //-------------------------
         for(Int_t ir=0; ir< kRho; ir++){
            name = (ic==0) ? Form("fHJetSpecGenMB%s",bgtype[ir].Data()) : 
                             Form("fHJetSpecGen%d%d%s",TMath::Nint(fCentralityBins[ic-1]),
                                TMath::Nint(fCentralityBins[ic]),bgtype[ir].Data()); 
 

            fHJetSpecGen[ic][ir] = (TH2D*)fHJetSpec[ic][ir]->Clone(name.Data());
            fOutput->Add((TH2D*) fHJetSpecGen[ic][ir]);
            //fHJetSpecGen[ic][ir] = (THnSparseF*)fHJetSpec[ic][ir]->Clone(name.Data());
            //fOutput->Add((THnSparseF*) fHJetSpecGen[ic][ir]);
         }
        
         //-------------------------
         name = (ic==0) ? Form("fhJetPhiGenMB") : 
                          Form("fhJetPhiGen%d%d",TMath::Nint(fCentralityBins[ic-1]),
                            TMath::Nint(fCentralityBins[ic])); 

         fhJetPhiGen[ic] = (TH2F*)  fhJetPhi[ic]->Clone(name.Data());
         fOutput->Add((TH2F*) fhJetPhiGen[ic]);
         //-------------------------
         name = (ic==0) ? Form("fhJetEtaGenMB") : 
                          Form("fhJetEtaGen%d%d",TMath::Nint(fCentralityBins[ic-1]),
                            TMath::Nint(fCentralityBins[ic])); 


         fhJetEtaGen[ic]   = (TH2F*) fhJetEta[ic]->Clone(name.Data());
         fOutput->Add((TH2F*) fhJetEtaGen[ic]);
         //-------------------------
         name = (ic==0) ? Form("fhTrackPtGenMB") : 
                          Form("fhTrackPtGen%d%d",TMath::Nint(fCentralityBins[ic-1]),
                            TMath::Nint(fCentralityBins[ic])); 

         fhTrackPtGen[ic] = (TH1F*) fhTrackPt[ic]->Clone(name.Data());
         fOutput->Add((TH1F*) fhTrackPtGen[ic]);
         //------------------------- Rongrong's analysis
         for(Int_t ir=0; ir< kRho; ir++){
 
            name = Form("%sGen",fhDphiTriggerJet[ic][ir]->GetName());
            fhDphiTriggerJetGen[ic][ir] = (TH2F*) fhDphiTriggerJet[ic][ir]->Clone(name.Data());
            fOutput->Add((TH2F*) fhDphiTriggerJetGen[ic][ir]);
         } 
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
Double_t AliAnalysisTaskHJetSpectra::EstimateBgCone(AliJetContainer *jetCont, AliParticleContainer *trkCont, AliVParticle* triggerHadron, Bool_t isGen){
   //Estimate background rho by means of integrating track pT outside TT jet + recoil jet region 
   //if TT exists find jet that containts TT and exclude  range +- phiCut around the TT/TTjet in azimuth

   if(!trkCont) return 0.0; 
   
   AliEmcalJet*  jet        = NULL;
   AliVParticle* track      = NULL;
   Double_t phiTT = fRandom->Rndm()*TMath::TwoPi(); //in case of no TT make random dice 
   Double_t etaTT = -fTrackEtaWindow  + fRandom->Rndm()*2*fTrackEtaWindow;
   Bool_t bTTJetFound = kFALSE;

   if(triggerHadron){

      phiTT = triggerHadron->Phi();
      etaTT = triggerHadron->Eta();

      if(jetCont){
         //find ANY jet that contains TT if it exists
         jetCont->ResetCurrentID();
         while((jet = jetCont->GetNextAcceptJet())){ //loop over reconstructed jets
            if(!jet) continue;
            if(jet->Pt() < triggerHadron->Pt()*0.5) continue;
            //CUT ON JET ACCEPTANCE IS NOT NEEDED, WE SEARCH FOR ANY JET THAT CONTAINS TT
      
            for(Int_t iq=0; iq < jet->GetNumberOfTracks(); iq++){
               track = (AliVParticle*) (jet->TrackAt(iq,trkCont->GetArray())); //matched rec and emb tracks
               if(!track) continue;
               if(track != triggerHadron) continue;
      
               phiTT = jet->Phi(); // used  phi,eta coordinates of the jet to exclude the TT jet
               etaTT = jet->Eta();
               bTTJetFound = kTRUE;
               break;  
            }
            if(bTTJetFound) break; //skip the rest of jets when the jet with TT is found
         }
      } 
   }

   phiTT = RelativePhi(phiTT,0.); //convert phi TT to (-pi,pi)

   if(TMath::Abs(etaTT) > fTrackEtaWindow){
        etaTT = (etaTT<0) ? -fTrackEtaWindow : fTrackEtaWindow;
   }
   //Sum pT outside TT+recoil jet region  
   Double_t sumPt = 0.;

   trkCont->ResetCurrentID();
   while((track = trkCont->GetNextAcceptParticle())){
      if(!track) continue;
      if(!IsTrackInAcceptance(track, isGen)) continue;
      
      if(TMath::Abs(RelativePhi(phiTT+TMath::Pi(),track->Phi())) < fCutPhi)  continue; //exclude recoil region of TT 
      if(GetDeltaR(phiTT, track->Phi(), etaTT, track->Eta()) < fCutPhi)  continue; //exclude region around TT 
 
      if(fTypeOfAnal == kEmb ||  fTypeOfAnal == kEmbSingl){
         if(TMath::Abs(track->GetLabel()) == 99999) continue;//reject embedded stuff 9999 set in AddTaskHJetSpectra.C
      }
 
      sumPt += track->Pt();
   }
  //Calculate area
   Double_t area = 2*fTrackEtaWindow*2*TMath::Pi();
   Double_t alpha; 

   area -= 2*fTrackEtaWindow*2*fCutPhi;  // subtract area of the recoil region

   if(TMath::Abs(etaTT) < fTrackEtaWindow - fCutPhi){ //TT cicle fully inside acceptance
      area -= fCutPhi*fCutPhi*TMath::Pi();  // subtract area of the trigger region 
   }else{  //TT circle partly around acceptance
      alpha = TMath::ACos((fTrackEtaWindow - TMath::Abs(etaTT))/fCutPhi);
      area -= (fCutPhi*fCutPhi*(TMath::Pi()-alpha) + fCutPhi*TMath::Sin(alpha)*(fTrackEtaWindow - TMath::Abs(etaTT)));  
   }

   return  sumPt/area;

}
//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::EstimateBgKT(AliJetContainer *jetCont, AliParticleContainer *trkCont, AliVParticle* triggerHadron){
   //Estimate rho from KT jet median. Ignore jet that contains TT
   Double_t rhoKT    = 0.0;
 
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

      if(triggerHadron){ //identify the KT jet which contains TT 
         if(jet->Pt() > triggerHadron->Pt()*0.5){ //jet containing TT has pT larger than pT of TT
            for(Int_t iq=0; iq < jet->GetNumberOfTracks(); iq++) {
               constTrack = (AliVParticle*) (jet->TrackAt(iq,trkCont->GetArray())); //matched rec and emb tracks
               if(!constTrack) continue;
               if(constTrack == triggerHadron){
                  bKTJetCloseToTT = kTRUE; 
                  break; 
               }
            }
         }
      }
      if(bKTJetCloseToTT) continue; //skip the jet that contains TT 

      sumEmbPt = 0.;//sum pt of embedded tracks in jet which is to be subtracted
      if(fTypeOfAnal == kEmb ||  fTypeOfAnal == kEmbSingl){
         for(Int_t iq=0; iq < jet->GetNumberOfTracks(); iq++) {
            constTrack = (AliVParticle*) (jet->TrackAt(iq,trkCont->GetArray())); 
            if(!constTrack) continue;
            if(TMath::Abs(constTrack->GetLabel()) == 99999){
               sumEmbPt += constTrack->Pt();
            }
         }
      }
 
      jetpt = jet->Pt()- sumEmbPt; //subtract embedded pt
      if(triggerHadron) fhKTAreaPt->Fill(jetpt,jet->Area());

      if(jetpt <0.005) jetpt = 0.; //set pt of ghost jets identical to zero
      frhovec[nJetAcc] = jetpt/jet->Area();
      nJetAcc++;
   }

   if(nJetAcc>0){
      rhoKT = TMath::Median(nJetAcc, frhovec);
   }
 
  return rhoKT; 
}
//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::EstimateBgKTcms(AliJetContainer *jetCont, AliParticleContainer *trkCont, AliVParticle* triggerHadron){
   //Estimate rho from KT jet median ala CMS. Ignore jet that contains TT
   Double_t rhoKTcms    = 0.0;
 
   if(!jetCont) return rhoKTcms;   

   AliEmcalJet*  jet        = NULL;
   AliVParticle* constTrack = NULL;
   Bool_t bKTJetCloseToTT = kFALSE;
   Int_t nJetAcc = 0;
   Double_t areaPhysJets = 0.0;
   Double_t areaAllJets  = 0.0;
   Double_t jetpt;
   Double_t sumEmbPt = 0.;

   jetCont->ResetCurrentID();
   while((jet = jetCont->GetNextAcceptJet())){ //loop over KT jets
      if(!jet) continue;
      if(!IsSignalJetInAcceptance(jet,kFALSE)) continue;

      bKTJetCloseToTT = kFALSE;

      if(triggerHadron){ //identify the KT jet which contains TT 
         if(jet->Pt() > triggerHadron->Pt()*0.5){ //jet containing TT has pT larger than pT of TT
            for(Int_t iq=0; iq < jet->GetNumberOfTracks(); iq++) {
               constTrack = (AliVParticle*) (jet->TrackAt(iq,trkCont->GetArray())); //matched rec and emb tracks
               if(!constTrack) continue;
               if(constTrack != triggerHadron) continue;
               bKTJetCloseToTT = kTRUE; 
               break;  
            }
         }
      }
      if(bKTJetCloseToTT) continue; //skip the jet that contains TT 

      sumEmbPt = 0.;
      if(fTypeOfAnal == kEmb ||  fTypeOfAnal == kEmbSingl){
         for(Int_t iq=0; iq < jet->GetNumberOfTracks(); iq++) {
            constTrack = (AliVParticle*) (jet->TrackAt(iq,trkCont->GetArray())); //matched rec and emb tracks
            if(!constTrack) continue;
            if(TMath::Abs(constTrack->GetLabel()) == 99999){
               sumEmbPt += constTrack->Pt();
            }
         }
      }
 

      areaAllJets += jet->Area();

      jetpt = jet->Pt()- sumEmbPt; //subtract pt of embedded tracks

      if(jetpt > 0.1){
         areaPhysJets += jet->Area();
         frhovec[nJetAcc] = jetpt/jet->Area();
         nJetAcc++;
      }
   }

   if(nJetAcc>0){
      rhoKTcms = TMath::Median(nJetAcc, frhovec)*(areaPhysJets/areaAllJets);
   }
 
  return rhoKTcms; 
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
         pgen = (AliVParticle*) (jGen->TrackAt(ig, jconGen->GetParticleContainer()->GetArray()));
         if(!pgen) continue;
      
         for(Int_t ir=0; ir< jRec->GetNumberOfTracks(); ir++){
            prec = (AliVParticle*) (jRec->TrackAt(ir, jconRec->GetParticleContainer()->GetArray()));
            if(!prec) continue;

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
