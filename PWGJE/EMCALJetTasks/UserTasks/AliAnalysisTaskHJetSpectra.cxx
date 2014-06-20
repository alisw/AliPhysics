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

#include "AliAnalysisTaskHJetSpectra.h"
using std::min;

// ANALYSIS OF HIGH PT HADRON TRIGGER ASSOCIATED SPECTRUM OF RECOIL JETS IN P+PB
// Author Filip Krizek   (17.May. 2014)

//TODO: Not accessing the particles when using MC
//TODO: FillHistogram can be done better with virtual TH1(?)
ClassImp(AliAnalysisTaskHJetSpectra)
//________________________________________________________________________________________

AliAnalysisTaskHJetSpectra::AliAnalysisTaskHJetSpectra(): 
AliAnalysisTaskSE(), fOutputList(0), fAnalyzePythia(0), fAnalyzeHijing(0),  fIsKinematics(0), fUseDefaultVertexCut(1), fUsePileUpCut(1),  
fJetArray(0), fTrackArray(0), fBackgroundJetArray(0), fJetArrayName(0), fTrackArrayName(0), fBackgroundJetArrayName(0),  fRhoTaskName(), 
fRandConeRadius(0.4),fRandConeRadiusSquared(fRandConeRadius*fRandConeRadius), fSignalJetRadius(0.4), fBackgroundJetRadius(0.3), fBackgroundJetPtMin(15.0),
fSignalJetEtaWindow(0.5), fBackgroundJetEtaWindow(0.5), fTrackEtaWindow(0.9), fMinTrackPt(0.150), fMinJetArea(0.5), fNumberOfCentralityBins(20), fCentralityType("V0A"),  
fCrossSection(0.0), fTrials(0.0), fImpParam(-1.0), fRandom(0), fHelperClass(0), fInitialized(0),
fTTlow(8.0), fTThigh(9.0), fTTtype(0), fDphiCut(TMath::Pi()-0.6), fUseDoubleBinPrecision(0),
fHistEvtSelection(0x0), fh2Ntriggers(0x0), fHJetSpec(0x0), fHJetSpecSubUeMedian(0x0), fHJetSpecSubUeCone(0x0), fHJetSpecSubUeCMS(0x0),
fhRhoCellMedian(0x0), fhRhoCone(0x0), fhRhoCMS(0x0), 
fhRhoCellMedianIncl(0x0), fhRhoConeIncl(0x0), fhRhoCMSIncl(0x0), 
fARhoCellMedian(0x0), fARhoCone(0x0), fARhoCMS(0x0), 
fhDeltaPtMedian(0x0), fhDeltaPtCone(0x0), fhDeltaPtCMS(0x0),
fhDeltaPtMedianIncl(0x0), fhDeltaPtConeIncl(0x0), fhDeltaPtCMSIncl(0x0),
fhDeltaPtMedianNearSide(0x0), fhDeltaPtMedianAwaySide(0x0), fhDeltaPtCMSNearSide(0x0), fhDeltaPtCMSAwaySide(0x0),
fhDeltaPtMedianExclTrigCone(0x0),fhDeltaPtCMSExclTrigCone(0x0), fhDeltaPtMedianExclAwayJet(0x0), fhDeltaPtCMSExclAwayJet(0x0),
fhJetPhi(0x0), fhTrackPhi(0x0), fhJetEta(0x0), fhTrackEta(0x0), fhTrackCentVsPt(0x0), fhVertexZ(0x0), fhVertexZAccept(0x0),
fhDphiTriggerJetMinBias(0x0),fhDphiTriggerJetCent20(0x0), fhDphiTriggerJetAccept(0x0),
fhCentrality(0x0), fhCentralityV0M(0x0), fhCentralityV0A(0x0), fhCentralityV0C(0x0), fhCentralityZNA(0x0),
fNofRndTrials(2000), fJetFreeAreaFrac(0.8), fnEta(2), fnPhi(11), fEtaSize(0.9), fPhiSize(2*TMath::Pi()/fnPhi), fCellArea(fPhiSize*fEtaSize),
fh1Xsec(0x0), fh1Trials(0x0), fh1PtHard(0x0), fhImpactParameter(0x0), fhImpactParameterTT(0x0),
fNofRandomCones(1),
fRConesR(0.1),fRConesRSquared(fRConesR*fRConesR),fnRCones(16)
{
    //default constructor
    for(Int_t k=0; k<50; k++){
       fRConePhi[k] = 0.0; 
       fRConeEta[k] = 0.0; 
    } 
}

//________________________________________________________________________
AliAnalysisTaskHJetSpectra::AliAnalysisTaskHJetSpectra(const char *name, const char* trackArrayName, const char* jetArrayName, const char* backgroundJetArrayName) : 
AliAnalysisTaskSE(name), fOutputList(0), fAnalyzePythia(0), fAnalyzeHijing(0), fIsKinematics(0), fUseDefaultVertexCut(1), fUsePileUpCut(1),
 fJetArray(0), fTrackArray(0), fBackgroundJetArray(0), fJetArrayName(0), fTrackArrayName(0), fBackgroundJetArrayName(0), fRhoTaskName(), 
fRandConeRadius(0.4),fRandConeRadiusSquared(fRandConeRadius*fRandConeRadius), fSignalJetRadius(0.4), fBackgroundJetRadius(0.3), fBackgroundJetPtMin(15.0),
fSignalJetEtaWindow(0.5), fBackgroundJetEtaWindow(0.5), fTrackEtaWindow(0.9), fMinTrackPt(0.150), fMinJetArea(0.5),  fNumberOfCentralityBins(20), fCentralityType("V0A"),  
fCrossSection(0.0), fTrials(0.0), fImpParam(-1.0), fRandom(0), fHelperClass(0), fInitialized(0), 
fTTlow(8.0), fTThigh(9.0), fTTtype(0), fDphiCut(TMath::Pi()-0.6), fUseDoubleBinPrecision(0),
fHistEvtSelection(0x0), fh2Ntriggers(0x0), fHJetSpec(0x0), fHJetSpecSubUeMedian(0x0), fHJetSpecSubUeCone(0x0), fHJetSpecSubUeCMS(0x0),
fhRhoCellMedian(0x0), fhRhoCone(0x0), fhRhoCMS(0x0), 
fhRhoCellMedianIncl(0x0), fhRhoConeIncl(0x0), fhRhoCMSIncl(0x0), 
fARhoCellMedian(0x0), fARhoCone(0x0), fARhoCMS(0x0), 
fhDeltaPtMedian(0x0), fhDeltaPtCone(0x0), fhDeltaPtCMS(0x0),
fhDeltaPtMedianIncl(0x0), fhDeltaPtConeIncl(0x0), fhDeltaPtCMSIncl(0x0),
fhDeltaPtMedianNearSide(0x0), fhDeltaPtMedianAwaySide(0x0), fhDeltaPtCMSNearSide(0x0), fhDeltaPtCMSAwaySide(0x0),
fhDeltaPtMedianExclTrigCone(0x0),fhDeltaPtCMSExclTrigCone(0x0), fhDeltaPtMedianExclAwayJet(0x0), fhDeltaPtCMSExclAwayJet(0x0),
fhJetPhi(0x0), fhTrackPhi(0x0), fhJetEta(0x0), fhTrackEta(0x0), fhTrackCentVsPt(0x0), fhVertexZ(0x0), fhVertexZAccept(0x0),
fhDphiTriggerJetMinBias(0x0), fhDphiTriggerJetCent20(0x0), fhDphiTriggerJetAccept(0x0),
fhCentrality(0x0), fhCentralityV0M(0x0), fhCentralityV0A(0x0), fhCentralityV0C(0x0), fhCentralityZNA(0x0),
fNofRndTrials(2000), fJetFreeAreaFrac(0.8), fnEta(2), fnPhi(11), fEtaSize(0.9), fPhiSize(2*TMath::Pi()/fnPhi), fCellArea(fPhiSize*fEtaSize),
fh1Xsec(0x0), fh1Trials(0x0), fh1PtHard(0x0), fhImpactParameter(0x0), fhImpactParameterTT(0x0),
fNofRandomCones(1),
fRConesR(0.1),fRConesRSquared(fRConesR*fRConesR), fnRCones(16)
{
   //constructor that is called 
   //LIST OF TRACKS 
   fTrackArrayName = new TString(trackArrayName);
   if((fTrackArrayName->Contains("MC") && fTrackArrayName->Contains("Particles")) || 
      (fTrackArrayName->Contains("mc") && fTrackArrayName->Contains("particles"))){
      fIsKinematics = kTRUE;
   }

   //LIST of JETS 
   fJetArrayName = new TString(jetArrayName);
   if(strcmp(fJetArrayName->Data(),"") == 0){
      AliError(Form("%s: Jet branch missing !", GetName())); 
   }
     
   //LIST OF JETS TO BE IGNORED WHILE RHO ESTIMATE
   fBackgroundJetArrayName = new TString(backgroundJetArrayName); //jets to be removed from cell median rho estimate
   if(strcmp(fBackgroundJetArrayName->Data(),"") == 0){
      AliError(Form("%s: Bg Jet branch missing !", GetName())); 
   }
 
   for(Int_t k=0; k<50; k++){
      fRConePhi[k] = 0.0; 
      fRConeEta[k] = 0.0; 
   } 

   DefineOutput(1, TList::Class());
}
  
//________________________________________________________________________
void  AliAnalysisTaskHJetSpectra::SetAnalyzeMC(Int_t val){
   if(val==1){
      fAnalyzePythia = kTRUE;
      fAnalyzeHijing = kFALSE;
      return; 
   }
   if(val==2){
      fAnalyzeHijing = kTRUE;
      fAnalyzePythia = kFALSE;
      return;
   }
 
   fAnalyzeHijing = kFALSE;
   fAnalyzePythia = kFALSE;
   return;
}
//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetConePt(Double_t eta, Double_t phi, Double_t radius){
   //sum up pt inside a cone
   Double_t tmpConePt = 0.0;
   Double_t dphi      = 0.0;
   Double_t deta      = 0.0;
   Double_t radiussquared = radius*radius;

   for(Int_t i = 0; i < fTrackArray->GetEntries(); i++){
      AliVTrack* tmpTrack = static_cast<AliVTrack*>(fTrackArray->At(i));
      if(!tmpTrack) continue; 
      if(IsTrackInAcceptance(tmpTrack)){
         dphi = RelativePhi(tmpTrack->Phi(),phi);
         deta = tmpTrack->Eta() - eta;
         if( dphi*dphi + deta*deta < radiussquared ){
            tmpConePt = tmpConePt + tmpTrack->Pt();
         }
      }
   }
   return tmpConePt;
}


//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetPtHard(){
   //get pt hard from pythia
   AliGenPythiaEventHeader* pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
   if(MCEvent()){ 
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
            return kFALSE;
         }
         Float_t xsection = 0.0;
         Float_t trials  = 1.0;

         AliAnalysisHelperJetTasks::PythiaInfoFromFile(curfile->GetName(),xsection,trials);

         fCrossSection = (Double_t) xsection;//pythiaHeader->GetXsection();

         if(fCrossSection>0.){ //save cross-section and the number of trials
            fTrials = (Double_t) trials; //pythiaHeader->Trials();
            fh1Xsec->Fill("<#sigma>", fCrossSection);
            fh1Trials->Fill("#sum{ntrials}",fTrials);
         }
      }
      return pythiaHeader->GetPtHard();
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
   if(MCEvent()){
      if(fAnalyzePythia){ 
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

      if(fAnalyzeHijing){ 
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
   if(mcHeader){
      
      TArrayF pyVtx(3);
      mcHeader->PrimaryVertex(pyVtx);
      return (Double_t) (pyVtx[2]);
   }
   AliWarning(Form("In task %s: Pythia Vertex failed!", GetName()));
   return 9999.0;
}



//________________________________________________________________________
/*Double_t AliAnalysisTaskHJetSpectra::GetPythiaTrials()
{
  #ifdef DEBUGMODE
    AliInfo("Starting GetPythiaTrials.");
  #endif
  AliGenPythiaEventHeader* pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());
  if (MCEvent()) 
    if (!pythiaHeader)
    {
      // Check if AOD
      AliAODMCHeader* aodMCH = dynamic_cast<AliAODMCHeader*>(InputEvent()->FindListObject(AliAODMCHeader::StdBranchName()));

      if (aodMCH)
      {
        for(UInt_t i = 0;i<aodMCH->GetNCocktailHeaders();i++)
        {
          pythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(aodMCH->GetCocktailHeader(i));
          if (pythiaHeader) break;
        }
      }
    }

  #ifdef DEBUGMODE
    AliInfo("Ending GetPythiaTrials.");
  #endif
  if (pythiaHeader)
    return pythiaHeader->Trials();

  AliWarning(Form("In task %s: GetPythiaTrials() failed!", GetName()));
  return -1.0;
}
*/

//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::GetExternalRho(){
   // Get rho from event using CMS approach

   AliRhoParameter *rho = 0;
   if(!fRhoTaskName.IsNull()) {
      rho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRhoTaskName.Data()));
      if(!rho){
         AliWarning(Form("%s: Could not retrieve rho with name %s!", GetName(), fRhoTaskName.Data())); 
         return 0.0;
      }
   }else return 0.0;

   return (rho->GetVal());
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHJetSpectra::IsEventInAcceptance(AliVEvent* event){
   //EVENT SELECTION


   if(!event) return kFALSE;

   //___________________________________________________

   if(fAnalyzePythia || fAnalyzeHijing){ //PURE MC
      if(!MCEvent()) return kFALSE;

       //BEFORE VERTEX CUT
      Double_t vtxMC = GetSimPrimaryVertex();
      fhVertexZ->Fill(vtxMC);

      if(TMath::Abs(vtxMC) > 10.0){
         fHistEvtSelection->Fill(3); //count events rejected by vertex cut 
         return kFALSE;
      }
      fhVertexZAccept->Fill(vtxMC);

      return kTRUE;
   }
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
      if(TMath::Abs(event->GetPrimaryVertex()->GetZ()) > 10.0){
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
Bool_t AliAnalysisTaskHJetSpectra::IsTrackInAcceptance(AliVParticle* track){
   // Check if the track pt and eta range 
   if(track != 0){
      if(fIsKinematics){
         // TODO: Only working for AOD MC
         if((!track->Charge()) || (!(static_cast<AliAODMCParticle*>(track))->IsPhysicalPrimary()) )
            return kFALSE;
      }
      if(TMath::Abs(track->Eta()) <= fTrackEtaWindow){ //APPLY TRACK ETA CUT
         if(track->Pt() >= fMinTrackPt){   //APPLY TRACK CUT
            return kTRUE;
         }
      }
   }
   return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHJetSpectra::IsBackgroundJetInAcceptance(AliEmcalJet *jet){   
   //find jets to be removed from bg calculation 
   if(jet != 0){
      if(TMath::Abs(jet->Eta()) <= fBackgroundJetEtaWindow){
         if(jet->Pt() >= fBackgroundJetPtMin){ //accept only hard jets
            return kTRUE;
         }
      }
   }
   return kFALSE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHJetSpectra::IsSignalJetInAcceptance(AliEmcalJet *jet){   
   //select jets in acceptance 
   if(jet == 0) return kFALSE;
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
void AliAnalysisTaskHJetSpectra::ExecOnce(){
   //Read arrays of jets and tracks
 

   fInitialized = kTRUE; //change flag to skip this function next time when processing UserExec


   fnRCones = TMath::Nint(fRandConeRadiusSquared/fRConesRSquared); //the number of small R=0.1 random cones

   // Check for track array
   if(strcmp(fTrackArrayName->Data(), "") != 0){
      fTrackArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTrackArrayName->Data()));
      if(!fTrackArray){
         AliWarning(Form("%s: Could not retrieve tracks %s!", GetName(), fTrackArrayName->Data())); 
      }else{
         TClass *cl = fTrackArray->GetClass();
         if(!cl->GetBaseClass("AliVParticle")){
      	    AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fTrackArrayName->Data())); 
      	    fTrackArray = 0;
         }
      }
   }

   // Check for jet array
   if(strcmp(fJetArrayName->Data(), "") != 0){
      fJetArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetArrayName->Data()));

      if(!fJetArray){
         AliWarning(Form("%s: Could not retrieve jets %s!", GetName(), fJetArrayName->Data())); 
      }else{
         if(!fJetArray->GetClass()->GetBaseClass("AliEmcalJet")){
            AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fJetArrayName->Data())); 
            fJetArray = 0;
         }
      }
   }

   // Check for list of jets to be removed from background
   if(strcmp(fBackgroundJetArrayName->Data(), "") != 0){
      fBackgroundJetArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fBackgroundJetArrayName->Data()));
      if(!fBackgroundJetArray){
         AliInfo(Form("%s: Could not retrieve background jets %s!", GetName(), fBackgroundJetArrayName->Data())); 
      }else{
         if(!fBackgroundJetArray->GetClass()->GetBaseClass("AliEmcalJet")){
            AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fBackgroundJetArrayName->Data())); 
            fBackgroundJetArray = 0;
         }
      }
   }

   // Look, if initialization is OK
  
   // Initialize helper class (for vertex selection & pile up correction)
   fHelperClass = new AliAnalysisUtils();
   fHelperClass->SetCutOnZVertexSPD(kFALSE); // kFALSE: no cut; kTRUE: |zvtx-SPD - zvtx-TPC|<0.5cm

   return;
}


//________________________________________________________________________
void  AliAnalysisTaskHJetSpectra::GetDeltaPt(Double_t rho1, Double_t &dpt1, Double_t rho2, Double_t &dpt2, 
                                           Double_t rho3, Double_t &dpt3, 
                                           Double_t &rcPhi, Double_t &rcEta,
                                           Double_t leadingJetExclusionProbability){

   //delta pt = random cone - rho

   // Define an invalid delta pt
   dpt1 = -10000.0;
   dpt2 = -10000.0;
   dpt3 = -10000.0;

   // Define eta range
   Double_t etaMin, etaMax;
   etaMin = -(fTrackEtaWindow-fRandConeRadius);
   etaMax = +(fTrackEtaWindow-fRandConeRadius);
 
   // Define random cone Eta+Phi
   Bool_t coneValid = kTRUE;
   Double_t tmpRandConeEta = etaMin + fRandom->Rndm()*(etaMax-etaMin);
   Double_t tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();
 
   // if there is a jet, check for overlap if demanded
   if(leadingJetExclusionProbability){
      AliEmcalJet* tmpLeading = NULL; 
      Double_t lpt = -1.0; 
      // Get leading jet (regardless of pT)
      for(Int_t i = 0; i<fJetArray->GetEntries(); i++){
         AliEmcalJet* tmpJet = static_cast<AliEmcalJet*>(fJetArray->At(i));
         if(!tmpJet) continue;
         if((TMath::Abs(tmpJet->Eta()) <= fSignalJetEtaWindow) && (tmpJet->Area() >= fMinJetArea)){
            if(tmpJet->Pt() > lpt){
               tmpLeading = tmpJet;
               lpt =  tmpJet->Pt();
            }
         }
      }
      if(tmpLeading){
         Double_t excludedJetPhi = tmpLeading->Phi();
         Double_t tmpDeltaPhi    = RelativePhi(tmpRandConePhi, excludedJetPhi);
         Double_t excludedJetEta = tmpLeading->Eta()-tmpRandConeEta;
 
         // Check, if cone has overlap with jet
         if(tmpDeltaPhi*tmpDeltaPhi + excludedJetEta*excludedJetEta <= fRandConeRadiusSquared){
            // Define probability to exclude the RC
            Double_t probability = leadingJetExclusionProbability;
 
            // Only exclude cone with a given probability
            if(fRandom->Rndm()<=probability)  coneValid = kFALSE;
         }
      }
   }
 
   rcPhi = 9999.0; 
   rcEta = 9999.0; 
 
   // Get the cones' pt and calculate delta pt
   if(coneValid){
      rcPhi = tmpRandConePhi;
      rcEta = tmpRandConeEta;
      Double_t conePt = GetConePt(tmpRandConeEta,tmpRandConePhi,fRandConeRadius);
      dpt1 =  conePt - (rho1*fRandConeRadiusSquared*TMath::Pi());
      dpt2 =  conePt - (rho2*fRandConeRadiusSquared*TMath::Pi());
      dpt3 =  conePt - (rho3*fRandConeRadiusSquared*TMath::Pi());
   }
 
}


//________________________________________________________________________
void AliAnalysisTaskHJetSpectra::Calculate(AliVEvent* event){
   //Analyze the event and Fill histograms

   if(fAnalyzePythia){
      fh1PtHard->Fill(GetPtHard());
   }

   if(fAnalyzeHijing){
      fImpParam = GetImpactParameter(); 
      fhImpactParameter->Fill(fImpParam);
   }
   //_________________________________________________________________
   //  FILL EVENT STATISTICS
   fHistEvtSelection->Fill(1); //Count input event

   if(!IsEventInAcceptance(event)) return; //post data is in UserExec
   

   // Get centrality
   AliCentrality* tmpCentrality = event->GetCentrality();
   if(!tmpCentrality){
      fHistEvtSelection->Fill(4);
      return; //post data is in UserExec
   }
   Double_t centralityPercentile    = -1.0;
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

   if((centralityPercentile < 0.0) || (centralityPercentile > 100.0)){
      AliWarning(Form("Centrality value not valid (c=%E)",centralityPercentile)); 
      fHistEvtSelection->Fill(4);
      return;
   }
   fhCentrality->Fill(centralityPercentile);
   fhCentralityV0M->Fill(centralityPercentileV0M); 
   fhCentralityV0A->Fill(centralityPercentileV0A);
   fhCentralityV0C->Fill(centralityPercentileV0C); 
   fhCentralityZNA->Fill(centralityPercentileZNA);
 
   fHistEvtSelection->Fill(0); //Count input event

   // END EVENT SELECTION
   //___________________________________________________________

   //LOOP OVER TRACKS  SEARCH FOR TRIGGER
   std::vector<Int_t> trigTracks; //list pf trigger particle indices
   //Bool_t bContainesHighPtTrack = kFALSE;

   Int_t nTracks = fTrackArray->GetEntries();

   for(Int_t i = 0; i < nTracks; i++){
      AliVTrack* track = static_cast<AliVTrack*>(fTrackArray->At(i));

      if(!track) continue;

      if(IsTrackInAcceptance(track)){
         fhTrackPhi->Fill(track->Pt(), RelativePhi(track->Phi(),0.0)); // phi = -pi az pi
         fhTrackEta->Fill(track->Pt(), track->Eta());
         fhTrackCentVsPt->Fill(track->Pt(), centralityPercentile);
                
         if(fTTlow <= track->Pt() && track->Pt() < fTThigh){
            trigTracks.push_back(i);  //trigger candidates
         }

         //if(track->Pt()>=8.0) bContainesHighPtTrack = kTRUE;
      }
   }

   Int_t ntriggers = (Int_t) trigTracks.size();
   Int_t indexSingleRndTrig = -1; //index of single random trigger
   Double_t areaJet,  pTJet; 
   Double_t tmpArray[3];              
   Double_t rhoFromCellMedian = 0.0; //UE density cell median
   Double_t rhoCone           = 0.0; //UE density perp cone
   Double_t rhoCMS            = 0.0; //UE density ala CMS
   Double_t deltaptCellMedian, deltaptCone, deltaptCMS, randConePhi, randConeEta;       
   Double_t distanceFromTrigger; 

   if(ntriggers>0){
      if(fTTtype==0){ //select single inclusive trigger
         indexSingleRndTrig = fRandom->Integer(ntriggers); //Integer 0 ... ntriggers-1
      }
   }

   rhoFromCellMedian = EstimateBgRhoMedian();
   rhoCone           = EstimateBgCone();
   rhoCMS            = GetExternalRho();

   fhRhoCellMedianIncl->Fill((Float_t) rhoFromCellMedian,(Float_t) centralityPercentile);
   fhRhoConeIncl->Fill(      (Float_t) rhoCone,          (Float_t) centralityPercentile); 
   fhRhoCMSIncl->Fill(       (Float_t) rhoCMS,           (Float_t) centralityPercentile); 

   for(Int_t irc=0; irc<fNofRandomCones; irc++){ //generate 4 random cones per event
      GetDeltaPt(rhoFromCellMedian, deltaptCellMedian,rhoCone, deltaptCone, rhoCMS, deltaptCMS, randConePhi, randConeEta, 0);
   
      fhDeltaPtMedianIncl->Fill(deltaptCellMedian, (Double_t) centralityPercentile); 
      fhDeltaPtConeIncl->Fill( deltaptCone,        (Double_t) centralityPercentile); 
      fhDeltaPtCMSIncl->Fill( deltaptCMS,          (Double_t) centralityPercentile); 
  
      if(ntriggers>0){
         //fill delta pt histograms near side + away side
         fhDeltaPtMedian->Fill( deltaptCellMedian, (Double_t) centralityPercentile); 
         fhDeltaPtCone->Fill( deltaptCone,         (Double_t) centralityPercentile); 
         fhDeltaPtCMS->Fill(  deltaptCMS,          (Double_t) centralityPercentile);

         if(indexSingleRndTrig>-1){
            AliVTrack* triggHad = static_cast<AliVTrack*>(fTrackArray->At(trigTracks[indexSingleRndTrig]));
            Double_t dphiTrigRC =  RelativePhi(triggHad->Phi(), randConePhi); 
            Double_t detaTrigRC =  triggHad->Eta()- randConeEta; 
            if(TMath::Abs(dphiTrigRC)< TMath::Pi()/2){ //near side
               fhDeltaPtMedianNearSide->Fill( deltaptCellMedian, (Double_t) centralityPercentile); 
               fhDeltaPtCMSNearSide->Fill(  deltaptCMS,          (Double_t) centralityPercentile);
            }else{ //away side
               fhDeltaPtMedianAwaySide->Fill( deltaptCellMedian, (Double_t) centralityPercentile); 
               fhDeltaPtCMSAwaySide->Fill(  deltaptCMS,          (Double_t) centralityPercentile);
            }

            distanceFromTrigger = sqrt(dphiTrigRC*dphiTrigRC+detaTrigRC*detaTrigRC);
            while(distanceFromTrigger<0.5 + fRandConeRadius){
               GetDeltaPt(rhoFromCellMedian, deltaptCellMedian,rhoCone, deltaptCone, rhoCMS, deltaptCMS, randConePhi, randConeEta, 0);
               dphiTrigRC =  RelativePhi(triggHad->Phi(), randConePhi); 
               detaTrigRC =  triggHad->Eta()- randConeEta; 
               distanceFromTrigger = sqrt(dphiTrigRC*dphiTrigRC+detaTrigRC*detaTrigRC);
            }
            if(distanceFromTrigger>0.5 + fRandConeRadius){
               fhDeltaPtMedianExclTrigCone->Fill( deltaptCellMedian, (Double_t) centralityPercentile); 
               fhDeltaPtCMSExclTrigCone->Fill(  deltaptCMS,          (Double_t) centralityPercentile);
            }
         } 
      }
   }
   
   //_______________________________________
   Int_t    idxLeadingJetAwaySide = -1;
   Double_t ptLeadingJetAwaySide  = -1.0;
   Double_t phiTrigger            = -1.0; //-pi,pi

   if(ntriggers>0){
      //Estimate UE density
      //Fill once per event
      fhRhoCellMedian->Fill((Float_t) rhoFromCellMedian,(Float_t) centralityPercentile);
      fhRhoCone->Fill(      (Float_t) rhoCone,          (Float_t) centralityPercentile); 
      fhRhoCMS->Fill(       (Float_t) rhoCMS,           (Float_t) centralityPercentile); 


      //TRIGGER PARTICLE LOOP 
      for(Int_t it=0; it<ntriggers; it++){ //loop over trigger configurations
     
         if(fTTtype==0){
            if(it != indexSingleRndTrig) continue;
         }

         AliVTrack* triggerHadron = static_cast<AliVTrack*>(fTrackArray->At(trigTracks[it]));
         if(!triggerHadron) continue;
 
         fh2Ntriggers->Fill((Float_t) centralityPercentile, (Float_t) triggerHadron->Pt()); //trigger p 

         if(fAnalyzeHijing){ //impact parameter for triggered events
            fhImpactParameterTT->Fill(fImpParam);
         }

         phiTrigger = RelativePhi(triggerHadron->Phi(),0.0); //-pi,pi

         //JET LOOP
         for(Int_t ij = 0; ij < fJetArray->GetEntries(); ij++){
            AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJetArray->At(ij));
            if(!jet){
               AliError(Form("%s: Could not receive jet %d", GetName(), ij));
               continue;
            }
            if(!IsSignalJetInAcceptance(jet)) continue;

            areaJet = jet->Area();
            pTJet   = jet->Pt();

            if(it==0 || it == indexSingleRndTrig){
               fhJetPhi->Fill( pTJet, RelativePhi(jet->Phi(),0.0));
               fhJetEta->Fill( pTJet, jet->Eta());
            }

            Double_t dphi = RelativePhi(triggerHadron->Phi(), jet->Phi());
       
            Double_t dfi = dphi; //-0.5*pi to 1.5*Pi
            if(dfi<-0.5*TMath::Pi()) dfi += 2*TMath::Pi();
            if(dfi> 1.5*TMath::Pi()) dfi -= 2*TMath::Pi();
            fhDphiTriggerJetMinBias->Fill((Float_t) jet->Pt(),(Float_t) dfi); 
            if(centralityPercentile<20.) fhDphiTriggerJetCent20->Fill((Float_t) jet->Pt(),(Float_t) dfi); 
            //-------------------------
 
            if(TMath::Abs(dphi) < fDphiCut) continue;  //Dphi cut between trigger and assoc
            fhDphiTriggerJetAccept->Fill(dfi); //Accepted

            if(pTJet > ptLeadingJetAwaySide){ //search for the leading away side jet
               idxLeadingJetAwaySide = ij; 
               ptLeadingJetAwaySide  = pTJet;
            }

           //Centrality, A, pTjet
            tmpArray[0] =  centralityPercentile;
            tmpArray[1] =  areaJet; 
            tmpArray[2] =  pTJet;
            fHJetSpec->Fill(tmpArray);

            //Subtract cell median
            tmpArray[2] = pTJet - areaJet*rhoFromCellMedian;
            fHJetSpecSubUeMedian->Fill(tmpArray);
            fARhoCellMedian->Fill((Float_t) (areaJet*rhoFromCellMedian));

            //Subtract perp cone 
            tmpArray[2] = pTJet - areaJet*rhoCone;
            fHJetSpecSubUeCone->Fill(tmpArray);
            fARhoCone->Fill((Float_t) (areaJet*rhoCone));

            //Subtract CMS bg 
            tmpArray[2] = pTJet - areaJet*rhoCMS;
            fHJetSpecSubUeCMS->Fill(tmpArray);
            fARhoCMS->Fill((Float_t) (areaJet*rhoCMS));

         }//JET LOOP
      }
   }

   //_______________________________________
   // Get delta phi from small R=0.1 cones
   if(ntriggers>0 && indexSingleRndTrig>-1){

      AliEmcalJet* jet = NULL;
      Double_t phiExclJet =0., etaExclJet = 9999., rExclJet = 0.0;
 
      if(idxLeadingJetAwaySide>-1){
         jet = static_cast<AliEmcalJet*>(fJetArray->At(idxLeadingJetAwaySide));
         if(!jet){
            AliError(Form("%s: Could not receive leading jet %d", GetName(), idxLeadingJetAwaySide));
         }else{
             phiExclJet = jet->Phi();
             etaExclJet = jet->Eta();
             rExclJet   = TMath::Sqrt(jet->Area()/TMath::Pi());
         }
      } 

      Int_t countPlacedRcones = 0;
      Int_t inwhile=0;
      Double_t dphiMaxFromPiForRcones = TMath::Pi() - fDphiCut + fRandConeRadius - fRConesR;
      Double_t detaMaxForRcones = fTrackEtaWindow - fRConesR;
      Double_t rcphi,rceta;
      Bool_t goodrc;

      while(countPlacedRcones < fnRCones){ //generate fnRCones cones of radius R=0.1
         inwhile++;
         if(inwhile>500){
            AliError(Form("%s: Small space where to put another random cone %d", GetName(), idxLeadingJetAwaySide));
            break;
         }
         rcphi = RelativePhi(phiTrigger + TMath::Pi() + dphiMaxFromPiForRcones*fRandom->Uniform(-1.0,1.0), 0.0); //-pi,pi
         rceta = detaMaxForRcones*fRandom->Uniform(-1.0,1.0);
         if(jet){//do not merge random cones with the leading jet  
            if(!DistantCones(rcphi,rceta,fRConesR, phiExclJet, etaExclJet, rExclJet)) continue;
         }

         goodrc = kTRUE; //generate disjoint random cones
         for(int k=0; k<countPlacedRcones; k++){
            if(!DistantCones(rcphi, rceta, fRConesR, (Double_t) fRConePhi[k], (Double_t) fRConeEta[k], fRConesR)){
               goodrc = kFALSE;  
               break;
            }
         }//end loop over already placed cones 

         if(!goodrc) continue;
         fRConePhi[countPlacedRcones] = rcphi;
         fRConeEta[countPlacedRcones] = rceta;
         countPlacedRcones++;
      }//end of  loop generating small R=0.1 random cones
      //sum track pT in the random cones
      Double_t sumPtofRandomCones = 0.0;
 
      for(Int_t itrk = 0; itrk < nTracks; itrk++){
         AliVTrack* track = static_cast<AliVTrack*>(fTrackArray->At(itrk));
         if(!track) continue;

         if(IsTrackInAcceptance(track)){
            for(int k=0; k<countPlacedRcones; k++){
               
               rcphi = RelativePhi(track->Phi(), fRConePhi[k]); // phi = -pi az pi
               rceta = track->Eta() - fRConeEta[k];
               if(rcphi*rcphi + rceta*rceta < fRConesRSquared){
                  sumPtofRandomCones += track->Pt(); //track is in the cone
                  break;
               } 
            }//loop over cones
         }
      }//loop over tracks
      Double_t totarea = countPlacedRcones*TMath::Pi()*fRConesRSquared;
      fhDeltaPtMedianExclAwayJet->Fill( sumPtofRandomCones - rhoFromCellMedian*totarea, (Double_t) centralityPercentile );
      fhDeltaPtCMSExclAwayJet->Fill(    sumPtofRandomCones - rhoCMS*totarea , (Double_t) centralityPercentile );

   }//end delta phi from small R=0.1 random cones

   return;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskHJetSpectra::UserNotify(){
   // Implemented Notify() to read the cross sections
   // and number of trials from pyxsec.root
   /*
   if(fAnalyzePythia){
   
     TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
     TFile *currFile = tree->GetCurrentFile();
 
     TString file(currFile->GetName());
 
     if(file.Contains("root_archive.zip#")){
       Ssiz_t pos1 = file.Index("root_archive",12,TString::kExact);
       Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
       file.Replace(pos+1,20,"");
     }
     else {
       // not an archive take the basename....
       file.ReplaceAll(gSystem->BaseName(file.Data()),"");
     }
    
     TFile *fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec.root")); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
     if(!fxsec){
       // next trial fetch the histgram file
       fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
       if(!fxsec){
           // not a severe condition but inciate that we have no information
         return kFALSE;
       }
       else{
         // find the tlist we want to be independtent of the name so use the Tkey
         TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0); 
         if(!key){
           fxsec->Close();
           return kFALSE;
         }
         TList *list = dynamic_cast<TList*>(key->ReadObj());
         if(!list){
           fxsec->Close();
           return kFALSE;
         }
         fCrossSection = ((TProfile*)list->FindObject("h1Xsec"))->GetBinContent(1);
         fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
         fxsec->Close();
       }
     } // no tree pyxsec.root
     else {
       TTree *xtree = (TTree*)fxsec->Get("Xsection");
       if(!xtree){
         fxsec->Close();
         return kFALSE;
       }
       UInt_t   ntrials  = 0;
       Double_t  xsection  = 0;
       xtree->SetBranchAddress("xsection",&xsection);
       xtree->SetBranchAddress("ntrials",&ntrials);
       xtree->GetEntry(0);
       fTrials = ntrials;
       fCrossSection = xsection;
       fxsec->Close();
     }


     fh1Xsec->Fill("<#sigma>", fCrossSection);
     fh1Trials->Fill("#sum{ntrials}",fTrials);
     
   }
   */
   return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskHJetSpectra::Terminate(Option_t *){
   //Treminate 
   PostData(1, fOutputList);

   // Mandatory
   fOutputList = dynamic_cast<TList*> (GetOutputData(1)); // '1' refers to the output slot
   if(!fOutputList) {
      printf("ERROR: Output list not available\n");
      return;
   }
}

//________________________________________________________________________
AliAnalysisTaskHJetSpectra::~AliAnalysisTaskHJetSpectra(){
   // Destructor. Clean-up the output list, but not the histograms that are put inside
   // (the list is owner and will clean-up these histograms). Protect in PROOF case.
   if(fOutputList && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutputList;
   }
   delete fRandom;
   delete fTrackArrayName;
   delete fJetArrayName;
   delete fBackgroundJetArrayName;
   delete fHelperClass;
 
} 

//________________________________________________________________________
void AliAnalysisTaskHJetSpectra::UserCreateOutputObjects(){
  // called once to create user defined output objects like histograms, plots etc. 
  // and to put it on the output list.
  // Note: Saving to file with e.g. OpenFile(0) is must be before creating other objects.

   fRandom = new TRandom3(0);

   fOutputList = new TList();
   fOutputList->SetOwner(); // otherwise it produces leaks in merging
   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);

   //__________________________________________________________
   // Event statistics
   fHistEvtSelection = new TH1I("fHistEvtSelection", "event selection", 6, -0.5, 5.5);
   fHistEvtSelection->GetXaxis()->SetBinLabel(1,"ACCEPTED");
   fHistEvtSelection->GetXaxis()->SetBinLabel(2,"events IN");
   fHistEvtSelection->GetXaxis()->SetBinLabel(3,"pile up (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(4,"vertex cut (rejected)");
   fHistEvtSelection->GetXaxis()->SetBinLabel(5,"centrality (rejected)");

   fOutputList->Add(fHistEvtSelection);
   //___________________________________________________________
   // Hard trigger counter
   fh2Ntriggers = new TH2F("fh2Ntriggers","# of triggers",
                            fNumberOfCentralityBins,0.0,100.0,50,0.0,50.0);
   fOutputList->Add(fh2Ntriggers);
   //___________________________________________________________
   // trigger associated jet spectra (jet pT not corrected for UE)
   Int_t bw = (fUseDoubleBinPrecision==0) ? 1 : 2; //make larger bin width

   //jet associated to given TT 
   //Centrality, A, pTjet  
   const Int_t    dimSpec   = 3;
   const Int_t    nBinsSpec[dimSpec]  = {fNumberOfCentralityBins,  50, bw*110};
   const Double_t lowBinSpec[dimSpec] = {0.0,                     0.0,  -20.0};
   const Double_t hiBinSpec[dimSpec]  = {100.0,                   1.5,  200.0};
   fHJetSpec = new THnSparseF("fHJetSpec",
                   "Recoil jet spectrum [cent,A,pTjet]",
                   dimSpec,nBinsSpec,lowBinSpec,hiBinSpec);
   fOutputList->Add(fHJetSpec);
   //___________________________________________________________
   //jet associated to given TT (jet pT corrected with rho from cell median)
   fHJetSpecSubUeMedian = (THnSparseF*) fHJetSpec->Clone("fHJetSpecSubUeMedian");
   fHJetSpecSubUeMedian->SetTitle("Recoil jet spectrum [cent,A,pTjet-pTUe]");
   fOutputList->Add(fHJetSpecSubUeMedian);
   //___________________________________________________________
   //jet associated to given TT (jet pT corrected with rho from perp cone)
   fHJetSpecSubUeCone = (THnSparseF*) fHJetSpec->Clone("fHJetSpecSubUeCone");
   fHJetSpecSubUeCone->SetTitle("Recoil jet spectrum [cent,A,pTjet-pTUe]");
   fOutputList->Add(fHJetSpecSubUeCone);
   //___________________________________________________________
   //jet associated to given TT (jet pT corrected with rho from CMS approach)
   fHJetSpecSubUeCMS = (THnSparseF*) fHJetSpec->Clone("fHJetSpecSubUeCMS");
   fHJetSpecSubUeCMS->SetTitle("Recoil jet spectrum [cent,A,pTjet-pTUe]");
   fOutputList->Add(fHJetSpecSubUeCMS);

   //____________________________________________________________________
   //UE from cell median  [Centrality, rho, pTUe ]

   fhRhoCellMedian = new TH2F("fhRhoCellMedian","Rho",40, 0.0, 20.0, fNumberOfCentralityBins, 0.0, 100.);
   fOutputList->Add(fhRhoCellMedian);

   fhRhoCone = (TH2F*) fhRhoCellMedian->Clone("fhRhoCone");
   fOutputList->Add(fhRhoCone);

   fhRhoCMS = (TH2F*) fhRhoCellMedian->Clone("fhRhoCMS");
   fOutputList->Add(fhRhoCMS);

   fhRhoCellMedianIncl = (TH2F*) fhRhoCellMedian->Clone("fhRhoCellMedianIncl");
   fOutputList->Add(fhRhoCellMedianIncl);

   fhRhoConeIncl = (TH2F*) fhRhoCellMedian->Clone("fhRhoConeIncl");
   fOutputList->Add(fhRhoConeIncl);

   fhRhoCMSIncl = (TH2F*) fhRhoCellMedian->Clone("fhRhoCMSIncl");
   fOutputList->Add(fhRhoCMSIncl);
 
   //_______________________________________________________________________
   // rho times area 
   fARhoCellMedian = new TH1F("fARhoCellMedian","Area times rho",40, 0.0, 20.0);
   fOutputList->Add(fARhoCellMedian);

   fARhoCone = (TH1F*) fARhoCellMedian->Clone("fARhoCone");
   fOutputList->Add(fARhoCone);

   fARhoCMS  = (TH1F*) fARhoCellMedian->Clone("fARhoCMS");
   fOutputList->Add(fARhoCMS);

   //_______________________________________________________________________
   // Delta pt distributions   
   fhDeltaPtMedian = new TH2D("fhDeltaPtMedian","DeltaPt", bw*110, -20, 200, fNumberOfCentralityBins,0.0,100.0);
   fOutputList->Add(fhDeltaPtMedian);

   fhDeltaPtCone = (TH2D*) fhDeltaPtMedian->Clone("fhDeltaPtCone");
   fOutputList->Add(fhDeltaPtCone);

   fhDeltaPtCMS = (TH2D*) fhDeltaPtMedian->Clone("fhDeltaPtCMS");
   fOutputList->Add(fhDeltaPtCMS);

   fhDeltaPtMedianIncl = (TH2D*) fhDeltaPtMedian->Clone("fhDeltaPtMedianIncl");
   fOutputList->Add(fhDeltaPtMedianIncl);
 
   fhDeltaPtConeIncl = (TH2D*) fhDeltaPtMedian->Clone("fhDeltaPtConeIncl");
   fOutputList->Add(fhDeltaPtConeIncl);

   fhDeltaPtCMSIncl = (TH2D*) fhDeltaPtMedian->Clone("fhDeltaPtCMSIncl");
   fOutputList->Add(fhDeltaPtCMSIncl);

   fhDeltaPtMedianNearSide= (TH2D*) fhDeltaPtMedian->Clone("fhDeltaPtMedianNearSide");
   fOutputList->Add(fhDeltaPtMedianNearSide);

   fhDeltaPtMedianAwaySide= (TH2D*) fhDeltaPtMedian->Clone("fhDeltaPtMedianAwaySide");
   fOutputList->Add(fhDeltaPtMedianAwaySide);

   fhDeltaPtCMSNearSide= (TH2D*) fhDeltaPtMedian->Clone("fhDeltaPtCMSNearSide");
   fOutputList->Add(fhDeltaPtCMSNearSide);

   fhDeltaPtCMSAwaySide= (TH2D*) fhDeltaPtMedian->Clone("fhDeltaPtCMSAwaySide");
   fOutputList->Add(fhDeltaPtCMSAwaySide);

   fhDeltaPtMedianExclTrigCone= (TH2D*) fhDeltaPtMedian->Clone("fhDeltaPtMedianExclTrigCone");
   fOutputList->Add(fhDeltaPtMedianExclTrigCone);

   fhDeltaPtCMSExclTrigCone= (TH2D*) fhDeltaPtMedian->Clone("fhDeltaPtCMSExclTrigCone");
   fOutputList->Add(fhDeltaPtCMSExclTrigCone);

   fhDeltaPtMedianExclAwayJet = (TH2D*) fhDeltaPtMedian->Clone("fhDeltaPtMedianExclAwayJet");
   fOutputList->Add(fhDeltaPtMedianExclAwayJet);

   fhDeltaPtCMSExclAwayJet = (TH2D*) fhDeltaPtMedian->Clone("fhDeltaPtCMSExclAwayJet");
   fOutputList->Add(fhDeltaPtCMSExclAwayJet);

   //_______________________________________________________________________
   //inclusive azimuthal and pseudorapidity histograms
   fhJetPhi   = new TH2F("fhJetPhi","Azim dist jets vs pTjet", 50, 0, 100, 50,-TMath::Pi(),TMath::Pi());
   fOutputList->Add(fhJetPhi);
   //-------------------------
   fhTrackPhi = new TH2F("fhTrackPhi","azim dist trig had vs pT,trk", 50, 0, 50, 50,-TMath::Pi(),TMath::Pi());
   fOutputList->Add(fhTrackPhi);
   //-------------------------
   fhJetEta   = new TH2F("fhJetEta","Eta dist jets vs pTjet", 50,0, 100, 40,-0.9,0.9);
   fOutputList->Add(fhJetEta);
   //-------------------------
   fhTrackEta = new TH2F("fhTrackEta","Eta dist trig had vs pT,trk", 50, 0, 50, 40,-0.9,0.9);
   fOutputList->Add(fhTrackEta);
   //-------------------------
   fhTrackCentVsPt = new TH2F("fhTrackCentVsPt","pT,trk vs centrality", 50, 0, 50, fNumberOfCentralityBins,0,100);
   fOutputList->Add(fhTrackCentVsPt);
   //-------------------------
   fhVertexZ = new TH1F("fhVertexZ","z vertex",40,-20,20);
   fOutputList->Add(fhVertexZ);
   //-------------------------
   fhVertexZAccept = new TH1F("fhVertexZAccept","z vertex after cut",40,-20,20);
   fOutputList->Add(fhVertexZAccept);
   //-------------------------
   //fhContribVtx = new TH1F("fhContribVtx","contrib to vtx",200,0,200);
   //fOutputList->Add(fhContribVtx);
   //-------------------------
   //fhContribVtxAccept = new TH1F("fhContribVtxAccept","contrib to vtx after cut",200,0,200);
   //fOutputList->Add(fhContribVtxAccept);
   //-------------------------
   fhDphiTriggerJetMinBias = new TH2F("fhDphiTriggerJetMinBias","Deltaphi trig-jet",50,0,100, 100, -0.5*TMath::Pi(),1.5*TMath::Pi());
   fOutputList->Add(fhDphiTriggerJetMinBias);

   fhDphiTriggerJetCent20 = (TH2F*) fhDphiTriggerJetMinBias->Clone("fhDphiTriggerJetCent20");
   fOutputList->Add(fhDphiTriggerJetCent20);
   //-------------------------

   fhDphiTriggerJetAccept = new TH1F("fhDphiTriggerJetAccept","Deltaphi trig-jet after cut",50, -0.5*TMath::Pi(),1.5*TMath::Pi());
   fOutputList->Add(fhDphiTriggerJetAccept);
   //-------------------------
   fhCentrality = new TH1F("fhCentrality","Centrality",100,0,100);
   fOutputList->Add(fhCentrality);
   //-------------------------
   fhCentralityV0M = new TH1F("hCentralityV0M","hCentralityV0M",100,0,100);
   fOutputList->Add(fhCentralityV0M); 
   //-------------------------
   fhCentralityV0A = new TH1F("hCentralityV0A","hCentralityV0A",100,0,100);
   fOutputList->Add(fhCentralityV0A); 
   //-------------------------
   fhCentralityV0C = new TH1F("hCentralityV0C","hCentralityV0C",100,0,100);
   fOutputList->Add(fhCentralityV0C);
   //-------------------------
   fhCentralityZNA = new TH1F("hCentralityZNA","hCentralityZNA",100,0,100);
   fOutputList->Add(fhCentralityZNA);

   fh1Xsec = new TProfile("fh1Xsec","xsec from pyxsec.root",1,0,1);
   fh1Xsec->GetXaxis()->SetBinLabel(1,"<#sigma>");
   fOutputList->Add(fh1Xsec);
 
   fh1Trials = new TH1F("fh1Trials","trials root file",1,0,1);
   fh1Trials->GetXaxis()->SetBinLabel(1,"#sum{ntrials}");
   fOutputList->Add(fh1Trials);

   fh1PtHard = new TH1F("fh1PtHard","PYTHIA Pt hard;p_{T,hard}",300,0,300);
   fOutputList->Add(fh1PtHard);

   fhImpactParameter = new TH1D("fhImpactParameter","impact parameter distribution from HIJING",50,0,10);
   fOutputList->Add(fhImpactParameter);

   fhImpactParameterTT = new TH1D("fhImpactParameterTT","b versus TT",50,0,10);
   fOutputList->Add(fhImpactParameterTT);
   // =========== Switch on Sumw2 for all histos ===========
   for(Int_t i=0; i<fOutputList->GetEntries(); i++){
      TH1 *h1 = dynamic_cast<TH1*>(fOutputList->At(i));
      if(h1){
         h1->Sumw2();
         continue;
      }
      THnSparse *hn = dynamic_cast<THnSparse*>(fOutputList->At(i));
      if(hn){
         hn->Sumw2();
      }
   }
   TH1::AddDirectory(oldStatus);


   PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskHJetSpectra::UserExec(Option_t *){
   //executed in each event 

   if(!InputEvent()){
      AliError("??? Event pointer == 0 ???");
      return;
   }

   //Execute only once:  Get tracks, jets, background from arrays if not already given 
   if(!fInitialized) ExecOnce(); 
   if(fJetArray && fTrackArray && fBackgroundJetArray){ 
      Calculate(InputEvent());
   }
   PostData(1, fOutputList);
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
Double_t AliAnalysisTaskHJetSpectra::EstimateBgRhoMedian(){
   //Estimate background rho by means of integrating track pT outside identified jet cones
   Double_t rhoMedian = 0.0;

   //phi,eta and R2 of jets to be removed
   std::vector<Double_t>  jphi;
   std::vector<Double_t>  jeta;
   std::vector<Double_t> jRsquared;

   if(!fBackgroundJetArray) return 0.0;
   if(!fTrackArray)         return 0.0;

   for(Int_t i = 0; i < fBackgroundJetArray->GetEntries(); i++){
      AliEmcalJet* backgroundJet = static_cast<AliEmcalJet*>(fBackgroundJetArray->At(i));

      if(!backgroundJet){
         AliError(Form("%s: Could not receive jet %d", GetName(), i));
         continue;
      }
      if(!IsBackgroundJetInAcceptance(backgroundJet)) continue; //apply minimum pT cut on jet to be removed from bg
      jphi.push_back(RelativePhi(backgroundJet->Phi(),0.0)); //-pi,pi
      jeta.push_back(backgroundJet->Eta());
      jRsquared.push_back(1.78*backgroundJet->Area()/TMath::Pi()); //1.78 = JetArea_R04/JetArea_R03
   }


   static Double_t nOutCone[10][4];
   static Double_t sumPtOutOfCone[10][4];
   Double_t rndphi, rndeta;
   Double_t rndphishift, rndetashift;
   Double_t dphi, deta;
   Bool_t   bIsInCone;


   for(Int_t ie=0; ie < fnEta; ie++){
      for(Int_t ip=0; ip < fnPhi; ip++){
         nOutCone[ip][ie]       = 0.0;     //initialize counter
         sumPtOutOfCone[ip][ie] = 0.0;
      }
   }

   //get area in cells out of identified jet cones
   if(jphi.size()==0){ //no jet to be removed from the bg => all areas have their nominal area
      for(Int_t ie=0; ie < fnEta; ie++){
         for(Int_t ip=0; ip < fnPhi; ip++){
            nOutCone[ip][ie] = fNofRndTrials; 
         }
      } 
   }else{
      for(Int_t it=0; it<fNofRndTrials; it++){

         rndphi = fRandom->Uniform(0, fPhiSize);
         rndeta = fRandom->Uniform(0, fEtaSize);

         for(Int_t ip=0; ip<fnPhi; ip++){  //move radom position to each cell
            rndphishift = rndphi + ip*fPhiSize - TMath::Pi();
            for(Int_t ie=0; ie<fnEta; ie++){
               rndetashift = rndeta + ie*fEtaSize - fEtaSize;

               bIsInCone = 0; //tag if trial is in the jet cone
               for(Int_t ij=0; ij< (Int_t) jRsquared.size(); ij++){
                  deta = jeta[ij] - rndetashift;
                  dphi = RelativePhi(rndphishift,jphi[ij]);
                  if((dphi*dphi + deta*deta) < jRsquared[ij]){
                     bIsInCone = 1;
                     break;
                  }
               }
               if(!bIsInCone) nOutCone[ip][ie]++;
            }
         }
      }
   }

   Int_t phicell,etacell;
   for(Int_t ip=0; ip < fTrackArray->GetEntries(); ip++){
      AliVTrack* part = static_cast<AliVTrack*>(fTrackArray->At(ip));
      if(!part) continue;
      if(!IsTrackInAcceptance((AliVParticle*) part)) continue; 

      bIsInCone = 0; //init
      for(Int_t ij=0; ij<(Int_t) jRsquared.size(); ij++){
         dphi = RelativePhi(jphi[ij], part->Phi());
         deta = jeta[ij] - part->Eta();
         if((dphi*dphi + deta*deta) < jRsquared[ij]){
            bIsInCone = 1;
            break;
         }
      }
      if(!bIsInCone){
         phicell = TMath::Nint(TMath::Floor((RelativePhi(part->Phi(),0.0) + TMath::Pi())/fPhiSize));
         etacell = TMath::Nint(TMath::Floor((part->Eta()+fEtaSize)/fEtaSize));
         sumPtOutOfCone[phicell][etacell]+= part->Pt();
      }
   }
   // Calculate rho
   static Double_t rhoInCells[20];
   Double_t  relativeArea;
   Int_t  nCells=0;
   Double_t bufferArea=0.0, bufferPt=0.0; //sum cells where A< fJetFreeAreaFrac
   for(Int_t ip=0; ip<fnPhi; ip++){
      for(Int_t ie=0; ie<fnEta; ie++){
         relativeArea = nOutCone[ip][ie]/fNofRndTrials;

         bufferArea += relativeArea;
         bufferPt   += sumPtOutOfCone[ip][ie];
         if(bufferArea > fJetFreeAreaFrac){
            rhoInCells[nCells] = bufferPt/(bufferArea*fCellArea);

            bufferArea = 0.0;
            bufferPt   = 0.0;
            nCells++;
         }
      }
   }

   if(nCells>0){
      rhoMedian = TMath::Median(nCells, rhoInCells);
   }  
   return rhoMedian;
}

//________________________________________________________________________
Double_t AliAnalysisTaskHJetSpectra::EstimateBgCone(){
   //Estimate background rho by means of integrating track pT outside identified jet cones
   Double_t rhoPerpCone = 0.0;
   
   Double_t pTleading  = -1.0;
   Double_t phiLeading = 1000.;
   Double_t etaLeading = 1000.;

   if(!fJetArray) return 0.0;

   for(Int_t ij = 0; ij < fJetArray->GetEntries(); ij++){
      AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJetArray->At(ij));
      if(!jet){
         AliError(Form("%s: Could not receive jet %d", GetName(), ij));
         continue;
      }
      if(!IsSignalJetInAcceptance(jet)) continue;

      if(pTleading < jet->Pt()){
         pTleading  = jet->Pt();
         phiLeading = jet->Phi();
         etaLeading = jet->Eta();
      }
   } 
   if(pTleading < 0.0) return 0.0;

   Double_t phileftcone  = phiLeading + TMath::Pi()/2;
   Double_t phirightcone = phiLeading - TMath::Pi()/2;

   /* Double_t dp, de;

   for(Int_t ip=0; ip < fTrackArray->GetEntries(); ip++){
 
      AliVTrack* part = static_cast<AliVTrack*>(fTrackArray->At(ip));
      if(!part) continue;
      if(!IsTrackInAcceptance((AliVParticle*) part)) continue; 


      dp = RelativePhi(phileftcone, part->Phi());
      de = etaLeading - part->Eta();
      if( dp*dp + de*de < fRandConeRadiusSquared ) rhoPerpCone += part->Pt();

      dp = RelativePhi(phirightcone, part->Phi());
      if( dp*dp + de*de < fRandConeRadiusSquared) rhoPerpCone += part->Pt();

   }*/

   rhoPerpCone +=  GetConePt(etaLeading, phileftcone,  fRandConeRadius);
   rhoPerpCone +=  GetConePt(etaLeading, phirightcone, fRandConeRadius);

   //normalize total pT by two times cone are 
   rhoPerpCone = rhoPerpCone/(2*TMath::Pi()*fRandConeRadiusSquared);
 


   return rhoPerpCone;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskHJetSpectra::DistantCones(Double_t phi1, Double_t eta1, Double_t r1, Double_t phi2, Double_t eta2, Double_t r2){
   //checks if the two cones are farther away than the sum of their radii

   Double_t dphi = RelativePhi(phi1,phi2);
   Double_t deta = eta1-eta2;
   Double_t d = r1+r2;
   if( dphi*dphi + deta*deta < d*d ) return kFALSE;

   return kTRUE;
}
