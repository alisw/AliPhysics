/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *									  *
 * Author: Martin Wilde, Daniel Lohner, Friederike Bock	       		  *
 * Version 1.0								  *
 *									  *
 * based on: on older version (see aliroot up to v5-04-42-AN)             *
 *           AliAnalysisTaskGammaConversion.cxx                           *
 *           Authors: Kathrin Koch, Kenneth Aamodt, Ana Marin             *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its	  *
 * documentation strictly for non-commercial purposes is hereby granted	  *
 * without fee, provided that the above copyright notice appears in all	  *
 * copies and that both the copyright notice and this permission notice	  *
 * appear in the supporting documentation. The authors make no claims	  *
 * about the suitability of this software for any purpose. It is	  *
 * provided "as is" without express or implied warranty.       		  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------

// Class used to do analysis on conversion pairs
//---------------------------------------------
///////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliESDVZERO.h"
#include "AliESDpid.h"
#include "AliAnalysisTaskGammaConvV1.h"
#include "AliVParticle.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliKFVertex.h"
#include "AliV0ReaderV1.h"
#include "AliGenCocktailEventHeader.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliEventplane.h"

ClassImp(AliAnalysisTaskGammaConvV1)

//________________________________________________________________________
AliAnalysisTaskGammaConvV1::AliAnalysisTaskGammaConvV1(): AliAnalysisTaskSE(),
   fV0Reader(NULL),
   fBGHandler(NULL),
   fBGHandlerRP(NULL),
   fInputEvent(NULL),
   fMCEvent(NULL),
   fMCStack(NULL),
   fCutFolder(NULL),
   fESDList(NULL),
   fBackList(NULL),
   fMotherList(NULL),
   fPhotonDCAList(NULL),
   fMesonDCAList(NULL),        
   fTrueList(NULL),
   fMCList(NULL),
   fHeaderNameList(NULL),
   fOutputContainer(0),
   fReaderGammas(NULL),
   fGammaCandidates(NULL),
   fCutArray(NULL),
   fConversionCuts(NULL),
   fMesonCutArray(NULL),
   fMesonCuts(NULL),
   hESDConvGammaPt(NULL),
   hESDConvGammaR(NULL),
   hESDConvGammaEta(NULL),
   tESDConvGammaPtDcazCat(NULL),
   fPtGamma(0),
   fDCAzPhoton(0),
   fRConvPhoton(0),
   fEtaPhoton(0),
   iCatPhoton(0),
   iPhotonMCInfo(0),
   hESDMotherInvMassPt(NULL),
   sESDMotherInvMassPtZM(NULL),
   hESDMotherBackInvMassPt(NULL),
   sESDMotherBackInvMassPtZM(NULL),
   hESDMotherInvMassEalpha(NULL),
   hESDMotherPi0PtY(NULL),
   hESDMotherEtaPtY(NULL),
   hESDMotherPi0PtAlpha(NULL),
   hESDMotherEtaPtAlpha(NULL),
   hESDMotherPi0PtOpenAngle(NULL),
   hESDMotherEtaPtOpenAngle(NULL),
   hMCHeaders(NULL),
   hMCAllGammaPt(NULL),
   hMCDecayGammaPi0Pt(NULL),
   hMCDecayGammaRhoPt(NULL),
   hMCDecayGammaEtaPt(NULL),
   hMCDecayGammaOmegaPt(NULL),
   hMCDecayGammaEtapPt(NULL),
   hMCDecayGammaPhiPt(NULL),
   hMCDecayGammaSigmaPt(NULL),
   hMCConvGammaPt(NULL),
   hMCConvGammaR(NULL),
   hMCConvGammaEta(NULL),
   hMCPi0Pt(NULL),
   hMCPi0WOWeightPt(NULL),
   hMCEtaPt(NULL),
   hMCEtaWOWeightPt(NULL),
   hMCPi0InAccPt(NULL),
   hMCEtaInAccPt(NULL),
   hMCPi0PtY(NULL),
   hMCEtaPtY(NULL),
   hMCK0sPt(NULL),
   hMCK0sWOWeightPt(NULL),
   hMCK0sPtY(NULL),
   hESDTrueMotherInvMassPt(NULL),
   hESDTruePrimaryMotherInvMassPt(NULL),
   hESDTruePrimaryMotherW0WeightingInvMassPt(NULL),
   pESDTruePrimaryMotherWeightsInvMassPt(NULL),
   hESDTruePrimaryPi0MCPtResolPt(NULL),
   hESDTruePrimaryEtaMCPtResolPt(NULL),
   hESDTrueSecondaryMotherInvMassPt(NULL),
   hESDTrueSecondaryMotherFromK0sInvMassPt(NULL),
   hESDTrueK0sWithPi0DaughterMCPt(NULL),
   hESDTrueSecondaryMotherFromEtaInvMassPt(NULL),
   hESDTrueEtaWithPi0DaughterMCPt(NULL),
   hESDTrueBckGGInvMassPt(NULL),
   hESDTrueBckContInvMassPt(NULL),
   hESDTruePi0PtY(NULL),
   hESDTrueEtaPtY(NULL),
   hESDTruePi0PtAlpha(NULL),
   hESDTrueEtaPtAlpha(NULL),
   hESDTruePi0PtOpenAngle(NULL),
   hESDTrueEtaPtOpenAngle(NULL),
   hESDTrueMotherDalitzInvMassPt(NULL),
   hESDTrueConvGammaPt(NULL),
   hESDTrueConvGammaEta(NULL),
   hESDCombinatorialPt(NULL),
   hESDTruePrimaryConvGammaPt(NULL),
   hESDTruePrimaryConvGammaESDPtMCPt(NULL),
   hESDTrueSecondaryConvGammaPt(NULL),
   hESDTrueSecondaryConvGammaFromXFromK0sPt(NULL),
   hESDTrueSecondaryConvGammaFromXFromLambdaPt(NULL),
   hESDTrueDalitzPsiPairDeltaPhi(NULL),
   hESDTrueGammaPsiPairDeltaPhi(NULL),
   hNEvents(NULL),
   hNGoodESDTracks(NULL),
   hNGammaCandidates(NULL),
   hNV0Tracks(NULL),
   hEtaShift(NULL),
   tESDMesonsInvMassPtDcazMinDcazMaxFlag(NULL),
   fInvMass(0),
   fPt(0),
   fDCAzGammaMin(0),
   fDCAzGammaMax(0),
   iFlag(0),
   iMesonMCInfo(0),
   fEventPlaneAngle(-100),
   fRandom(0),
   fnGammaCandidates(0),
   fUnsmearedPx(NULL),
   fUnsmearedPy(NULL),
   fUnsmearedPz(NULL),
   fUnsmearedE(NULL),
   fMCStackPos(NULL),
   fMCStackNeg(NULL),
   fESDArrayPos(NULL),
   fESDArrayNeg(NULL),
   fnCuts(0),
   fiCut(0),
   fMoveParticleAccordingToVertex(kTRUE),
   fIsHeavyIon(0),
   fDoMesonAnalysis(kTRUE),
   fDoMesonQA(0),
   fDoPhotonQA(0),
   fIsFromMBHeader(kTRUE),
   fIsMC(kFALSE)
{

}

//________________________________________________________________________
AliAnalysisTaskGammaConvV1::AliAnalysisTaskGammaConvV1(const char *name):
   AliAnalysisTaskSE(name),
   fV0Reader(NULL),
   fBGHandler(NULL),
   fBGHandlerRP(NULL),
   fInputEvent(NULL),
   fMCEvent(NULL),
   fMCStack(NULL),
   fCutFolder(NULL),
   fESDList(NULL),
   fBackList(NULL),
   fMotherList(NULL),
   fPhotonDCAList(NULL),
   fMesonDCAList(NULL),
   fTrueList(NULL),
   fMCList(NULL),
   fHeaderNameList(NULL),
   fOutputContainer(0),
   fReaderGammas(NULL),
   fGammaCandidates(NULL),
   fCutArray(NULL),
   fConversionCuts(NULL),
   fMesonCutArray(NULL),
   fMesonCuts(NULL),
   hESDConvGammaPt(NULL),
   hESDConvGammaR(NULL),
   hESDConvGammaEta(NULL),
   tESDConvGammaPtDcazCat(NULL),
   fPtGamma(0),
   fDCAzPhoton(0),
   fRConvPhoton(0),
   fEtaPhoton(0),
   iCatPhoton(0),
   iPhotonMCInfo(0),
   hESDMotherInvMassPt(NULL),
   sESDMotherInvMassPtZM(NULL),
   hESDMotherBackInvMassPt(NULL),
   sESDMotherBackInvMassPtZM(NULL),
   hESDMotherInvMassEalpha(NULL),
   hESDMotherPi0PtY(NULL),
   hESDMotherEtaPtY(NULL),
   hESDMotherPi0PtAlpha(NULL),
   hESDMotherEtaPtAlpha(NULL),
   hESDMotherPi0PtOpenAngle(NULL),
   hESDMotherEtaPtOpenAngle(NULL),
   hMCHeaders(NULL),
   hMCAllGammaPt(NULL),
   hMCDecayGammaPi0Pt(NULL),
   hMCDecayGammaRhoPt(NULL),
   hMCDecayGammaEtaPt(NULL),
   hMCDecayGammaOmegaPt(NULL),
   hMCDecayGammaEtapPt(NULL),
   hMCDecayGammaPhiPt(NULL),
   hMCDecayGammaSigmaPt(NULL),
   hMCConvGammaPt(NULL),
   hMCConvGammaR(NULL),
   hMCConvGammaEta(NULL),
   hMCPi0Pt(NULL),
   hMCPi0WOWeightPt(NULL),
   hMCEtaPt(NULL),
   hMCEtaWOWeightPt(NULL),
   hMCPi0InAccPt(NULL),
   hMCEtaInAccPt(NULL),
   hMCPi0PtY(NULL),
   hMCEtaPtY(NULL),
   hMCK0sPt(NULL),
   hMCK0sWOWeightPt(NULL),
   hMCK0sPtY(NULL),
   hESDTrueMotherInvMassPt(NULL),
   hESDTruePrimaryMotherInvMassPt(NULL),
   hESDTruePrimaryMotherW0WeightingInvMassPt(NULL),
   pESDTruePrimaryMotherWeightsInvMassPt(NULL),
   hESDTruePrimaryPi0MCPtResolPt(NULL),
   hESDTruePrimaryEtaMCPtResolPt(NULL),
   hESDTrueSecondaryMotherInvMassPt(NULL),
   hESDTrueSecondaryMotherFromK0sInvMassPt(NULL),
   hESDTrueK0sWithPi0DaughterMCPt(NULL),
   hESDTrueSecondaryMotherFromEtaInvMassPt(NULL),
   hESDTrueEtaWithPi0DaughterMCPt(NULL),
   hESDTrueBckGGInvMassPt(NULL),
   hESDTrueBckContInvMassPt(NULL),
   hESDTruePi0PtY(NULL),
   hESDTrueEtaPtY(NULL),
   hESDTruePi0PtAlpha(NULL),
   hESDTrueEtaPtAlpha(NULL),
   hESDTruePi0PtOpenAngle(NULL),
   hESDTrueEtaPtOpenAngle(NULL),
   hESDTrueMotherDalitzInvMassPt(NULL),
   hESDTrueConvGammaPt(NULL),
   hESDTrueConvGammaEta(NULL),
   hESDCombinatorialPt(NULL),
   hESDTruePrimaryConvGammaPt(NULL),
   hESDTruePrimaryConvGammaESDPtMCPt(NULL),
   hESDTrueSecondaryConvGammaPt(NULL),
   hESDTrueSecondaryConvGammaFromXFromK0sPt(NULL),
   hESDTrueSecondaryConvGammaFromXFromLambdaPt(NULL),
   hESDTrueDalitzPsiPairDeltaPhi(NULL),
   hESDTrueGammaPsiPairDeltaPhi(NULL),
   hNEvents(NULL),
   hNGoodESDTracks(NULL),
   hNGammaCandidates(NULL),
   hNV0Tracks(NULL),
   hEtaShift(NULL),
   tESDMesonsInvMassPtDcazMinDcazMaxFlag(NULL),
   fInvMass(0),
   fPt(0),
   fDCAzGammaMin(0),
   fDCAzGammaMax(0),
   iFlag(0),
   iMesonMCInfo(0),
   fEventPlaneAngle(-100),
   fRandom(0),
   fnGammaCandidates(0),
   fUnsmearedPx(NULL),
   fUnsmearedPy(NULL),
   fUnsmearedPz(NULL),
   fUnsmearedE(NULL),
   fMCStackPos(NULL),
   fMCStackNeg(NULL),
   fESDArrayPos(NULL),
   fESDArrayNeg(NULL),
   fnCuts(0),
   fiCut(0),
   fMoveParticleAccordingToVertex(kTRUE),
   fIsHeavyIon(0),
   fDoMesonAnalysis(kTRUE),
   fDoMesonQA(0),
   fDoPhotonQA(0),
   fIsFromMBHeader(kTRUE),
   fIsMC(kFALSE)
{
   // Define output slots here
   DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaConvV1::~AliAnalysisTaskGammaConvV1()
{
   if(fGammaCandidates){
      delete fGammaCandidates;
      fGammaCandidates = 0x0;
   }
   if(fBGHandler){
      delete[] fBGHandler;
      fBGHandler = 0x0;
   }
   if(fBGHandlerRP){
      delete[] fBGHandlerRP;
      fBGHandlerRP = 0x0;
   }
}
//___________________________________________________________
void AliAnalysisTaskGammaConvV1::InitBack(){

   const Int_t nDim = 4;
   Int_t nBins[nDim] = {800,250,7,4};
   Double_t xMin[nDim] = {0,0, 0,0};
   Double_t xMax[nDim] = {0.8,25,7,4};
   
   sESDMotherInvMassPtZM = new THnSparseF*[fnCuts];
   sESDMotherBackInvMassPtZM = new THnSparseF*[fnCuts];

   fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];
   fBGHandlerRP = new AliConversionAODBGHandlerRP*[fnCuts];
   for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      if (((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
         TString cutstring = ((AliConversionCuts*)fCutArray->At(iCut))->GetCutNumber();
         TString cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
         
         Int_t collisionSystem = atoi((TString)(((AliConversionCuts*)fCutArray->At(iCut))->GetCutNumber())(0,1));
         Int_t centMin = atoi((TString)(((AliConversionCuts*)fCutArray->At(iCut))->GetCutNumber())(1,1));
         Int_t centMax = atoi((TString)(((AliConversionCuts*)fCutArray->At(iCut))->GetCutNumber())(2,1));
         
         if(collisionSystem == 1 || collisionSystem == 2 ||
            collisionSystem == 5 || collisionSystem == 8 ||
            collisionSystem == 9){
            centMin = centMin*10;
            centMax = centMax*10; 
            if(centMax ==0 && centMax!=centMin) centMax=100;
         }
         else if(collisionSystem == 3 || collisionSystem == 6){
            centMin = centMin*5;
            centMax = centMax*5;
         }
         else if(collisionSystem == 4 || collisionSystem == 7){
            centMin = ((centMin*5)+45);
            centMax = ((centMax*5)+45);
         }
         
         fBackList[iCut] = new TList();
         fBackList[iCut]->SetName(Form("%s_%s Back histograms",cutstring.Data(),cutstringMeson.Data()));
         fBackList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fBackList[iCut]);

         sESDMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m","Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
         fBackList[iCut]->Add(sESDMotherBackInvMassPtZM[iCut]);

         fMotherList[iCut] = new TList();
         fMotherList[iCut]->SetName(Form("%s_%s Mother histograms",cutstring.Data(),cutstringMeson.Data()));
         fMotherList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fMotherList[iCut]);

         sESDMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m","Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
         fMotherList[iCut]->Add(sESDMotherInvMassPtZM[iCut]);

         if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
            fBGHandler[iCut] = new AliGammaConversionAODBGHandler(
                                                                  collisionSystem,centMin,centMax,
                                                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents(),
                                                                  ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity());
            fBGHandlerRP[iCut] = NULL;
         }
         else{
            fBGHandlerRP[iCut] = new AliConversionAODBGHandlerRP(
                                                                 ((AliConversionCuts*)fCutArray->At(fiCut))->IsHeavyIon(),
                                                                 ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity(),
                                                                 ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents());
            fBGHandler[iCut] = NULL;
         }
      }
   }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::UserCreateOutputObjects()
{

   // Create histograms
   if(fOutputContainer != NULL){
      delete fOutputContainer;
      fOutputContainer = NULL;
   }
   if(fOutputContainer == NULL){
      fOutputContainer = new TList();
      fOutputContainer->SetOwner(kTRUE);
   }

   // Array of current cut's gammas
   fGammaCandidates = new TList();

   fCutFolder = new TList*[fnCuts];
   fESDList = new TList*[fnCuts];
   fBackList = new TList*[fnCuts];
   fMotherList = new TList*[fnCuts];
   hNEvents = new TH1I*[fnCuts];
   hNGoodESDTracks = new TH1I*[fnCuts];
   hNGammaCandidates = new TH1I*[fnCuts];
   hNV0Tracks = new TH1I*[fnCuts];
   hEtaShift = new TProfile*[fnCuts];
   hESDConvGammaPt = new TH1F*[fnCuts];

   if (fDoPhotonQA == 2){
      fPhotonDCAList = new TList*[fnCuts];
      tESDConvGammaPtDcazCat = new TTree*[fnCuts];
   }
   if (fDoPhotonQA > 0){
      hESDConvGammaR = new TH1F*[fnCuts];
      hESDConvGammaEta = new TH1F*[fnCuts];
   } 
   
   if(fDoMesonAnalysis){
      hESDMotherInvMassPt = new TH2F*[fnCuts];
      hESDMotherBackInvMassPt = new TH2F*[fnCuts];
      hESDMotherInvMassEalpha = new TH2F*[fnCuts];
      if (fDoMesonQA == 2){
         fMesonDCAList = new TList*[fnCuts];
         tESDMesonsInvMassPtDcazMinDcazMaxFlag = new TTree*[fnCuts];
      }
      if (fDoMesonQA > 0){
         hESDMotherPi0PtY =  new TH2F*[fnCuts];
         hESDMotherEtaPtY =  new TH2F*[fnCuts];
         hESDMotherPi0PtAlpha =  new TH2F*[fnCuts];
         hESDMotherEtaPtAlpha =  new TH2F*[fnCuts];
         hESDMotherPi0PtOpenAngle =  new TH2F*[fnCuts];
         hESDMotherEtaPtOpenAngle =  new TH2F*[fnCuts];
      }   
   }

   for(Int_t iCut = 0; iCut<fnCuts;iCut++){

      TString cutstring = ((AliConversionCuts*)fCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson = "NoMesonCut";
      if(fDoMesonAnalysis)cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

      fCutFolder[iCut] = new TList();
      fCutFolder[iCut]->SetName(Form("Cut Number %s_%s",cutstring.Data(),cutstringMeson.Data()));
      fCutFolder[iCut]->SetOwner(kTRUE);
      fOutputContainer->Add(fCutFolder[iCut]);
      fESDList[iCut] = new TList();
      fESDList[iCut]->SetName(Form("%s_%s ESD histograms",cutstring.Data(),cutstringMeson.Data()));
      fESDList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fESDList[iCut]);

      hNEvents[iCut] = new TH1I("NEvents","NEvents",9,-0.5,8.5);
      hNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Missing MC");
      if (((AliConversionCuts*)fCutArray->At(iCut))->IsSpecialTrigger() == 4 ){
         TString TriggerNames = "Not Trigger: ";
         TriggerNames = TriggerNames+ ( (AliConversionCuts*)fCutArray->At(iCut))->GetSpecialTriggerName();
         hNEvents[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
      } else {
         hNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
      }
      hNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
      fESDList[iCut]->Add(hNEvents[iCut]);
      
      if(fIsHeavyIon == 1) hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",4000,0,4000);
      else if(fIsHeavyIon == 2) hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",400,0,400);
      else hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
      fESDList[iCut]->Add(hNGoodESDTracks[iCut]);
      if(fIsHeavyIon == 1) hNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",100,0,100);
      else if(fIsHeavyIon == 2) hNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",50,0,50);
      else hNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",50,0,50);
      fESDList[iCut]->Add(hNGammaCandidates[iCut]);
      if(fIsHeavyIon == 1) hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",30000,0,30000);
      else if(fIsHeavyIon == 2) hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",2500,0,2500);
      else hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",1500,0,1500);
      fESDList[iCut]->Add(hNV0Tracks[iCut]);
      hEtaShift[iCut] = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
      fESDList[iCut]->Add(hEtaShift[iCut]);
      hESDConvGammaPt[iCut] = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",250,0,25);
      fESDList[iCut]->Add(hESDConvGammaPt[iCut]);

      if (fDoPhotonQA == 2){
         fPhotonDCAList[iCut] = new TList();
         fPhotonDCAList[iCut]->SetName(Form("%s_%s Photon DCA tree",cutstring.Data(),cutstringMeson.Data()));
         fPhotonDCAList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fPhotonDCAList[iCut]);
            
        tESDConvGammaPtDcazCat[iCut] = new TTree("ESD_ConvGamma_Pt_Dcaz_R_Eta","ESD_ConvGamma_Pt_Dcaz_R_Eta_Cat");   
         tESDConvGammaPtDcazCat[iCut]->Branch("Pt",&fPtGamma,"fPtGamma/F");
         tESDConvGammaPtDcazCat[iCut]->Branch("DcaZPhoton",&fDCAzPhoton,"fDCAzPhoton/F");
//          tESDConvGammaPtDcazCat[iCut]->Branch("R",&fRConvPhoton,"fRConvPhoton/F");
//          tESDConvGammaPtDcazCat[iCut]->Branch("Eta",&fEtaPhoton,"fEtaPhoton/F");
         
         tESDConvGammaPtDcazCat[iCut]->Branch("cat",&iCatPhoton,"iCatPhoton/b");
         if(fIsMC){
            tESDConvGammaPtDcazCat[iCut]->Branch("photonMCInfo",&iPhotonMCInfo,"iPhotonMCInfo/b");
         }
         fPhotonDCAList[iCut]->Add(tESDConvGammaPtDcazCat[iCut]);
      }

      if (fDoPhotonQA > 0){
         hESDConvGammaR[iCut] = new TH1F("ESD_ConvGamma_R","ESD_ConvGamma_R",800,0,200);
         fESDList[iCut]->Add(hESDConvGammaR[iCut]);
         hESDConvGammaEta[iCut] = new TH1F("ESD_ConvGamma_Eta","ESD_ConvGamma_Eta",2000,-2,2);
         fESDList[iCut]->Add(hESDConvGammaEta[iCut]);
      }
      
      if(fDoMesonAnalysis){
         hESDMotherInvMassPt[iCut] = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",800,0,0.8,250,0,25);
         fESDList[iCut]->Add(hESDMotherInvMassPt[iCut]);
         hESDMotherBackInvMassPt[iCut] = new TH2F("ESD_Background_InvMass_Pt","ESD_Background_InvMass_Pt",800,0,0.8,250,0,25);
         fESDList[iCut]->Add(hESDMotherBackInvMassPt[iCut]);
         hESDMotherInvMassEalpha[iCut] = new TH2F("ESD_Mother_InvMass_vs_E_alpha","ESD_Mother_InvMass_vs_E_alpha",800,0,0.8,250,0,25);
         fESDList[iCut]->Add(hESDMotherInvMassEalpha[iCut]);
         if (fDoMesonQA == 2){    
            fMesonDCAList[iCut] = new TList();
            fMesonDCAList[iCut]->SetName(Form("%s_%s Meson DCA tree",cutstring.Data(),cutstringMeson.Data()));
            fMesonDCAList[iCut]->SetOwner(kTRUE);
            fCutFolder[iCut]->Add(fMesonDCAList[iCut]);
            
            tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut] = new TTree("ESD_Mesons_InvMass_Pt_DcazMin_DcazMax_Flag","ESD_Mesons_InvMass_Pt_DcazMin_DcazMax_Flag");   
            tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("InvMass",&fInvMass,"fInvMass/F");
            tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("Pt",&fPt,"fPt/F");
            tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("DcaZMin",&fDCAzGammaMin,"fDCAzGammaMin/F");
            tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("DcaZMax",&fDCAzGammaMax,"fDCAzGammaMax/F");
            tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("kind",&iFlag,"iFlag/b");
            if(fIsMC){
               tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]->Branch("mesonMCInfo",&iMesonMCInfo,"iMesonMCInfo/b");
            }
            fMesonDCAList[iCut]->Add(tESDMesonsInvMassPtDcazMinDcazMaxFlag[iCut]);
              
         }
         if (fDoMesonQA > 0 ){
            hESDMotherPi0PtY[iCut] = new TH2F("ESD_MotherPi0_Pt_Y","ESD_MotherPi0_Pt_Y",150,0.03,15.,150,-1.5,1.5);            
            SetLogBinningXTH2(hESDMotherPi0PtY[iCut]);
            fESDList[iCut]->Add(hESDMotherPi0PtY[iCut]);
            hESDMotherEtaPtY[iCut] = new TH2F("ESD_MotherEta_Pt_Y","ESD_MotherEta_Pt_Y",150,0.03,15.,150,-1.5,1.5);
            SetLogBinningXTH2(hESDMotherEtaPtY[iCut]);
            fESDList[iCut]->Add(hESDMotherEtaPtY[iCut]);
            hESDMotherPi0PtAlpha[iCut] = new TH2F("ESD_MotherPi0_Pt_Alpha","ESD_MotherPi0_Pt_Alpha",150,0.03,15.,100,0,1);            
            SetLogBinningXTH2(hESDMotherPi0PtAlpha[iCut]);
            fESDList[iCut]->Add(hESDMotherPi0PtAlpha[iCut]);
            hESDMotherEtaPtAlpha[iCut] = new TH2F("ESD_MotherEta_Pt_Alpha","ESD_MotherEta_Pt_Alpha",150,0.03,15.,100,0,1);
            SetLogBinningXTH2(hESDMotherEtaPtAlpha[iCut]);
            fESDList[iCut]->Add(hESDMotherEtaPtAlpha[iCut]);
            hESDMotherPi0PtOpenAngle[iCut] = new TH2F("ESD_MotherPi0_Pt_OpenAngle","ESD_MotherPi0_Pt_OpenAngle",150,0.03,15.,200,0,2*TMath::Pi());            
            SetLogBinningXTH2(hESDMotherPi0PtOpenAngle[iCut]);
            fESDList[iCut]->Add(hESDMotherPi0PtOpenAngle[iCut]);
            hESDMotherEtaPtOpenAngle[iCut] = new TH2F("ESD_MotherEta_Pt_OpenAngle","ESD_MotherEta_Pt_OpenAngle",150,0.03,15.,200,0,2*TMath::Pi());            
            SetLogBinningXTH2(hESDMotherEtaPtOpenAngle[iCut]);
            fESDList[iCut]->Add(hESDMotherEtaPtOpenAngle[iCut]);
         }
         
            
      }


   }
   if(fDoMesonAnalysis){
      InitBack(); // Init Background Handler
   }

   if(fIsMC){
      // MC Histogramms
      fMCList = new TList*[fnCuts];
      // True Histogramms
      fTrueList = new TList*[fnCuts];
      // Selected Header List
      fHeaderNameList = new TList*[fnCuts];
      hMCHeaders = new TH1I*[fnCuts];
      hMCAllGammaPt = new TH1F*[fnCuts];
      hMCDecayGammaPi0Pt = new TH1F*[fnCuts];
      hMCDecayGammaRhoPt = new TH1F*[fnCuts];
      hMCDecayGammaEtaPt = new TH1F*[fnCuts];
      hMCDecayGammaOmegaPt = new TH1F*[fnCuts];
      hMCDecayGammaEtapPt = new TH1F*[fnCuts];
      hMCDecayGammaPhiPt = new TH1F*[fnCuts];
      hMCDecayGammaSigmaPt = new TH1F*[fnCuts];
      hMCConvGammaPt = new TH1F*[fnCuts];
      hESDTrueConvGammaPt = new TH1F*[fnCuts];

      hESDCombinatorialPt = new TH2F*[fnCuts];
      hESDTruePrimaryConvGammaPt = new TH1F*[fnCuts];
      hESDTruePrimaryConvGammaESDPtMCPt = new TH2F*[fnCuts];
      hESDTrueSecondaryConvGammaPt = new TH1F*[fnCuts];
      hESDTrueSecondaryConvGammaFromXFromK0sPt = new TH1F*[fnCuts];
      hESDTrueSecondaryConvGammaFromXFromLambdaPt = new TH1F*[fnCuts];

      hESDTrueDalitzPsiPairDeltaPhi= new TH2F*[fnCuts];
      hESDTrueGammaPsiPairDeltaPhi= new TH2F*[fnCuts];

      if (fDoPhotonQA > 0){
         hMCConvGammaR = new TH1F*[fnCuts];
         hMCConvGammaEta = new TH1F*[fnCuts];
	 hESDTrueConvGammaEta = new TH1F*[fnCuts];
      }

      if(fDoMesonAnalysis){
         hMCPi0Pt = new TH1F*[fnCuts];
         hMCPi0WOWeightPt = new TH1F*[fnCuts];
         hMCEtaPt = new TH1F*[fnCuts];
         hMCEtaWOWeightPt = new TH1F*[fnCuts];
         hMCPi0InAccPt = new TH1F*[fnCuts];
         hMCEtaInAccPt = new TH1F*[fnCuts];

         hESDTrueMotherInvMassPt = new TH2F*[fnCuts];
         hESDTruePrimaryMotherInvMassPt = new TH2F*[fnCuts];
         hESDTruePrimaryMotherW0WeightingInvMassPt = new TH2F*[fnCuts];
         pESDTruePrimaryMotherWeightsInvMassPt = new TProfile2D*[fnCuts];
         hESDTrueSecondaryMotherInvMassPt = new TH2F*[fnCuts];
         hESDTrueSecondaryMotherFromK0sInvMassPt = new TH2F*[fnCuts];
         hESDTrueSecondaryMotherFromEtaInvMassPt = new TH2F*[fnCuts];
         hESDTrueMotherDalitzInvMassPt = new TH2F*[fnCuts];
         if (fDoMesonQA > 0){
            hMCPi0PtY = new TH2F*[fnCuts];
            hMCEtaPtY = new TH2F*[fnCuts];
            hMCK0sPt = new TH1F*[fnCuts];
            hMCK0sWOWeightPt = new TH1F*[fnCuts];
            hMCK0sPtY = new TH2F*[fnCuts];
            hESDTruePrimaryPi0MCPtResolPt = new TH2F*[fnCuts];
            hESDTruePrimaryEtaMCPtResolPt = new TH2F*[fnCuts];
            hESDTrueK0sWithPi0DaughterMCPt = new TH1F*[fnCuts];
            hESDTrueEtaWithPi0DaughterMCPt = new TH1F*[fnCuts];
            hESDTrueBckGGInvMassPt = new TH2F*[fnCuts];
            hESDTrueBckContInvMassPt = new TH2F*[fnCuts];
            hESDTruePi0PtY = new TH2F*[fnCuts];
            hESDTrueEtaPtY = new TH2F*[fnCuts];
            hESDTruePi0PtAlpha = new TH2F*[fnCuts];
            hESDTrueEtaPtAlpha = new TH2F*[fnCuts];
            hESDTruePi0PtOpenAngle = new TH2F*[fnCuts];
            hESDTrueEtaPtOpenAngle = new TH2F*[fnCuts];
         }
      }

      for(Int_t iCut = 0; iCut<fnCuts;iCut++){
         TString cutstring = ((AliConversionCuts*)fCutArray->At(iCut))->GetCutNumber();
         TString cutstringMeson = "NoMesonCut";
         if(fDoMesonAnalysis)cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

         fMCList[iCut] = new TList();
         fMCList[iCut]->SetName(Form("%s_%s MC histograms",cutstring.Data(),cutstringMeson.Data()));
         fMCList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fMCList[iCut]);
         hMCHeaders[iCut] = new TH1I("MC_Headers","MC_Headers",20,0,20);
         fMCList[iCut]->Add(hMCHeaders[iCut]);
         hMCAllGammaPt[iCut] = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",250,0,25);
         fMCList[iCut]->Add(hMCAllGammaPt[iCut]);
         hMCDecayGammaPi0Pt[iCut] = new TH1F("MC_DecayGammaPi0_Pt","MC_DecayGammaPi0_Pt",250,0,25);
         fMCList[iCut]->Add(hMCDecayGammaPi0Pt[iCut]);
         hMCDecayGammaRhoPt[iCut] = new TH1F("MC_DecayGammaRho_Pt","MC_DecayGammaRho_Pt",250,0,25);
         fMCList[iCut]->Add(hMCDecayGammaRhoPt[iCut]);
         hMCDecayGammaEtaPt[iCut] = new TH1F("MC_DecayGammaEta_Pt","MC_DecayGammaEta_Pt",250,0,25);
         fMCList[iCut]->Add(hMCDecayGammaEtaPt[iCut]);
         hMCDecayGammaOmegaPt[iCut] = new TH1F("MC_DecayGammaOmega_Pt","MC_DecayGammaOmmega_Pt",250,0,25);
         fMCList[iCut]->Add(hMCDecayGammaOmegaPt[iCut]);
         hMCDecayGammaEtapPt[iCut] = new TH1F("MC_DecayGammaEtap_Pt","MC_DecayGammaEtap_Pt",250,0,25);
         fMCList[iCut]->Add(hMCDecayGammaEtapPt[iCut]);
         hMCDecayGammaPhiPt[iCut] = new TH1F("MC_DecayGammaPhi_Pt","MC_DecayGammaPhi_Pt",250,0,25);
         fMCList[iCut]->Add(hMCDecayGammaPhiPt[iCut]);
         hMCDecayGammaSigmaPt[iCut] = new TH1F("MC_DecayGammaSigma_Pt","MC_DecayGammaSigma_Pt",250,0,25);
         fMCList[iCut]->Add(hMCDecayGammaSigmaPt[iCut]);
         hMCConvGammaPt[iCut] = new TH1F("MC_ConvGamma_Pt","MC_ConvGamma_Pt",250,0,25);
         fMCList[iCut]->Add(hMCConvGammaPt[iCut]);
      
         if (fDoPhotonQA > 0){
            hMCConvGammaR[iCut] = new TH1F("MC_ConvGamma_R","MC_ConvGamma_R",800,0,200);
            fMCList[iCut]->Add(hMCConvGammaR[iCut]);
            hMCConvGammaEta[iCut] = new TH1F("MC_ConvGamma_Eta","MC_ConvGamma_Eta",2000,-2,2);
            fMCList[iCut]->Add(hMCConvGammaEta[iCut]);
         }

         if(fDoMesonAnalysis){
            hMCPi0Pt[iCut] = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",250,0,25);
            hMCPi0Pt[iCut]->Sumw2();
            fMCList[iCut]->Add(hMCPi0Pt[iCut]);
            hMCPi0WOWeightPt[iCut] = new TH1F("MC_Pi0_WOWeights_Pt","MC_Pi0_WOWeights_Pt",250,0,25);
            hMCPi0WOWeightPt[iCut]->Sumw2();
            fMCList[iCut]->Add(hMCPi0WOWeightPt[iCut]);
            
            hMCEtaPt[iCut] = new TH1F("MC_Eta_Pt","MC_Eta_Pt",250,0,25);
            hMCEtaPt[iCut]->Sumw2();
            fMCList[iCut]->Add(hMCEtaPt[iCut]);
            hMCEtaWOWeightPt[iCut] = new TH1F("MC_Eta_WOWeights_Pt","MC_Eta_WOWeights_Pt",250,0,25);
            hMCEtaWOWeightPt[iCut]->Sumw2();
            fMCList[iCut]->Add(hMCEtaWOWeightPt[iCut]);
            hMCPi0InAccPt[iCut] = new TH1F("MC_Pi0InAcc_Pt","MC_Pi0InAcc_Pt",250,0,25);
            hMCPi0InAccPt[iCut]->Sumw2();
            fMCList[iCut]->Add(hMCPi0InAccPt[iCut]);
            hMCEtaInAccPt[iCut] = new TH1F("MC_EtaInAcc_Pt","MC_EtaInAcc_Pt",250,0,25);
            hMCEtaInAccPt[iCut]->Sumw2();
            fMCList[iCut]->Add(hMCEtaInAccPt[iCut]);
            if (fDoMesonQA > 0){
               hMCPi0PtY[iCut] = new TH2F("MC_Pi0_Pt_Y","MC_Pi0_Pt_Y",150,0.03,15.,150,-1.5,1.5);
               hMCPi0PtY[iCut]->Sumw2();
               SetLogBinningXTH2(hMCPi0PtY[iCut]);
               fMCList[iCut]->Add(hMCPi0PtY[iCut]);
               hMCEtaPtY[iCut] = new TH2F("MC_Eta_Pt_Y","MC_Eta_Pt_Y",150,0.03,15.,150,-1.5,1.5);
               hMCEtaPtY[iCut]->Sumw2();
               SetLogBinningXTH2(hMCEtaPtY[iCut]);
               fMCList[iCut]->Add(hMCEtaPtY[iCut]);
               hMCK0sPt[iCut] = new TH1F("MC_K0s_Pt","MC_K0s_Pt",150,0,15);
               hMCK0sPt[iCut]->Sumw2();
               fMCList[iCut]->Add(hMCK0sPt[iCut]);
               hMCK0sWOWeightPt[iCut] = new TH1F("MC_K0s_WOWeights_Pt","MC_K0s_WOWeights_Pt",150,0,15);
               hMCK0sWOWeightPt[iCut]->Sumw2();
               fMCList[iCut]->Add(hMCK0sWOWeightPt[iCut]);
               hMCK0sPtY[iCut] = new TH2F("MC_K0s_Pt_Y","MC_K0s_Pt_Y",150,0.03,15.,150,-1.5,1.5);
               hMCK0sPtY[iCut]->Sumw2();
               SetLogBinningXTH2(hMCK0sPtY[iCut]);
               fMCList[iCut]->Add(hMCK0sPtY[iCut]);
               
            }

         }
         fTrueList[iCut] = new TList();
         fTrueList[iCut]->SetName(Form("%s_%s True histograms",cutstring.Data(),cutstringMeson.Data()));
         fTrueList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fTrueList[iCut]);

         hESDTrueConvGammaPt[iCut] = new TH1F("ESD_TrueConvGamma_Pt","ESD_TrueConvGamma_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueConvGammaPt[iCut]);

         hESDCombinatorialPt[iCut] = new TH2F("ESD_TrueCombinatorial_Pt","ESD_TrueCombinatorial_Pt",250,0,25,16,-0.5,15.5);
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 1,"Elec+Elec");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 2,"Elec+Pion");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 3,"Elec+Kaon");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 4,"Elec+Proton");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 5,"Elec+Muon");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 6,"Pion+Pion");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 7,"Pion+Kaon");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 8,"Pion+Proton");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel( 9,"Pion+Muon");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(10,"Kaon+Kaon");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(11,"Kaon+Proton");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(12,"Kaon+Muon");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(13,"Proton+Proton");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(14,"Proton+Muon");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(15,"Muon+Muon");
         hESDCombinatorialPt[iCut]->GetYaxis()->SetBinLabel(16,"Rest");
         fTrueList[iCut]->Add(hESDCombinatorialPt[iCut]);
         hESDTruePrimaryConvGammaPt[iCut] = new TH1F("ESD_TruePrimaryConvGamma_Pt","ESD_TruePrimaryConvGamma_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTruePrimaryConvGammaPt[iCut]);
         hESDTrueSecondaryConvGammaPt[iCut] = new TH1F("ESD_TrueSecondaryConvGamma_Pt","ESD_TrueSecondaryConvGamma_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueSecondaryConvGammaPt[iCut]);

         hESDTrueSecondaryConvGammaFromXFromK0sPt[iCut]
            = new TH1F("ESD_TrueSecondaryConvGammaFromXFromK0s_Pt", "ESD_TrueSecondaryConvGammaFromXFromK0s_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueSecondaryConvGammaFromXFromK0sPt[iCut]);
         hESDTrueSecondaryConvGammaFromXFromLambdaPt[iCut]
            = new TH1F("ESD_TrueSecondaryConvGammaFromXFromLambda_Pt", "ESD_TrueSecondaryConvGammaFromXFromLambda_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueSecondaryConvGammaFromXFromLambdaPt[iCut]);

         hESDTrueDalitzPsiPairDeltaPhi[iCut]
            = new TH2F("ESD_TrueDalitzPsiPairDeltaPhi_Pt", "ESD_TrueDalitzPsiPairDeltaPhi_Pt",100,-0.5,2,100,-0.5,0.5);
         fTrueList[iCut]->Add(hESDTrueDalitzPsiPairDeltaPhi[iCut]);
         
         hESDTrueGammaPsiPairDeltaPhi[iCut]
            = new TH2F("ESD_TrueGammaPsiPairDeltaPhi_Pt", "ESD_TrueGammaPsiPairDeltaPhi_Pt",100,-0.5,2,100,-0.5,0.5);
         fTrueList[iCut]->Add(hESDTrueGammaPsiPairDeltaPhi[iCut]);

         
         hESDTruePrimaryConvGammaESDPtMCPt[iCut] = new TH2F("ESD_TruePrimaryConvGammaESD_PtMCPt", "ESD_TruePrimaryConvGammaESD_PtMCPt",250,0,25,250,0,25);
         fTrueList[iCut]->Add(hESDTruePrimaryConvGammaESDPtMCPt[iCut]);
         
         if(fDoMesonAnalysis){
            hESDTrueMotherInvMassPt[iCut] = new TH2F("ESD_TrueMother_InvMass_Pt","ESD_TrueMother_InvMass_Pt",800,0,0.8,250,0,25);
            fTrueList[iCut]->Add(hESDTrueMotherInvMassPt[iCut]);
            hESDTruePrimaryMotherInvMassPt[iCut]
               = new TH2F("ESD_TruePrimaryMother_InvMass_Pt", "ESD_TruePrimaryMother_InvMass_Pt", 800,0,0.8,250,0,25);
            hESDTruePrimaryMotherInvMassPt[iCut]->Sumw2();
            fTrueList[iCut]->Add(hESDTruePrimaryMotherInvMassPt[iCut]);
            hESDTruePrimaryMotherW0WeightingInvMassPt[iCut]
               = new TH2F("ESD_TruePrimaryMotherW0Weights_InvMass_Pt", "ESD_TruePrimaryMotherW0Weights_InvMass_Pt", 800,0,0.8,250,0,25);
            hESDTruePrimaryMotherW0WeightingInvMassPt[iCut]->Sumw2();
            fTrueList[iCut]->Add(hESDTruePrimaryMotherW0WeightingInvMassPt[iCut]);
            pESDTruePrimaryMotherWeightsInvMassPt[iCut]
               = new TProfile2D("ESD_TruePrimaryMotherWeights_InvMass_Pt", "ESD_TruePrimaryMotherWeights_InvMass_Pt", 800,0,0.8,250,0,25);
            pESDTruePrimaryMotherWeightsInvMassPt[iCut]->Sumw2();
            fTrueList[iCut]->Add(pESDTruePrimaryMotherWeightsInvMassPt[iCut]);
            hESDTrueSecondaryMotherInvMassPt[iCut]
               = new TH2F("ESD_TrueSecondaryMother_InvMass_Pt", "ESD_TrueSecondaryMother_InvMass_Pt", 800,0,0.8,250,0,25);
            hESDTrueSecondaryMotherInvMassPt[iCut]->Sumw2();
            fTrueList[iCut]->Add(hESDTrueSecondaryMotherInvMassPt[iCut]);
            hESDTrueSecondaryMotherFromK0sInvMassPt[iCut]
               = new TH2F("ESD_TrueSecondaryMotherFromK0s_InvMass_Pt","ESD_TrueSecondaryMotherFromK0s_InvMass_Pt",800,0,0.8,250,0,25);
            hESDTrueSecondaryMotherFromK0sInvMassPt[iCut]->Sumw2();
            fTrueList[iCut]->Add(hESDTrueSecondaryMotherFromK0sInvMassPt[iCut]);
            hESDTrueSecondaryMotherFromEtaInvMassPt[iCut]
               = new TH2F("ESD_TrueSecondaryMotherFromEta_InvMass_Pt","ESD_TrueSecondaryMotherFromEta_InvMass_Pt",800,0,0.8,250,0,25);
            fTrueList[iCut]->Add(hESDTrueSecondaryMotherFromEtaInvMassPt[iCut]);
            hESDTrueMotherDalitzInvMassPt[iCut] = new TH2F("ESD_TrueDalitz_InvMass_Pt","ESD_TrueDalitz_InvMass_Pt",800,0,0.8,250,0,25);
            fTrueList[iCut]->Add(hESDTrueMotherDalitzInvMassPt[iCut]);         
            if (fDoMesonQA > 0){
               hESDTruePrimaryPi0MCPtResolPt[iCut] = new TH2F("ESD_TruePrimaryPi0_MCPt_ResolPt","ESD_TruePrimaryPi0_ResolPt_MCPt",500,0.03,25,1000,-1.,1.);
               hESDTruePrimaryPi0MCPtResolPt[iCut]->Sumw2();
               SetLogBinningXTH2(hESDTruePrimaryPi0MCPtResolPt[iCut]);
               fTrueList[iCut]->Add(hESDTruePrimaryPi0MCPtResolPt[iCut]);
               hESDTruePrimaryEtaMCPtResolPt[iCut]  = new TH2F("ESD_TruePrimaryEta_MCPt_ResolPt","ESD_TruePrimaryEta_ResolPt_MCPt",500,0.03,25,1000,-1.,1.);
               hESDTruePrimaryEtaMCPtResolPt[iCut]->Sumw2();
               SetLogBinningXTH2(hESDTruePrimaryEtaMCPtResolPt[iCut]);
               fTrueList[iCut]->Add(hESDTruePrimaryEtaMCPtResolPt[iCut]);
               hESDTrueBckGGInvMassPt[iCut] = new TH2F("ESD_TrueBckGG_InvMass_Pt","ESD_TrueBckGG_InvMass_Pt",800,0,0.8,250,0,25);
               fTrueList[iCut]->Add(hESDTrueBckGGInvMassPt[iCut]);
               hESDTrueBckContInvMassPt[iCut] = new TH2F("ESD_TrueBckCont_InvMass_Pt","ESD_TrueBckCont_InvMass_Pt",800,0,0.8,250,0,25);
               fTrueList[iCut]->Add(hESDTrueBckContInvMassPt[iCut]);
               hESDTrueK0sWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueK0sWithPi0Daughter_MCPt","ESD_TrueK0sWithPi0Daughter_MCPt",250,0,25);
               fTrueList[iCut]->Add(hESDTrueK0sWithPi0DaughterMCPt[iCut]);
               hESDTrueEtaWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueEtaWithPi0Daughter_MCPt","ESD_TrueEtaWithPi0Daughter_MCPt",250,0,25);
               fTrueList[iCut]->Add(hESDTrueEtaWithPi0DaughterMCPt[iCut]);
               hESDTruePi0PtY[iCut] = new TH2F("ESD_TruePi0_Pt_Y","ESD_TruePi0_Pt_Y",150,0.03,15.,150,-1.5,1.5);
               SetLogBinningXTH2(hESDTruePi0PtY[iCut]);
               fTrueList[iCut]->Add(hESDTruePi0PtY[iCut]);
               hESDTrueEtaPtY[iCut] = new TH2F("ESD_TrueEta_Pt_Y","ESD_TrueEta_Pt_Y",150,0.03,15.,150,-1.5,1.5);
               SetLogBinningXTH2(hESDTrueEtaPtY[iCut]);
               fTrueList[iCut]->Add(hESDTrueEtaPtY[iCut]);
               hESDTruePi0PtAlpha[iCut] = new TH2F("ESD_TruePi0_Pt_Alpha","ESD_TruePi0_Pt_Alpha",150,0.03,15.,100,0,1);
               SetLogBinningXTH2(hESDTruePi0PtAlpha[iCut]);
               fTrueList[iCut]->Add(hESDTruePi0PtAlpha[iCut]);
               hESDTrueEtaPtAlpha[iCut] = new TH2F("ESD_TrueEta_Pt_Alpha","ESD_TrueEta_Pt_Alpha",150,0.03,15.,100,0,1);
               SetLogBinningXTH2(hESDTrueEtaPtAlpha[iCut]);
               fTrueList[iCut]->Add(hESDTrueEtaPtAlpha[iCut]);
               
               hESDTruePi0PtOpenAngle[iCut] = new TH2F("ESD_TruePi0_Pt_OpenAngle","ESD_TruePi0_Pt_OpenAngle",150,0.03,15.,200,0,2*TMath::Pi());            
               SetLogBinningXTH2(hESDTruePi0PtOpenAngle[iCut]);
               fTrueList[iCut]->Add(hESDTruePi0PtOpenAngle[iCut]);
               hESDTrueEtaPtOpenAngle[iCut] = new TH2F("ESD_TrueEta_Pt_OpenAngle","ESD_TrueEta_Pt_OpenAngle",150,0.03,15.,200,0,2*TMath::Pi());            
               SetLogBinningXTH2(hESDTrueEtaPtOpenAngle[iCut]);
               fTrueList[iCut]->Add(hESDTrueEtaPtOpenAngle[iCut]);

               hESDTrueConvGammaEta[iCut] = new TH1F("ESD_TrueConvGamma_Eta","ESD_TrueConvGamma_Eta",2000,-2,2);
	       fTrueList[iCut]->Add(hESDTrueConvGammaEta[iCut]);

            }
         }
      }
   }

   fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
   if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

   if(fV0Reader)
      if((AliConversionCuts*)fV0Reader->GetConversionCuts())
         if(((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
            fOutputContainer->Add(((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());

   for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      if(!((AliConversionCuts*)fCutArray->At(iCut))) continue;
      if(((AliConversionCuts*)fCutArray->At(iCut))->GetCutHistograms()){
         fCutFolder[iCut]->Add(((AliConversionCuts*)fCutArray->At(iCut))->GetCutHistograms());
      }
      if(fDoMesonAnalysis){
         if(!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))) continue;
         if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
            fCutFolder[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
         }
      }
   }
   PostData(1, fOutputContainer);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvV1::Notify()
{
   for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      if(!((AliConversionCuts*)fCutArray->At(iCut))->GetDoEtaShift()){
         hEtaShift[iCut]->Fill(0.,(((AliConversionCuts*)fCutArray->At(iCut))->GetEtaShift()));
         continue; // No Eta Shift requested, continue
      }
      if(((AliConversionCuts*)fCutArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
         ((AliConversionCuts*)fCutArray->At(iCut))->GetCorrectEtaShiftFromPeriod(fV0Reader->GetPeriodName());
          hEtaShift[iCut]->Fill(0.,(((AliConversionCuts*)fCutArray->At(iCut))->GetEtaShift()));
         ((AliConversionCuts*)fCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once   
         continue;
      }
      else{
         printf(" Gamma Conversion Task %s :: Eta Shift Manually Set to %f \n\n",
                (((AliConversionCuts*)fCutArray->At(iCut))->GetCutNumber()).Data(),((AliConversionCuts*)fCutArray->At(iCut))->GetEtaShift());
          hEtaShift[iCut]->Fill(0.,(((AliConversionCuts*)fCutArray->At(iCut))->GetEtaShift()));
         ((AliConversionCuts*)fCutArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once   
      }
   }
   
   return kTRUE;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaConvV1::UserExec(Option_t *)
{
   //
   // Called for each event
   //
   Int_t eventQuality = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetEventQuality();
   if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1
      for(Int_t iCut = 0; iCut<fnCuts; iCut++){
         hNEvents[iCut]->Fill(eventQuality);
      }
      return;
   }

   if(fIsMC) fMCEvent = MCEvent();
   if(fMCEvent == NULL) fIsMC = kFALSE;
   
   fInputEvent = InputEvent();

   if(fIsMC && fInputEvent->IsA()==AliESDEvent::Class()){
      fMCStack = fMCEvent->Stack();
      if(fMCStack == NULL) fIsMC = kFALSE;
   }

   fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
   
   // ------------------- BeginEvent ----------------------------

   AliEventplane *EventPlane = fInputEvent->GetEventplane();
   if(fIsHeavyIon ==1)fEventPlaneAngle = EventPlane->GetEventplane("V0",fInputEvent,2);
   else fEventPlaneAngle=0.0;
   
   if(fIsMC && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
      RelabelAODPhotonCandidates(kTRUE);    // In case of AODMC relabeling MC
      fV0Reader->RelabelAODs(kTRUE);
   }
   for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      fiCut = iCut;
      Int_t eventNotAccepted =
         ((AliConversionCuts*)fCutArray->At(iCut))
         ->IsEventAcceptedByConversionCut(fV0Reader->GetConversionCuts(),fInputEvent,fMCEvent,fIsHeavyIon);
      if(eventNotAccepted){
         // cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
         hNEvents[iCut]->Fill(eventNotAccepted); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
         continue;
      }

      if(eventQuality != 0){// Event Not Accepted
         // cout << "event rejected due to: " <<eventQuality << endl;
         hNEvents[iCut]->Fill(eventQuality);
         continue;
      }

      hNEvents[iCut]->Fill(eventQuality); // Should be 0 here
      hNGoodESDTracks[iCut]->Fill(fV0Reader->GetNumberOfPrimaryTracks());
      if(((AliConversionCuts*)fCutArray->At(iCut))->IsHeavyIon() == 2)	hNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A());
      else hNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());

      if(fIsMC){
         // Process MC Particle
         if(((AliConversionCuts*)fCutArray->At(iCut))->GetSignalRejection() != 0){
            if(fInputEvent->IsA()==AliESDEvent::Class()){
               ((AliConversionCuts*)fCutArray->At(iCut))->GetNotRejectedParticles(((AliConversionCuts*)fCutArray->At(iCut))->GetSignalRejection(),
                                                                                  ((AliConversionCuts*)fCutArray->At(iCut))->GetAcceptedHeader(),
                                                                                  fMCEvent);
            }
            else if(fInputEvent->IsA()==AliAODEvent::Class()){
               ((AliConversionCuts*)fCutArray->At(iCut))->GetNotRejectedParticles(((AliConversionCuts*)fCutArray->At(iCut))->GetSignalRejection(),
                                                                                  ((AliConversionCuts*)fCutArray->At(iCut))->GetAcceptedHeader(),
                                                                                  fInputEvent);
            }

            if(((AliConversionCuts*)fCutArray->At(iCut))->GetAcceptedHeader()){
               for(Int_t i = 0;i<(((AliConversionCuts*)fCutArray->At(iCut))->GetAcceptedHeader())->GetEntries();i++){
                  TString nameBin= hMCHeaders[iCut]->GetXaxis()->GetBinLabel(i+1);
                  if (nameBin.CompareTo("")== 0){
                     TString nameHeader = ((TObjString*)((TList*)((AliConversionCuts*)fCutArray->At(iCut))
                                                         ->GetAcceptedHeader())->At(i))->GetString();
                     hMCHeaders[iCut]->GetXaxis()->SetBinLabel(i+1,nameHeader.Data());
                  }
               }
            }
         }
      }
      if(fIsMC){
         if(fInputEvent->IsA()==AliESDEvent::Class())
            ProcessMCParticles();
         if(fInputEvent->IsA()==AliAODEvent::Class())
            ProcessAODMCParticles();
      }

      ProcessPhotonCandidates(); // Process this cuts gammas

      hNGammaCandidates[iCut]->Fill(fGammaCandidates->GetEntries());
      if(fDoMesonAnalysis){ // Meson Analysis
         if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fIsMC){
            fUnsmearedPx = new Double_t[fGammaCandidates->GetEntries()]; // Store unsmeared Momenta
            fUnsmearedPy = new Double_t[fGammaCandidates->GetEntries()];
            fUnsmearedPz = new Double_t[fGammaCandidates->GetEntries()];
            fUnsmearedE =  new Double_t[fGammaCandidates->GetEntries()];

            for(Int_t gamma=0;gamma<fGammaCandidates->GetEntries();gamma++){ // Smear the AODPhotons in MC
               fUnsmearedPx[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Px();
               fUnsmearedPy[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Py();
               fUnsmearedPz[gamma] = ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->Pz();
               fUnsmearedE[gamma] =  ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->E();
               ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->SmearParticle(dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(gamma)));
            }
         }

         CalculatePi0Candidates(); // Combine Gammas
         if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
            if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->BackgroundHandlerType() == 0){
               CalculateBackground(); // Combinatorial Background
               UpdateEventByEventData(); // Store Event for mixed Events
            }
            else{
               CalculateBackgroundRP(); // Combinatorial Background
               fBGHandlerRP[iCut]->AddEvent(fGammaCandidates,fInputEvent); // Store Event for mixed Events
            }
         }
         if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fIsMC){
            for(Int_t gamma=0;gamma<fGammaCandidates->GetEntries();gamma++){ // Smear the AODPhotons in MC
               ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPx(fUnsmearedPx[gamma]); // Reset Unsmeared Momenta
               ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPy(fUnsmearedPy[gamma]);
               ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetPz(fUnsmearedPz[gamma]);
               ((AliAODConversionPhoton*)fGammaCandidates->At(gamma))->SetE(fUnsmearedE[gamma]);
            }
            delete[] fUnsmearedPx; fUnsmearedPx = 0x0;
            delete[] fUnsmearedPy; fUnsmearedPy = 0x0;
            delete[] fUnsmearedPz; fUnsmearedPz = 0x0;
            delete[] fUnsmearedE;  fUnsmearedE  = 0x0;
         }
      }
      fGammaCandidates->Clear(); // delete this cuts good gammas
   }

   if(fIsMC && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
      RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
      fV0Reader->RelabelAODs(kFALSE);
   }

   PostData(1, fOutputContainer);
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessPhotonCandidates()
{
   Int_t nV0 = 0;
   TList *GammaCandidatesStepOne = new TList();
   TList *GammaCandidatesStepTwo = new TList();
   // Loop over Photon Candidates allocated by ReaderV1
   for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
      AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
      if(!PhotonCandidate) continue;
      fIsFromMBHeader = kTRUE;
      if(fIsMC && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 0){
         Int_t isPosFromMBHeader
            = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack, fInputEvent);
         if(isPosFromMBHeader == 0 && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
         Int_t isNegFromMBHeader
            = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack, fInputEvent);
         if(isNegFromMBHeader == 0 && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 3) continue;

         if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
      }

      if(!((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
      if(!((AliConversionCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(PhotonCandidate->GetPhotonPhi(),fEventPlaneAngle)) continue;
      if(!((AliConversionCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
         !((AliConversionCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
         fGammaCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas

         if(fIsFromMBHeader){
            hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
            if (fDoPhotonQA > 0){
               hESDConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
               hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
            }   
         }
         if(fIsMC){
            if(fInputEvent->IsA()==AliESDEvent::Class())
               ProcessTruePhotonCandidates(PhotonCandidate);
            if(fInputEvent->IsA()==AliAODEvent::Class())
               ProcessTruePhotonCandidatesAOD(PhotonCandidate);
         }
         if (fIsFromMBHeader && fDoPhotonQA == 2){
            if (fIsHeavyIon == 1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
               fPtGamma = PhotonCandidate->Pt();
               fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
               fRConvPhoton = PhotonCandidate->GetConversionRadius();
               fEtaPhoton = PhotonCandidate->GetPhotonEta();
               iCatPhoton = PhotonCandidate->GetPhotonQuality();
               tESDConvGammaPtDcazCat[fiCut]->Fill();
            } else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
               fPtGamma = PhotonCandidate->Pt();
               fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
               fRConvPhoton = PhotonCandidate->GetConversionRadius();
               fEtaPhoton = PhotonCandidate->GetPhotonEta();
               iCatPhoton = PhotonCandidate->GetPhotonQuality();
               tESDConvGammaPtDcazCat[fiCut]->Fill();
            }   
         }   
      } else if(((AliConversionCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
         ((AliConversionCuts*)fCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
         nV0++;
         GammaCandidatesStepOne->Add(PhotonCandidate);
      } else if(!((AliConversionCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
              ((AliConversionCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
         GammaCandidatesStepTwo->Add(PhotonCandidate);
      }
   }
   if(((AliConversionCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){
      for(Int_t i = 0;i<GammaCandidatesStepOne->GetEntries();i++){
         AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GammaCandidatesStepOne->At(i);
         if(!PhotonCandidate) continue;
         fIsFromMBHeader = kTRUE;
         if(fMCStack && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 0){
            Int_t isPosFromMBHeader
               = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack, fInputEvent);
            Int_t isNegFromMBHeader
               = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack, fInputEvent);
            if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
         }
         if(!((AliConversionCuts*)fCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GammaCandidatesStepOne->GetEntries())) continue;
         if(!((AliConversionCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
            fGammaCandidates->Add(PhotonCandidate);
            if(fIsFromMBHeader){
               hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
               if (fDoPhotonQA > 0){
                  hESDConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
                  hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
               }
            }
         }
         if(fIsMC){
            if(fInputEvent->IsA()==AliESDEvent::Class())
               ProcessTruePhotonCandidates(PhotonCandidate);
            if(fInputEvent->IsA()==AliAODEvent::Class())
               ProcessTruePhotonCandidatesAOD(PhotonCandidate);
         } else GammaCandidatesStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
         
           if (fIsFromMBHeader && fDoPhotonQA == 2){
            if (fIsHeavyIon ==1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
               fPtGamma = PhotonCandidate->Pt();
               fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
               fRConvPhoton = PhotonCandidate->GetConversionRadius();
               fEtaPhoton = PhotonCandidate->GetPhotonEta();
               iCatPhoton = PhotonCandidate->GetPhotonQuality();
               tESDConvGammaPtDcazCat[fiCut]->Fill();
            } else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
               fPtGamma = PhotonCandidate->Pt();
               fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
               fRConvPhoton = PhotonCandidate->GetConversionRadius();
               fEtaPhoton = PhotonCandidate->GetPhotonEta();
               iCatPhoton = PhotonCandidate->GetPhotonQuality();
               tESDConvGammaPtDcazCat[fiCut]->Fill();
            }
         }
      }
   }
   if(((AliConversionCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
      for(Int_t i = 0;i<GammaCandidatesStepTwo->GetEntries();i++){
         AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GammaCandidatesStepTwo->At(i);
         if(!PhotonCandidate) continue;
         fIsFromMBHeader = kTRUE;
         if(fMCStack && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 0){
            Int_t isPosFromMBHeader
               = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack, fInputEvent);
            Int_t isNegFromMBHeader
               = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack, fInputEvent);
            if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
         }
         if(!((AliConversionCuts*)fCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GammaCandidatesStepTwo,i)) continue;
         fGammaCandidates->Add(PhotonCandidate); // Add gamma to current cut TList
         if(fIsFromMBHeader){
            hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
            if (fDoPhotonQA > 0){
               hESDConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
               hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
            }
         }
         if(fIsMC){
            if(fInputEvent->IsA()==AliESDEvent::Class())
               ProcessTruePhotonCandidates(PhotonCandidate);
            if(fInputEvent->IsA()==AliAODEvent::Class())
               ProcessTruePhotonCandidatesAOD(PhotonCandidate);
         }
         if (fIsFromMBHeader){
           if (fIsHeavyIon == 1 && PhotonCandidate->Pt() > 0.399 && PhotonCandidate->Pt() < 12.){
               fPtGamma = PhotonCandidate->Pt();
               fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
               fRConvPhoton = PhotonCandidate->GetConversionRadius();
               fEtaPhoton = PhotonCandidate->GetPhotonEta();
               iCatPhoton = PhotonCandidate->GetPhotonQuality();
               tESDConvGammaPtDcazCat[fiCut]->Fill();
            } else if ( PhotonCandidate->Pt() > 0.299 && PhotonCandidate->Pt() < 16.){
               fPtGamma = PhotonCandidate->Pt();
               fDCAzPhoton = PhotonCandidate->GetDCAzToPrimVtx();
               fRConvPhoton = PhotonCandidate->GetConversionRadius();
               fEtaPhoton = PhotonCandidate->GetPhotonEta();
               iCatPhoton = PhotonCandidate->GetPhotonQuality();
               tESDConvGammaPtDcazCat[fiCut]->Fill();
            }
         }
      }
   }

   delete GammaCandidatesStepOne;
   GammaCandidatesStepOne = 0x0;
   delete GammaCandidatesStepTwo;
   GammaCandidatesStepTwo = 0x0;

}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessTruePhotonCandidatesAOD(AliAODConversionPhoton *TruePhotonCandidate)
{

   Double_t magField = fInputEvent->GetMagneticField();
   if( magField  < 0.0 ){
      magField =  1.0;
   }
   else {
      magField =  -1.0;
   }

   TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
   AliAODMCParticle *posDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
   AliAODMCParticle *negDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());
   iPhotonMCInfo = 0;
   
   if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
   Int_t pdgCode[2] = {abs(posDaughter->GetPdgCode()),abs(negDaughter->GetPdgCode())};

   if(posDaughter->GetMother() != negDaughter->GetMother()){
      FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode);
      iPhotonMCInfo = 1;
      return;
   }
   else if(posDaughter->GetMother() == -1){
      FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode);
      iPhotonMCInfo = 1;
      return;
   }

   if(pdgCode[0]!=11 || pdgCode[1]!=11){
      iPhotonMCInfo = 1;
      return; //One Particle is not a electron
   }
   if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()){
      iPhotonMCInfo = 1; 
      return; // Same Charge
   }
   
   AliAODMCParticle *Photon = (AliAODMCParticle*) AODMCTrackArray->At(posDaughter->GetMother());
   AliVTrack * electronCandidate = ((AliConversionCuts*)fCutArray->At(fiCut))->GetTrack(fInputEvent,TruePhotonCandidate->GetTrackLabelNegative() );
   AliVTrack * positronCandidate = ((AliConversionCuts*)fCutArray->At(fiCut))->GetTrack(fInputEvent,TruePhotonCandidate->GetTrackLabelPositive() );
   Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronCandidate->Phi()-positronCandidate->Phi());

   if(Photon->GetPdgCode() != 22){
      hESDTrueDalitzPsiPairDeltaPhi[fiCut]->Fill(deltaPhi,TruePhotonCandidate->GetPsiPair());
      iPhotonMCInfo = 1;
      return; // Mother is no Photon
   }
   
   if(((posDaughter->GetMCProcessCode())) != 5 || ((negDaughter->GetMCProcessCode())) != 5){
      iPhotonMCInfo = 1;
      return;// check if the daughters come from a conversion 
   }   
   // STILL A BUG IN ALIROOT >>8 HAS TPO BE REMOVED AFTER FIX
   
   
   
   // True Photon
   if(fIsFromMBHeader){
     hESDTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
     if (fDoPhotonQA > 0) hESDTrueConvGammaEta[fiCut]->Fill(TruePhotonCandidate->Eta());
   }
   hESDTrueGammaPsiPairDeltaPhi[fiCut]->Fill(deltaPhi,TruePhotonCandidate->GetPsiPair());  
   if(Photon->IsPrimary()){ 
      // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
      if(fIsFromMBHeader){
         iPhotonMCInfo = 6;
         hESDTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
         hESDTruePrimaryConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
      }
      // (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
   }
   else{
      if(fIsFromMBHeader){
         hESDTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
         iPhotonMCInfo = 2;
         if(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother() > -1 &&
            ((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 3122){
            iPhotonMCInfo = 5;
            hESDTrueSecondaryConvGammaFromXFromLambdaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
         }
         if(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother() > -1 &&
            ((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 310){
            iPhotonMCInfo = 4;
            hESDTrueSecondaryConvGammaFromXFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt());
         }
         if(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother() > -1 &&
            ((AliAODMCParticle*)AODMCTrackArray->At(((AliAODMCParticle*)AODMCTrackArray->At(Photon->GetMother()))->GetMother()))->GetPdgCode() == 221){
            iPhotonMCInfo = 3;
         }
      }
   }
   
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessTruePhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
   
  Double_t magField = fInputEvent->GetMagneticField();
   if( magField  < 0.0 ){
      magField =  1.0;
   }
   else {
      magField =  -1.0;
   }

   // Process True Photons
   TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCStack);
   TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCStack);

   iPhotonMCInfo = 0;
   
   if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
   Int_t pdgCode[2] = {abs(posDaughter->GetPdgCode()),abs(negDaughter->GetPdgCode())};
   iPhotonMCInfo = 1;
   if(posDaughter->GetMother(0) != negDaughter->GetMother(0)){
      FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode);
      return;
   }
   else if(posDaughter->GetMother(0) == -1){
      FillPhotonCombinatorialBackgroundHist(TruePhotonCandidate, pdgCode);
      return;
   }

   if(pdgCode[0]!=11 || pdgCode[1]!=11) return; //One Particle is not a electron
   
   if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge

   TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCStack);
   AliVTrack * electronCandidate = ((AliConversionCuts*)fCutArray->At(fiCut))->GetTrack(fInputEvent,TruePhotonCandidate->GetTrackLabelNegative() );
   AliVTrack * positronCandidate = ((AliConversionCuts*)fCutArray->At(fiCut))->GetTrack(fInputEvent,TruePhotonCandidate->GetTrackLabelPositive() );
   Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronCandidate->Phi()-positronCandidate->Phi());
   
   if(Photon->GetPdgCode() != 22){
      hESDTrueDalitzPsiPairDeltaPhi[fiCut]->Fill(deltaPhi,TruePhotonCandidate->GetPsiPair());
      return; // Mother is no Photon
   }
   
   if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion

   
   
   // True Photon
    if(fIsFromMBHeader){
     hESDTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
     if (fDoPhotonQA > 0) hESDTrueConvGammaEta[fiCut]->Fill(TruePhotonCandidate->Eta());
   }
   hESDTrueGammaPsiPairDeltaPhi[fiCut]->Fill(deltaPhi,TruePhotonCandidate->GetPsiPair());  
   if(posDaughter->GetMother(0) <= fMCStack->GetNprimary()){
      // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
      if(fIsFromMBHeader){
         iPhotonMCInfo = 6;
         hESDTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());  
         hESDTruePrimaryConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled

      }
      // (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
   }
   else{
      if(fIsFromMBHeader){
         iPhotonMCInfo = 2;
         hESDTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
         if(fMCStack->Particle(Photon->GetMother(0))->GetMother(0) > -1 &&
            fMCStack->Particle(fMCStack->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 3122){
            hESDTrueSecondaryConvGammaFromXFromLambdaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
            iPhotonMCInfo = 5;
         }
         if(fMCStack->Particle(Photon->GetMother(0))->GetMother(0) > -1 &&
            fMCStack->Particle(fMCStack->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
            hESDTrueSecondaryConvGammaFromXFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt());
            iPhotonMCInfo = 4;
         }
         if(fMCStack->Particle(Photon->GetMother(0))->GetMother(0) > -1 &&
            fMCStack->Particle(fMCStack->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 221){
            iPhotonMCInfo = 3;
         }
      }
   }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessAODMCParticles()
{

   TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));

   // Loop over all primary MC particle
   for(Int_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {

      AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
      if (!particle) continue;
      if (!particle->IsPrimary()) continue;

      Int_t isMCFromMBHeader = -1;
      if(((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 0){
         isMCFromMBHeader
            = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent);
         if(isMCFromMBHeader == 0 && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }
      
      if(!((AliConversionCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
      if(((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kFALSE)){
         hMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
         if(particle->GetMother() >-1){ // Meson Decay Gamma
            switch((static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetMother())))->GetPdgCode()){
            case 111: // Pi0
               hMCDecayGammaPi0Pt[fiCut]->Fill(particle->Pt());
               break;
            case 113: // Rho0
               hMCDecayGammaRhoPt[fiCut]->Fill(particle->Pt());
               break;
            case 221: // Eta
               hMCDecayGammaEtaPt[fiCut]->Fill(particle->Pt());
               break;
            case 223: // Omega
               hMCDecayGammaOmegaPt[fiCut]->Fill(particle->Pt());
               break;
            case 331: // Eta'
               hMCDecayGammaEtapPt[fiCut]->Fill(particle->Pt());
               break;
            case 333: // Phi
               hMCDecayGammaPhiPt[fiCut]->Fill(particle->Pt());
               break;
            case 3212: // Sigma
               hMCDecayGammaSigmaPt[fiCut]->Fill(particle->Pt());
               break;
            }
         }
      }
      if(((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(particle,AODMCTrackArray,kTRUE)){
         Double_t rConv = 0;
         for(Int_t daughterIndex=particle->GetDaughter(0);daughterIndex<=particle->GetDaughter(1);daughterIndex++){
            AliAODMCParticle *tmpDaughter = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(daughterIndex));
            if(!tmpDaughter) continue;
            if(abs(tmpDaughter->GetPdgCode()) == 11){
               rConv = sqrt( (tmpDaughter->Xv()*tmpDaughter->Xv()) + (tmpDaughter->Yv()*tmpDaughter->Yv()) );
            }
         }
         hMCConvGammaPt[fiCut]->Fill(particle->Pt());
         if (fDoPhotonQA > 0){
            hMCConvGammaR[fiCut]->Fill(rConv);
            hMCConvGammaEta[fiCut]->Fill(particle->Eta());
         }
      }
      // Converted MC Gamma
      if(fDoMesonAnalysis){
         if(particle->GetPdgCode() == 310 && fDoMesonQA > 0){
            Double_t mesonY = 10.;
            if(particle->E() - particle->Pz() == 0 || particle->E() + particle->Pz() == 0){
               mesonY=10.-((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift();
            } else{
               mesonY = 0.5*(TMath::Log((particle->E()+particle->Pz()) / (particle->E()-particle->Pz())))
                  -((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift();
            }
            Float_t weightedK0s= 1;
            if(((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
               if (particle->Pt()>0.005){
                  weightedK0s= ((AliConversionCuts*)fCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, 0x0, fInputEvent);
                  //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
               }
            }
            hMCK0sPt[fiCut]->Fill(particle->Pt(),weightedK0s);
            hMCK0sWOWeightPt[fiCut]->Fill(particle->Pt());
            hMCK0sPtY[fiCut]->Fill(particle->Pt(),mesonY,weightedK0s);
         }         
         if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
            ->MesonIsSelectedAODMC(particle,AODMCTrackArray,((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift())){
            AliAODMCParticle* daughter0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughter(0)));
            AliAODMCParticle* daughter1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(particle->GetDaughter(1)));
            Float_t weighted= 1;
            if(((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
               if (particle->Pt()>0.005){
                  weighted= ((AliConversionCuts*)fCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, 0x0, fInputEvent);
//                   if(particle->GetPdgCode() == 221){
//                      cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
//                   }
               }
            }
            Double_t mesonY = 10.;
            if(particle->E() - particle->Pz() == 0 || particle->E() + particle->Pz() == 0){
               mesonY=10.-((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift();
            } else{
               mesonY = 0.5*(TMath::Log((particle->E()+particle->Pz()) / (particle->E()-particle->Pz())))
                  -((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift();
            }

            if(particle->GetPdgCode() == 111){
               hMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted); // All MC Pi0
               hMCPi0WOWeightPt[fiCut]->Fill(particle->Pt());
               if (fDoMesonQA > 0) hMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
            } else if(particle->GetPdgCode() == 221){
               hMCEtaPt[fiCut]->Fill(particle->Pt(),weighted); // All MC Eta
               hMCEtaWOWeightPt[fiCut]->Fill(particle->Pt());
               if (fDoMesonQA > 0) hMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
            }
            
            // Check the acceptance for both gammas
            if(((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter0,AODMCTrackArray,kFALSE) &&
               ((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelectedAODMC(daughter1,AODMCTrackArray,kFALSE)  &&
               ((AliConversionCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
               ((AliConversionCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){

               if(particle->GetPdgCode() == 111){
                  hMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted); // MC Pi0 with gamma in acc
               } else if(particle->GetPdgCode() == 221){
                  hMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted); // MC Eta with gamma in acc
               }
            }
         }
      }
   }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessMCParticles()
{
   // Loop over all primary MC particle
   for(Int_t i = 0; i < fMCStack->GetNprimary(); i++) {
      TParticle* particle = (TParticle *)fMCStack->Particle(i);
      if (!particle) continue;

      Int_t isMCFromMBHeader = -1;
      if(((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 0){
         isMCFromMBHeader
            = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent);
         if(isMCFromMBHeader == 0 && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
      }

      if(!((AliConversionCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(particle->Phi(),fEventPlaneAngle,kFALSE)) continue;
      if(((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kFALSE)){
         hMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
         if(particle->GetMother(0) >-1){ // Meson Decay Gamma
            switch(fMCStack->Particle(particle->GetMother(0))->GetPdgCode()){
            case 111: // Pi0
               hMCDecayGammaPi0Pt[fiCut]->Fill(particle->Pt());
               break;
            case 113: // Rho0
               hMCDecayGammaRhoPt[fiCut]->Fill(particle->Pt());
               break;
            case 221: // Eta
               hMCDecayGammaEtaPt[fiCut]->Fill(particle->Pt());
               break;
            case 223: // Omega
               hMCDecayGammaOmegaPt[fiCut]->Fill(particle->Pt());
               break;
            case 331: // Eta'
               hMCDecayGammaEtapPt[fiCut]->Fill(particle->Pt());
               break;
            case 333: // Phi
               hMCDecayGammaPhiPt[fiCut]->Fill(particle->Pt());
               break;
            case 3212: // Sigma
               hMCDecayGammaSigmaPt[fiCut]->Fill(particle->Pt());
               break;
            }
         }
      }
      if(((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kTRUE)){
         hMCConvGammaPt[fiCut]->Fill(particle->Pt());
         if (fDoPhotonQA > 0){
            hMCConvGammaR[fiCut]->Fill(((TParticle*)fMCStack->Particle(particle->GetFirstDaughter()))->R());
            hMCConvGammaEta[fiCut]->Fill(particle->Eta());
         }
      } // Converted MC Gamma
      if(fDoMesonAnalysis){
         if(particle->GetPdgCode() == 310 && fDoMesonQA > 0){
            Double_t mesonY = 10.;
            if(particle->Energy() - particle->Pz() == 0 || particle->Energy() + particle->Pz() == 0){
               mesonY=10.-((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift();
            } else{
               mesonY = 0.5*(TMath::Log((particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz())))
                  -((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift();
            }
            Float_t weightedK0s= 1;
            if(((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
               if (particle->Pt()>0.005){
                  weightedK0s= ((AliConversionCuts*)fCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, fMCStack, fInputEvent);
                  //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
               }
            }
            hMCK0sPt[fiCut]->Fill(particle->Pt(),weightedK0s);
            hMCK0sWOWeightPt[fiCut]->Fill(particle->Pt());
            hMCK0sPtY[fiCut]->Fill(particle->Pt(),mesonY,weightedK0s);
         }
         if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
            ->MesonIsSelectedMC(particle,fMCStack,((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift())){
            TParticle* daughter0 = (TParticle*)fMCStack->Particle(particle->GetFirstDaughter());
            TParticle* daughter1 = (TParticle*)fMCStack->Particle(particle->GetLastDaughter());

            Float_t weighted= 1;
            if(((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
               if (particle->Pt()>0.005){
                  weighted= ((AliConversionCuts*)fCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, fMCStack, fInputEvent);
//                   if(particle->GetPdgCode() == 221){
//                      cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
//                   }
               }
            }
            Double_t mesonY = 10.;
            if(particle->Energy() - particle->Pz() == 0 || particle->Energy() + particle->Pz() == 0){
               mesonY=10.-((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift();
            } else{
               mesonY = 0.5*(TMath::Log((particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz())))
                  -((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift();
            }

            if(particle->GetPdgCode() == 111){
               hMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted); // All MC Pi0
               hMCPi0WOWeightPt[fiCut]->Fill(particle->Pt());
               if (fDoMesonQA > 0) hMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
            } else if(particle->GetPdgCode() == 221){
               hMCEtaPt[fiCut]->Fill(particle->Pt(),weighted); // All MC Eta
               hMCEtaWOWeightPt[fiCut]->Fill(particle->Pt());
               if (fDoMesonQA > 0) hMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
            } 

            // Check the acceptance for both gammas
            if(((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter0,fMCStack,kFALSE) &&
               ((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter1,fMCStack,kFALSE)  &&
               ((AliConversionCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter0->Phi(),fEventPlaneAngle,kFALSE) &&
               ((AliConversionCuts*)fCutArray->At(fiCut))->InPlaneOutOfPlaneCut(daughter1->Phi(),fEventPlaneAngle,kFALSE)){

               if(particle->GetPdgCode() == 111){
                  hMCPi0InAccPt[fiCut]->Fill(particle->Pt(),weighted); // MC Pi0 with gamma in acc
               } else if(particle->GetPdgCode() == 221){
                  hMCEtaInAccPt[fiCut]->Fill(particle->Pt(),weighted); // MC Eta with gamma in acc
               }
            }
         }
      }
   }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::CalculatePi0Candidates(){

   // Conversion Gammas
   if(fGammaCandidates->GetEntries()>1){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries()-1;firstGammaIndex++){
         AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
         if (gamma0==NULL) continue;
         for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGammaCandidates->GetEntries();secondGammaIndex++){
            AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(secondGammaIndex));
            //Check for same Electron ID
            if (gamma1==NULL) continue;
            if(gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelPositive() ||
               gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelNegative() ||
               gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelPositive() ||
               gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelNegative() ) continue;

            AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
            pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);

            pi0cand->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
            if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
                ->MesonIsSelected(pi0cand,kTRUE,((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift()))){
               hESDMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
               if(pi0cand->GetAlpha()<0.1)
                  hESDMotherInvMassEalpha[fiCut]->Fill(pi0cand->M(),pi0cand->E());
               
               if (fDoMesonQA > 0){
                  if ( pi0cand->M() > 0.05 && pi0cand->M() < 0.17){
                     hESDMotherPi0PtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift());  
                     hESDMotherPi0PtAlpha[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetAlpha());  
                     hESDMotherPi0PtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle()); 
                     
                  } 
                  if ( pi0cand->M() > 0.45 && pi0cand->M() < 0.65){
                     hESDMotherEtaPtY[fiCut]->Fill(pi0cand->Pt(),pi0cand->Rapidity()-((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift());  
                     hESDMotherEtaPtAlpha[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetAlpha());  
                     hESDMotherEtaPtOpenAngle[fiCut]->Fill(pi0cand->Pt(),pi0cand->GetOpeningAngle());       
                  }
               }   
               if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoBGCalculation()){
                  Int_t zbin = 0;
                  Int_t mbin = 0;

                  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() == 0){
                     zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
                     if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                        mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
                     } else {
                        mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
                     }
                  }
                  else{
                     zbin = fBGHandlerRP[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
                     if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                        mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
                     } else {
                        mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
                     }
                  }
                  Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),(Double_t)zbin,(Double_t)mbin};
                  sESDMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
               }
               

               if(fIsMC){
                  if(fInputEvent->IsA()==AliESDEvent::Class())
                     ProcessTrueMesonCandidates(pi0cand,gamma0,gamma1);
                  if(fInputEvent->IsA()==AliAODEvent::Class())
                     ProcessTrueMesonCandidatesAOD(pi0cand,gamma0,gamma1);
               }
               if (fDoMesonQA == 2){
                  fInvMass = pi0cand->M();
                  fPt  = pi0cand->Pt();
                  if (abs(gamma0->GetDCAzToPrimVtx()) < abs(gamma1->GetDCAzToPrimVtx())){
                     fDCAzGammaMin = gamma0->GetDCAzToPrimVtx();
                     fDCAzGammaMax = gamma1->GetDCAzToPrimVtx();
                  } else {
                     fDCAzGammaMin = gamma1->GetDCAzToPrimVtx();
                     fDCAzGammaMax = gamma0->GetDCAzToPrimVtx();
                  }
                  iFlag = pi0cand->GetMesonQuality();
//                   cout << "gamma 0: " << gamma0->GetV0Index()<< "\t" << gamma0->GetPx() << "\t" << gamma0->GetPy() << "\t" <<  gamma0->GetPz() << "\t" << endl; 
//                   cout << "gamma 1: " << gamma1->GetV0Index()<< "\t"<< gamma1->GetPx() << "\t" << gamma1->GetPy() << "\t" <<  gamma1->GetPz() << "\t" << endl; 
//                    cout << "pi0: "<<fInvMass << "\t" << fPt <<"\t" << fDCAzGammaMin << "\t" << fDCAzGammaMax << "\t" << (Int_t)iFlag << "\t" << (Int_t)iMesonMCInfo <<endl;
                  if (fIsHeavyIon == 1 && fPt > 0.399 && fPt < 20. ) {
                     if (fInvMass > 0.08 && fInvMass < 0.2) tESDMesonsInvMassPtDcazMinDcazMaxFlag[fiCut]->Fill();
                     if ((fInvMass > 0.45 && fInvMass < 0.6) &&  (fPt > 0.999 && fPt < 20.) )tESDMesonsInvMassPtDcazMinDcazMaxFlag[fiCut]->Fill();
                  } else if (fPt > 0.299 && fPt < 20. )  {
                     if ( (fInvMass > 0.08 && fInvMass < 0.2) || (fInvMass > 0.45 && fInvMass < 0.6)) tESDMesonsInvMassPtDcazMinDcazMaxFlag[fiCut]->Fill();
                  }   
               }
            }
            delete pi0cand;
            pi0cand=0x0;
         }
      }
   }
}
//______________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{
   // Process True Mesons
   AliStack *MCStack = fMCEvent->Stack();
   iMesonMCInfo = 0;
   if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
      Bool_t isTruePi0 = kFALSE;
      Bool_t isTrueEta = kFALSE;
      Bool_t isTruePi0Dalitz = kFALSE;
      Bool_t isTrueEtaDalitz = kFALSE;
      Bool_t gamma0DalitzCand = kFALSE;
      Bool_t gamma1DalitzCand = kFALSE;
      Int_t gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(MCStack);
      Int_t gamma0MotherLabel = -1;
      if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
         // Daughters Gamma 0
         TParticle * negativeMC = (TParticle*)TrueGammaCandidate0->GetNegativeMCDaughter(MCStack);
         TParticle * positiveMC = (TParticle*)TrueGammaCandidate0->GetPositiveMCDaughter(MCStack);
         TParticle * gammaMC0 = (TParticle*)MCStack->Particle(gamma0MCLabel);
         if(abs(negativeMC->GetPdgCode())==11 && abs(positiveMC->GetPdgCode())==11){  // Electrons ...
            if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
               if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
                  gamma0MotherLabel=gammaMC0->GetFirstMother();
               }
            }
            if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
               gamma0DalitzCand = kTRUE;
               gamma0MotherLabel=-111;
            }
            if(gammaMC0->GetPdgCode() ==221){ // Dalitz candidate
               gamma0DalitzCand = kTRUE;
               gamma0MotherLabel=-221;
            }
         }
      }
      if(TrueGammaCandidate1->GetV0Index()<fInputEvent->GetNumberOfV0s()){
         Int_t gamma1MCLabel = TrueGammaCandidate1->GetMCParticleLabel(MCStack);
         Int_t gamma1MotherLabel = -1;
         if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
            // Daughters Gamma 1
            TParticle * negativeMC = (TParticle*)TrueGammaCandidate1->GetNegativeMCDaughter(MCStack);
            TParticle * positiveMC = (TParticle*)TrueGammaCandidate1->GetPositiveMCDaughter(MCStack);
            TParticle * gammaMC1 = (TParticle*)MCStack->Particle(gamma1MCLabel);
            if(abs(negativeMC->GetPdgCode())==11 && abs(positiveMC->GetPdgCode())==11){  // Electrons ...
               if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
                  if(gammaMC1->GetPdgCode() == 22){ // ... with Gamma Mother
                     gamma1MotherLabel=gammaMC1->GetFirstMother();
                  }
               }
               if(gammaMC1->GetPdgCode() ==111 ){ // Dalitz candidate
                 gamma1DalitzCand = kTRUE;
                 gamma1MotherLabel=-111;
               }
               if(gammaMC1->GetPdgCode() ==221){ // Dalitz candidate
                 gamma1DalitzCand = kTRUE;
                 gamma1MotherLabel=-221;
               }
            }
         }
         if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
            if(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 111){
               isTruePi0=kTRUE;
            }
            if(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetPdgCode() == 221){
               isTrueEta=kTRUE;
            }
         }
         
         //Identify Dalitz candidate
         if (gamma1DalitzCand || gamma0DalitzCand){
            if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
               if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
               if (gamma0MotherLabel == -221) isTrueEtaDalitz = kTRUE;   
            }   
            if (gamma1DalitzCand && gamma1MCLabel >=0 && gamma1MCLabel==gamma0MotherLabel){
               if (gamma1MotherLabel == -111) isTruePi0Dalitz = kTRUE;
               if (gamma1MotherLabel == -221) isTrueEtaDalitz = kTRUE;   
            }
         }
         
         
         if(isTruePi0 || isTrueEta){// True Pion or Eta
            hESDTrueMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt()); 
            if (fDoMesonQA > 0){
               if (isTruePi0){
                  if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
                     hESDTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift()); 
                     hESDTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetAlpha()); 
                     hESDTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle()); 
                  }
               } else if (isTrueEta){   
                  if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
                     hESDTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift()); 
                     hESDTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetAlpha()); 
                     hESDTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle()); 
                  }
               }
            }   
            if(gamma0MotherLabel >= MCStack->GetNprimary()){ // Secondary Meson
               Int_t secMotherLabel = ((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetMother(0);
               Float_t weightedSec= 1;
               if(((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, fMCStack, fInputEvent) && MCStack->Particle(secMotherLabel)->GetPdgCode()==310){
                  weightedSec= ((AliConversionCuts*)fCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),secMotherLabel, fMCStack, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
                  //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
               }
               hESDTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
               iMesonMCInfo = 2;
               if (secMotherLabel >-1){
                  if(MCStack->Particle(secMotherLabel)->GetPdgCode()==310){
                     iMesonMCInfo = 4;
                     hESDTrueSecondaryMotherFromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
                     if (fDoMesonQA > 0)hESDTrueK0sWithPi0DaughterMCPt[fiCut]
                                       ->Fill(MCStack->Particle(secMotherLabel)->Pt());
                  }
                  if(MCStack->Particle(secMotherLabel)->GetPdgCode()==221){
                     iMesonMCInfo = 3;
                     hESDTrueSecondaryMotherFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
                     if (fDoMesonQA > 0)hESDTrueEtaWithPi0DaughterMCPt[fiCut]
                                       ->Fill(MCStack->Particle(secMotherLabel)->Pt());
                  }
               }
            }else{ // Only primary pi0 for efficiency calculation
               iMesonMCInfo = 6;
               Float_t weighted= 1;
               if(((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, fMCStack, fInputEvent)){
                  if (((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt()>0.005){
                     weighted= ((AliConversionCuts*)fCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),gamma1MotherLabel, fMCStack, fInputEvent);
//                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
                  }
               }
               hESDTruePrimaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
               hESDTruePrimaryMotherW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
               pESDTruePrimaryMotherWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
              
               
               if (fDoMesonQA > 0){
                  if(isTruePi0){ // Only primary pi0 for resolution
                     hESDTruePrimaryPi0MCPtResolPt[fiCut]->Fill(((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt())/((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),weighted);
                  }
                  if (isTrueEta){ // Only primary eta for resolution
                     hESDTruePrimaryEtaMCPtResolPt[fiCut]->Fill(((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt())/((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),weighted);
                  }
               }
            }
         } else if(!isTruePi0 && !isTrueEta){ // Background
            if (fDoMesonQA > 0){
               if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
                  hESDTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
                  iMesonMCInfo = 1;
               } else { // No photon or without mother
                  hESDTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
               }
            }
            if( isTruePi0Dalitz || isTrueEtaDalitz ){
               // Dalitz
               iMesonMCInfo = 5;
               hESDTrueMotherDalitzInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            } else if (gamma0DalitzCand || gamma1DalitzCand){
               if (fDoMesonQA > 0)hESDTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            }   
         }
      }
   }
}
//______________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1)
{

   // Process True Mesons
   TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
   Bool_t isTruePi0 = kFALSE;
   Bool_t isTrueEta = kFALSE;
   Bool_t isTruePi0Dalitz = kFALSE;
   Bool_t isTrueEtaDalitz = kFALSE;
   Bool_t gamma0DalitzCand = kFALSE;
   Bool_t gamma1DalitzCand = kFALSE;
   
   AliAODMCParticle *positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelPositive()));
   AliAODMCParticle *negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate0->GetMCLabelNegative()));

   iMesonMCInfo = 0;
   Int_t gamma0MCLabel = -1;
   Int_t gamma0MotherLabel = -1;
   if(!positiveMC||!negativeMC)
      return;
   
   if(positiveMC->GetMother()>-1&&(negativeMC->GetMother() == positiveMC->GetMother())){
      gamma0MCLabel = positiveMC->GetMother();
   }

   if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 0
      AliAODMCParticle * gammaMC0 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MCLabel));
      if(abs(negativeMC->GetPdgCode())==11 && abs(positiveMC->GetPdgCode())==11){  // Electrons ...
         if(((positiveMC->GetMCProcessCode())) == 5 && ((negativeMC->GetMCProcessCode())) == 5){ // ... From Conversion ...     
            if(gammaMC0->GetPdgCode() == 22){ // ... with Gamma Mother
               gamma0MotherLabel=gammaMC0->GetMother();
            }
         }
         if(gammaMC0->GetPdgCode() ==111){ // Dalitz candidate
            gamma0DalitzCand = kTRUE;
            gamma0MotherLabel=-111;
         }
         if(gammaMC0->GetPdgCode() ==221){ // Dalitz candidate
            gamma0DalitzCand = kTRUE;
            gamma0MotherLabel=-221;
         }
      }
   }
   positiveMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelPositive()));
   negativeMC = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(TrueGammaCandidate1->GetMCLabelNegative()));
   
   Int_t gamma1MCLabel = -1;
   Int_t gamma1MotherLabel = -1;
   if(!positiveMC||!negativeMC)
      return;
   
   if(positiveMC->GetMother()>-1&&(negativeMC->GetMother() == positiveMC->GetMother())){
      gamma1MCLabel = positiveMC->GetMother();
   }
   if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
      // Daughters Gamma 1
      AliAODMCParticle * gammaMC1 = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MCLabel));
      if(abs(negativeMC->GetPdgCode())==11 && abs(positiveMC->GetPdgCode())==11){  // Electrons ...
         if(((positiveMC->GetMCProcessCode())) == 5 && ((negativeMC->GetMCProcessCode())) == 5){ // ... From Conversion ...     
            if(gammaMC1->GetPdgCode() == 22){ // ... with Gamma Mother
               gamma1MotherLabel=gammaMC1->GetMother();
            }
         }
         if(gammaMC1->GetPdgCode() ==111 ){ // Dalitz candidate
                 gamma1DalitzCand = kTRUE;
                 gamma1MotherLabel=-111;
         }
         if(gammaMC1->GetPdgCode() ==221){ // Dalitz candidate
            gamma1DalitzCand = kTRUE;
            gamma1MotherLabel=-221;
         }
      }
   }
   if(gamma0MotherLabel>=0 && gamma0MotherLabel==gamma1MotherLabel){
      if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 111){
         isTruePi0=kTRUE;
      }
      if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetPdgCode() == 221){
         isTrueEta=kTRUE;
      }
   }
   
   //Identify Dalitz candidate
   if (gamma1DalitzCand || gamma0DalitzCand){
      if (gamma0DalitzCand && gamma0MCLabel >=0 && gamma0MCLabel==gamma1MotherLabel){
         if (gamma0MotherLabel == -111) isTruePi0Dalitz = kTRUE;
         if (gamma0MotherLabel == -221) isTrueEtaDalitz = kTRUE;   
      }   
      if (gamma1DalitzCand && gamma1MCLabel >=0 && gamma1MCLabel==gamma0MotherLabel){
         if (gamma1MotherLabel == -111) isTruePi0Dalitz = kTRUE;
         if (gamma1MotherLabel == -221) isTrueEtaDalitz = kTRUE;   
      }
   }
            
   if(isTruePi0 || isTrueEta){// True Pion or Eta
      hESDTrueMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      if (fDoMesonQA > 0){
         if (isTruePi0){
            if ( Pi0Candidate->M() > 0.05 && Pi0Candidate->M() < 0.17){
               hESDTruePi0PtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift());
               hESDTruePi0PtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetAlpha()); 
               hESDTruePi0PtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle()); 
            }
         } else if (isTrueEta){   
            if ( Pi0Candidate->M() > 0.45 && Pi0Candidate->M() < 0.65){
               hESDTrueEtaPtY[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->Rapidity()-((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift()); 
               hESDTrueEtaPtAlpha[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetAlpha()); 
               hESDTrueEtaPtOpenAngle[fiCut]->Fill(Pi0Candidate->Pt(),Pi0Candidate->GetOpeningAngle()); 
            }
         }
      }
      if(!(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma0MotherLabel))->IsPrimary())){ // Secondary Meson
         Int_t secMotherLabel = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->GetMother();
         Float_t weightedSec= 1;
         if(((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(secMotherLabel, 0x0, fInputEvent) && static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
            weightedSec= ((AliConversionCuts*)fCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),secMotherLabel, 0x0, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
            //cout << "MC input \t"<<i << "\t" <<  particle->Pt()<<"\t"<<weighted << endl;
         }
         hESDTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
         iMesonMCInfo = 2;   
         if (secMotherLabel >-1){
            if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==310){
               iMesonMCInfo = 4;
               hESDTrueSecondaryMotherFromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
               if (fDoMesonQA > 0)hESDTrueK0sWithPi0DaughterMCPt[fiCut]
                                 ->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt());
            }
            if(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->GetPdgCode()==221){
               iMesonMCInfo = 3;
               hESDTrueSecondaryMotherFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
               if (fDoMesonQA > 0)hESDTrueEtaWithPi0DaughterMCPt[fiCut]
                                 ->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(secMotherLabel))->Pt());
            }
         }
      }else{ // Only primary pi0 for efficiency calculation
         Float_t weighted= 1;
         iMesonMCInfo = 6;
         if(((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, 0x0, fInputEvent)){
            if (static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt()>0.005){
               weighted= ((AliConversionCuts*)fCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),gamma1MotherLabel, 0x0, fInputEvent);
               //                      cout << "rec \t " <<gamma1MotherLabel << "\t" <<  weighted << endl;
            }
         }
         hESDTruePrimaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
         hESDTruePrimaryMotherW0WeightingInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
         pESDTruePrimaryMotherWeightsInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
              
         if (fDoMesonQA > 0){
            if(isTruePi0){ // Only primary pi0 for resolution
               hESDTruePrimaryPi0MCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                                                          (Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted);
               
            }
            if (isTrueEta){ // Only primary eta for resolution
               hESDTruePrimaryEtaMCPtResolPt[fiCut]->Fill(static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),
                                                          (Pi0Candidate->Pt()-static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt())/static_cast<AliAODMCParticle*>(AODMCTrackArray->At(gamma1MotherLabel))->Pt(),weighted);
            }
         }
      }
   } else if(!isTruePi0 && !isTrueEta) { // Background
      if (fDoMesonQA > 0){
         if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
            hESDTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            iMesonMCInfo = 1;
         } else { // No photon or without mother
            hESDTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
         }
      }
      if( isTruePi0Dalitz || isTrueEtaDalitz ){
         // Dalitz
         iMesonMCInfo = 5;
         hESDTrueMotherDalitzInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      } else if (gamma0DalitzCand || gamma1DalitzCand){
         if (fDoMesonQA > 0)hESDTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
      }
   }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::CalculateBackground(){

   Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
   Int_t mbin = 0;

   if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
   } else {
      mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
   }

   if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseRotationMethod()){

      for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
         AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
         for(Int_t iCurrent2=iCurrent+1;iCurrent2<fGammaCandidates->GetEntries();iCurrent2++){
            for(Int_t nRandom=0;nRandom<((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents();nRandom++){
               AliAODConversionPhoton currentEventGoodV02 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent2));

               if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoBGProbability()){
                  AliAODConversionMother *backgroundCandidateProb = new AliAODConversionMother(&currentEventGoodV0,&currentEventGoodV02);
                  Double_t massBGprob = backgroundCandidateProb->M();
                  if(massBGprob>0.1 && massBGprob<0.14){
                     if(fRandom.Rndm()>fBGHandler[fiCut]->GetBGProb(zbin,mbin)){
                        delete backgroundCandidateProb;
                        continue;
                     }
                  }
                  delete backgroundCandidateProb;
                  backgroundCandidateProb = 0x0;
               }

               RotateParticle(&currentEventGoodV02);
               AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&currentEventGoodV02);
               backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
               if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
                   ->MesonIsSelected(backgroundCandidate,kFALSE,((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift()))){
                  hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
                  Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                  sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
               }
               delete backgroundCandidate;
               backgroundCandidate = 0x0;
            }
         }
      }
   }else{
      AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertex = NULL;

      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
         for(Int_t nEventsInBG=0;nEventsInBG<fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
            AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
            if(fMoveParticleAccordingToVertex == kTRUE){
               bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
            }

            for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
               AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
               for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
                  AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
                  if(fMoveParticleAccordingToVertex == kTRUE){
                     MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
                  }
                  if(((AliConversionCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
                     RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
                  }

                  AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
                  backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
                  if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
                      ->MesonIsSelected(backgroundCandidate,kFALSE,((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift()))){
                     hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
                     Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                     sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
                  }
                  delete backgroundCandidate;
                  backgroundCandidate = 0x0;
               }
            }
         }
      }
      else{
         for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
            AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
            if(previousEventV0s){
               if(fMoveParticleAccordingToVertex == kTRUE){
                  bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
               }
               for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){
                  AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));
                  for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){

                     AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));

                     if(fMoveParticleAccordingToVertex == kTRUE){
                        MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
                     }
                     if(((AliConversionCuts*)fCutArray->At(fiCut))->GetInPlaneOutOfPlaneCut() != 0){
                        RotateParticleAccordingToEP(&previousGoodV0,bgEventVertex->fEP,fEventPlaneAngle);
                     }


                     AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
                     backgroundCandidate->CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
                     if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
                         ->MesonIsSelected(backgroundCandidate,kFALSE,((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift()))){
                        hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
                        Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                        sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
                     }
                     delete backgroundCandidate;
                     backgroundCandidate = 0x0;
                  }
               }
            }
         }
      }
   }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::CalculateBackgroundRP(){

   Int_t zbin= fBGHandlerRP[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
   Int_t mbin = 0;
   if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fV0Reader->GetNumberOfPrimaryTracks());
   } else {
      mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
   }


   //Rotation Method
   if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseRotationMethod()){
      // Correct for the number of rotations
      // BG is for rotation the same, except for factor NRotations
      Double_t weight=1./Double_t(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents());

      for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries();firstGammaIndex++){

         AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
         if (gamma0==NULL) continue;
         for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGammaCandidates->GetEntries();secondGammaIndex++){
            AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(secondGammaIndex));
            if (gamma1 == NULL) continue;
            if(!((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelected(gamma1,fInputEvent))continue;
            for(Int_t nRandom=0;nRandom<((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->GetNumberOfBGEvents();nRandom++){

               RotateParticle(gamma1);

               AliAODConversionMother backgroundCandidate(gamma0,gamma1);
               backgroundCandidate.CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
               if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
                  ->MesonIsSelected(&backgroundCandidate,kFALSE,((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift())){
                  hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt());
                  Double_t sparesFill[4] = {backgroundCandidate.M(),backgroundCandidate.Pt(),(Double_t)zbin,(Double_t)mbin};
                  sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,weight);
               }
            }
         }
      }
   }
   else{
      // Do Event Mixing
      for(Int_t nEventsInBG=0;nEventsInBG <fBGHandlerRP[fiCut]->GetNBGEvents(fGammaCandidates,fInputEvent);nEventsInBG++){

         AliGammaConversionPhotonVector *previousEventGammas = fBGHandlerRP[fiCut]->GetBGGoodGammas(fGammaCandidates,fInputEvent,nEventsInBG);

         if(previousEventGammas){
            // test weighted background
            Double_t weight=1.0;
            // Correct for the number of eventmixing:
            // N gammas -> (N-1) + (N-2) +(N-3) ...+ (N-(N-1))  using sum formula sum(i)=N*(N-1)/2  -> N*(N-1)/2
            // real combinations (since you cannot combine a photon with its own)
            // but BG leads to N_{a}*N_{b} combinations
            weight*=0.5*(Double_t(fGammaCandidates->GetEntries()-1))/Double_t(previousEventGammas->size());

            for(Int_t iCurrent=0;iCurrent<fGammaCandidates->GetEntries();iCurrent++){

               AliAODConversionPhoton *gamma0 = (AliAODConversionPhoton*)(fGammaCandidates->At(iCurrent));

               for(UInt_t iPrevious=0;iPrevious<previousEventGammas->size();iPrevious++){

                  AliAODConversionPhoton *gamma1 = (AliAODConversionPhoton*)(previousEventGammas->at(iPrevious));

                  AliAODConversionMother backgroundCandidate(gamma0,gamma1);
                  backgroundCandidate.CalculateDistanceOfClossetApproachToPrimVtx(fInputEvent->GetPrimaryVertex());
                  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))
                     ->MesonIsSelected(&backgroundCandidate,kFALSE,((AliConversionCuts*)fCutArray->At(fiCut))->GetEtaShift())){
                     hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate.M(),backgroundCandidate.Pt());
                     Double_t sparesFill[4] = {backgroundCandidate.M(),backgroundCandidate.Pt(),(Double_t)zbin,(Double_t)mbin};
                     sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,weight);
                  }
               }
            }
         }
      }
   }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::RotateParticle(AliAODConversionPhoton *gamma){
   Int_t fNDegreesPMBackground= ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->NDegreesRotation();
   Double_t nRadiansPM = fNDegreesPMBackground*TMath::Pi()/180;
   Double_t rotationValue = fRandom.Rndm()*2*nRadiansPM + TMath::Pi()-nRadiansPM;
   gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP){

   previousEventEP=previousEventEP+TMath::Pi();
   thisEventEP=thisEventEP+TMath::Pi();
   Double_t rotationValue= thisEventEP-previousEventEP;
   gamma->RotateZ(rotationValue);
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
   //see header file for documentation

   Double_t dx = vertex->fX - fInputEvent->GetPrimaryVertex()->GetX();
   Double_t dy = vertex->fY - fInputEvent->GetPrimaryVertex()->GetY();
   Double_t dz = vertex->fZ - fInputEvent->GetPrimaryVertex()->GetZ();

   Double_t movedPlace[3] = {particle->GetConversionX() - dx,particle->GetConversionY() - dy,particle->GetConversionZ() - dz};
   particle->SetConversionPoint(movedPlace);
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::UpdateEventByEventData(){
   //see header file for documentation
   if(fGammaCandidates->GetEntries() >0 ){
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
         fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fV0Reader->GetNumberOfPrimaryTracks(),fEventPlaneAngle);
      }
      else{ // means we use #V0s for multiplicity
         fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries(),fEventPlaneAngle);
      }
   }
}


//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::FillPhotonCombinatorialBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode[])
{
   // Combinatorial Bck = 0 ee, 1 ep,i 2 ek, 3 ep, 4 emu, 5 pipi, 6 pik, 7 pip, 8 pimu, 9 kk, 10 kp, 11 kmu, 12 pp, 13 pmu, 14 mumu, 15 Rest
   if(pdgCode[0]==11   && pdgCode[1]==11){if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),0);}
   else if( (pdgCode[0]==11   && pdgCode[1]==211) || (pdgCode[0]==211  && pdgCode[1]==11) )
      {if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),1);}
   else if( (pdgCode[0]==11   && pdgCode[1]==321) || (pdgCode[0]==321  && pdgCode[1]==11) )
      {if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),2);}
   else if( (pdgCode[0]==11   && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==11) )
      {if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),3);}
   else if( (pdgCode[0]==11   && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==11) )
      {if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),4);}
   else if(  pdgCode[0]==211  && pdgCode[1]==211 ){if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),5);}
   else if( (pdgCode[0]==211  && pdgCode[1]==321) || (pdgCode[0]==321  && pdgCode[1]==211) )
      {if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),6);}
   else if( (pdgCode[0]==211  && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==211) )
      {if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),7);}
   else if( (pdgCode[0]==211  && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==211) )
      {if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),8);}
   else if(  pdgCode[0]==321  && pdgCode[1]==321 ){if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),9);}
   else if( (pdgCode[0]==321  && pdgCode[1]==2212) || (pdgCode[0]==2212 && pdgCode[1]==321) )
      {if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),10);}
   else if( (pdgCode[0]==321  && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==321) )
      {if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),11);}
   else if(  pdgCode[0]==2212   && pdgCode[1]==2212  ){if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),12);}
   else if( (pdgCode[0]==2212  && pdgCode[1]==13) || (pdgCode[0]==13   && pdgCode[1]==2212) )
      {if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),13);}
   else if(  pdgCode[0]==13   && pdgCode[1]==13  ){if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),14);}
   else {if(fIsFromMBHeader)hESDCombinatorialPt[fiCut]->Fill(TruePhotonCandidate->Pt(),15);}
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::RelabelAODPhotonCandidates(Bool_t mode){

   // Relabeling For AOD Event
   // ESDiD -> AODiD
   // MCLabel -> AODMCLabel
   
   if(mode){
      fMCStackPos = new Int_t[fReaderGammas->GetEntries()];
      fMCStackNeg = new Int_t[fReaderGammas->GetEntries()];
      fESDArrayPos = new Int_t[fReaderGammas->GetEntries()];
      fESDArrayNeg = new Int_t[fReaderGammas->GetEntries()];
   }
   
   for(Int_t iGamma = 0;iGamma<fReaderGammas->GetEntries();iGamma++){
      AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(iGamma);
      if(!PhotonCandidate) continue;
      if(!mode){// Back to ESD Labels
         PhotonCandidate->SetMCLabelPositive(fMCStackPos[iGamma]);
         PhotonCandidate->SetMCLabelNegative(fMCStackNeg[iGamma]);
         PhotonCandidate->SetLabelPositive(fESDArrayPos[iGamma]);
         PhotonCandidate->SetLabelNegative(fESDArrayNeg[iGamma]);
         continue;
      }
      fMCStackPos[iGamma] =  PhotonCandidate->GetMCLabelPositive();
      fMCStackNeg[iGamma] =  PhotonCandidate->GetMCLabelNegative();
      fESDArrayPos[iGamma] = PhotonCandidate->GetTrackLabelPositive();
      fESDArrayNeg[iGamma] = PhotonCandidate->GetTrackLabelNegative();

      Bool_t AODLabelPos = kFALSE;
      Bool_t AODLabelNeg = kFALSE;

      for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
         AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
         if(!AODLabelPos){
            if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive() ){
               PhotonCandidate->SetMCLabelPositive(abs(tempDaughter->GetLabel()));
               PhotonCandidate->SetLabelPositive(i);
               AODLabelPos = kTRUE;
            }
         }
         if(!AODLabelNeg){
            if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
               PhotonCandidate->SetMCLabelNegative(abs(tempDaughter->GetLabel()));
               PhotonCandidate->SetLabelNegative(i);
               AODLabelNeg = kTRUE;
            }
         }
         if(AODLabelNeg && AODLabelPos){
            break;
         }
      }
      if(!AODLabelPos || !AODLabelNeg){
         cout<<"WARNING!!! AOD TRACKS NOT FOUND FOR"<<endl;
      }
   }
   
   
   if(!mode){
      delete[] fMCStackPos;
      delete[] fMCStackNeg;
      delete[] fESDArrayPos;
      delete[] fESDArrayNeg;
   }
}

void AliAnalysisTaskGammaConvV1::SetLogBinningXTH2(TH2* histoRebin){
   TAxis *axisafter = histoRebin->GetXaxis(); 
   Int_t bins = axisafter->GetNbins();
   Double_t from = axisafter->GetXmin();
   Double_t to = axisafter->GetXmax();
   Double_t *newbins = new Double_t[bins+1];
   newbins[0] = from;
   Double_t factor = TMath::Power(to/from, 1./bins);
   for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
   axisafter->Set(bins, newbins);
   delete [] newbins;

}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::Terminate(const Option_t *)
{

   //fOutputContainer->Print(); // Will crash on GRID
}
