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
#include "AliESDtrackCuts.h"
#include "AliKFVertex.h"
#include "AliV0ReaderV1.h"
#include "AliGenCocktailEventHeader.h"
#include "AliConversionAODBGHandlerRP.h"

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
   fMotherRapList(NULL),
   fTrueList(NULL),
   fTrueMotherRapList(NULL),
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
   hESDMotherInvMassPt(NULL),
   sESDMotherInvMassPtZM(NULL),
   sESDMotherInvMassPtY(NULL),
   hESDMotherBackInvMassPt(NULL),
   sESDMotherBackInvMassPtZM(NULL),
   hESDMotherInvMassEalpha(NULL),
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
   hMCConvGammaRSPt(NULL),
   hMCConvGammaRSR(NULL),
   hMCConvGammaRSEta(NULL),
   hMCPi0Pt(NULL),
   hMCEtaPt(NULL),
   hMCPi0InAccPt(NULL),
   hMCEtaInAccPt(NULL),
   hMCPi0PtY(NULL),
   hMCEtaPtY(NULL),
   hESDTrueMotherInvMassPt(NULL),
   hESDTruePrimaryMotherInvMassPt(NULL),
   hESDTruePrimaryPi0MCPtResolPt(NULL),
   hESDTruePrimaryEtaMCPtResolPt(NULL),
   sESDTruePrimaryMotherInvMassPtY(NULL),
   hESDTrueSecondaryMotherInvMassPt(NULL),
   hESDTrueSecondaryMotherFromK0sInvMassPt(NULL),
   hESDTrueK0sWithPi0DaughterMCPt(NULL),
   hESDTrueSecondaryMotherFromEtaInvMassPt(NULL),
   hESDTrueEtaWithPi0DaughterMCPt(NULL),
   hESDTrueBckGGInvMassPt(NULL),
   hESDTrueBckContInvMassPt(NULL),
   hESDTrueMotherDalitzInvMassPt(NULL),
   hESDTrueConvGammaPt(NULL),
   hESDCombinatorialPt(NULL),
   hESDTruePrimaryConvGammaPt(NULL),
   hESDTruePrimaryConvGammaR(NULL),
   hESDTruePrimaryConvGammaEta(NULL),
   hESDTruePrimaryConvGammaESDPtMCPt(NULL),
   hESDTruePrimaryConvGammaRSESDPtMCPt(NULL),
   hESDTrueSecondaryConvGammaPt(NULL),
   hESDTrueSecondaryConvGammaR(NULL),
   hESDTrueSecondaryConvGammaFromXFromK0sPt(NULL),
   hNEvents(NULL),
   hNGoodESDTracks(NULL),
   hNGammaCandidates(NULL),
   hNV0Tracks(NULL),
   fRandom(0),
   fnGammaCandidates(0),
   fUnsmearedPx(NULL),
   fUnsmearedPy(NULL),
   fUnsmearedPz(NULL),
   fUnsmearedE(NULL),
   fnCuts(0),
   fiCut(0),
   fNumberOfESDTracks(0),
   fMoveParticleAccordingToVertex(kTRUE),
   fIsHeavyIon(kFALSE),
   fDoMesonAnalysis(kTRUE),
   fDoMesonQA(kFALSE),
   fDoPhotonQA(kFALSE),
   fIsFromMBHeader(kTRUE)
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
   fMotherRapList(NULL),
   fTrueList(NULL),
   fTrueMotherRapList(NULL),
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
   hESDMotherInvMassPt(NULL),
   sESDMotherInvMassPtZM(NULL),
   sESDMotherInvMassPtY(NULL),
   hESDMotherBackInvMassPt(NULL),
   sESDMotherBackInvMassPtZM(NULL),
   hESDMotherInvMassEalpha(NULL),
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
   hMCConvGammaRSPt(NULL),
   hMCConvGammaRSR(NULL),
   hMCConvGammaRSEta(NULL),
   hMCPi0Pt(NULL),
   hMCEtaPt(NULL),
   hMCPi0InAccPt(NULL),
   hMCEtaInAccPt(NULL),
   hMCPi0PtY(NULL),
   hMCEtaPtY(NULL),
   hESDTrueMotherInvMassPt(NULL),
   hESDTruePrimaryMotherInvMassPt(NULL),
   hESDTruePrimaryPi0MCPtResolPt(NULL),
   hESDTruePrimaryEtaMCPtResolPt(NULL),
   sESDTruePrimaryMotherInvMassPtY(NULL),
   hESDTrueSecondaryMotherInvMassPt(NULL),
   hESDTrueSecondaryMotherFromK0sInvMassPt(NULL),
   hESDTrueK0sWithPi0DaughterMCPt(NULL),
   hESDTrueSecondaryMotherFromEtaInvMassPt(NULL),
   hESDTrueEtaWithPi0DaughterMCPt(NULL),
   hESDTrueBckGGInvMassPt(NULL),
   hESDTrueBckContInvMassPt(NULL),
   hESDTrueMotherDalitzInvMassPt(NULL),
   hESDTrueConvGammaPt(NULL),
   hESDCombinatorialPt(NULL),
   hESDTruePrimaryConvGammaPt(NULL),
   hESDTruePrimaryConvGammaR(NULL),
   hESDTruePrimaryConvGammaEta(NULL),
   hESDTruePrimaryConvGammaESDPtMCPt(NULL),
   hESDTruePrimaryConvGammaRSESDPtMCPt(NULL),
   hESDTrueSecondaryConvGammaPt(NULL),
   hESDTrueSecondaryConvGammaR(NULL),
   hESDTrueSecondaryConvGammaFromXFromK0sPt(NULL),
   hNEvents(NULL),
   hNGoodESDTracks(NULL),
   hNGammaCandidates(NULL),
   hNV0Tracks(NULL),
   fRandom(0),
   fnGammaCandidates(0),
   fUnsmearedPx(NULL),
   fUnsmearedPy(NULL),
   fUnsmearedPz(NULL),
   fUnsmearedE(NULL),
   fnCuts(0),
   fiCut(0),
   fNumberOfESDTracks(0),
   fMoveParticleAccordingToVertex(kTRUE),
   fIsHeavyIon(kFALSE),
   fDoMesonAnalysis(kTRUE),
   fDoMesonQA(kFALSE),
   fDoPhotonQA(kFALSE),
   fIsFromMBHeader(kTRUE)
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

   Double_t *zBinLimitsArray= new Double_t[9] ;
   zBinLimitsArray[0] = -50.00;
   zBinLimitsArray[1] = -3.375;
   zBinLimitsArray[2] = -1.605;
   zBinLimitsArray[3] = -0.225;
   zBinLimitsArray[4] = 1.065;
   zBinLimitsArray[5] = 2.445;
   zBinLimitsArray[6] = 4.245;
   zBinLimitsArray[7] = 50.00;
   zBinLimitsArray[8] = 1000.00;

   Double_t *multiplicityBinLimitsArrayTracks= new Double_t[6];
   multiplicityBinLimitsArrayTracks[0] = 0;
   multiplicityBinLimitsArrayTracks[1] = 8.5;
   multiplicityBinLimitsArrayTracks[2] = 16.5;
   multiplicityBinLimitsArrayTracks[3] = 27.5;
   multiplicityBinLimitsArrayTracks[4] = 41.5;
   multiplicityBinLimitsArrayTracks[5] = 200.;
   if(fIsHeavyIon){
      multiplicityBinLimitsArrayTracks[0] = 0;
      multiplicityBinLimitsArrayTracks[1] = 200.;
      multiplicityBinLimitsArrayTracks[2] = 500.;
      multiplicityBinLimitsArrayTracks[3] = 1000.;
      multiplicityBinLimitsArrayTracks[4] = 1500.;
      multiplicityBinLimitsArrayTracks[5] = 5000.;
   }

   Double_t *multiplicityBinLimitsArrayV0s= new Double_t[5];
   multiplicityBinLimitsArrayV0s[0] = 2;
   multiplicityBinLimitsArrayV0s[1] = 3;
   multiplicityBinLimitsArrayV0s[2] = 4;
   multiplicityBinLimitsArrayV0s[3] = 5;
   multiplicityBinLimitsArrayV0s[4] = 9999;
   if(fIsHeavyIon){
      multiplicityBinLimitsArrayV0s[0] = 2;
      multiplicityBinLimitsArrayV0s[1] = 10;
      multiplicityBinLimitsArrayV0s[2] = 30;
      multiplicityBinLimitsArrayV0s[3] = 50;
      multiplicityBinLimitsArrayV0s[4] = 9999;
   }

   const Int_t nDim = 4;
   Int_t nBins[nDim] = {800,250,8,5};
   Double_t xMin[nDim] = {0,0, 0,0};
   Double_t xMax[nDim] = {0.8,25,8,5};

   sESDMotherInvMassPtZM = new THnSparseF*[fnCuts];
   sESDMotherBackInvMassPtZM = new THnSparseF*[fnCuts];

   fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];
   fBGHandlerRP = new AliConversionAODBGHandlerRP*[fnCuts];
   for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      if (((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->DoBGCalculation()){
         TString cutstring = ((AliConversionCuts*)fCutArray->At(iCut))->GetCutNumber();
         TString cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
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
            if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity()){
               fBGHandler[iCut] = new AliGammaConversionAODBGHandler(9,6,((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents());
               fBGHandler[iCut]->Initialize(zBinLimitsArray, multiplicityBinLimitsArrayTracks);
               fBGHandlerRP[iCut] = NULL;
            }
            else{
               fBGHandler[iCut] = new AliGammaConversionAODBGHandler(9,5,((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetNumberOfBGEvents());
               fBGHandler[iCut]->Initialize(zBinLimitsArray, multiplicityBinLimitsArrayV0s);
               fBGHandlerRP[iCut] = NULL;
            }
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
   hESDConvGammaPt = new TH1F*[fnCuts];
   if (fDoPhotonQA){
      hESDConvGammaR = new TH1F*[fnCuts];
   }
   const Int_t nDim = 3;
   Int_t nBins[nDim] = {800,250,40};
   Double_t xMin[nDim] = {0,0, -1};
   Double_t xMax[nDim] = {0.8,25,1};

   if(fDoMesonAnalysis){
      hESDMotherInvMassPt = new TH2F*[fnCuts];
      hESDMotherBackInvMassPt = new TH2F*[fnCuts];
      hESDMotherInvMassEalpha = new TH2F*[fnCuts];
      if (fDoMesonQA){
         fMotherRapList = new TList*[fnCuts];
         sESDMotherInvMassPtY = new THnSparseF*[fnCuts];
      }
   }

   for(Int_t iCut = 0; iCut<fnCuts;iCut++){

      TString cutstring = ((AliConversionCuts*)fCutArray->At(iCut))->GetCutNumber();
      TString cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

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
      hNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
      fESDList[iCut]->Add(hNEvents[iCut]);
      if(fIsHeavyIon) hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",3000,0,3000);
      else hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
      fESDList[iCut]->Add(hNGoodESDTracks[iCut]);
      if(fIsHeavyIon) hNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",100,0,100);
      else hNGammaCandidates[iCut] = new TH1I("GammaCandidates","GammaCandidates",50,0,50);
      fESDList[iCut]->Add(hNGammaCandidates[iCut]);
      if(fIsHeavyIon) hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",30000,0,30000);
      else hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",2000,0,2000);
      fESDList[iCut]->Add(hNV0Tracks[iCut]);
      hESDConvGammaPt[iCut] = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",250,0,25);
      fESDList[iCut]->Add(hESDConvGammaPt[iCut]);

      if (fDoPhotonQA){
         hESDConvGammaR[iCut] = new TH1F("ESD_ConvGamma_R","ESD_ConvGamma_R",800,0,200);
         fESDList[iCut]->Add(hESDConvGammaR[iCut]);
      }

      if(fDoMesonAnalysis){
         hESDMotherInvMassPt[iCut] = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",800,0,0.8,250,0,25);
         fESDList[iCut]->Add(hESDMotherInvMassPt[iCut]);
         hESDMotherBackInvMassPt[iCut] = new TH2F("ESD_Background_InvMass_Pt","ESD_Background_InvMass_Pt",800,0,0.8,250,0,25);
         fESDList[iCut]->Add(hESDMotherBackInvMassPt[iCut]);
         hESDMotherInvMassEalpha[iCut] = new TH2F("ESD_Mother_InvMass_vs_E_alpha","ESD_Mother_InvMass_vs_E_alpha",800,0,0.8,250,0,25);
         fESDList[iCut]->Add(hESDMotherInvMassEalpha[iCut]);
         if (fDoMesonQA){
            fMotherRapList[iCut] = new TList();
            fMotherRapList[iCut]->SetName(Form("%s_%s Mother Y histograms",cutstring.Data(),cutstringMeson.Data()));
            fMotherRapList[iCut]->SetOwner(kTRUE);
            fCutFolder[iCut]->Add(fMotherRapList[iCut]);
            sESDMotherInvMassPtY[iCut] = new THnSparseF("Mother_InvMass_Pt_Y","Mother_InvMass_Pt_Y",nDim,nBins,xMin,xMax);
            fMotherRapList[iCut]->Add(sESDMotherInvMassPtY[iCut]);
         }
      }


   }
   if(fDoMesonAnalysis){
      InitBack(); // Init Background Handler
   }

   if(MCEvent()){
      // MC Histogramms
      fMCList = new TList*[fnCuts];
      // True Histogramms
      fTrueList = new TList*[fnCuts];
      // Selected Header List
      fHeaderNameList = new TList*[fnCuts];

      hMCAllGammaPt = new TH1F*[fnCuts];
      hMCDecayGammaPi0Pt = new TH1F*[fnCuts];
      hMCDecayGammaRhoPt = new TH1F*[fnCuts];
      hMCDecayGammaEtaPt = new TH1F*[fnCuts];
      hMCDecayGammaOmegaPt = new TH1F*[fnCuts];
      hMCDecayGammaEtapPt = new TH1F*[fnCuts];
      hMCDecayGammaPhiPt = new TH1F*[fnCuts];
      hMCDecayGammaSigmaPt = new TH1F*[fnCuts];
      hMCConvGammaPt = new TH1F*[fnCuts];
      hMCConvGammaRSPt = new TH1F*[fnCuts];
      hESDTrueConvGammaPt = new TH1F*[fnCuts];

      hESDCombinatorialPt = new TH2F*[fnCuts];
      hESDTruePrimaryConvGammaPt = new TH1F*[fnCuts];
      hESDTruePrimaryConvGammaESDPtMCPt = new TH2F*[fnCuts];
      hESDTruePrimaryConvGammaRSESDPtMCPt = new TH2F*[fnCuts];
      hESDTrueSecondaryConvGammaPt = new TH1F*[fnCuts];

      hESDTrueSecondaryConvGammaFromXFromK0sPt = new TH1F*[fnCuts];

      if (fDoPhotonQA){
         hMCConvGammaR = new TH1F*[fnCuts];
         hMCConvGammaEta = new TH1F*[fnCuts];
         hMCConvGammaRSR = new TH1F*[fnCuts];
         hMCConvGammaRSEta = new TH1F*[fnCuts];
         hESDTruePrimaryConvGammaR = new TH1F*[fnCuts];
         hESDTruePrimaryConvGammaEta = new TH1F*[fnCuts];
         hESDTrueSecondaryConvGammaR = new TH1F*[fnCuts];
      }

      if(fDoMesonAnalysis){
         hMCPi0Pt = new TH1F*[fnCuts];
         hMCEtaPt = new TH1F*[fnCuts];
         hMCPi0InAccPt = new TH1F*[fnCuts];
         hMCEtaInAccPt = new TH1F*[fnCuts];

         hESDTrueMotherInvMassPt = new TH2F*[fnCuts];
         hESDTruePrimaryMotherInvMassPt = new TH2F*[fnCuts];
         hESDTrueSecondaryMotherInvMassPt = new TH2F*[fnCuts];
         hESDTrueSecondaryMotherFromK0sInvMassPt = new TH2F*[fnCuts];
         hESDTrueSecondaryMotherFromEtaInvMassPt = new TH2F*[fnCuts];
         if (fDoMesonQA){
            hMCPi0PtY = new TH2F*[fnCuts];
            hMCEtaPtY = new TH2F*[fnCuts];
            hESDTruePrimaryPi0MCPtResolPt = new TH2F*[fnCuts];
            hESDTruePrimaryEtaMCPtResolPt = new TH2F*[fnCuts];
            hESDTrueK0sWithPi0DaughterMCPt = new TH1F*[fnCuts];
            hESDTrueEtaWithPi0DaughterMCPt = new TH1F*[fnCuts];
            hESDTrueBckGGInvMassPt = new TH2F*[fnCuts];
            hESDTrueBckContInvMassPt = new TH2F*[fnCuts];
            hESDTrueMotherDalitzInvMassPt = new TH2F*[fnCuts];
            fTrueMotherRapList = new TList*[fnCuts];
            sESDTruePrimaryMotherInvMassPtY = new THnSparseF*[fnCuts];
         }
      }

      for(Int_t iCut = 0; iCut<fnCuts;iCut++){
         TString cutstring = ((AliConversionCuts*)fCutArray->At(iCut))->GetCutNumber();
         TString cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();

         fMCList[iCut] = new TList();
         fMCList[iCut]->SetName(Form("%s_%s MC histograms",cutstring.Data(),cutstringMeson.Data()));
         fMCList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fMCList[iCut]);

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
         hMCConvGammaRSPt[iCut] = new TH1F("MC_ConvGamma_RS_Pt","MC_ConvGamma_RS_Pt",250,0,25);
         fMCList[iCut]->Add(hMCConvGammaRSPt[iCut]);

         if (fDoPhotonQA){
            hMCConvGammaR[iCut] = new TH1F("MC_ConvGamma_R","MC_ConvGamma_R",800,0,200);
            fMCList[iCut]->Add(hMCConvGammaR[iCut]);
            hMCConvGammaEta[iCut] = new TH1F("MC_ConvGamma_Eta","MC_ConvGamma_Eta",100,-4,4);
            fMCList[iCut]->Add(hMCConvGammaEta[iCut]);
            hMCConvGammaRSR[iCut] = new TH1F("MC_ConvGamma_RS_R","MC_ConvGamma_RS_R",800,0,200);
            fMCList[iCut]->Add(hMCConvGammaRSR[iCut]);
            hMCConvGammaRSEta[iCut] = new TH1F("MC_ConvGamma_RS_Eta","MC_ConvGamma_RS_Eta",100,-4,4);
            fMCList[iCut]->Add(hMCConvGammaRSEta[iCut]);
         }

         if(fDoMesonAnalysis){
            hMCPi0Pt[iCut] = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",250,0,25);
            hMCPi0Pt[iCut]->Sumw2();
            fMCList[iCut]->Add(hMCPi0Pt[iCut]);
            hMCEtaPt[iCut] = new TH1F("MC_Eta_Pt","MC_Eta_Pt",250,0,25);
            hMCEtaPt[iCut]->Sumw2();
            fMCList[iCut]->Add(hMCEtaPt[iCut]);
            hMCPi0InAccPt[iCut] = new TH1F("MC_Pi0InAcc_Pt","MC_Pi0InAcc_Pt",250,0,25);
            hMCPi0InAccPt[iCut]->Sumw2();
            fMCList[iCut]->Add(hMCPi0InAccPt[iCut]);
            hMCEtaInAccPt[iCut] = new TH1F("MC_EtaInAcc_Pt","MC_EtaInAcc_Pt",250,0,25);
            hMCEtaInAccPt[iCut]->Sumw2();
            fMCList[iCut]->Add(hMCEtaInAccPt[iCut]);
            if (fDoMesonQA){
               hMCPi0PtY[iCut] = new TH2F("MC_Pi0_Pt_Y","MC_Pi0_Pt_Y",250,0,25,20,-1,1);
               hMCPi0PtY[iCut]->Sumw2();
               fMCList[iCut]->Add(hMCPi0PtY[iCut]);
               hMCEtaPtY[iCut] = new TH2F("MC_Eta_Pt_Y","MC_Eta_Pt_Y",250,0,25,20,-1,1);
               hMCEtaPtY[iCut]->Sumw2();
               fMCList[iCut]->Add(hMCEtaPtY[iCut]);
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
         hESDTruePrimaryConvGammaESDPtMCPt[iCut] = new TH2F("ESD_TruePrimaryConvGammaESD_PtMCPt", "ESD_TruePrimaryConvGammaESD_PtMCPt",250,0,25,250,0,25);
         fTrueList[iCut]->Add(hESDTruePrimaryConvGammaESDPtMCPt[iCut]);
         hESDTruePrimaryConvGammaRSESDPtMCPt[iCut]
            = new TH2F("ESD_TruePrimaryConvGammaESD_RS_PtMCPt", "ESD_TruePrimaryConvGammaESD_RS_PtMCPt",250,0,25,250,0,25);
         fTrueList[iCut]->Add(hESDTruePrimaryConvGammaRSESDPtMCPt[iCut]);

         if (fDoPhotonQA){
            hESDTruePrimaryConvGammaR[iCut] = new TH1F("ESD_TruePrimaryConvGamma_R","ESD_TruePrimaryConvGamma_R",800,0,200);
            fTrueList[iCut]->Add(hESDTruePrimaryConvGammaR[iCut]);
            hESDTrueSecondaryConvGammaR[iCut] = new TH1F("ESD_TrueSecondaryConvGamma_R","ESD_TrueSecondaryConvGamma_R",800,0,200);
            fTrueList[iCut]->Add(hESDTrueSecondaryConvGammaR[iCut]);
            hESDTruePrimaryConvGammaEta[iCut] = new TH1F("ESD_TruePrimaryConvGamma_Eta","ESD_TruePrimaryConvGamma_Eta",100,-4,4);
            fTrueList[iCut]->Add(hESDTruePrimaryConvGammaEta[iCut]);
         }

         if(fDoMesonAnalysis){
            hESDTrueMotherInvMassPt[iCut] = new TH2F("ESD_TrueMother_InvMass_Pt","ESD_TrueMother_InvMass_Pt",800,0,0.8,250,0,25);
            fTrueList[iCut]->Add(hESDTrueMotherInvMassPt[iCut]);
            hESDTruePrimaryMotherInvMassPt[iCut]
               = new TH2F("ESD_TruePrimaryMother_InvMass_Pt", "ESD_TruePrimaryMother_InvMass_Pt", 800,0,0.8,250,0,25);
            hESDTruePrimaryMotherInvMassPt[iCut]->Sumw2();
            fTrueList[iCut]->Add(hESDTruePrimaryMotherInvMassPt[iCut]);
            hESDTrueSecondaryMotherInvMassPt[iCut]
               = new TH2F("ESD_TrueSecondaryMother_InvMass_Pt", "ESD_TrueSecondaryMother_InvMass_Pt", 800,0,0.8,250,0,25);
            fTrueList[iCut]->Add(hESDTrueSecondaryMotherInvMassPt[iCut]);
            hESDTrueSecondaryMotherFromK0sInvMassPt[iCut]
               = new TH2F("ESD_TrueSecondaryMotherFromK0s_InvMass_Pt","ESD_TrueSecondaryMotherFromK0s_InvMass_Pt",800,0,0.8,250,0,25);
            fTrueList[iCut]->Add(hESDTrueSecondaryMotherFromK0sInvMassPt[iCut]);
            hESDTrueSecondaryMotherFromEtaInvMassPt[iCut]
               = new TH2F("ESD_TrueSecondaryMotherFromEta_InvMass_Pt","ESD_TrueSecondaryMotherFromEta_InvMass_Pt",800,0,0.8,250,0,25);
            fTrueList[iCut]->Add(hESDTrueSecondaryMotherFromEtaInvMassPt[iCut]);

            if (fDoMesonQA){
               hESDTruePrimaryPi0MCPtResolPt[iCut] = new TH2F("ESD_TruePrimaryPi0_MCPt_ResolPt","ESD_TruePrimaryPi0_ResolPt_MCPt",500,0,25,1000,-1.,1.);
               hESDTruePrimaryPi0MCPtResolPt[iCut]->Sumw2();
               fTrueList[iCut]->Add(hESDTruePrimaryPi0MCPtResolPt[iCut]);
               hESDTruePrimaryEtaMCPtResolPt[iCut]  = new TH2F("ESD_TruePrimaryEta_MCPt_ResolPt","ESD_TruePrimaryEta_ResolPt_MCPt",500,0,25,1000,-1.,1.);
               hESDTruePrimaryEtaMCPtResolPt[iCut]->Sumw2();
               fTrueList[iCut]->Add(hESDTruePrimaryEtaMCPtResolPt[iCut]);
               hESDTrueBckGGInvMassPt[iCut] = new TH2F("ESD_TrueBckGG_InvMass_Pt","ESD_TrueBckGG_InvMass_Pt",800,0,0.8,250,0,25);
               fTrueList[iCut]->Add(hESDTrueBckGGInvMassPt[iCut]);
               hESDTrueBckContInvMassPt[iCut] = new TH2F("ESD_TrueBckCont_InvMass_Pt","ESD_TrueBckCont_InvMass_Pt",800,0,0.8,250,0,25);
               fTrueList[iCut]->Add(hESDTrueBckContInvMassPt[iCut]);
               hESDTrueMotherDalitzInvMassPt[iCut] = new TH2F("ESD_TrueDalitz_InvMass_Pt","ESD_TrueDalitz_InvMass_Pt",800,0,0.8,250,0,25);
               fTrueList[iCut]->Add(hESDTrueMotherDalitzInvMassPt[iCut]);
               hESDTrueK0sWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueK0sWithPi0Daughter_MCPt","ESD_TrueK0sWithPi0Daughter_MCPt",250,0,25);
               fTrueList[iCut]->Add(hESDTrueK0sWithPi0DaughterMCPt[iCut]);
               hESDTrueEtaWithPi0DaughterMCPt[iCut] = new TH1F("ESD_TrueEtaWithPi0Daughter_MCPt","ESD_TrueEtaWithPi0Daughter_MCPt",250,0,25);
               fTrueList[iCut]->Add(hESDTrueEtaWithPi0DaughterMCPt[iCut]);

               fTrueMotherRapList[iCut] = new TList();
               fTrueMotherRapList[iCut]->SetName(Form("%s_%s True Mother Y histograms",cutstring.Data(),cutstringMeson.Data()));
               fTrueMotherRapList[iCut]->SetOwner(kTRUE);
               fCutFolder[iCut]->Add(fTrueMotherRapList[iCut]);
               sESDTruePrimaryMotherInvMassPtY[iCut] = new THnSparseF("TruePrimaryMother_InvMass_Pt_Y","TruePrimaryMother_InvMass_Pt_Y",nDim,nBins,xMin,xMax);
               fTrueMotherRapList[iCut]->Add(sESDTruePrimaryMotherInvMassPtY[iCut]);

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
      if(!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))) continue;
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
         fCutFolder[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
      }
   }

   PostData(1, fOutputContainer);
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

   fMCEvent = MCEvent();
   if(fMCEvent){
      fMCStack = fMCEvent->Stack();
   }
   fInputEvent = InputEvent();

   fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
   CountESDTracks(); // Estimate Event Multiplicity

   // ------------------- BeginEvent ----------------------------

   for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      fiCut = iCut;

      Int_t eventNotAccepted =
         ((AliConversionCuts*)fCutArray->At(iCut))
         ->IsEventAcceptedByConversionCut(fV0Reader->GetConversionCuts(),fInputEvent,fMCEvent,fIsHeavyIon);
      if(eventNotAccepted){
         // 			cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
         hNEvents[iCut]->Fill(eventNotAccepted); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
         continue;
      }

      if(eventQuality != 0){// Event Not Accepted
         //          cout << "event rejected due to: " <<eventQuality << endl;
         hNEvents[iCut]->Fill(eventQuality);
         continue;
      }

      hNEvents[iCut]->Fill(eventQuality); // Should be 0 here
      hNGoodESDTracks[iCut]->Fill(fNumberOfESDTracks);
      hNV0Tracks[iCut]->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());

      if(fMCEvent){ // Process MC Particle
         if(((AliConversionCuts*)fCutArray->At(iCut))->GetSignalRejection() != 0){
            ((AliConversionCuts*)fCutArray->At(iCut))->GetNotRejectedParticles(((AliConversionCuts*)fCutArray->At(iCut))->GetSignalRejection(),
                                                                               ((AliConversionCuts*)fCutArray->At(iCut))->GetAcceptedHeader(),
                                                                               fMCEvent);
         }
         ProcessMCParticles();
      }

      ProcessPhotonCandidates(); // Process this cuts gammas

      hNGammaCandidates[iCut]->Fill(fGammaCandidates->GetEntries());
      if(fDoMesonAnalysis){ // Meson Analysis
         if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fMCEvent){
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
         if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fMCEvent){
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
      if(fMCEvent && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 0){
         Int_t isPosFromMBHeader
            = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack);
         if(isPosFromMBHeader == 0 && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
         Int_t isNegFromMBHeader
            = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack);
         if(isNegFromMBHeader == 0 && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 3) continue;

         if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
      }

      if(!((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fInputEvent)) continue;
      if(!((AliConversionCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
         !((AliConversionCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
         fGammaCandidates->Add(PhotonCandidate); // if no second loop is required add to events good gammas

         if(fIsFromMBHeader){
            hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
            if (fDoPhotonQA)hESDConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
         }
         if(fMCEvent){
            ProcessTruePhotonCandidates(PhotonCandidate);
         }
      }
      else if(((AliConversionCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
         ((AliConversionCuts*)fCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
         nV0++;
         GammaCandidatesStepOne->Add(PhotonCandidate);
      }
      else if(!((AliConversionCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
              ((AliConversionCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
         GammaCandidatesStepTwo->Add(PhotonCandidate);
      }
   }
   if(((AliConversionCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){
      for(Int_t i = 0;i<GammaCandidatesStepOne->GetEntries();i++){
         AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GammaCandidatesStepOne->At(i);
         if(!PhotonCandidate) continue;
         fIsFromMBHeader = kTRUE;
         if(fMCEvent && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 0){
            Int_t isPosFromMBHeader
               = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack);
            Int_t isNegFromMBHeader
               = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack);
            if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
         }
         if(!((AliConversionCuts*)fCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GammaCandidatesStepOne->GetEntries())) continue;
         if(!((AliConversionCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
            fGammaCandidates->Add(PhotonCandidate);
            if(fIsFromMBHeader){
               hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
               if (fDoPhotonQA)hESDConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
            }
            if(fMCEvent){
               ProcessTruePhotonCandidates(PhotonCandidate);
            }
         }
         else GammaCandidatesStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
      }
   }
   if(((AliConversionCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
      for(Int_t i = 0;i<GammaCandidatesStepTwo->GetEntries();i++){
         AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GammaCandidatesStepTwo->At(i);
         if(!PhotonCandidate) continue;
         fIsFromMBHeader = kTRUE;
         if(fMCEvent && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 0){
            Int_t isPosFromMBHeader
               = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack);
            Int_t isNegFromMBHeader
               = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack);
            if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
         }
         if(!((AliConversionCuts*)fCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GammaCandidatesStepTwo,i)) continue;
         fGammaCandidates->Add(PhotonCandidate); // Add gamma to current cut TList
         if(fIsFromMBHeader){
            hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
            if (fDoPhotonQA)hESDConvGammaR[fiCut]->Fill(PhotonCandidate->GetConversionRadius());
         }
         if(fMCEvent){
            ProcessTruePhotonCandidates(PhotonCandidate);
         }
      }
   }

   delete GammaCandidatesStepOne;
   GammaCandidatesStepOne = 0x0;
   delete GammaCandidatesStepTwo;
   GammaCandidatesStepTwo = 0x0;

}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessTruePhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
   // Process True Photons
   AliStack *MCStack = fMCEvent->Stack();
   TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(MCStack);
   TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(MCStack);

   if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist

   Int_t pdgCode[2] = {abs(posDaughter->GetPdgCode()),abs(negDaughter->GetPdgCode())};

   if(posDaughter->GetMother(0) != negDaughter->GetMother(0)){
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
      return;
   }
   else if(posDaughter->GetMother(0) == -1){
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
      return;
   }

   if(pdgCode[0]!=11 || pdgCode[1]!=11) return; //One Particle is not a electron

   if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge

   if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion

   TParticle *Photon = TruePhotonCandidate->GetMCParticle(MCStack);
   if(Photon->GetPdgCode() != 22) return; // Mother is no Photon

   // True Photon
   if(fIsFromMBHeader)hESDTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());

   if(posDaughter->GetMother(0) <= MCStack->GetNprimary()){
      // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
      if(fIsFromMBHeader){
         hESDTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
         if (fDoPhotonQA){
            hESDTruePrimaryConvGammaR[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius());
            hESDTruePrimaryConvGammaEta[fiCut]->Fill(TruePhotonCandidate->Eta());
         }
         hESDTruePrimaryConvGammaRSESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
      }
      hESDTruePrimaryConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt()); // Allways Filled
      // (Not Filled for i6, Extra Signal Gamma (parambox) are secondary)
   }
   else{
      if(fIsFromMBHeader){
         hESDTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
         if (fDoPhotonQA) hESDTrueSecondaryConvGammaR[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius());
         if(MCStack->Particle(Photon->GetMother(0))->GetMother(0) > -1 &&
            MCStack->Particle(MCStack->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
            hESDTrueSecondaryConvGammaFromXFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt());
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

      Bool_t mcIsFromMB = kTRUE;
      Int_t isMCFromMBHeader = -1;
      if(((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 0){
         isMCFromMBHeader
            = ((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack);
         if(isMCFromMBHeader == 0 && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 3) continue;
         if(isMCFromMBHeader != 2) mcIsFromMB = kFALSE;
      }

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
         if (fDoPhotonQA){
            hMCConvGammaR[fiCut]->Fill(((TParticle*)fMCStack->Particle(particle->GetFirstDaughter()))->R());
            hMCConvGammaEta[fiCut]->Fill(particle->Eta());
         }
         if(mcIsFromMB){
            hMCConvGammaRSPt[fiCut]->Fill(particle->Pt());
            if (fDoPhotonQA){
               hMCConvGammaRSR[fiCut]->Fill(((TParticle*)fMCStack->Particle(particle->GetFirstDaughter()))->R());
               hMCConvGammaRSEta[fiCut]->Fill(particle->Eta());
            }
         }
      } // Converted MC Gamma
      if(fDoMesonAnalysis){
         if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMC(particle,fMCStack)){
            TParticle* daughter0 = (TParticle*)fMCStack->Particle(particle->GetFirstDaughter());
            TParticle* daughter1 = (TParticle*)fMCStack->Particle(particle->GetLastDaughter());

            Float_t weighted= 1;
            if(((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack)){
               if (particle->Pt()>0.005){
                  weighted= ((AliConversionCuts*)fCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, fMCStack);
               }
            }
            Double_t mesonY = 10.;
            if(particle->Energy() - particle->Pz() == 0 || particle->Energy() + particle->Pz() == 0){
               mesonY=10.;
            } else{
               mesonY = 0.5*(TMath::Log((particle->Energy()+particle->Pz()) / (particle->Energy()-particle->Pz())));
            }

            if(particle->GetPdgCode() == 111){
               hMCPi0Pt[fiCut]->Fill(particle->Pt(),weighted); // All MC Pi0
               if (fDoMesonQA) hMCPi0PtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
            } else if(particle->GetPdgCode() == 221){
               hMCEtaPt[fiCut]->Fill(particle->Pt(),weighted); // All MC Eta
               if (fDoMesonQA) hMCEtaPtY[fiCut]->Fill(particle->Pt(),mesonY,weighted); // All MC Pi0
            }

            // Check the acceptance for both gammas
            if(((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter0,fMCStack,kFALSE) &&
               ((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter1,fMCStack,kFALSE) ){

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

            if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE))){
               hESDMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());

               if(pi0cand->GetAlpha()<0.1)
                  hESDMotherInvMassEalpha[fiCut]->Fill(pi0cand->M(),pi0cand->E());
               if (fDoMesonQA){
                  Double_t sparesFill2[3] = {pi0cand->M(),pi0cand->Pt(),pi0cand->Rapidity()};
                  sESDMotherInvMassPtY[fiCut]->Fill(sparesFill2,1);
               }
               if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->DoBGCalculation()){
                  Int_t zbin = 0;
                  Int_t mbin = 0;

                  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->BackgroundHandlerType() == 0){

                     zbin = fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
                     if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                        mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
                     } else {
                        mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
                     }
                  }
                  else{
                     zbin = fBGHandlerRP[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
                     if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                        mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
                     } else {
                        mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fGammaCandidates->GetEntries());
                     }
                  }
                  Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),(Double_t)zbin,(Double_t)mbin};
                  sESDMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
               }
               if(fMCEvent){
                  ProcessTrueMesonCandidates(pi0cand,gamma0,gamma1);
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

   if(TrueGammaCandidate0->GetV0Index()<fInputEvent->GetNumberOfV0s()){
      Bool_t isTruePi0 = kFALSE;
      Bool_t isTrueEta = kFALSE;
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
            if(gammaMC0->GetPdgCode() ==111){ // Conversion but Pi0 Mother
               gamma0MotherLabel=-111;
            }
            if(gammaMC0->GetPdgCode() ==221){ // Conversion but Eta Mother
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
               if(gammaMC1->GetPdgCode() ==111){ // Conversion but Pi0 Mother
                  gamma1MotherLabel=-111;
               }
               if(gammaMC1->GetPdgCode() ==221){ // Conversion but Eta Mother
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
         if(isTruePi0 || isTrueEta){// True Pion or Eta
            hESDTrueMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());

            if(gamma0MotherLabel >= MCStack->GetNprimary()){ // Secondary Meson
               hESDTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
               if (((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetMother(0) >-1){
                  if(MCStack->Particle(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetMother(0))->GetPdgCode()==kK0Short){
                     hESDTrueSecondaryMotherFromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
                     if (fDoMesonQA)hESDTrueK0sWithPi0DaughterMCPt[fiCut]
                                       ->Fill(MCStack->Particle(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetMother(0))->Pt());
                  }
                  if(MCStack->Particle(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetMother(0))->GetPdgCode()==221){
                     hESDTrueSecondaryMotherFromEtaInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
                     if (fDoMesonQA)hESDTrueEtaWithPi0DaughterMCPt[fiCut]
                                       ->Fill(MCStack->Particle(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetMother(0))->Pt());
                  }
               }
            }else{ // Only primary pi0 for efficiency calculation
               Float_t weighted= 1;
               if(((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(gamma1MotherLabel, fMCStack)){
                  if (((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt()>0.005){
                     weighted= ((AliConversionCuts*)fCutArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),gamma1MotherLabel, fMCStack);
                  }
               }
               hESDTruePrimaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);

               if (fDoMesonQA){
                  Double_t sparesFill[3] = {Pi0Candidate->M(),Pi0Candidate->Pt(),Pi0Candidate->Rapidity()};
                  sESDTruePrimaryMotherInvMassPtY[fiCut]->Fill(sparesFill,1);
                  if(isTruePi0){ // Only primary pi0 for resolution
                     hESDTruePrimaryPi0MCPtResolPt[fiCut]->Fill(((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt())/((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),weighted);
                  }
                  if (isTrueEta){ // Only primary eta for resolution
                     hESDTruePrimaryEtaMCPtResolPt[fiCut]->Fill(((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),(Pi0Candidate->Pt()-((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt())/((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt(),weighted);
                  }
               }
            }
         }
         else if(!isTruePi0 && !isTrueEta && fDoMesonQA){ // Background
            if(gamma0MotherLabel>-1 && gamma1MotherLabel>-1){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
               hESDTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            } else { // No photon or without mother
               hESDTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            }
            if((gamma0MotherLabel==-111 || gamma1MotherLabel==-111 || gamma0MotherLabel==-221 || gamma1MotherLabel==-221) ){
               // Dalitz
               hESDTrueMotherDalitzInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            }
         }
      }
   }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::CalculateBackground(){

   Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fInputEvent->GetPrimaryVertex()->GetZ());
   Int_t mbin = 0;

   if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
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
               if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE))){
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

                  AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
                  if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE))){
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

                     AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);

                     if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE))){
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
      mbin = fBGHandlerRP[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
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

               if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(&backgroundCandidate,kFALSE)){
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

                  if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(&backgroundCandidate,kFALSE)){
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
         fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fNumberOfESDTracks);
      }
      else{ // means we use #V0s for multiplicity
         fBGHandler[fiCut]->AddEvent(fGammaCandidates,fInputEvent->GetPrimaryVertex()->GetX(),fInputEvent->GetPrimaryVertex()->GetY(),fInputEvent->GetPrimaryVertex()->GetZ(),fGammaCandidates->GetEntries());
      }
   }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::CountESDTracks(){

   // Using standard function for setting Cuts
   Bool_t selectPrimaries=kTRUE;
   AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
   EsdTrackCuts->SetMaxDCAToVertexZ(2);
   EsdTrackCuts->SetEtaRange(-0.8, 0.8);
   EsdTrackCuts->SetPtRange(0.15);

   fNumberOfESDTracks = 0;

   for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
      AliESDtrack* curTrack = (AliESDtrack*) fInputEvent->GetTrack(iTracks);
      if(!curTrack) continue;
      // if(fMCEvent && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 0){
      //    if(!((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(abs(curTrack->GetLabel()), fMCStack)) continue;
      // }
      if(EsdTrackCuts->AcceptTrack(curTrack) ) fNumberOfESDTracks++;
   }
   delete EsdTrackCuts;
   EsdTrackCuts=0x0;

   return;
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::Terminate(const Option_t *)
{

   // Not Executed by GRID on SubJobLevel
   for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      if(!((AliConversionCuts*)fCutArray->At(iCut))) continue;
      if(((AliConversionCuts*)fCutArray->At(iCut))->GetSignalRejection() == 2 && fMCEvent){
         fHeaderNameList[iCut] = new TList();
         TString HeaderNames = "Header:";
         for(Int_t i = 0;i<(((AliConversionCuts*)fCutArray->At(iCut))->GetAcceptedHeader())->GetEntries();i++){
            HeaderNames = HeaderNames+"_"+ ((TObjString*)((TList*) ( (AliConversionCuts*)fCutArray->At(iCut))
                                                          ->GetAcceptedHeader())->At(i))->GetString();
         }
         fHeaderNameList[iCut]->SetName(HeaderNames);
         fHeaderNameList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fHeaderNameList[iCut]);
      }
      else if(((AliConversionCuts*)fCutArray->At(iCut))->GetSignalRejection() == 0 &&
              (((AliConversionCuts*)fCutArray->At(iCut))->GetFoundHeader()) && fMCEvent){
         fHeaderNameList[iCut] = new TList();
         TString HeaderNames = (((AliConversionCuts*)fCutArray->At(iCut))->GetFoundHeader())[0];
         fHeaderNameList[iCut]->SetName(HeaderNames);
         fHeaderNameList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fHeaderNameList[iCut]);
      }
   }

   //fOutputContainer->Print(); // Will crash on GRID
}
