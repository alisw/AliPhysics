/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *									  *
 * Author: Martin Wilde, Daniel Lohner					  *
 * Version 1.0								  *
 *									  *
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
#include "AliGammaConversionAODBGHandler.h"
#include "AliGenCocktailEventHeader.h"

ClassImp(AliAnalysisTaskGammaConvV1)

//________________________________________________________________________
AliAnalysisTaskGammaConvV1::AliAnalysisTaskGammaConvV1(): AliAnalysisTaskSE(),
   fV0Reader(NULL),
   fBGHandler(NULL),
   fESDEvent(NULL),
   fMCEvent(NULL),
   fMCStack(NULL),
   fCutFolder(NULL),
   fESDList(NULL),
   fBackList(NULL),
   fMotherList(NULL),
   fTrueList(NULL),
   fMCList(NULL),
   fHeaderNameList(NULL),
   fOutputContainer(0),
   fReaderGammas(NULL),
   fGoodGammas(NULL),
   fCutArray(NULL),
   fConversionCuts(NULL),
   fMesonCutArray(NULL),
   fMesonCuts(NULL),
   hESDConvGammaPt(NULL),
   hESDMotherInvMassPt(NULL),
   sESDMotherInvMassPtZM(NULL),
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
   hMCConvGammaPt(NULL),
   hMCConvGammaR(NULL),
   hMCConvGammaEta(NULL),
   hMCPi0Pt(NULL),
   hMCEtaPt(NULL),
   hMCPi0InAccPt(NULL),
   hMCEtaInAccPt(NULL),
   hESDTrueMotherInvMassPt(NULL),
   hESDTruePrimaryMotherInvMassMCPt(NULL),
   hESDTruePrimaryPi0ESDPtMCPt(NULL),
   hESDTrueSecondaryMotherInvMassPt(NULL),
   hESDTrueSecondaryMotherFromK0sInvMassPt(NULL),
   hESDTrueBckGGInvMassPt(NULL),
   hESDTrueBckContInvMassPt(NULL),
   hESDTrueMotherDalitzInvMassPt(NULL),
   hESDTrueConvGammaPt(NULL),
   hESDTrueTwoElecCombPt(NULL),
   hESDTrueTwoPionCombPt(NULL),
   hESDTrueElecPionCombPt(NULL),
   hESDTrueCombPt(NULL),
   hESDTruePrimaryConvGammaPt(NULL),
   hESDTruePrimaryConvGammaR(NULL),
   hESDTruePrimaryConvGammaEta(NULL),
   hESDTruePrimaryConvGammaESDPtMCPt(NULL),
   hESDTrueSecondaryConvGammaPt(NULL),
   hESDTrueSecondaryConvGammaFromK0sPt(NULL),
   hESDTrueSecondaryConvGammaFromXFromK0sPt(NULL),
   hNEvents(NULL),
   hNGoodESDTracks(NULL),
   hNV0Tracks(NULL),
   fRandom(0),
   fUnsmearedPx(NULL),
   fUnsmearedPy(NULL),
   fUnsmearedPz(NULL),
   fUnsmearedE(NULL),
   fIsHeavyIon(kFALSE),
   fDoMesonAnalysis(kTRUE)

{
   // default Constructor
   DefineInput(0, TChain::Class());
   DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskGammaConvV1::AliAnalysisTaskGammaConvV1(const char *name):
   AliAnalysisTaskSE(name),
   fV0Reader(NULL),
   fBGHandler(NULL),
   fESDEvent(NULL),
   fMCEvent(NULL),
   fMCStack(NULL),
   fCutFolder(NULL),
   fESDList(NULL),
   fBackList(NULL),
   fMotherList(NULL),
   fTrueList(NULL),
   fMCList(NULL),
   fHeaderNameList(NULL),
   fOutputContainer(0),
   fReaderGammas(NULL),
   fGoodGammas(NULL),
   fCutArray(NULL),
   fConversionCuts(NULL),
   fMesonCutArray(NULL),
   fMesonCuts(NULL),
   hESDConvGammaPt(NULL),
   hESDMotherInvMassPt(NULL),
   sESDMotherInvMassPtZM(NULL),
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
   hMCConvGammaPt(NULL),
   hMCConvGammaR(NULL),
   hMCConvGammaEta(NULL),
   hMCPi0Pt(NULL),
   hMCEtaPt(NULL),
   hMCPi0InAccPt(NULL),
   hMCEtaInAccPt(NULL),
   hESDTrueMotherInvMassPt(NULL),
   hESDTruePrimaryMotherInvMassMCPt(NULL),
   hESDTruePrimaryPi0ESDPtMCPt(NULL),
   hESDTrueSecondaryMotherInvMassPt(NULL),
   hESDTrueSecondaryMotherFromK0sInvMassPt(NULL),
   hESDTrueBckGGInvMassPt(NULL),
   hESDTrueBckContInvMassPt(NULL),
   hESDTrueMotherDalitzInvMassPt(NULL),
   hESDTrueConvGammaPt(NULL),
   hESDTrueTwoElecCombPt(NULL),
   hESDTrueTwoPionCombPt(NULL),
   hESDTrueElecPionCombPt(NULL),
   hESDTrueCombPt(NULL),
   hESDTruePrimaryConvGammaPt(NULL),
   hESDTruePrimaryConvGammaR(NULL),
   hESDTruePrimaryConvGammaEta(NULL),
   hESDTruePrimaryConvGammaESDPtMCPt(NULL),
   hESDTrueSecondaryConvGammaPt(NULL),
   hESDTrueSecondaryConvGammaFromK0sPt(NULL),
   hESDTrueSecondaryConvGammaFromXFromK0sPt(NULL),
   hNEvents(NULL),
   hNGoodESDTracks(NULL),
   hNV0Tracks(NULL),
   fRandom(0),
   fUnsmearedPx(NULL),
   fUnsmearedPy(NULL),
   fUnsmearedPz(NULL),
   fUnsmearedE(NULL),
   fIsHeavyIon(kFALSE),
   fDoMesonAnalysis(kTRUE)

{
   // Define input and output slots here
   DefineInput(0, TChain::Class());
   DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaConvV1::~AliAnalysisTaskGammaConvV1()
{
   if(fGoodGammas){
      delete fGoodGammas;
      fGoodGammas = 0x0;
   }
   if(fBGHandler){
      delete[] fBGHandler;
      fBGHandler = 0x0;
   }
}
//___________________________________________________________
void AliAnalysisTaskGammaConvV1::InitBack(){

   Double_t *zBinLimitsArray = new Double_t[9];
   zBinLimitsArray[0] = -50.00;
   zBinLimitsArray[1] = -3.375;
   zBinLimitsArray[2] = -1.605;
   zBinLimitsArray[3] = -0.225;
   zBinLimitsArray[4] = 1.065;
   zBinLimitsArray[5] = 2.445;
   zBinLimitsArray[6] = 4.245;
   zBinLimitsArray[7] = 50.00;
   zBinLimitsArray[8] = 1000.00;

   Double_t *multiplicityBinLimitsArrayTracks = new Double_t[6];
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

   Double_t *multiplicityBinLimitsArrayV0s = new Double_t[5];
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
   Int_t nBins[nDim] = {1000,250,8,5};
   Double_t xMin[nDim] = {0,0, 0,0};
   Double_t xMax[nDim] = {1,25,8,5};

   sESDMotherInvMassPtZM = new THnSparseF*[fnCuts];
   sESDMotherBackInvMassPtZM = new THnSparseF*[fnCuts];

   fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];
   for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      TString cutstring = ((AliConversionCuts*)fCutArray->At(iCut))->GetCutNumber();
		TString cutstringMeson = ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutNumber();
      fBackList[iCut] = new TList();
      fBackList[iCut]->SetName(Form("%s_%s Back histograms",cutstring.Data(),cutstringMeson.Data()));
      fBackList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fBackList[iCut]);

      sESDMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m","Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
      sESDMotherBackInvMassPtZM[iCut]->Sumw2();
      fBackList[iCut]->Add(sESDMotherBackInvMassPtZM[iCut]);
      
      fMotherList[iCut] = new TList();
      fMotherList[iCut]->SetName(Form("%s_%s Mother histograms",cutstring.Data(),cutstringMeson.Data()));
      fMotherList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fMotherList[iCut]);

      sESDMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m","Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
      sESDMotherInvMassPtZM[iCut]->Sumw2();
      fMotherList[iCut]->Add(sESDMotherInvMassPtZM[iCut]);

      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseTrackMultiplicity()){
         fBGHandler[iCut] = new AliGammaConversionAODBGHandler(9,6,((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->NumberOfBGEvents());
         fBGHandler[iCut]->Initialize(zBinLimitsArray, multiplicityBinLimitsArrayTracks);
      }
      else{
         fBGHandler[iCut] = new AliGammaConversionAODBGHandler(9,5,((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->NumberOfBGEvents());
         fBGHandler[iCut]->Initialize(zBinLimitsArray, multiplicityBinLimitsArrayV0s);
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
   fGoodGammas = new TList();

   fCutFolder = new TList*[fnCuts];
   fESDList = new TList*[fnCuts];
   fBackList = new TList*[fnCuts];
   fMotherList = new TList*[fnCuts];
   hESDConvGammaPt = new TH1F*[fnCuts];
   hNEvents = new TH1I*[fnCuts];
   hNGoodESDTracks = new TH1I*[fnCuts];
   hNV0Tracks = new TH1I*[fnCuts];

   if(fDoMesonAnalysis){
      hESDMotherInvMassPt = new TH2F*[fnCuts];
      hESDMotherBackInvMassPt = new TH2F*[fnCuts];
      hESDMotherInvMassEalpha = new TH2F*[fnCuts];
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

      hNEvents[iCut] = new TH1I("NEvents","NEvents",7,-0.5,6.5);
      fESDList[iCut]->Add(hNEvents[iCut]);
      if(fIsHeavyIon) hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",3000,0,3000);
      else hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
      fESDList[iCut]->Add(hNGoodESDTracks[iCut]);
      if(fIsHeavyIon) hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",25000,0,25000);
      else hNV0Tracks[iCut] = new TH1I("V0 Multiplicity","V0 Multiplicity",2000,0,2000);
      fESDList[iCut]->Add(hNV0Tracks[iCut]);

      hESDConvGammaPt[iCut] = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",250,0,25);
      fESDList[iCut]->Add(hESDConvGammaPt[iCut]);

      if(fDoMesonAnalysis){
         hESDMotherInvMassPt[iCut] = new TH2F("ESD_Mother_InvMass_Pt","ESD_Mother_InvMass_Pt",1000,0,1,250,0,25);
         fESDList[iCut]->Add(hESDMotherInvMassPt[iCut]);
         hESDMotherBackInvMassPt[iCut] = new TH2F("ESD_Background_InvMass_Pt","ESD_Background_InvMass_Pt",1000,0,1,250,0,25);
         fESDList[iCut]->Add(hESDMotherBackInvMassPt[iCut]);
         hESDMotherInvMassEalpha[iCut] = new TH2F("ESD_Mother_InvMass_vs_E_alpha","ESD_Mother_InvMass_vs_E_alpha",1000,0,1,250,0,25);
         fESDList[iCut]->Add(hESDMotherInvMassEalpha[iCut]);
      }

      fCutFolder[iCut]->Add(fESDList[iCut]);
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
      hMCConvGammaPt = new TH1F*[fnCuts];
      hMCConvGammaR = new TH1F*[fnCuts];
      hMCConvGammaEta = new TH1F*[fnCuts];
      hESDTrueConvGammaPt = new TH1F*[fnCuts];
      hESDTrueTwoElecCombPt = new TH1F*[fnCuts];
      hESDTrueTwoPionCombPt = new TH1F*[fnCuts];
      hESDTrueElecPionCombPt = new TH1F*[fnCuts];
      hESDTrueCombPt = new TH1F*[fnCuts];
      hESDTruePrimaryConvGammaPt = new TH1F*[fnCuts];
      hESDTruePrimaryConvGammaR = new TH1F*[fnCuts];
      hESDTruePrimaryConvGammaEta = new TH1F*[fnCuts];
      hESDTruePrimaryConvGammaESDPtMCPt = new TH2F*[fnCuts];
      hESDTrueSecondaryConvGammaPt = new TH1F*[fnCuts];
      hESDTrueSecondaryConvGammaFromK0sPt = new TH1F*[fnCuts];
      hESDTrueSecondaryConvGammaFromXFromK0sPt = new TH1F*[fnCuts];

      if(fDoMesonAnalysis){
         hMCPi0Pt = new TH1F*[fnCuts];
         hMCEtaPt = new TH1F*[fnCuts];
         hMCPi0InAccPt = new TH1F*[fnCuts];
         hMCEtaInAccPt = new TH1F*[fnCuts];

         hESDTrueMotherInvMassPt = new TH2F*[fnCuts];
         hESDTruePrimaryPi0ESDPtMCPt = new TH2F*[fnCuts];
         hESDTruePrimaryMotherInvMassMCPt = new TH2F*[fnCuts];
         hESDTrueSecondaryMotherInvMassPt = new TH2F*[fnCuts];
         hESDTrueSecondaryMotherFromK0sInvMassPt = new TH2F*[fnCuts];
         hESDTrueBckGGInvMassPt = new TH2F*[fnCuts];
         hESDTrueBckContInvMassPt = new TH2F*[fnCuts];
         hESDTrueMotherDalitzInvMassPt = new TH2F*[fnCuts];
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
         hMCConvGammaPt[iCut] = new TH1F("MC_ConvGamma_Pt","MC_ConvGamma_Pt",250,0,25);
         fMCList[iCut]->Add(hMCConvGammaPt[iCut]);
         hMCConvGammaR[iCut] = new TH1F("MC_ConvGamma_R","MC_ConvGamma_R",1000,0,250);
         fMCList[iCut]->Add(hMCConvGammaR[iCut]);
         hMCConvGammaEta[iCut] = new TH1F("MC_ConvGamma_Eta","MC_ConvGamma_Eta",100,-4,4);
         fMCList[iCut]->Add(hMCConvGammaEta[iCut]);
         if(fDoMesonAnalysis){
            hMCPi0Pt[iCut] = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",250,0,25);
            fMCList[iCut]->Add(hMCPi0Pt[iCut]);
            hMCEtaPt[iCut] = new TH1F("MC_Eta_Pt","MC_Eta_Pt",250,0,25);
            fMCList[iCut]->Add(hMCEtaPt[iCut]);
            hMCPi0InAccPt[iCut] = new TH1F("MC_Pi0InAcc_Pt","MC_Pi0InAcc_Pt",250,0,25);
            fMCList[iCut]->Add(hMCPi0InAccPt[iCut]);
            hMCEtaInAccPt[iCut] = new TH1F("MC_EtaInAcc_Pt","MC_EtaInAcc_Pt",250,0,25);
            fMCList[iCut]->Add(hMCEtaInAccPt[iCut]);
         }
         fTrueList[iCut] = new TList();
         fTrueList[iCut]->SetName(Form("%s_%s True histograms",cutstring.Data(),cutstringMeson.Data()));
         fTrueList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fTrueList[iCut]);

         hESDTrueConvGammaPt[iCut] = new TH1F("ESD_TrueConvGamma_Pt","ESD_TrueConvGamma_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueConvGammaPt[iCut]);
         hESDTrueTwoElecCombPt[iCut] = new TH1F("ESD_TrueTwoElecComb_Pt","ESD_TrueTwoElecComb_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueTwoElecCombPt[iCut]);
         hESDTrueTwoPionCombPt[iCut] = new TH1F("ESD_TrueTwoPionComb_Pt","ESD_TrueTwoPionComb_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueTwoPionCombPt[iCut]);
         hESDTrueElecPionCombPt[iCut] = new TH1F("ESD_TrueElecPionComb_Pt","ESD_TrueElecPionComb_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueElecPionCombPt[iCut]);
         hESDTrueCombPt[iCut] = new TH1F("ESD_TrueComb_Pt","ESD_TrueComb_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueCombPt[iCut]);
         hESDTruePrimaryConvGammaPt[iCut] = new TH1F("ESD_TruePrimaryConvGamma_Pt","ESD_TruePrimaryConvGamma_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTruePrimaryConvGammaPt[iCut]);
         hESDTruePrimaryConvGammaR[iCut] = new TH1F("ESD_TruePrimaryConvGamma_R","ESD_TruePrimaryConvGamma_R",1000,0,250);
         fTrueList[iCut]->Add(hESDTruePrimaryConvGammaR[iCut]);
         hESDTruePrimaryConvGammaEta[iCut] = new TH1F("ESD_TruePrimaryConvGamma_Eta","ESD_TruePrimaryConvGamma_Eta",100,-4,4);
         fTrueList[iCut]->Add(hESDTruePrimaryConvGammaEta[iCut]);
         hESDTrueSecondaryConvGammaPt[iCut] = new TH1F("ESD_TrueSecondaryConvGamma_Pt","ESD_TrueSecondaryConvGamma_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueSecondaryConvGammaPt[iCut]);
         hESDTrueSecondaryConvGammaFromK0sPt[iCut] = new TH1F("ESD_TrueSecondaryConvGammaFromK0s_Pt","ESD_TrueSecondaryConvGammaFromK0s_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueSecondaryConvGammaFromK0sPt[iCut]);
         hESDTrueSecondaryConvGammaFromXFromK0sPt[iCut] = new TH1F("ESD_TrueSecondaryConvGammaFromXFromK0s_Pt", "ESD_TrueSecondaryConvGammaFromXFromK0s_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueSecondaryConvGammaFromXFromK0sPt[iCut]);
         hESDTruePrimaryConvGammaESDPtMCPt[iCut] = new TH2F("ESD_TruePrimaryConvGammaESD_PtMCPt", "ESD_TruePrimaryConvGammaESD_PtMCPt",250,0,25,250,0,25);
         fTrueList[iCut]->Add(hESDTruePrimaryConvGammaESDPtMCPt[iCut]);

         if(fDoMesonAnalysis){
            hESDTrueMotherInvMassPt[iCut] = new TH2F("ESD_TrueMother_InvMass_Pt","ESD_TrueMother_InvMass_Pt",1000,0,1,250,0,25);
            fTrueList[iCut]->Add(hESDTrueMotherInvMassPt[iCut]);
            hESDTruePrimaryPi0ESDPtMCPt[iCut] = new TH2F("ESD_TruePrimaryPi0_ESDPt_MCPt","ESD_TruePrimaryPi0_ESDPt_MCPt",250,0,25,250,0,25);
            fTrueList[iCut]->Add(hESDTruePrimaryPi0ESDPtMCPt[iCut]);
            hESDTruePrimaryMotherInvMassMCPt[iCut] = new TH2F("ESD_TruePrimaryMother_InvMass_MCPt", "ESD_TruePrimaryMother_InvMass_MCPt", 1000,0,1,250,0,25);
            fTrueList[iCut]->Add(hESDTruePrimaryMotherInvMassMCPt[iCut]);
            hESDTrueSecondaryMotherInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryMother_InvMass_Pt", "ESD_TrueSecondaryMother_InvMass_Pt", 1000,0,1,250,0,25);
            fTrueList[iCut]->Add(hESDTrueSecondaryMotherInvMassPt[iCut]);
            hESDTrueSecondaryMotherFromK0sInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryMotherFromK0s_InvMass_Pt","ESD_TrueSecondaryMotherFromK0s_InvMass_Pt",1000,0,1,250,0,25);
            fTrueList[iCut]->Add(hESDTrueSecondaryMotherFromK0sInvMassPt[iCut]);
            hESDTrueBckGGInvMassPt[iCut] = new TH2F("ESD_TrueBckGG_InvMass_Pt","ESD_TrueBckGG_InvMass_Pt",1000,0,1,250,0,25);
            fTrueList[iCut]->Add(hESDTrueBckGGInvMassPt[iCut]);
            hESDTrueBckContInvMassPt[iCut] = new TH2F("ESD_TrueBckCont_InvMass_Pt","ESD_TrueBckCont_InvMass_Pt",1000,0,1,250,0,25);
            fTrueList[iCut]->Add(hESDTrueBckContInvMassPt[iCut]);
            hESDTrueMotherDalitzInvMassPt[iCut] = new TH2F("ESD_TrueDalitz_InvMass_Pt","ESD_TrueDalitz_InvMass_Pt",1000,0,1,250,0,25);
            fTrueList[iCut]->Add(hESDTrueMotherDalitzInvMassPt[iCut]);
         }
         
         if(((AliConversionCuts*)fCutArray->At(iCut))->GetSignalRejection() == 2){
            fHeaderNameList[iCut] = new TList();
            TString HeaderNames = "Header";
            for(Int_t i = 0;i<(((AliConversionCuts*)fCutArray->At(iCut))->GetAcceptedHeader())->GetEntries();i++){
               HeaderNames = HeaderNames+"_"+ ((TObjString*)((TList*) ( (AliConversionCuts*)fCutArray->At(iCut)) ->GetAcceptedHeader())->At(i))->GetString();
            }
            fHeaderNameList[iCut]->SetName(HeaderNames);
            fHeaderNameList[iCut]->SetOwner(kTRUE);
            fCutFolder[iCut]->Add(fHeaderNameList[iCut]);
         }
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
   fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
   if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

   Int_t eventQuality = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetEventQuality();
   if(eventQuality != 0){// Event Not Accepted
      for(Int_t iCut = 0; iCut<fnCuts; iCut++){
         hNEvents[iCut]->Fill(eventQuality);
      }
      return;
   }

   fMCEvent = MCEvent();
   fESDEvent = (AliESDEvent*) InputEvent();
   fReaderGammas = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
   CountESDTracks(); // Estimate Event Multiplicity

   for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      fiCut = iCut;
      if(fIsHeavyIon && !((AliConversionCuts*)fCutArray->At(iCut))->IsCentralitySelected(fESDEvent)){
         hNEvents[iCut]->Fill(1); // Check Centrality --> Not Accepted => eventQuality = 1
         continue;
      }
      hNEvents[iCut]->Fill(eventQuality);

      hNGoodESDTracks[iCut]->Fill(fNumberOfESDTracks);
      hNV0Tracks[iCut]->Fill(fESDEvent->GetVZEROData()->GetMTotV0A()+fESDEvent->GetVZEROData()->GetMTotV0C());
      if(fMCEvent){ // Process MC Particle
         fMCStack = fMCEvent->Stack();
         if(((AliConversionCuts*)fCutArray->At(iCut))->GetSignalRejection() != 0){
            ((AliConversionCuts*)fCutArray->At(iCut))->GetNotRejectedParticles(((AliConversionCuts*)fCutArray->At(iCut))->GetSignalRejection(),
                                                                               ((AliConversionCuts*)fCutArray->At(iCut))->GetAcceptedHeader(),
                                                                               fMCEvent);
         }
         ProcessMCParticles();
      }
      
      ProcessPhotonCandidates(); // Process this cuts gammas

      if(fDoMesonAnalysis){ // Meson Analysis
         if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fMCEvent){
            fUnsmearedPx = new Double_t[fGoodGammas->GetEntries()]; // Store unsmeared Momenta
            fUnsmearedPy = new Double_t[fGoodGammas->GetEntries()];
            fUnsmearedPz = new Double_t[fGoodGammas->GetEntries()];
            fUnsmearedE =  new Double_t[fGoodGammas->GetEntries()];

            for(Int_t gamma=0;gamma<fGoodGammas->GetEntries();gamma++){ // Smear the AODPhotons in MC
               fUnsmearedPx[gamma] = ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->Px();
               fUnsmearedPy[gamma] = ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->Py();
               fUnsmearedPz[gamma] = ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->Pz();
               fUnsmearedE[gamma] =  ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->E();
               ((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->SmearParticle(dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(gamma)));
            }
         }

         CalculatePi0Candidates(); // Combine Gammas
         CalculateBackground(); // Combinatorial Background
         UpdateEventByEventData(); // Store Event for mixed Events

         if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->UseMCPSmearing() && fMCEvent){
            for(Int_t gamma=0;gamma<fGoodGammas->GetEntries();gamma++){ // Smear the AODPhotons in MC
               ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetPx(fUnsmearedPx[gamma]); // Reset Unsmeared Momenta
               ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetPy(fUnsmearedPy[gamma]);
               ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetPz(fUnsmearedPz[gamma]);
               ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetE(fUnsmearedE[gamma]);
            }
            delete[] fUnsmearedPx; fUnsmearedPx = 0x0;
            delete[] fUnsmearedPy; fUnsmearedPy = 0x0;
            delete[] fUnsmearedPz; fUnsmearedPz = 0x0;
            delete[] fUnsmearedE;  fUnsmearedE  = 0x0;
         }
      }
      fGoodGammas->Clear(); // delete this cuts good gammas

   }
   
   PostData(1, fOutputContainer);
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessPhotonCandidates()
{
   Int_t nV0 = 0;
   TList *GoodGammasStepOne = new TList();
   TList *GoodGammasStepTwo = new TList();
   // Loop over Photon Candidates allocated by ReaderV1
   for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
      AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
      if(!PhotonCandidate) continue;

      if(fMCEvent && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 0){
         if(!((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack)) continue;
         if(!((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack)) continue;
      }

      if(!((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fESDEvent)) continue;

      if(!((AliConversionCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
         !((AliConversionCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // if no post reader loop is required add to events good gammas
         fGoodGammas->Add(PhotonCandidate);
         hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
         if(fMCEvent){
            ProcessTruePhotonCandidates(PhotonCandidate);
         }
      }
      else if(((AliConversionCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
         ((AliConversionCuts*)fCutArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
         nV0++;
         GoodGammasStepOne->Add(PhotonCandidate);
      }
      else if(!((AliConversionCuts*)fCutArray->At(fiCut))->UseElecSharingCut() &&
              ((AliConversionCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
         GoodGammasStepTwo->Add(PhotonCandidate);
      }
   }
   if(((AliConversionCuts*)fCutArray->At(fiCut))->UseElecSharingCut()){
      for(Int_t i = 0;i<GoodGammasStepOne->GetEntries();i++){
         AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GoodGammasStepOne->At(i);
         if(!PhotonCandidate) continue;
         if(!((AliConversionCuts*)fCutArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GoodGammasStepOne->GetEntries())) continue;
         if(!((AliConversionCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
            fGoodGammas->Add(PhotonCandidate);
            hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
            if(fMCEvent){
               ProcessTruePhotonCandidates(PhotonCandidate);
            }
         }
         else GoodGammasStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
      }
   }
   if(((AliConversionCuts*)fCutArray->At(fiCut))->UseToCloseV0sCut()){
      for(Int_t i = 0;i<GoodGammasStepTwo->GetEntries();i++){
         AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GoodGammasStepTwo->At(i);
         if(!PhotonCandidate) continue;
         if(!((AliConversionCuts*)fCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GoodGammasStepTwo,i)) continue;
         fGoodGammas->Add(PhotonCandidate); // Add gamma to current cut TList
         hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt()); // Differences to old V0Reader in p_t due to conversion KF->TLorentzVector
         if(fMCEvent){
            ProcessTruePhotonCandidates(PhotonCandidate);
         }
      }
   }

   delete GoodGammasStepOne;
   GoodGammasStepOne = 0x0;
   delete GoodGammasStepTwo;
   GoodGammasStepTwo = 0x0;

}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::ProcessTruePhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
   // Process True Photons
   AliStack *MCStack = fMCEvent->Stack();
   TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(MCStack);
   TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(MCStack);

   if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
   if(posDaughter->GetMother(0) != negDaughter->GetMother(0)){  // Not Same Mother == Combinatorial Bck
      if(TMath::Abs(posDaughter->GetPdgCode())==11 && TMath::Abs(negDaughter->GetPdgCode())==11)
         hESDTrueTwoElecCombPt[fiCut]->Fill(TruePhotonCandidate->Pt()); //Electron Combinatorial
      else if(TMath::Abs(posDaughter->GetPdgCode())==211 && TMath::Abs(negDaughter->GetPdgCode())==211)
         hESDTrueTwoPionCombPt[fiCut]->Fill(TruePhotonCandidate->Pt()); //At least on Pion Combinatorial
      else if( (TMath::Abs(posDaughter->GetPdgCode())==11 && TMath::Abs(negDaughter->GetPdgCode())==211) ||
               (TMath::Abs(posDaughter->GetPdgCode())==211 && TMath::Abs(negDaughter->GetPdgCode())==11) )
         hESDTrueElecPionCombPt[fiCut]->Fill(TruePhotonCandidate->Pt()); //At least on Pion Combinatorial
      else hESDTrueCombPt[fiCut]->Fill(TruePhotonCandidate->Pt()); //At least on Pion Combinatorial
      return;
   }
   else if(posDaughter->GetMother(0) == -1){
      if(TMath::Abs(posDaughter->GetPdgCode())==11 && TMath::Abs(negDaughter->GetPdgCode())==11)
         hESDTrueTwoElecCombPt[fiCut]->Fill(TruePhotonCandidate->Pt()); //Electron Combinatorial
      else if(TMath::Abs(posDaughter->GetPdgCode())==211 && TMath::Abs(negDaughter->GetPdgCode())==211)
         hESDTrueTwoPionCombPt[fiCut]->Fill(TruePhotonCandidate->Pt()); //At least on Pion Combinatorial
      else if( (TMath::Abs(posDaughter->GetPdgCode())==11 && TMath::Abs(negDaughter->GetPdgCode())==211) ||
               (TMath::Abs(posDaughter->GetPdgCode())==211 && TMath::Abs(negDaughter->GetPdgCode())==11) )
         hESDTrueElecPionCombPt[fiCut]->Fill(TruePhotonCandidate->Pt()); //At least on Pion Combinatorial
      else hESDTrueCombPt[fiCut]->Fill(TruePhotonCandidate->Pt()); //At least on Pion Combinatorial
      return;
   }
   if(TMath::Abs(posDaughter->GetPdgCode())!=11 || TMath::Abs(negDaughter->GetPdgCode())!=11) return; //One Particle is not electron
   if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge
   if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion

   TParticle *Photon = TruePhotonCandidate->GetMCParticle(MCStack);
   if(Photon->GetPdgCode() != 22) return; // Mother is no Photon

   // True Photon
   hESDTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());

   if(posDaughter->GetMother(0) <= MCStack->GetNprimary()){
      // Count just primary MC Gammas as true --> For Ratio esdtruegamma / mcconvgamma
      hESDTruePrimaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
      hESDTruePrimaryConvGammaR[fiCut]->Fill(TruePhotonCandidate->GetConversionRadius());
      hESDTruePrimaryConvGammaEta[fiCut]->Fill(TruePhotonCandidate->Eta());
      hESDTruePrimaryConvGammaESDPtMCPt[fiCut]->Fill(TruePhotonCandidate->Pt(),Photon->Pt());
   }
   else{
      hESDTrueSecondaryConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
      if(MCStack->Particle(Photon->GetMother(0))->GetPdgCode() == 310){
         hESDTrueSecondaryConvGammaFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt());
      }
      if(MCStack->Particle(Photon->GetMother(0))->GetMother(0) > -1 &&
         MCStack->Particle(MCStack->Particle(Photon->GetMother(0))->GetMother(0))->GetPdgCode() == 310){
         hESDTrueSecondaryConvGammaFromXFromK0sPt[fiCut]->Fill(TruePhotonCandidate->Pt());
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

      if(((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 0){
         if(!((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack)) continue;
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
            }
         }
      }
      if(((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kTRUE)){
         hMCConvGammaPt[fiCut]->Fill(particle->Pt());
         hMCConvGammaR[fiCut]->Fill(((TParticle*)fMCStack->Particle(particle->GetFirstDaughter()))->R());
         hMCConvGammaEta[fiCut]->Fill(particle->Eta());
      } // Converted MC Gamma
      if(fDoMesonAnalysis){
         if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelectedMC(particle,fMCStack,kFALSE)){
            if(particle->GetPdgCode() == 111)hMCPi0Pt[fiCut]->Fill(particle->Pt()); // All MC Pi0
            if(particle->GetPdgCode() == 221)hMCEtaPt[fiCut]->Fill(particle->Pt()); // All MC Eta
            // Check the acceptance for both gammas
            if(particle->GetNDaughters() == 2){
               TParticle* daughter0 = (TParticle*)fMCStack->Particle(particle->GetFirstDaughter());
               TParticle* daughter1 = (TParticle*)fMCStack->Particle(particle->GetLastDaughter());
               if(((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter0,fMCStack,kFALSE) &&
                  ((AliConversionCuts*)fCutArray->At(fiCut))->PhotonIsSelectedMC(daughter1,fMCStack,kFALSE) ){
                  if(particle->GetPdgCode() == 111)hMCPi0InAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0 with gamma in acc
                  if(particle->GetPdgCode() == 221)hMCEtaInAccPt[fiCut]->Fill(particle->Pt()); // MC Eta with gamma in acc
               }
            }
         }
      }
   }
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::CalculatePi0Candidates(){

   // Conversion Gammas
   if(fGoodGammas->GetEntries()>1){
      for(Int_t firstGammaIndex=0;firstGammaIndex<fGoodGammas->GetEntries()-1;firstGammaIndex++){
         AliAODConversionPhoton *gamma0=dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(firstGammaIndex));
         for(Int_t secondGammaIndex=firstGammaIndex+1;secondGammaIndex<fGoodGammas->GetEntries();secondGammaIndex++){
            AliAODConversionPhoton *gamma1=dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(secondGammaIndex));
            //Check for same Electron ID
            if(gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelPositive() ||
               gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelNegative() ||
               gamma0->GetTrackLabelNegative() == gamma1->GetTrackLabelPositive() ||
               gamma0->GetTrackLabelPositive() == gamma1->GetTrackLabelNegative() ) continue;

            AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma0,gamma1);
            pi0cand->SetLabels(firstGammaIndex,secondGammaIndex);

            if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE))){
               hESDMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
               if(pi0cand->GetAlpha()<0.1){
                  hESDMotherInvMassEalpha[fiCut]->Fill(pi0cand->M(),pi0cand->E());
               }
               Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fESDEvent->GetPrimaryVertex()->GetZ());
               Int_t mbin = 0;
               if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
                  mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
               } else {
                  mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGoodGammas->GetEntries());
               }
               Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),zbin,mbin};
               sESDMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
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

   if(TrueGammaCandidate0->GetV0Index()<fESDEvent->GetNumberOfV0s()){
      Bool_t isTruePi0 = kFALSE;
      Bool_t isTrueEta = kFALSE;
      Int_t gamma0MCLabel = TrueGammaCandidate0->GetMCParticleLabel(MCStack);
      Int_t gamma0MotherLabel = -1;
      if(gamma0MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
         // Daughters Gamma 0
         TParticle * negativeMC = (TParticle*)TrueGammaCandidate0->GetNegativeMCDaughter(MCStack);
         TParticle * positiveMC = (TParticle*)TrueGammaCandidate0->GetPositiveMCDaughter(MCStack);
         TParticle * gammaMC0 = (TParticle*)MCStack->Particle(gamma0MCLabel);
         if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
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
      if(TrueGammaCandidate1->GetV0Index()<fESDEvent->GetNumberOfV0s()){
         Int_t gamma1MCLabel = TrueGammaCandidate1->GetMCParticleLabel(MCStack);
         Int_t gamma1MotherLabel = -1;
         if(gamma1MCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
            // Daughters Gamma 1
            TParticle * negativeMC = (TParticle*)TrueGammaCandidate1->GetNegativeMCDaughter(MCStack);
            TParticle * positiveMC = (TParticle*)TrueGammaCandidate1->GetPositiveMCDaughter(MCStack);
            TParticle * gammaMC1 = (TParticle*)MCStack->Particle(gamma1MCLabel);
            if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...
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
         if(isTruePi0 || isTrueEta){ // True Pion or Eta
            hESDTrueMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            if(gamma0MotherLabel > MCStack->GetNprimary()){ // Secondary Meson
               hESDTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
               if (((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetMother(0) >-1){
                  if(MCStack->Particle(((TParticle*)MCStack->Particle(gamma1MotherLabel))->GetMother(0))->GetPdgCode()==kK0Short){
                     hESDTrueSecondaryMotherFromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
                  }
               }
            }
            if(gamma0MotherLabel <= MCStack->GetNprimary()){ // Only primary pi0 for efficiency calculation
               hESDTruePrimaryMotherInvMassMCPt[fiCut]->Fill(Pi0Candidate->M(),((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt());
               if(isTruePi0){ // Only primaries for unfolding
                  hESDTruePrimaryPi0ESDPtMCPt[fiCut]->Fill(Pi0Candidate->Pt(),((TParticle*)MCStack->Particle(gamma1MotherLabel))->Pt());
               }
            }
         }
         if(!isTruePi0 && !isTrueEta){ // Background
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

   Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fESDEvent->GetPrimaryVertex()->GetZ());
   Int_t mbin = 0;

   if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
      mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
   } else {
      mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGoodGammas->GetEntries());
   }

   if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseRotationMethod()){

      for(Int_t iCurrent=0;iCurrent<fGoodGammas->GetEntries();iCurrent++){
         AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGoodGammas->At(iCurrent));
         for(Int_t iCurrent2=iCurrent+1;iCurrent2<fGoodGammas->GetEntries();iCurrent2++){
            for(Int_t nRandom=0;nRandom<((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->NumberOfBGEvents();nRandom++){
               AliAODConversionPhoton currentEventGoodV02 = *(AliAODConversionPhoton*)(fGoodGammas->At(iCurrent2));

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
                  Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),zbin,mbin};
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

            for(Int_t iCurrent=0;iCurrent<fGoodGammas->GetEntries();iCurrent++){
               AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGoodGammas->At(iCurrent));
               for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
                  AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));
                  if(fMoveParticleAccordingToVertex == kTRUE){
                     MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
                  }

                  AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
                  if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE))){
                     hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
                     Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),zbin,mbin};
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
               for(Int_t iCurrent=0;iCurrent<fGoodGammas->GetEntries();iCurrent++){
                  AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGoodGammas->At(iCurrent));
                  for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){

                     AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));

                     if(fMoveParticleAccordingToVertex == kTRUE){
                        MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
                     }

                     AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);

                     if((((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE))){
                        hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
                        Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),zbin,mbin};
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
void AliAnalysisTaskGammaConvV1::RotateParticle(AliAODConversionPhoton *gamma){
   Int_t fNDegreesPMBackground= ((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->NDegreesRotation();
   Double_t nRadiansPM = fNDegreesPMBackground*TMath::Pi()/180;
   Double_t rotationValue = fRandom.Rndm()*2*nRadiansPM + TMath::Pi()-nRadiansPM;
   gamma->RotateZ(rotationValue);
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
   //see header file for documentation

   Double_t dx = vertex->fX - fESDEvent->GetPrimaryVertex()->GetX();
   Double_t dy = vertex->fY - fESDEvent->GetPrimaryVertex()->GetY();
   Double_t dz = vertex->fZ - fESDEvent->GetPrimaryVertex()->GetZ();

   Double_t movedPlace[3] = {particle->GetConversionX() - dx,particle->GetConversionY() - dy,particle->GetConversionZ() - dz};
   particle->SetConversionPoint(movedPlace);
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::UpdateEventByEventData(){
   //see header file for documentation
   if(fGoodGammas->GetEntries() >0 ){
      if(((AliConversionMesonCuts*)fMesonCutArray->At(fiCut))->UseTrackMultiplicity()){
         fBGHandler[fiCut]->AddEvent(fGoodGammas,fESDEvent->GetPrimaryVertex()->GetX(),fESDEvent->GetPrimaryVertex()->GetY(),fESDEvent->GetPrimaryVertex()->GetZ(),fNumberOfESDTracks);
      }
      else{ // means we use #V0s for multiplicity
         fBGHandler[fiCut]->AddEvent(fGoodGammas,fESDEvent->GetPrimaryVertex()->GetX(),fESDEvent->GetPrimaryVertex()->GetY(),fESDEvent->GetPrimaryVertex()->GetZ(),fGoodGammas->GetEntries());
      }
   }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::CountESDTracks(){

   AliESDtrackCuts *EsdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
   // Using standard function for setting Cuts
   Bool_t selectPrimaries=kTRUE;
   EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
   EsdTrackCuts->SetMaxDCAToVertexZ(2);
   EsdTrackCuts->SetEtaRange(-0.8, 0.8);
   EsdTrackCuts->SetPtRange(0.15);

   fNumberOfESDTracks = 0;
   for(Int_t iTracks = 0; iTracks < fESDEvent->GetNumberOfTracks(); iTracks++){
      AliESDtrack* curTrack = fESDEvent->GetTrack(iTracks);
      if(!curTrack) continue;
      if(EsdTrackCuts->AcceptTrack(curTrack) ) fNumberOfESDTracks++;
   }
   delete EsdTrackCuts;
   EsdTrackCuts=0x0;

   return;
}
//________________________________________________________________________
void AliAnalysisTaskGammaConvV1::Terminate(const Option_t *)
{
   fOutputContainer->Add(((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());

   for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      if(((AliConversionCuts*)fCutArray->At(iCut))->GetCutHistograms()){
         fCutFolder[iCut]->Add(((AliConversionCuts*)fCutArray->At(iCut))->GetCutHistograms());
      }
      if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
         fCutFolder[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
      }
   }
   fOutputContainer->Print();
}
