/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *																								*
 * Authors: Friederike Bock															*
 * Version 1.0																				*
 *																								*
 * Permission to use, copy, modify and distribute this software and its	 *
 * documentation strictly for non-commercial purposes is hereby granted	 *
 * without fee, provided that the above copyright notice appears in all	 *
 * copies and that both the copyright notice and this permission notice	 *
 * appear in the supporting documentation. The authors make no claims	 *
 * about the suitability of this software for any purpose. It is			*
 * provided "as is" without express or implied warranty.						*
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// QA Task for V0 Reader V1
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskMaterialHistos.h"
#include "TChain.h"
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "TVectorF.h"
#include "AliPIDResponse.h"
#include "AliESDtrackCuts.h"
#include "TFile.h"

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskMaterialHistos)

AliAnalysisTaskMaterialHistos::AliAnalysisTaskMaterialHistos() : AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fConversionGammas(NULL),
  fGammaCandidates(NULL),
  fConversionCutArray(NULL),
  fEventCutArray(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
	fOutputList(NULL),
	fAllMCGammaList(NULL),
	fAllMCConvGammaList(NULL),
	fPrimVtxZ(0.),
	fNContrVtx(0),
	fNESDtracksEta09(0),
	fNESDtracksEta0914(0),
	fNESDtracksEta14(0),
	fGammaMCPt(0.),
	fGammaMCTheta(0.),
	fGammaMCConvPt(0.),
	fGammaMCConvTheta(0.),
	fGammaPt(0.),
	fGammaTheta(0.),
	fGammaChi2NDF(0.),
	fKind(0),      
	fIsHeavyIon(0),
	fIsMC(0),
	fInputEvent(NULL),
	fMCEvent(NULL),
  fMCStack(NULL),
	fnCuts(0),
	fiCut(0),
  hNEvents(NULL),
	hNGoodESDTracks09(NULL),
	hNGoodESDTracks14(NULL),
	hNGoodESDTracks09_14(NULL),       
	hESDConversionMappingRPhi(NULL),
	hESDConversionMappingRZ(NULL),
	hESDConversionR(NULL),
        hESDConversionRInBins(NULL),
        hESDConversionPhiInBins(NULL),
	hESDConversionAsymP(NULL),
	hESDConversionMidPtR(NULL),
	hESDConversionEta(NULL),
	hESDConversionMidPtEta(NULL),
	hESDConversionPt(NULL),
	hESDConversionPt5cm(NULL),
	hESDConversionDCA(NULL),            
	hESDConversionMidPtDCA(NULL),            
 	hESDConversionPsiPair(NULL),       
 	hESDConversionChi2(NULL),         
 	hESDConversionMass(NULL),         
	hMCConversionMappingRPhi(NULL),
	hMCConversionR(NULL),
	hMCConversionMidPtR(NULL),
	hMCConversionEta(NULL),
	hMCConversionMidPtEta(NULL),
	hMCConversionPt(NULL),
	hMCAllGammaPt(NULL),
	hMCTrueConversionMappingRPhi(NULL),
	hMCTrueConversionMappingRZ(NULL),
	hMCTrueConversionR(NULL),
	hMCTrueConversionMidPtR(NULL),
	hMCTrueConversionEta(NULL),
	hMCTrueConversionMidPtEta(NULL),
	hMCTrueConversionPt(NULL),
	hMCTrueConversionPt5cm(NULL),
	hMCTrueConversionAsymP(NULL),
	hMCTrueConversionDCA(NULL),
	hMCTrueConversionMidPtDCA(NULL),
	hMCTrueConversionPsiPair(NULL),
	hMCTrueConversionChi2(NULL),
	hMCTrueConversionMass(NULL),
	hMCTruePi0DalConversionR(NULL),
	hMCTruePi0DalConversionMidPtR(NULL),
	hMCTruePi0DalConversionEta(NULL),
	hMCTruePi0DalConversionMidPtEta(NULL),
	hMCTruePi0DalConversionPt(NULL),
	hMCTrueEtaDalConversionR(NULL),
	hMCTrueEtaDalConversionMidPtR(NULL),
	hMCTrueEtaDalConversionEta(NULL),
	hMCTrueEtaDalConversionMidPtEta(NULL),
	hMCTrueEtaDalConversionPt(NULL),
	hMCTrueCombConversionR(NULL),
	hMCTrueCombConversionMidPtR(NULL),
	hMCTrueCombConversionEta(NULL),
	hMCTrueCombConversionMidPtEta(NULL),
	hMCTrueCombConversionPt(NULL)
{
   
}


//________________________________________________________________________
AliAnalysisTaskMaterialHistos::AliAnalysisTaskMaterialHistos(const char *name) : AliAnalysisTaskSE(name),
fV0Reader(NULL),
fV0ReaderName("V0ReaderV1"),
fConversionGammas(NULL),
fGammaCandidates(NULL),
fConversionCutArray(NULL),
	fEventCutArray(NULL),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
	fOutputList(NULL),
	fAllMCGammaList(NULL),
	fAllMCConvGammaList(NULL),
	fPrimVtxZ(0.),
	fNContrVtx(0),
	fNESDtracksEta09(0),
	fNESDtracksEta0914(0),
	fNESDtracksEta14(0),
	fGammaMCPt(0.),
	fGammaMCTheta(0.),
	fGammaMCConvPt(0.),
	fGammaMCConvTheta(0.),
	fGammaPt(0.),
	fGammaTheta(0.),
	fGammaChi2NDF(0.),
	fKind(0),      
	fIsHeavyIon(0),
	fIsMC(0),
	fInputEvent(NULL),
	fMCEvent(NULL),
        fMCStack(NULL),
  fnCuts(0),							
  fiCut(0),
  hNEvents(NULL),
	hNGoodESDTracks09(NULL),
	hNGoodESDTracks14(NULL),
	hNGoodESDTracks09_14(NULL),       
	hESDConversionMappingRPhi(NULL),
	hESDConversionMappingRZ(NULL),
	hESDConversionR(NULL),
        hESDConversionRInBins(NULL),
        hESDConversionPhiInBins(NULL),
	hESDConversionAsymP(NULL),
	hESDConversionMidPtR(NULL),
	hESDConversionEta(NULL),
	hESDConversionMidPtEta(NULL),
	hESDConversionPt(NULL),
	hESDConversionPt5cm(NULL),
	hESDConversionDCA(NULL),            
	hESDConversionMidPtDCA(NULL),            
 	hESDConversionPsiPair(NULL),       
 	hESDConversionChi2(NULL),         
 	hESDConversionMass(NULL),         
  hMCConversionMappingRPhi(NULL),
	hMCConversionR(NULL),
	hMCConversionMidPtR(NULL),
	hMCConversionEta(NULL),
	hMCConversionMidPtEta(NULL),
	hMCConversionPt(NULL),
	hMCAllGammaPt(NULL),
  hMCTrueConversionMappingRPhi(NULL),
  hMCTrueConversionMappingRZ(NULL),
	hMCTrueConversionR(NULL),
	hMCTrueConversionMidPtR(NULL),
	hMCTrueConversionEta(NULL),
	hMCTrueConversionMidPtEta(NULL),
	hMCTrueConversionPt(NULL),
	hMCTrueConversionPt5cm(NULL),
	hMCTrueConversionAsymP(NULL),
	hMCTrueConversionDCA(NULL),
	hMCTrueConversionMidPtDCA(NULL),
	hMCTrueConversionPsiPair(NULL),
	hMCTrueConversionChi2(NULL),
	hMCTrueConversionMass(NULL),
	hMCTruePi0DalConversionR(NULL),
	hMCTruePi0DalConversionMidPtR(NULL),
	hMCTruePi0DalConversionEta(NULL),
	hMCTruePi0DalConversionMidPtEta(NULL),
	hMCTruePi0DalConversionPt(NULL),
  hMCTrueEtaDalConversionR(NULL),
	hMCTrueEtaDalConversionMidPtR(NULL),
	hMCTrueEtaDalConversionEta(NULL),
	hMCTrueEtaDalConversionMidPtEta(NULL),
	hMCTrueEtaDalConversionPt(NULL),
  hMCTrueCombConversionR(NULL),
	hMCTrueCombConversionMidPtR(NULL),
	hMCTrueCombConversionEta(NULL),
	hMCTrueCombConversionMidPtEta(NULL),
	hMCTrueCombConversionPt(NULL)
{
	// Default constructor
	
	DefineInput(0, TChain::Class());
	DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskMaterialHistos::~AliAnalysisTaskMaterialHistos()
{
	// default deconstructor
  if(fGammaCandidates){
    delete fGammaCandidates;
    fGammaCandidates = 0x0;
  }
}
//________________________________________________________________________
void AliAnalysisTaskMaterialHistos::UserCreateOutputObjects()
{
	// Create User Output Objects

	if(fOutputList != NULL){
		delete fOutputList;
		fOutputList = NULL;
	}
	if(fOutputList == NULL){
		fOutputList = new TList();
		fOutputList->SetOwner(kTRUE);
	}

  // Array of current cut's gammas
	
	fGammaCandidates          = new TList();
	fCutFolder                = new TList*[fnCuts];
	fESDList                  = new TList*[fnCuts];
	fMCList                   = new TList*[fnCuts];
	fTrueList                 = new TList*[fnCuts];	


	hNEvents                       = new TH1F*[fnCuts];
	hNGoodESDTracks09              = new TH1F*[fnCuts];
	hNGoodESDTracks14              = new TH1F*[fnCuts];
	hNGoodESDTracks09_14           = new TH1F*[fnCuts];       
	hESDConversionMappingRPhi      = new TH2F*[fnCuts];
	hESDConversionMappingRZ        = new TH2F*[fnCuts];
	hESDConversionR                = new TH1F*[fnCuts];
        hESDConversionRInBins          = new TH3F*[fnCuts];
        hESDConversionPhiInBins        = new TH3F*[fnCuts];
	hESDConversionAsymP            = new TH2F*[fnCuts];
	hESDConversionMidPtR           = new TH1F*[fnCuts];
	hESDConversionEta              = new TH1F*[fnCuts];
	hESDConversionMidPtEta         = new TH1F*[fnCuts];
	hESDConversionPt               = new TH1F*[fnCuts];
	hESDConversionPt5cm            = new TH1F*[fnCuts];
	hESDConversionDCA              = new TH1F*[fnCuts];
	hESDConversionMidPtDCA         = new TH1F*[fnCuts];
	hESDConversionPsiPair          = new TH1F*[fnCuts];
	hESDConversionChi2             = new TH1F*[fnCuts];
	hESDConversionMass             = new TH1F*[fnCuts];
	
	hMCConversionMappingRPhi       = new TH2F*[fnCuts];
	hMCConversionR                 = new TH1F*[fnCuts];
	hMCConversionMidPtR            = new TH1F*[fnCuts];
	hMCConversionEta               = new TH1F*[fnCuts];
	hMCConversionMidPtEta          = new TH1F*[fnCuts];
	hMCConversionPt                = new TH1F*[fnCuts];
	hMCAllGammaPt                  = new TH1F*[fnCuts];
	
	hMCTrueConversionMappingRPhi   = new TH2F*[fnCuts];
	hMCTrueConversionMappingRZ     = new TH2F*[fnCuts];
	hMCTrueConversionR             = new TH1F*[fnCuts];
	hMCTrueConversionMidPtR        = new TH1F*[fnCuts];
	hMCTrueConversionEta           = new TH1F*[fnCuts];
	hMCTrueConversionMidPtEta      = new TH1F*[fnCuts];
	hMCTrueConversionPt            = new TH1F*[fnCuts];
	hMCTrueConversionPt5cm         = new TH1F*[fnCuts];
	
	hMCTrueConversionAsymP         = new TH2F*[fnCuts];
	hMCTrueConversionDCA           = new TH1F*[fnCuts];
	hMCTrueConversionMidPtDCA      = new TH1F*[fnCuts];
	hMCTrueConversionPsiPair       = new TH1F*[fnCuts];
	hMCTrueConversionChi2          = new TH1F*[fnCuts];
	hMCTrueConversionMass          = new TH1F*[fnCuts];
	

	hMCTruePi0DalConversionR                = new TH1F*[fnCuts];
	hMCTruePi0DalConversionMidPtR           = new TH1F*[fnCuts];
	hMCTruePi0DalConversionEta              = new TH1F*[fnCuts];
	hMCTruePi0DalConversionMidPtEta         = new TH1F*[fnCuts];
	hMCTruePi0DalConversionPt               = new TH1F*[fnCuts];
	
	
	hMCTrueEtaDalConversionEta              = new TH1F*[fnCuts];
	hMCTrueEtaDalConversionMidPtEta         = new TH1F*[fnCuts];
	hMCTrueEtaDalConversionPt               = new TH1F*[fnCuts];
	hMCTrueEtaDalConversionR                = new TH1F*[fnCuts];
	hMCTrueEtaDalConversionMidPtR           = new TH1F*[fnCuts];
	
	hMCTrueCombConversionR                = new TH1F*[fnCuts];
	hMCTrueCombConversionMidPtR           = new TH1F*[fnCuts];
	hMCTrueCombConversionEta                = new TH1F*[fnCuts];
	hMCTrueCombConversionMidPtEta           = new TH1F*[fnCuts];
	hMCTrueCombConversionPt                 = new TH1F*[fnCuts];
	


	for(Int_t iCut = 0; iCut<fnCuts;iCut++){

		
		TString cutstringEvent          = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
		TString cutstringPhoton         = ((AliConversionPhotonCuts*)fConversionCutArray->At(iCut))->GetCutNumber();
		fCutFolder[iCut]            = new TList();
    fCutFolder[iCut]->SetName(Form("Cut Number %s_%s",cutstringEvent.Data() ,cutstringPhoton.Data()));
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputList->Add(fCutFolder[iCut]);

    fESDList[iCut]              = new TList();
    fESDList[iCut]->SetName(Form("%s_%s ESD histograms",cutstringEvent.Data() ,cutstringPhoton.Data()));
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);
		
    Int_t nBinsR=500;
    Int_t nBinsX=2000;
    Int_t nBinsY=2000;
    Int_t nBinsZ=500;
    Int_t nBinsPhi=750;
		Int_t nBinsEta=2000;
		Int_t nBinsPt=500;
		


			hNEvents[iCut]            = new TH1F("NEvents","NEvents",12,-0.5,11.5);
			hNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
			hNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
			hNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");
			if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
				TString TriggerNames      = "Not Trigger: ";
				TriggerNames              = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
				hNEvents[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
			} else {
				hNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
			}
			hNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
			hNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
			hNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
			hNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
			hNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
			hNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL problem");
			hNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
			hNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
			fESDList[iCut]->Add(hNEvents[iCut]);  
			

			hNGoodESDTracks09[iCut]              = new TH1F("GoodESDTracks09","GoodESDTracks09",4000,0,4000);
			fESDList[iCut]->Add(hNGoodESDTracks09[iCut]);

			hNGoodESDTracks14[iCut]              = new TH1F("GoodESDTracks14","GoodESDTracks14",4000,0,4000);
			fESDList[iCut]->Add(hNGoodESDTracks14[iCut]);

			hNGoodESDTracks09_14[iCut]               = new TH1F("GoodESDTracks09_14","GoodESDTracks09_14",4000,0,4000);
			fESDList[iCut]->Add(hNGoodESDTracks09_14[iCut]);

			hESDConversionMappingRPhi[iCut]       = new TH2F("ESD_ConversionMapping_RPhi","ESD_ConversionMapping_RPhi",nBinsPhi,0.,2*TMath::Pi(),nBinsR,0.,200.);

			fESDList[iCut]->Add(hESDConversionMappingRPhi[iCut]);
			hESDConversionMappingRZ[iCut]         = new TH2F("ESD_ConversionMapping_RZ","ESD_ConversionMapping_RZ",nBinsPhi,-180.,180.,nBinsR,0.,200.);
			fESDList[iCut]->Add(hESDConversionMappingRZ[iCut]);
			hESDConversionR[iCut]               = new TH1F("ESD_ConversionMapping_R","ESD_ConversionMapping_R",nBinsR,0.,200.);
			fESDList[iCut]->Add(hESDConversionR[iCut]);

			hESDConversionRInBins[iCut]  = new TH3F("ESD_ConversionMapping_RInBins","ESD_ConversionMapping_RInBins",100,0.,20.,8,0.,4000,nBinsR,0.,200.);
                        fESDList[iCut]->Add(hESDConversionRInBins[iCut]);
		        hESDConversionPhiInBins[iCut]  = new TH3F("ESD_ConversionMapping_PhiInBins","ESD_ConversionMapping_PhiInBins",100,0.,200.,8,0.,4000,nBinsPhi,0.,2*TMath::Pi());
                        fESDList[iCut]->Add(hESDConversionPhiInBins[iCut]);

			hESDConversionAsymP[iCut]               = new TH2F("ESD_ConversionMapping_AsymP","ESD_ConversionMapping_AsymP",nBinsPt,0.01,20.,500,0.,1.);
			fESDList[iCut]->Add(hESDConversionAsymP[iCut]);

			hESDConversionMidPtR[iCut]          = new TH1F("ESD_ConversionMappingMidPt_R","ESD_ConversionMappingMidPt_R",nBinsR,0.,200.);
			fESDList[iCut]->Add(hESDConversionMidPtR[iCut]);	

			hESDConversionEta[iCut]               = new TH1F("ESD_ConversionMapping_Eta","ESD_ConversionMapping_Eta",nBinsEta,-2.,2.);
			fESDList[iCut]->Add(hESDConversionEta[iCut]);
			hESDConversionMidPtEta[iCut]          = new TH1F("ESD_ConversionMappingMidPt_Eta","ESD_ConversionMappingMidPt_Eta",nBinsEta,-2.,2.);
			fESDList[iCut]->Add(hESDConversionMidPtEta[iCut]);	
			hESDConversionPt[iCut]                = new TH1F("ESD_ConversionMapping_Pt","ESD_ConversionMapping_Pt",nBinsPt,0.,20.);
			fESDList[iCut]->Add(hESDConversionPt[iCut]);

			hESDConversionPt5cm[iCut]                = new TH1F("ESD_ConversionMapping_Pt5cm","ESD_ConversionMapping_Pt5cm",nBinsPt,0.,20.);
			fESDList[iCut]->Add(hESDConversionPt5cm[iCut]);

			hESDConversionDCA[iCut]                = new TH1F("ESD_ConversionMapping_DCA","ESD_ConversionMapping_DCA",nBinsPt,0.,5.);
			fESDList[iCut]->Add(hESDConversionDCA[iCut]);

			hESDConversionMidPtDCA[iCut]           = new TH1F("ESD_ConversionMappingMidPt_DCA","ESD_ConversionMappingMidPt_DCA",nBinsPt,0.,5.);
			fESDList[iCut]->Add(hESDConversionMidPtDCA[iCut]);


			hESDConversionPsiPair[iCut]                = new TH1F("ESD_ConversionMapping_PsiPair","ESD_ConversionMapping_PsiPair",nBinsPt,0.,5.);
			fESDList[iCut]->Add(hESDConversionPsiPair[iCut]);

			hESDConversionChi2[iCut]                = new TH1F("ESD_ConversionMapping_Chi2","ESD_ConversionMapping_Chi2",nBinsPt,0.,50.);
			fESDList[iCut]->Add(hESDConversionChi2[iCut]);

			hESDConversionMass[iCut]                = new TH1F("ESD_ConversionMapping_Mass","ESD_ConversionMapping_Mass",nBinsPt,0.,1.);
			fESDList[iCut]->Add(hESDConversionMass[iCut]);




		TAxis *AxisAfter = hESDConversionAsymP[iCut]->GetXaxis();
		Int_t bins = AxisAfter->GetNbins();
		Double_t from = AxisAfter->GetXmin();
		Double_t to = AxisAfter->GetXmax();
		Double_t *newBins = new Double_t[bins+1];
		newBins[0] = from;
		Double_t factor = TMath::Power(to/from, 1./bins);
		for(Int_t i=1; i<=bins; ++i) newBins[i] = factor * newBins[i-1];
		AxisAfter->Set(bins, newBins);



		if (fIsMC>0) {		
			fMCList[iCut]                   = new TList();
			fMCList[iCut]->SetName(Form("%s_%s MC histograms",cutstringEvent.Data() ,cutstringPhoton.Data()));
			fMCList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fMCList[iCut]);

			fTrueList[iCut]                 = new TList();
			fTrueList[iCut]->SetName(Form("%s_%s True histograms",cutstringEvent.Data() ,cutstringPhoton.Data()));
			fTrueList[iCut]->SetOwner(kTRUE);
			fCutFolder[iCut]->Add(fTrueList[iCut]);
				



				hMCConversionMappingRPhi[iCut]       = new TH2F("MC_ConversionMapping_RPhi","MC_ConversionMapping_RPhi",nBinsPhi,0.,2*TMath::Pi(),nBinsR,0.,200.);
				fMCList[iCut]->Add(hMCConversionMappingRPhi[iCut]);
				hMCConversionR[iCut]               = new TH1F("MC_ConversionMapping_R","MC_ConversionMapping_R",nBinsR,0.,200.);
				fMCList[iCut]->Add(hMCConversionR[iCut]);
				hMCConversionMidPtR[iCut]          = new TH1F("MC_ConversionMappingMidPt_R","MC_ConversionMappingMidPt_R",nBinsR,0.,200.);
				fMCList[iCut]->Add(hMCConversionMidPtR[iCut]);	
				hMCConversionEta[iCut]               = new TH1F("MC_ConversionMapping_Eta","MC_ConversionMapping_Eta",nBinsEta,-2.,2.);
				fMCList[iCut]->Add(hMCConversionEta[iCut]);
				hMCConversionMidPtEta[iCut]          = new TH1F("MC_ConversionMappingMidPt_Eta","MC_ConversionMappingMidPt_Eta",nBinsEta,-2.,2.);
				fMCList[iCut]->Add(hMCConversionMidPtEta[iCut]);	
				hMCConversionPt[iCut]                = new TH1F("MC_ConversionMapping_Pt","MC_ConversionMapping_Pt",nBinsPt,0.,20.);
				fMCList[iCut]->Add(hMCConversionPt[iCut]);
				hMCAllGammaPt[iCut]                = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",nBinsPt,0.,20.);
				fMCList[iCut]->Add(hMCAllGammaPt[iCut]);
				


				hMCTrueConversionMappingRPhi[iCut]    = new TH2F("ESD_TrueConversionMapping_RPhi","ESD_TrueConversionMapping_RPhi",nBinsPhi,0.,2*TMath::Pi(),nBinsR,0.,200.);
				fTrueList[iCut]->Add(hMCTrueConversionMappingRPhi[iCut]);

				hMCTrueConversionMappingRZ[iCut]    = new TH2F("ESD_TrueConversionMapping_RZ","ESD_TrueConversionMapping_RZ",nBinsPhi,-180.,180.,nBinsR,0.,200.);
				fTrueList[iCut]->Add(hMCTrueConversionMappingRZ[iCut]);


				hMCTrueConversionR[iCut]               = new TH1F("ESD_TrueConversionMapping_R","ESD_TrueConversionMapping_R",nBinsR,0.,200.);
				fTrueList[iCut]->Add(hMCTrueConversionR[iCut]);
				hMCTrueConversionMidPtR[iCut]          = new TH1F("ESD_TrueConversionMappingMidPt_R","ESD_TrueConversionMappingMidPt_R",nBinsR,0.,200.);
				fTrueList[iCut]->Add(hMCTrueConversionMidPtR[iCut]);	

				hMCTrueConversionEta[iCut]              = new TH1F("ESD_TrueConversionMapping_Eta","ESD_TrueConversionMapping_Eta",nBinsEta,-2.,2.);
				fTrueList[iCut]->Add(hMCTrueConversionEta[iCut]);
				hMCTrueConversionMidPtEta[iCut]         = new TH1F("ESD_TrueConversionMappingMidPt_Eta","ESD_TrueConversionMappingMidPt_Eta",nBinsEta,-2.,2.);
				fTrueList[iCut]->Add(hMCTrueConversionMidPtEta[iCut]);
				hMCTrueConversionPt[iCut]               = new TH1F("ESD_TrueConversionMapping_Pt","ESD_TrueConversionMapping_Pt",nBinsPt,0.,20.);
				fTrueList[iCut]->Add(hMCTrueConversionPt[iCut]);

				hMCTrueConversionPt5cm[iCut]               = new TH1F("ESD_TrueConversionMapping_Pt5cm","ESD_TrueConversionMapping_Pt5cm",nBinsPt,0.,20.);
				fTrueList[iCut]->Add(hMCTrueConversionPt5cm[iCut]);
				
				hMCTrueConversionAsymP[iCut]               = new TH2F("ESD_TrueConversionMapping_AsymP","ESD_TrueConversionMapping_AsymP",nBinsPt,0.01,20.,500,0.,1.);
				fTrueList[iCut]->Add(hMCTrueConversionAsymP[iCut]);

				AxisAfter = hMCTrueConversionAsymP[iCut]->GetXaxis();
				AxisAfter->Set(bins, newBins);

				hMCTrueConversionDCA[iCut]                = new TH1F("ESD_TrueConversionMapping_DCA","ESD_TrueConversionMapping_DCA",nBinsPt,0.,5.);
				fTrueList[iCut]->Add(hMCTrueConversionDCA[iCut]);

				hMCTrueConversionMidPtDCA[iCut]           = new TH1F("ESD_TrueConversionMappingMidPt_DCA","ESD_TrueConversionMappingMidPt_DCA",nBinsPt,0.,5.);
				fTrueList[iCut]->Add(hMCTrueConversionMidPtDCA[iCut]);

 
				hMCTrueConversionPsiPair[iCut]                = new TH1F("ESD_TrueConversionMapping_PsiPair","ESD_TrueConversionMapping_PsiPair",nBinsPt,0.,5.);
				fTrueList[iCut]->Add(hMCTrueConversionPsiPair[iCut]);

				hMCTrueConversionChi2[iCut]                = new TH1F("ESD_TrueConversionMapping_Chi2","ESD_TrueConversionMapping_Chi2",nBinsPt,0.,50.);
				fTrueList[iCut]->Add(hMCTrueConversionChi2[iCut]);

				hMCTrueConversionMass[iCut]                = new TH1F("ESD_TrueConversionMapping_Mass","ESD_TrueConversionMapping_Mass",nBinsPt,0.,1.);
				fTrueList[iCut]->Add(hMCTrueConversionMass[iCut]);


				hMCTruePi0DalConversionR[iCut]               = new TH1F("ESD_TruePi0DalConversionMapping_R","ESD_TruePi0DalConversionMapping_R",nBinsR,0.,200.);
				fTrueList[iCut]->Add(hMCTruePi0DalConversionR[iCut]);
				hMCTruePi0DalConversionMidPtR[iCut]          = new TH1F("ESD_TruePi0DalConversionMappingMidPt_R","ESD_TruePi0DalConversionMappingMidPt_R",nBinsR,0.,200.);
				fTrueList[iCut]->Add(hMCTruePi0DalConversionMidPtR[iCut]);	

				hMCTruePi0DalConversionEta[iCut]              = new TH1F("ESD_TruePi0DalConversionMapping_Eta","ESD_TruePi0DalConversionMapping_Eta",nBinsEta,-2.,2.);
				fTrueList[iCut]->Add(hMCTruePi0DalConversionEta[iCut]);
				hMCTruePi0DalConversionMidPtEta[iCut]         = new TH1F("ESD_TruePi0DalConversionMappingMidPt_Eta","ESD_TruePi0DalConversionMappingMidPt_Eta",nBinsEta,-2.,2.);
				fTrueList[iCut]->Add(hMCTruePi0DalConversionMidPtEta[iCut]);
				hMCTruePi0DalConversionPt[iCut]               = new TH1F("ESD_TruePi0DalConversionMapping_Pt","ESD_TruePi0DalConversionMapping_Pt",nBinsPt,0.,20.);
				fTrueList[iCut]->Add(hMCTruePi0DalConversionPt[iCut]);
				
		 			
				hMCTrueEtaDalConversionR[iCut]               = new TH1F("ESD_TrueEtaDalConversionMapping_R","ESD_TrueEtaDalConversionMapping_R",nBinsR,0.,200.);
				fTrueList[iCut]->Add(hMCTrueEtaDalConversionR[iCut]);
				hMCTrueEtaDalConversionMidPtR[iCut]          = new TH1F("ESD_TrueEtaDalConversionMappingMidPt_R","ESD_TrueEtaDalConversionMappingMidPt_R",nBinsR,0.,200.);
				fTrueList[iCut]->Add(hMCTrueEtaDalConversionMidPtR[iCut]);	

				hMCTrueEtaDalConversionEta[iCut]              = new TH1F("ESD_TrueEtaDalConversionMapping_Eta","ESD_TrueEtaDalConversionMapping_Eta",nBinsEta,-2.,2.);
				fTrueList[iCut]->Add(hMCTrueEtaDalConversionEta[iCut]);
				hMCTrueEtaDalConversionMidPtEta[iCut]         = new TH1F("ESD_TrueEtaDalConversionMappingMidPt_Eta","ESD_TrueEtaDalConversionMappingMidPt_Eta",nBinsEta,-2.,2.);
				fTrueList[iCut]->Add(hMCTrueEtaDalConversionMidPtEta[iCut]);
				hMCTrueEtaDalConversionPt[iCut]               = new TH1F("ESD_TrueEtaDalConversionMapping_Pt","ESD_TrueEtaDalConversionMapping_Pt",nBinsPt,0.,20.);
				fTrueList[iCut]->Add(hMCTrueEtaDalConversionPt[iCut]);
				
				
				
				hMCTrueCombConversionR[iCut]               = new TH1F("ESD_TrueCombConversionMapping_R","ESD_TrueCombConversionMapping_R",nBinsR,0.,200.);
				fTrueList[iCut]->Add(hMCTrueCombConversionR[iCut]);
				hMCTrueCombConversionMidPtR[iCut]          = new TH1F("ESD_TrueCombConversionMappingMidPt_R","ESD_TrueCombConversionMappingMidPt_R",nBinsR,0.,200.);
				fTrueList[iCut]->Add(hMCTrueCombConversionMidPtR[iCut]);	
				hMCTrueCombConversionEta[iCut]              = new TH1F("ESD_TrueCombConversionMapping_Eta","ESD_TrueCombConversionMapping_Eta",nBinsEta,-2.,2.);
				fTrueList[iCut]->Add(hMCTrueCombConversionEta[iCut]);
				hMCTrueCombConversionMidPtEta[iCut]         = new TH1F("ESD_TrueCombConversionMappingMidPt_Eta","ESD_TrueCombConversionMappingMidPt_Eta",nBinsEta,-2.,2.);
				fTrueList[iCut]->Add(hMCTrueCombConversionMidPtEta[iCut]);
				hMCTrueCombConversionPt[iCut]               = new TH1F("ESD_TrueCombConversionMapping_Pt","ESD_TrueCombConversionMapping_Pt",nBinsPt,0.,20.);
				fTrueList[iCut]->Add(hMCTrueCombConversionPt[iCut]);
				

		}
	}
	
	
    fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
	if(fV0Reader && fV0Reader->GetProduceV0FindingEfficiency())
		if (fV0Reader->GetV0FindingEfficiencyHistograms())
			fOutputList->Add(fV0Reader->GetV0FindingEfficiencyHistograms());
	
  if(fV0Reader)
    if((AliConvEventCuts*)fV0Reader->GetEventCuts())
      if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
        fOutputList->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());

  if(fV0Reader)
    if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
      if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
        fOutputList->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());

   for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
    }
    if(!((AliConversionPhotonCuts*)fConversionCutArray->At(iCut))) continue;
    if(((AliConversionPhotonCuts*)fConversionCutArray->At(iCut))->GetCutHistograms()){
      fCutFolder[iCut]->Add(((AliConversionPhotonCuts*)fConversionCutArray->At(iCut))->GetCutHistograms());
    }
	 }

	PostData(1, fOutputList);
   
}

//________________________________________________________________________
void AliAnalysisTaskMaterialHistos::UserExec(Option_t *){
	
	if(fIsMC>0) fMCEvent = MCEvent();
	if(fMCEvent == NULL) fIsMC = 0;

	fInputEvent = InputEvent();
	if (fInputEvent==NULL) return;

  if(fIsMC>0 && fInputEvent->IsA()==AliESDEvent::Class()){
    fMCStack = fMCEvent->Stack();
    if(fMCStack == NULL) fIsMC = 0;
  }
	
	
	Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(fInputEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or because it is incomplere or  wrong trigger for V0ReaderV1

			for(Int_t iCut = 0; iCut<fnCuts; iCut++){
				hNEvents[iCut]->Fill(eventQuality);
			}

    return;
  }

	

	fConversionGammas=fV0Reader->GetReconstructedGammas();// Gammas from default Cut
	
  // ------------------- BeginEvent ----------------------------

  for(Int_t iCut = 0; iCut<fnCuts; iCut++){
    fiCut = iCut;
		Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
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


		if(fIsMC > 0){
			// Process MC Particle
		  if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection() != 0){
		    if(fInputEvent->IsA()==AliESDEvent::Class()){
		      ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
											     ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
											     fMCEvent);
		    }
		  }
		}

		if(fIsMC>0){
		  ProcessMCPhotons();
		}


		fNESDtracksEta09 = CountTracks09(); // Estimate Event Multiplicity
		fNESDtracksEta0914 = CountTracks0914(); // Estimate Event Multiplicity
		fNESDtracksEta14 = fNESDtracksEta09 + fNESDtracksEta0914;


			hNGoodESDTracks09[iCut]->Fill(fNESDtracksEta09);
			hNGoodESDTracks14[iCut]->Fill(fNESDtracksEta14);
			hNGoodESDTracks09_14[iCut]->Fill(fNESDtracksEta0914);


		if(fInputEvent){
			if(fInputEvent->GetPrimaryVertexTracks()->GetNContributors()>0) {
				fNContrVtx = fInputEvent->GetPrimaryVertexTracks()->GetNContributors();
			} else {
				fNContrVtx = 0;
			}	

		}
		//		fPrimVtxZ = fInputEvent->GetPrimaryVertex()->GetZ();

		ProcessPhotons();
		fGammaCandidates->Clear(); // delete this cuts good gammas
	}
	
	//cout<<" done with the event"<<endl;

	PostData(1, fOutputList);
}

// ///________________________________________________________________________
void AliAnalysisTaskMaterialHistos::FillMCTree(Int_t stackPos){
	//	AliStack *MCStack = fMCEvent->Stack();
	TParticle* candidate = (TParticle *)fMCStack->Particle(stackPos);
	
	if(((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->PhotonIsSelectedMC(candidate,fMCStack,kFALSE)){
		fGammaMCPt = candidate->Pt();
		fGammaMCTheta = candidate->Theta();

		hMCAllGammaPt[fiCut]->Fill(candidate->Pt());                          

	}
	Double_t minPt=0.4;
	Double_t maxPt=1.5;

	if(((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->PhotonIsSelectedMC(candidate,fMCStack,kTRUE)){
		fGammaMCConvPt = candidate->Pt();
		fGammaMCConvTheta = candidate->Theta();
		TParticle* daughter1 = (TParticle *)fMCStack->Particle(candidate->GetFirstDaughter()); 
		//		TParticle* daughter2 = (TParticle *)fMCStack->Particle(candidate->GetLastDaughter()); 

			//		hMCConversionMappingRZ[fiCut]->Fill(daughter1->Vz(),daughter1->R());         
			hMCConversionMappingRPhi[fiCut]->Fill(candidate->Phi(),daughter1->R());      
			hMCConversionEta[fiCut]->Fill(candidate->Eta());             
			hMCConversionPt[fiCut]->Fill(candidate->Pt());                          
			hMCConversionR[fiCut]->Fill(daughter1->R());      
			if (fGammaMCConvPt>minPt && fGammaMCConvPt<maxPt){
				hMCConversionMidPtR[fiCut]->Fill(daughter1->R());    
			}

		
	} // Converted MC Gamma		
}
// 
// ///________________________________________________________________________
void AliAnalysisTaskMaterialHistos::ProcessMCPhotons(){
	// Loop over all primary MC particle
	//	AliStack *ffMCStack = fMCEvent->Stack();
	for(Int_t i = 0; i < fMCStack->GetNprimary(); i++) {
		TParticle* particle = (TParticle *)fMCStack->Particle(i);
		if (!particle) continue;
		
		
		if(fMCStack && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
			Int_t isPosFromMBHeader
				= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent);
			Int_t isNegFromMBHeader
				= ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent);
			if( (isNegFromMBHeader < 1) || (isPosFromMBHeader < 1)) continue;
		}
		
		
		if (particle->GetPdgCode() == 111 && particle->GetFirstDaughter() >= fMCStack->GetNprimary()){
			// 			cout << "Undecayed pi0 found with mother: " << particle->GetMother(0) << endl;
			for (Int_t j = 0; j < 2 ; j++){
 				FillMCTree(particle->GetDaughter(j));
			}
		} else {
			FillMCTree(i);
		}	
		
	}	
}

///________________________________________________________________________
void AliAnalysisTaskMaterialHistos::ProcessPhotons(){
	
	// Fill Histograms for QA and MC
  TList *GammaCandidatesStepTwo = new TList();

  //not needed here already calculated
  //fNESDtracksEta09 = CountTracks09();
  for(Int_t firstGammaIndex=0;firstGammaIndex<fConversionGammas->GetEntriesFast();firstGammaIndex++){
    AliAODConversionPhoton *gamma= (AliAODConversionPhoton*)fConversionGammas->At(firstGammaIndex);
    
    if (gamma == NULL) continue;
    
    if(!((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->PhotonIsSelected(gamma,fInputEvent))continue;
    
    if( ! ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->UseToCloseV0sCut()){
      fGammaCandidates->Add(gamma); // if no second loop is required add to events good gammas
    }else if(((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->UseToCloseV0sCut()) { // shared electron is disabled, step one not needed -> step two
      GammaCandidatesStepTwo->Add(gamma);
    }
  }
  
  if(((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->UseToCloseV0sCut()){
    for(Int_t i = 0;i<GammaCandidatesStepTwo->GetEntries();i++){
      AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GammaCandidatesStepTwo->At(i);
      if(!PhotonCandidate) continue;
      if(!((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GammaCandidatesStepTwo,i)) continue;
      fGammaCandidates->Add(PhotonCandidate); // Add gamma to current cut TList
    }
  }
  

 

  for(Int_t firstGammaIndex=0;firstGammaIndex<fGammaCandidates->GetEntries()-1;firstGammaIndex++){
    AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fGammaCandidates->At(firstGammaIndex));
    if (gamma==NULL) continue;

    fGammaPt = gamma->GetPhotonPt();
    fGammaTheta = gamma->GetPhotonTheta();
    fGammaChi2NDF = gamma->GetChi2perNDF();
    
  
    AliVTrack * negTrack = ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetTrack(fInputEvent, gamma->GetTrackLabelNegative());
    AliVTrack * posTrack = ((AliConversionPhotonCuts*)fConversionCutArray->At(fiCut))->GetTrack(fInputEvent, gamma->GetTrackLabelPositive());
    
    
    Double_t asym=0.;
    if(gamma->GetPhotonP()!=0){
      asym=negTrack->P()/gamma->GetPhotonP();
    }
    fKind = 9;	
  
    if(fIsMC>0){
    // 			cout << "generating MC stack"<< endl;
      AliStack *fMCStack = fMCEvent->Stack();
      if (!fMCStack) continue;
      
      const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
      Double_t mcProdVtxX 	= primVtxMC->GetX();
      Double_t mcProdVtxY 	= primVtxMC->GetY();
      Double_t mcProdVtxZ 	= primVtxMC->GetZ();
      
      TParticle *posDaughter = gamma->GetPositiveMCDaughter(fMCStack);
      TParticle *negDaughter = gamma->GetNegativeMCDaughter(fMCStack);
      // 			cout << "generate Daughters: "<<posDaughter << "\t" << negDaughter << endl;
		
      if(fMCStack && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
	Int_t isPosFromMBHeader
	  = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma->GetMCLabelPositive(), fMCStack, fInputEvent);
	Int_t isNegFromMBHeader
	  = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(gamma->GetMCLabelNegative(), fMCStack, fInputEvent);
	if( (isNegFromMBHeader < 1) || (isPosFromMBHeader < 1)) continue;
      }
      
      if(posDaughter == NULL || negDaughter == NULL){ 
	fKind = 9; // garbage
	// 				cout << "one of the daughters not available" << endl;
      } else if(posDaughter->GetMother(0) != negDaughter->GetMother(0) || (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1)){ 
	// Not Same Mother == Combinatorial Bck
	fKind = 1;
	// 				cout << "not the same mother" << endl;
	Int_t pdgCodePos; 
	if (posDaughter->GetPdgCode()) pdgCodePos = posDaughter->GetPdgCode(); else continue;
	Int_t pdgCodeNeg; 
	if (negDaughter->GetPdgCode()) pdgCodeNeg = negDaughter->GetPdgCode(); else continue;
	// 				cout << "PDG codes daughters: " << pdgCodePos << "\t" << pdgCodeNeg << endl;
	if(TMath::Abs(pdgCodePos)==11 && TMath::Abs(pdgCodeNeg)==11)
	  fKind = 10; //Electron Combinatorial
	if(TMath::Abs(pdgCodePos)==11 && TMath::Abs(pdgCodeNeg)==11 && (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1))
	  fKind = 15; //direct Electron Combinatorial
	
	if(TMath::Abs(pdgCodePos)==211 && TMath::Abs(pdgCodeNeg)==211)
	  fKind = 11; //Pion Combinatorial
	if((TMath::Abs(pdgCodePos)==211 && TMath::Abs(pdgCodeNeg)==2212) ||
	   (TMath::Abs(pdgCodePos)==2212 && TMath::Abs(pdgCodeNeg)==211))
	  fKind = 12; //Pion, Proton Combinatorics
	if((TMath::Abs(pdgCodePos)==211 && TMath::Abs(pdgCodeNeg)==11) ||
	   (TMath::Abs(pdgCodePos)==11 && TMath::Abs(pdgCodeNeg)==211))
	  fKind = 13; //Pion, Electron Combinatorics
	if (TMath::Abs(pdgCodePos)==321 || TMath::Abs(pdgCodeNeg)==321)	
	  fKind = 14; //Kaon combinatorics
	
      } else {
	// 				cout << "same mother" << endl;
	Int_t pdgCodePos; 
	if (posDaughter->GetPdgCode()) pdgCodePos = posDaughter->GetPdgCode(); else continue;
	Int_t pdgCodeNeg; 
	if (negDaughter->GetPdgCode()) pdgCodeNeg = negDaughter->GetPdgCode(); else continue;
	// 				cout << "PDG codes daughters: " << pdgCodePos << "\t" << pdgCodeNeg << endl;
	Int_t pdgCode; 
	if (gamma->GetMCParticle(fMCStack)->GetPdgCode()) pdgCode = gamma->GetMCParticle(fMCStack)->GetPdgCode(); else continue;
	// 				cout << "PDG code: " << pdgCode << endl;
	if(TMath::Abs(pdgCodePos)!=11 || TMath::Abs(pdgCodeNeg)!=11)
	  fKind = 2; // combinatorics from hadronic decays
	else if ( !(pdgCodeNeg==pdgCodePos)){
	  TParticle *truePhotonCanditate = gamma->GetMCParticle(fMCStack);
	  Int_t motherLabelPhoton = truePhotonCanditate->GetMother(0);
	  Bool_t gammaIsPrimary = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsConversionPrimaryESD( fMCStack, posDaughter->GetMother(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
	  if(pdgCode == 111) 
	    fKind = 3; // pi0 Dalitz
	  else if (pdgCode == 221) 
	    fKind = 4; // eta Dalitz
	  else if (!(negDaughter->GetUniqueID() != 5 || posDaughter->GetUniqueID() !=5)){
	    if(pdgCode == 22 && gammaIsPrimary){
	      fKind = 0; // primary photons
	    } else if (pdgCode == 22){
	      fKind = 5; //secondary photons
	    }		
	  } else 	fKind = 9; //garbage
	} else fKind = 9; //garbage
      }					
    }
    
    
    
    
    Double_t minPt=0.4;
    Double_t maxPt=1.5;
    

    hESDConversionMappingRPhi[fiCut]->Fill(gamma->GetPhotonPhi(),gamma->GetConversionRadius());  
    hESDConversionMappingRZ[fiCut]->Fill(gamma->GetConversionZ(),gamma->GetConversionRadius());  
    hESDConversionEta[fiCut]->Fill(gamma->GetPhotonEta());              
    hESDConversionPt[fiCut]->Fill(gamma->GetPhotonPt());
    if(gamma->GetConversionRadius()>5){
      hESDConversionPt5cm[fiCut]->Fill(gamma->GetPhotonPt());    
    }
    hESDConversionR[fiCut]->Fill(gamma->GetConversionRadius());
    hESDConversionAsymP[fiCut]->Fill(gamma->GetPhotonP(),asym);
    
    if(fInputEvent->IsA()==AliESDEvent::Class()){
      AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
      if(esdEvent){
	AliESDv0 *v0 = esdEvent->GetV0(gamma->GetV0Index());
	hESDConversionDCA[fiCut]->Fill(v0->GetDcaV0Daughters());
	if(gamma->GetPhotonPt()>minPt && gamma->GetPhotonPt()<maxPt){
	  hESDConversionMidPtDCA[fiCut]->Fill(v0->GetDcaV0Daughters());         
	}
      }
    }
    
    
    hESDConversionPsiPair[fiCut]->Fill(gamma->GetPsiPair()); 
    hESDConversionChi2[fiCut]->Fill(gamma->GetChi2perNDF());    
    hESDConversionMass[fiCut]->Fill(gamma->GetInvMassPair());     
    
    
    if(gamma->GetPhotonPt()>minPt && gamma->GetPhotonPt()<maxPt){
      hESDConversionMidPtEta[fiCut]->Fill(gamma->GetPhotonEta());         
      hESDConversionMidPtR[fiCut]->Fill(gamma->GetConversionRadius());
    }
    
    Double_t convR = gamma->GetConversionRadius();
    hESDConversionRInBins[fiCut]->Fill(fGammaPt,fNESDtracksEta09,convR);
    hESDConversionPhiInBins[fiCut]->Fill(convR,fNESDtracksEta09,gamma->GetPhotonPhi());
    
    if(fIsMC>0){
      if(fKind==0 || fKind==5){
	hMCTrueConversionMappingRPhi[fiCut]->Fill(gamma->GetPhotonPhi(),gamma->GetConversionRadius());       
	hMCTrueConversionMappingRZ[fiCut]->Fill(gamma->GetConversionZ(),gamma->GetConversionRadius());       
	hMCTrueConversionEta[fiCut]->Fill(gamma->GetPhotonEta());                  
	hMCTrueConversionPt[fiCut]->Fill(gamma->GetPhotonPt()); 
	if(gamma->GetConversionRadius()>5){
	  hMCTrueConversionPt5cm[fiCut]->Fill(gamma->GetPhotonPt());    
	}
	hMCTrueConversionR[fiCut]->Fill(gamma->GetConversionRadius());  
	hMCTrueConversionAsymP[fiCut]->Fill(gamma->GetPhotonP(),asym);
	hMCTrueConversionPsiPair[fiCut]->Fill(gamma->GetPsiPair()); 
	hMCTrueConversionChi2[fiCut]->Fill(gamma->GetChi2perNDF());    
	hMCTrueConversionMass[fiCut]->Fill(gamma->GetInvMassPair());     
	if(fInputEvent->IsA()==AliESDEvent::Class()){
	  AliESDEvent *esdEvent = dynamic_cast<AliESDEvent*>(fInputEvent);
	  if(esdEvent){
	  AliESDv0 *v0 = esdEvent->GetV0(gamma->GetV0Index());
	  hMCTrueConversionDCA[fiCut]->Fill(v0->GetDcaV0Daughters());
	  if(gamma->GetPhotonPt()>minPt && gamma->GetPhotonPt()<maxPt){
	    hMCTrueConversionMidPtDCA[fiCut]->Fill(v0->GetDcaV0Daughters());         
	  }
	  }
	}
	
	
	if(gamma->GetPhotonPt()>minPt && gamma->GetPhotonPt()<maxPt){
	  hMCTrueConversionMidPtR[fiCut]->Fill(gamma->GetConversionRadius());      
	  hMCTrueConversionMidPtEta[fiCut]->Fill(gamma->GetPhotonEta());          
	}
      }
      
      if(fKind==3){
	hMCTruePi0DalConversionEta[fiCut]->Fill(gamma->GetPhotonEta());                    
	hMCTruePi0DalConversionPt[fiCut]->Fill(gamma->GetPhotonPt());               
	hMCTruePi0DalConversionR[fiCut]->Fill(gamma->GetConversionRadius());      
	if(gamma->GetPhotonPt()>minPt && gamma->GetPhotonPt()<maxPt){
	  hMCTruePi0DalConversionMidPtEta[fiCut]->Fill(gamma->GetPhotonEta());            
	  hMCTruePi0DalConversionMidPtR[fiCut]->Fill(gamma->GetConversionRadius());      
	}
      }
      if(fKind==4){
	hMCTrueEtaDalConversionEta[fiCut]->Fill(gamma->GetPhotonEta());                 
	hMCTrueEtaDalConversionPt[fiCut]->Fill(gamma->GetPhotonPt());
	hMCTrueEtaDalConversionR[fiCut]->Fill(gamma->GetConversionRadius());                     
	if(gamma->GetPhotonPt()>minPt && gamma->GetPhotonPt()<maxPt){
	  hMCTrueEtaDalConversionMidPtEta[fiCut]->Fill(gamma->GetPhotonEta());              
	  hMCTrueEtaDalConversionMidPtR[fiCut]->Fill(gamma->GetConversionRadius());      
	}
      }
      if(fKind==1 || fKind==10 || fKind==11 || fKind==12 || fKind==13 || fKind==14 || fKind==15){
	hMCTrueCombConversionEta[fiCut]->Fill(gamma->GetPhotonEta());                 
	hMCTrueCombConversionPt[fiCut]->Fill(gamma->GetPhotonPt()); 
	hMCTrueCombConversionR[fiCut]->Fill(gamma->GetConversionRadius());                     
	if(gamma->GetPhotonPt()>minPt && gamma->GetPhotonPt()<maxPt){
	  hMCTrueCombConversionMidPtEta[fiCut]->Fill(gamma->GetPhotonEta());                
	  hMCTrueCombConversionMidPtR[fiCut]->Fill(gamma->GetConversionRadius());      
	}
      }   
    }
  }

  delete GammaCandidatesStepTwo;
  GammaCandidatesStepTwo = 0x0;
  
  
}

//________________________________________________________________________
Int_t AliAnalysisTaskMaterialHistos::CountTracks09(){
	Int_t fNumberOfESDTracks = 0;
	if(fInputEvent->IsA()==AliESDEvent::Class()){
	// Using standard function for setting Cuts
		
		AliStack *fMCStack = NULL;
		if (fMCEvent){
			fMCStack= fMCEvent->Stack();
			if (!fMCStack) return 0;
		}	
				
		Bool_t selectPrimaries=kTRUE;
               AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
		EsdTrackCuts->SetMaxDCAToVertexZ(2);
		EsdTrackCuts->SetEtaRange(-0.9, 0.9);
		EsdTrackCuts->SetPtRange(0.15);
		
		for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
			AliESDtrack* curTrack = (AliESDtrack*) fInputEvent->GetTrack(iTracks);
			if(!curTrack) continue;
			if(EsdTrackCuts->AcceptTrack(curTrack) ){
				if (fMCEvent){
					if(fMCStack && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
						Int_t isFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(abs(curTrack->GetLabel()), fMCStack, fInputEvent);
						if( (isFromMBHeader < 1) ) continue;
					}					
				}	
				fNumberOfESDTracks++;
			}	
		}
		delete EsdTrackCuts;
		EsdTrackCuts=0x0;
	} 
	return fNumberOfESDTracks;
}

//________________________________________________________________________
Int_t AliAnalysisTaskMaterialHistos::CountTracks0914(){
	Int_t fNumberOfESDTracks = 0;
	if(fInputEvent->IsA()==AliESDEvent::Class()){
		// Using standard function for setting Cuts
		
		AliStack *fMCStack = NULL;
		if (fMCEvent){
			fMCStack= fMCEvent->Stack();
			if (!fMCStack) return 0;
		}	

		AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
		EsdTrackCuts->SetMaxDCAToVertexZ(5);
		EsdTrackCuts->SetEtaRange(0.9, 1.4);
		EsdTrackCuts->SetPtRange(0.15);
		
		for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
			AliESDtrack* curTrack = (AliESDtrack*) fInputEvent->GetTrack(iTracks);
			if(!curTrack) continue;
			if(EsdTrackCuts->AcceptTrack(curTrack) ){
				if (fMCEvent){
					if(fMCStack && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
						Int_t isFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(abs(curTrack->GetLabel()), fMCStack, fInputEvent);
						if( (isFromMBHeader < 1) ) continue;
					}					
				}	
				fNumberOfESDTracks++;
			}
		}
		EsdTrackCuts->SetEtaRange(-1.4, -0.9);
		for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
			AliESDtrack* curTrack =(AliESDtrack*) fInputEvent->GetTrack(iTracks);
			if(!curTrack) continue;
			if(EsdTrackCuts->AcceptTrack(curTrack) ){
				if (fMCEvent){
					if(fMCStack && ((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
						Int_t isFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(abs(curTrack->GetLabel()), fMCStack, fInputEvent);
						if( (isFromMBHeader < 1) ) continue;
					}					
				}	
				fNumberOfESDTracks++;
			}	
		}
		delete EsdTrackCuts;
		EsdTrackCuts=0x0;
	} 
	
	return fNumberOfESDTracks;
}

//________________________________________________________________________
void AliAnalysisTaskMaterialHistos::Terminate(Option_t *)
{
//    if (fStreamMaterial){
//       fStreamMaterial->GetFile()->Write();
//    }
}
