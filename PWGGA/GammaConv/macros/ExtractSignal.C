// provided by Gamma Conversion Group, $ALICE_ROOT/PWGGA/GammaConv ;https://twiki.cern.ch/twiki/bin/view/ALICE/PWG4GammaConversion

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TApplication.h>
#include <TPaveLabel.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TString.h>
#include "TGaxis.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TDatabasePDG.h"
#include "PlottingMeson.h"

TString prefix;
TString prefix2;
TString suffix;
TFile* Output;
TFile* Output1;
TFile* Output2;

TH1D*  fBckNorm;
TH1D*  fSignal;
TF1 * fReco;
TF1 * fGausExp;
TF1 * fLinearBck;

Double_t fYields;
Double_t fYieldsError;

Double_t fFWHMFunc;
Double_t fFWHMFuncError;
Double_t fYieldsFunc;
Double_t fYieldsFuncError;
Double_t fIntLinearBck;
Double_t fIntLinearBckError;

Int_t fIsMC=0;
Double_t fYMaxMeson=0.9;
TString fCutSelection;
Int_t fNRebinGlobal = 2; 

// common meson analysis variables (So it is easy to switch between the pi0 and the eta without changing the code)
Int_t iPt0 = 0;
Int_t column = 0; 
Int_t row = 0;	
Int_t NBinsPt = 0;	
Double_t * BinsPt = NULL; 
Double_t * BinsMt = NULL;
Int_t * fNRebin = NULL;
Double_t *BGFitRange = NULL;
Double_t *BGFitRangeLeft = NULL;
Double_t *MesonPlotRange = NULL;
Double_t *MesonIntRange = NULL;
Double_t *MesonIntRangeWide = NULL;
Double_t *MesonIntRangeNarrow = NULL;

Double_t *MesonMassRange = NULL;
Double_t *MesonFitRange = NULL;
Double_t MesonMassExp=0;   // expected meson mass
Double_t MesonId=0;
Double_t MesonWidthExp=0;
Double_t MesonLambdaTail=0;
Double_t *MesonWidthRange = NULL;
Double_t *MesonLambdaTailRange = NULL;
// end common meson analysis variables

///////////////////////////

Double_t nEvt;
TH1F* fNumberOfGoodESDTracksVtx;

TString histonameGG;
TString histonameBack;
TString histonameTrue;
TString histoNameEffi;

Double_t* GGYields;
Double_t* BckYields;
Double_t* MesonYields;
Double_t* MesonTrueYields;
Double_t* MesonYieldsFunc;
Double_t* MesonYieldsResidualBckFunc;
Double_t* MesonYieldsCorResidualBckFunc;
Double_t* MesonYieldsPerEvent;
Double_t* MesonMass;
Double_t* MesonWidth;
Double_t* MesonSB;
Double_t* MesonSign;
Double_t* MesonFWHM;

// Normalization at the left of the peak
Double_t* GGYieldsLeft;
Double_t* BckYieldsLeft;
Double_t* MesonYieldsLeft;
Double_t* MesonYieldsFuncLeft;
Double_t* MesonYieldsResidualBckFuncLeft;
Double_t* MesonYieldsCorResidualBckFuncLeft;
Double_t* MesonYieldsLeftPerEvent;
Double_t* MesonMassLeft;
Double_t* MesonWidthLeft;
Double_t* MesonSBLeft;
Double_t* MesonSignLeft;
Double_t* MesonFWHMLeft;

// Narrow Integration Window
Double_t* GGYieldsNarrow;
Double_t* BckYieldsNarrow;
Double_t* MesonYieldsNarrow;
Double_t* MesonTrueYieldsNarrow;
Double_t* MesonYieldsFuncNarrow;
Double_t* MesonYieldsResidualBckFuncNarrow;
Double_t* MesonYieldsCorResidualBckFuncNarrow;
Double_t* MesonYieldsPerEventNarrow;
Double_t* MesonSBNarrow;
Double_t* MesonSignNarrow;

Double_t* GGYieldsLeftNarrow;
Double_t* BckYieldsLeftNarrow;
Double_t* MesonYieldsLeftNarrow;
Double_t* MesonYieldsFuncLeftNarrow;
Double_t* MesonYieldsResidualBckFuncLeftNarrow;
Double_t* MesonYieldsCorResidualBckFuncLeftNarrow;
Double_t* MesonYieldsLeftPerEventNarrow;
Double_t* MesonSBLeftNarrow;
Double_t* MesonSignLeftNarrow;

// Wide Integration Window
Double_t* GGYieldsWide;
Double_t* BckYieldsWide;
Double_t* MesonYieldsWide;
Double_t* MesonTrueYieldsWide;
Double_t* MesonYieldsFuncWide;
Double_t* MesonYieldsResidualBckFuncWide;
Double_t* MesonYieldsCorResidualBckFuncWide;
Double_t* MesonYieldsPerEventWide;
Double_t* MesonSBWide;
Double_t* MesonSignWide;

Double_t* GGYieldsLeftWide;
Double_t* BckYieldsLeftWide;
Double_t* MesonYieldsLeftWide;
Double_t* MesonYieldsFuncLeftWide;
Double_t* MesonYieldsResidualBckFuncLeftWide;
Double_t* MesonYieldsCorResidualBckFuncLeftWide;
Double_t* MesonYieldsLeftPerEventWide;
Double_t* MesonSBLeftWide;
Double_t* MesonSignLeftWide;

Double_t* GGYieldsError;
Double_t* BckYieldsError;
Double_t* MesonYieldsError;
Double_t* MesonYieldsFuncError;
Double_t* MesonYieldsResidualBckFuncError;
Double_t* MesonYieldsCorResidualBckFuncError;
Double_t* MesonYieldsPerEventError;
Double_t* MesonMassError;
Double_t* MesonWidthError;
Double_t* MesonSBError;
Double_t* MesonSignError;
Double_t* MesonFWHMError;

Double_t* GGYieldsLeftError;
Double_t* BckYieldsLeftError;
Double_t* MesonYieldsLeftError;
Double_t* MesonYieldsFuncLeftError;
Double_t* MesonYieldsResidualBckFuncLeftError;
Double_t* MesonYieldsCorResidualBckFuncLeftError;
Double_t* MesonYieldsLeftPerEventError;
Double_t* MesonMassLeftError;
Double_t* MesonWidthLeftError;
Double_t* MesonSBLeftError;
Double_t* MesonSignLeftError;
Double_t* MesonFWHMLeftError;

// Narrow integration Window
Double_t* GGYieldsNarrowError;
Double_t* BckYieldsNarrowError;
Double_t* MesonYieldsNarrowError;
Double_t* MesonYieldsFuncNarrowError;
Double_t* MesonYieldsResidualBckFuncNarrowError;
Double_t* MesonYieldsCorResidualBckFuncNarrowError;
Double_t* MesonYieldsPerEventNarrowError;
Double_t* MesonSBNarrowError;
Double_t* MesonSignNarrowError;

Double_t* GGYieldsLeftNarrowError;
Double_t* BckYieldsLeftNarrowError;
Double_t* MesonYieldsLeftNarrowError;
Double_t* MesonYieldsFuncLeftNarrowError;
Double_t* MesonYieldsResidualBckFuncLeftNarrowError;
Double_t* MesonYieldsCorResidualBckFuncLeftNarrowError;
Double_t* MesonYieldsLeftPerEventNarrowError;
Double_t* MesonSBLeftNarrowError;
Double_t* MesonSignLeftNarrowError;

// Wide integration Window
Double_t* GGYieldsWideError;
Double_t* BckYieldsWideError;
Double_t* MesonYieldsWideError;
Double_t* MesonYieldsFuncWideError;
Double_t* MesonYieldsResidualBckFuncWideError;
Double_t* MesonYieldsCorResidualBckFuncWideError;
Double_t* MesonYieldsPerEventWideError;
Double_t* MesonSBWideError;
Double_t* MesonSignWideError;

Double_t* GGYieldsLeftWideError;
Double_t* BckYieldsLeftWideError;
Double_t* MesonYieldsLeftWideError;
Double_t* MesonYieldsFuncLeftWideError;
Double_t* MesonYieldsResidualBckFuncLeftWideError;
Double_t* MesonYieldsCorResidualBckFuncLeftWideError;
Double_t* MesonYieldsLeftPerEventWideError;
Double_t* MesonSBLeftWideError;
Double_t* MesonSignLeftWideError;

TH1D* Mapping_BackNorm_InvMass;
TH1D* Mapping_Signal_InvMass;   

TH1D** Mapping_TrueMeson_InvMass_PtBin;    // array of histos for pt slices
TH1D** Mapping_GG_InvMass_PtBin;    // array of histos for pt slices
TH1D** Mapping_Back_InvMass_PtBin;
TH1D** Mapping_BackNorm_InvMass_PtBin;
TH1D** Mapping_Signal_InvMass_PtBin;   
TF1** Signal_InvMassFit_PtBin;   
TF1** SignalPeak_InvMassFit_PtBin;   
TF1** SignalBck_InvMassFit_PtBin;   

TH1F *fHistoYieldMeson;
TH1F *fHistoYieldMesonPerEvent;
TH1F *fHistoYieldMesonMt;
TH1F *fHistoYieldMesonMtPerEvent;
TH1F *fHistoSignMeson;
TH1F *fHistoSBMeson;

TH1F *fHistoYieldMesonNarrow;
TH1F *fHistoYieldMesonPerEventNarrow;
TH1F *fHistoYieldMesonMtNarrow;
TH1F *fHistoYieldMesonMtPerEventNarrow;
TH1F *fHistoSignMesonNarrow;
TH1F *fHistoSBMesonNarrow;

TH1F *fHistoYieldMesonWide;
TH1F *fHistoYieldMesonPerEventWide;
TH1F *fHistoYieldMesonMtWide;
TH1F *fHistoYieldMesonMtPerEventWide;
TH1F *fHistoSignMesonWide;
TH1F *fHistoSBMesonWide;

TH1F *fHistoMassMeson;
TH1F *fHistoWidthMeson;
TH1F *fHistoFWHMMeson;
TH1F *fdeltaPt;
TH1F *fdeltaMt;

// Histograms for normalization on the left of the peak
TH1D** Mapping_BackNorm_InvMass_Left_PtBin;
TH1D** Mapping_Signal_InvMass_Left_PtBin;   

TF1** Signal_InvMassFit_Left_PtBin;   
TF1** SignalPeak_InvMassFit_Left_PtBin;   
TF1** SignalBck_InvMassFit_Left_PtBin;   

TH1F *fHistoYieldMesonLeft;
TH1F *fHistoYieldMesonLeftPerEvent;
TH1F *fHistoYieldMesonLeftMt;
TH1F *fHistoYieldMesonLeftMtPerEvent;
TH1F *fHistoSignMesonLeft;
TH1F *fHistoSBMesonLeft;

TH1F *fHistoYieldMesonLeftNarrow;
TH1F *fHistoYieldMesonLeftPerEventNarrow;
TH1F *fHistoYieldMesonLeftMtNarrow;
TH1F *fHistoYieldMesonLeftMtPerEventNarrow;
TH1F *fHistoSignMesonLeftNarrow;
TH1F *fHistoSBMesonLeftNarrow;

TH1F *fHistoYieldMesonLeftWide;
TH1F *fHistoYieldMesonLeftPerEventWide;
TH1F *fHistoYieldMesonLeftMtWide;
TH1F *fHistoYieldMesonLeftMtPerEventWide;
TH1F *fHistoSignMesonLeftWide;
TH1F *fHistoSBMesonLeftWide;

TH1F *fHistoMassMesonLeft;
TH1F *fHistoWidthMesonLeft;
TH1F *fHistoFWHMMesonLeft;

TH2F *fMC_Meson_Pt_Eta_withinAcceptance;
TH1D *fMC_Meson_Pt;
 
TH1D * fMC_Meson_Pt_WithinAcceptance;
TH1D * fMC_MesonWithinAccepPt; // Proper bins in Pt
TH1D * fMC_MesonPt; // Proper bins in Pt

TH1F * fHistoYieldTrueMeson;
TH1F * fHistoYieldTrueMesonWide;
TH1F * fHistoYieldTrueMesonNarrow;

TH1D * fMC_MesonAccepPt;
TH2F * fTrueMesonInvMassVSPt;

TH1D * fMCMesonEffiPt;
TH1D * fMCTrueMesonEffiPt;

TH1D * fMC_MesonEffiPt;
TH1D * fMC_MesonNarrowEffiPt;
TH1D * fMC_MesonWideEffiPt;

TH1D * fMC_MesonLeftEffiPt;
TH1D * fMC_MesonLeftNarrowEffiPt;
TH1D * fMC_MesonLeftWideEffiPt;

TH1D * fMC_TrueMesonEffiPt;
TH1D * fMC_TrueMesonNarrowEffiPt;
TH1D * fMC_TrueMesonWideEffiPt;

void ProcessEM(TH1D*,TH1D*,Double_t *);
void FillMassHistosArray(TH2F*,TH2F*);
void FillMassMCTrueMesonHistosArray(TH2F*);
void CreatePtHistos();
void FillPtHistos();
void PlotPtHistos(TString,TString);
//void PlotInvMassInPtBins(TH1D * [NBinsPt],TH1D * [NBinsPt], TString, TString, TString);
//void PlotWithFitSubtractedInvMassInPtBins(TH1D * [NBinsPt], TF1 * [NBinsPt], TString , TString , TString); // Plots and fits the spectra
//void FitAllSubtractedInvMassInPtBins(TH1D * [NBinsPt]); // Only fits the spectra
void PlotInvMassInPtBins(TH1D **,TH1D **, TString, TString, TString);
void PlotWithFitSubtractedInvMassInPtBins(TH1D **, TF1 **, TString , TString , TString); // Plots and fits the spectra
void FitAllSubtractedInvMassInPtBins(TH1D **); // Only fits the spectra
void FitSubtractedInvMassInPtBins(TH1D * ,Double_t *);  // Fits the Invariant Mass histos with a given function
void IntegrateHistoInvMass(TH1D * , Double_t * );
void IntegrateFitFunc(TF1 * , TH1D *, Double_t *);   
void FillHistosArrayMC(TH2F* , TH1D*, TH1F* );
void CalculateMesonAcceptance();
void CalculateMesonEfficiency(TH1F*, TString);
void SaveHistos(Int_t, TString, TString);
void SaveHistosPtBins(TH1D * ,TString , TString);
void SaveCorrectionHistos(Int_t , TString , TString);
void Initialize(TString setPi0);
void CalculateFWHM(TF1 * , Double_t *);

// Main Function
void ExtractSignal(TString meson, TString file, TString cutSelection, TString suffix, TString isMC){

	//gROOT->Reset();
	//gROOT->SetStyle("Plain");
	prefix=meson;
	if(meson.CompareTo("Pi0") == 0){// means we want to plot values for the pi0
		Initialize("Pi0");
	}
	else if(meson.CompareTo("Eta") == 0){
		Initialize("Eta");
	}
	else if(meson.CompareTo("Pi0EtaBinning") == 0) {
		Initialize("Pi0EtaBinning");
	}
	else	{
		cout<<"ERROR: First argument in the ExtractSignal(....) has to be either Pi0 or Eta"<<endl;
		return;
	}

	if(cutSelection.Length() == 0){
		cout<<"ERROR: Cut selection is not set, please do!"<<endl;
		return;
	}

	if(isMC.CompareTo("kTRUE") == 0){ 
		fIsMC = 1; 
		prefix2 = "MC";
	} else { 
		fIsMC = 0;
		prefix2 = "data";
	}
  
	StyleSettings();

	TFile f(file.Data());
	fCutSelection = cutSelection;

	TDirectory *pwg4dir =(TDirectory*)f.Get(Form("PWGGA_GammaConversion_%s",fCutSelection.Data()));
	if(pwg4dir == NULL){
		cout<<Form("ERROR: PWGGA_GammaConversion_%s",fCutSelection.Data())<<" is not found in the file"<<endl;
		return;
	}
	TList *fHistosGammaConversion = (TList*)pwg4dir->Get(Form("histogramsAliGammaConversion_%s",fCutSelection.Data()));
	TList *fESDContainer = (TList*)fHistosGammaConversion->FindObject("ESD histograms");
	TList *fBackgroundContainer = (TList*)fHistosGammaConversion->FindObject("Background histograms");


	fNumberOfGoodESDTracksVtx = (TH1F*)fESDContainer->FindObject("ESD_NumberOfGoodESDTracksVtx");
	TH1D *fBck = (TH1D*)fBackgroundContainer->FindObject("ESD_Background_InvMass");
	TH1D *fGammaGamma = (TH1D*)fESDContainer->FindObject("ESD_Mother_InvMass");

	TH2F *fGammaGammaInvMassVSPt= (TH2F*)fESDContainer->FindObject("ESD_Mother_InvMass_vs_Pt");
	TH2F *fBckInvMassVSPt = (TH2F*)fBackgroundContainer->FindObject("ESD_Background_InvMass_vs_Pt");

	
	if(fIsMC){
		TList *fMCContainer = (TList*)fHistosGammaConversion->FindObject("MC histograms");
		if( MesonId == 111){
			fMC_Meson_Pt_Eta_withinAcceptance = (TH2F*)fMCContainer->FindObject("MC_Pi0_Pt_Eta_withinAcceptance");   // we should move to rapidity,also the cut !!!
			fMC_Meson_Pt = (TH1D*)fMCContainer->FindObject("MC_Pi0_Pt");   // Not the best; better having a 2D Pt_vs_Rapid in case we change limits 
		}
		if( MesonId == 221){
			fMC_Meson_Pt_Eta_withinAcceptance = (TH2F*)fMCContainer->FindObject("MC_Eta_Pt_Eta_withinAcceptance");   // we should move to rapidity,also the cut !!!
			fMC_Meson_Pt = (TH1D*)fMCContainer->FindObject("MC_Eta_Pt");   // Not the best; better having a 2D Pt_vs_Rapid in case we change limits 
		}

		// True Meson 
		fTrueMesonInvMassVSPt = (TH2F*)fESDContainer->FindObject("ESD_TruePi0_InvMass_vs_Pt");
		FillMassMCTrueMesonHistosArray(fTrueMesonInvMassVSPt);
	}
	
	MesonMassExp = TDatabasePDG::Instance()->GetParticle(MesonId)->Mass();

	nEvt = fNumberOfGoodESDTracksVtx->GetEntries();
	
	cout<< "The mass of the meson is::"<< MesonMassExp<< " Events analysed"<< nEvt<< endl;

	// Process the 1D invariant mass histos
	fGammaGamma->SetTitle(Form("%s %s",fGammaGamma->GetTitle(),fCutSelection.Data()));
	fGammaGamma->Rebin(fNRebinGlobal);
	fBck->Rebin(fNRebinGlobal);
	ProcessEM( fGammaGamma , fBck, BGFitRange);   
	Mapping_BackNorm_InvMass = fBckNorm;
	Mapping_Signal_InvMass = fSignal;   

	fGammaGamma->DrawCopy();
	Mapping_BackNorm_InvMass->DrawCopy("same");
	Mapping_Signal_InvMass->DrawCopy("same");

	// Function to Project the 2D histos InvariantMass VS Pt into Invariant Mass spectrum 
	FillMassHistosArray(fGammaGammaInvMassVSPt,fBckInvMassVSPt);
	
	for(Int_t iPt=iPt0;iPt<NBinsPt;iPt++){ // BEGIN ANALYSIS for each Pt bin
			// Function to subtract GG minus Bck
		ProcessEM( Mapping_GG_InvMass_PtBin[iPt], Mapping_Back_InvMass_PtBin[iPt], BGFitRange);   
		Mapping_Signal_InvMass_PtBin[iPt] = fSignal;
		Mapping_BackNorm_InvMass_PtBin[iPt] = fBckNorm;

		// Integrate the 2g histo
		IntegrateHistoInvMass( Mapping_GG_InvMass_PtBin[iPt], MesonIntRange);
		GGYields[iPt] = fYields;
		GGYieldsError[iPt] = fYieldsError;
		IntegrateHistoInvMass( Mapping_GG_InvMass_PtBin[iPt], MesonIntRangeWide);
		GGYieldsWide[iPt] = fYields;
		GGYieldsWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( Mapping_GG_InvMass_PtBin[iPt], MesonIntRangeNarrow);
		GGYieldsNarrow[iPt] = fYields;
		GGYieldsNarrowError[iPt] = fYieldsError;

		// Integrate the bck histo
		IntegrateHistoInvMass( Mapping_BackNorm_InvMass_PtBin[iPt], MesonIntRange);
		BckYields[iPt] = fYields;
		BckYieldsError[iPt] = fYieldsError;
		IntegrateHistoInvMass( Mapping_BackNorm_InvMass_PtBin[iPt], MesonIntRangeWide);
		BckYieldsWide[iPt] = fYields;
		BckYieldsWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( Mapping_BackNorm_InvMass_PtBin[iPt], MesonIntRangeNarrow);
		BckYieldsNarrow[iPt] = fYields;
		BckYieldsNarrowError[iPt] = fYieldsError;

		// Integrate the signal histo
		IntegrateHistoInvMass( Mapping_Signal_InvMass_PtBin[iPt], MesonIntRange);
		MesonYields[iPt] = fYields;
		MesonYieldsError[iPt] = fYieldsError;
		IntegrateHistoInvMass( Mapping_Signal_InvMass_PtBin[iPt], MesonIntRangeWide);
		MesonYieldsWide[iPt] = fYields;
		MesonYieldsWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( Mapping_Signal_InvMass_PtBin[iPt], MesonIntRangeNarrow);
		MesonYieldsNarrow[iPt] = fYields;
		MesonYieldsNarrowError[iPt] = fYieldsError;
		if(fIsMC){
			IntegrateHistoInvMass( Mapping_TrueMeson_InvMass_PtBin[iPt], MesonIntRange);
			MesonTrueYields[iPt] = fYields;
			IntegrateHistoInvMass( Mapping_TrueMeson_InvMass_PtBin[iPt], MesonIntRangeWide);
			MesonTrueYieldsWide[iPt] = fYields;
			IntegrateHistoInvMass( Mapping_TrueMeson_InvMass_PtBin[iPt], MesonIntRangeNarrow);
			MesonTrueYieldsNarrow[iPt] = fYields;
		}

		cout<< "iPt"<< iPt<< " "<< "standard range"<<endl;
		// Fitting the subtracted spectra
		FitSubtractedInvMassInPtBins(Mapping_Signal_InvMass_PtBin[iPt], MesonIntRange);

		Signal_InvMassFit_PtBin[iPt] = fReco;   
		SignalPeak_InvMassFit_PtBin[iPt] = fGausExp;   
		SignalBck_InvMassFit_PtBin[iPt] = fLinearBck;   
		MesonMass[iPt] = Signal_InvMassFit_PtBin[iPt]->GetParameter(1);
		MesonMassError[iPt] = Signal_InvMassFit_PtBin[iPt]->GetParError(1);
		MesonWidth[iPt] = Signal_InvMassFit_PtBin[iPt]->GetParameter(2);
		MesonWidthError[iPt] = Signal_InvMassFit_PtBin[iPt]->GetParError(2);

		MesonYieldsResidualBckFunc[iPt] = fIntLinearBck;
		MesonYieldsResidualBckFuncError[iPt] = fIntLinearBckError;
		MesonYieldsCorResidualBckFunc[iPt] = MesonYields[iPt]- MesonYieldsResidualBckFunc[iPt];
		MesonYieldsCorResidualBckFuncError[iPt] = 
			pow((MesonYieldsError[iPt]*MesonYieldsError[iPt]+ 
			MesonYieldsResidualBckFuncError[iPt]*MesonYieldsResidualBckFuncError[iPt]),0.5);
		MesonYieldsPerEvent[iPt]= MesonYieldsCorResidualBckFunc[iPt]/nEvt;
		MesonYieldsPerEventError[iPt]= MesonYieldsCorResidualBckFuncError[iPt]/nEvt;

		if( BckYields[iPt]!=0){
			MesonSB[iPt] = MesonYieldsCorResidualBckFunc[iPt]/BckYields[iPt];
		}else{
			MesonSB[iPt] = 0.;
		}
		MesonSign[iPt] = MesonYieldsCorResidualBckFunc[iPt]/pow(BckYields[iPt],0.5);


		//Integrate Fit Function
		IntegrateFitFunc( SignalPeak_InvMassFit_PtBin[iPt], Mapping_Signal_InvMass_PtBin[iPt], MesonIntRange);
		MesonYieldsFunc[iPt]=fYieldsFunc;

		//GetFWHM
		CalculateFWHM( Signal_InvMassFit_PtBin[iPt],MesonIntRange);
		MesonFWHM[iPt] = fFWHMFunc;
		MesonFWHMError[iPt] = fFWHMFuncError;
		
		cout<< "iPt"<< iPt<< " "<< "FWHM done"<<endl;
		
		
		// Wide integration mass window
		cout<< "iPt"<< iPt<< " "<< "wide range"<<endl;
		FitSubtractedInvMassInPtBins(Mapping_Signal_InvMass_PtBin[iPt],MesonIntRangeWide);
		MesonYieldsResidualBckFuncWide[iPt] = fIntLinearBck;
		MesonYieldsResidualBckFuncWideError[iPt] = fIntLinearBckError;
		MesonYieldsCorResidualBckFuncWide[iPt] = MesonYieldsWide[iPt]- MesonYieldsResidualBckFuncWide[iPt];
		MesonYieldsCorResidualBckFuncWideError[iPt] = 
			pow((MesonYieldsWideError[iPt]*MesonYieldsWideError[iPt]+ 
			MesonYieldsResidualBckFuncWideError[iPt]*MesonYieldsResidualBckFuncWideError[iPt]),0.5);
		MesonYieldsPerEventWide[iPt]= MesonYieldsCorResidualBckFuncWide[iPt]/nEvt;
		MesonYieldsPerEventWideError[iPt]= MesonYieldsCorResidualBckFuncWideError[iPt]/nEvt;


		// Narrow integration mass window
		cout<< "iPt"<< iPt<< " "<< "narrow range"<<endl;
		FitSubtractedInvMassInPtBins(Mapping_Signal_InvMass_PtBin[iPt],MesonIntRangeNarrow);
		MesonYieldsResidualBckFuncNarrow[iPt] = fIntLinearBck;
		MesonYieldsResidualBckFuncNarrowError[iPt] = fIntLinearBckError;
		MesonYieldsCorResidualBckFuncNarrow[iPt] = MesonYieldsNarrow[iPt]- MesonYieldsResidualBckFuncNarrow[iPt];
		MesonYieldsCorResidualBckFuncNarrowError[iPt] = 
			pow((MesonYieldsNarrowError[iPt]*MesonYieldsNarrowError[iPt]+ 
			MesonYieldsResidualBckFuncNarrowError[iPt]*MesonYieldsResidualBckFuncNarrowError[iPt]),0.5);
		MesonYieldsPerEventNarrow[iPt]= MesonYieldsCorResidualBckFuncNarrow[iPt]/nEvt;
		MesonYieldsPerEventNarrowError[iPt]= MesonYieldsCorResidualBckFuncNarrowError[iPt]/nEvt;


		//////////////////////////////// Start Analysis with  Normalization at the left of the Meson Peak

		// Function to subtract GG minus Bck
		ProcessEM( Mapping_GG_InvMass_PtBin[iPt], Mapping_Back_InvMass_PtBin[iPt], BGFitRangeLeft);   
		Mapping_Signal_InvMass_Left_PtBin[iPt] = fSignal;
		Mapping_BackNorm_InvMass_Left_PtBin[iPt] = fBckNorm;

		// Integrate the 2g histo  , Not needed in the GG nothing changes


		// Integrate the bck histo
		IntegrateHistoInvMass( Mapping_BackNorm_InvMass_Left_PtBin[iPt], MesonIntRange);
		BckYieldsLeft[iPt] = fYields;
		BckYieldsLeftError[iPt] = fYieldsError;
		IntegrateHistoInvMass( Mapping_BackNorm_InvMass_Left_PtBin[iPt], MesonIntRangeWide);
		BckYieldsLeftWide[iPt] = fYields;
		BckYieldsLeftWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( Mapping_BackNorm_InvMass_Left_PtBin[iPt], MesonIntRangeNarrow);
		BckYieldsLeftNarrow[iPt] = fYields;
		BckYieldsLeftNarrowError[iPt] = fYieldsError;
		
		// Integrate the signal histo
		IntegrateHistoInvMass( Mapping_Signal_InvMass_Left_PtBin[iPt], MesonIntRange);
		MesonYieldsLeft[iPt] = fYields;
		MesonYieldsLeftError[iPt] = fYieldsError;
		IntegrateHistoInvMass( Mapping_Signal_InvMass_Left_PtBin[iPt], MesonIntRangeWide);
		MesonYieldsLeftWide[iPt] = fYields;
		MesonYieldsLeftWideError[iPt] = fYieldsError;
		IntegrateHistoInvMass( Mapping_Signal_InvMass_Left_PtBin[iPt], MesonIntRangeNarrow);
		MesonYieldsLeftNarrow[iPt] = fYields;
		MesonYieldsLeftNarrowError[iPt] = fYieldsError;

		cout<< "iPt"<< iPt<< " "<< "standard range"<<endl;
		// Fitting the subtracted spectra
		FitSubtractedInvMassInPtBins(Mapping_Signal_InvMass_Left_PtBin[iPt], MesonIntRange);

		Signal_InvMassFit_Left_PtBin[iPt] = fReco;   
		SignalPeak_InvMassFit_Left_PtBin[iPt] = fGausExp;   
		SignalBck_InvMassFit_Left_PtBin[iPt] = fLinearBck;   
		MesonMassLeft[iPt] = Signal_InvMassFit_Left_PtBin[iPt]->GetParameter(1);
		MesonMassLeftError[iPt] = Signal_InvMassFit_Left_PtBin[iPt]->GetParError(1);
		MesonWidthLeft[iPt] = Signal_InvMassFit_Left_PtBin[iPt]->GetParameter(2);
		MesonWidthLeftError[iPt] = Signal_InvMassFit_Left_PtBin[iPt]->GetParError(2);

		MesonYieldsResidualBckFuncLeft[iPt] = fIntLinearBck;
		MesonYieldsResidualBckFuncLeftError[iPt] = fIntLinearBckError;
		MesonYieldsCorResidualBckFuncLeft[iPt] = MesonYieldsLeft[iPt]- MesonYieldsResidualBckFuncLeft[iPt];
		MesonYieldsCorResidualBckFuncLeftError[iPt] = 
			pow((MesonYieldsLeftError[iPt]*MesonYieldsLeftError[iPt]+ 
			MesonYieldsResidualBckFuncLeftError[iPt]*MesonYieldsResidualBckFuncLeftError[iPt]),0.5);
		MesonYieldsLeftPerEvent[iPt]= MesonYieldsCorResidualBckFuncLeft[iPt]/nEvt;
		MesonYieldsLeftPerEventError[iPt]= MesonYieldsCorResidualBckFuncLeftError[iPt]/nEvt;

		if( BckYieldsLeft[iPt]!=0){
			MesonSBLeft[iPt] = MesonYieldsCorResidualBckFuncLeft[iPt]/BckYieldsLeft[iPt];
		}else{
			MesonSBLeft[iPt] = 0.;
		}
		MesonSignLeft[iPt] = MesonYieldsCorResidualBckFuncLeft[iPt]/pow(BckYieldsLeft[iPt],0.5);


		//Integrate Fit Function
		IntegrateFitFunc( SignalPeak_InvMassFit_Left_PtBin[iPt], Mapping_Signal_InvMass_Left_PtBin[iPt], MesonIntRange);
		MesonYieldsFuncLeft[iPt]=fYieldsFunc;

		//GetFWHM
		CalculateFWHM(Signal_InvMassFit_Left_PtBin[iPt],MesonIntRange);
		cout<< "iPt left "<< iPt<< " "<< "B"<<endl;
		cout << fFWHMFunc << "     " << fFWHMFuncError << endl;
		MesonFWHMLeft[iPt] = fFWHMFunc;
		cout<< "iPt left "<< iPt<< " "<< "C"<<endl;
		MesonFWHMLeftError[iPt] = fFWHMFuncError;
		
		
		// Wide integration mass window
		cout<< "iPt"<< iPt<< " "<< "wide range"<<endl;
		FitSubtractedInvMassInPtBins(Mapping_Signal_InvMass_Left_PtBin[iPt],MesonIntRangeWide);
		MesonYieldsResidualBckFuncLeftWide[iPt] = fIntLinearBck;
		MesonYieldsResidualBckFuncLeftWideError[iPt] = fIntLinearBckError;
		MesonYieldsCorResidualBckFuncLeftWide[iPt] = MesonYieldsLeftWide[iPt]- MesonYieldsResidualBckFuncLeftWide[iPt];
		MesonYieldsCorResidualBckFuncLeftWideError[iPt] = 
			pow((MesonYieldsLeftWideError[iPt]*MesonYieldsLeftWideError[iPt]+ 
			MesonYieldsResidualBckFuncLeftWideError[iPt]*MesonYieldsResidualBckFuncLeftWideError[iPt]),0.5);
		MesonYieldsLeftPerEventWide[iPt]= MesonYieldsCorResidualBckFuncLeftWide[iPt]/nEvt;
		MesonYieldsLeftPerEventWideError[iPt]= MesonYieldsCorResidualBckFuncLeftWideError[iPt]/nEvt;


		// Narrow integration mass window
		cout<< "iPt"<< iPt<< " "<< "narrow range"<<endl;
		FitSubtractedInvMassInPtBins(Mapping_Signal_InvMass_Left_PtBin[iPt],MesonIntRangeNarrow);
		MesonYieldsResidualBckFuncLeftNarrow[iPt] = fIntLinearBck;
		MesonYieldsResidualBckFuncLeftNarrowError[iPt] = fIntLinearBckError;
		MesonYieldsCorResidualBckFuncLeftNarrow[iPt] = MesonYieldsLeftNarrow[iPt]- MesonYieldsResidualBckFuncLeftNarrow[iPt];
		MesonYieldsCorResidualBckFuncLeftNarrowError[iPt] = 
			pow((MesonYieldsLeftNarrowError[iPt]*MesonYieldsLeftNarrowError[iPt]+ 
			MesonYieldsResidualBckFuncLeftNarrowError[iPt]*MesonYieldsResidualBckFuncLeftNarrowError[iPt]),0.5);
		MesonYieldsLeftPerEventNarrow[iPt]= MesonYieldsCorResidualBckFuncLeftNarrow[iPt]/nEvt;
		MesonYieldsLeftPerEventNarrowError[iPt]= MesonYieldsCorResidualBckFuncLeftNarrowError[iPt]/nEvt;
	}
	//******************** Data OUTPUTFILE ***************************************************
	const char* SysErrDatname = Form("%s_%s_SystematicErrorYieldExtraction_RAWDATA_%s.dat",prefix.Data(),prefix2.Data(),cutSelection.Data());
	fstream SysErrDat;
	SysErrDat.open(SysErrDatname, ios::out);
	SysErrDat << "Calculation of the systematic error due to the yield extraction RAWDATA " << endl;
	SysErrDat <<  endl;		
	SysErrDat << "GGYields" << endl;	
	SysErrDat << "Bin \t Right \t Right Wide \t Right Narr" << endl;
		for(Int_t iPt=iPt0;iPt<NBinsPt;iPt++){
			SysErrDat << iPt << "\t" << GGYields[iPt] << "+-" << GGYieldsError[iPt] << "\t" <<
			GGYieldsWide[iPt] << "+-" << GGYieldsWideError[iPt] << "\t" <<
			GGYieldsNarrow[iPt] << "+-" << GGYieldsNarrowError[iPt] << endl;
			
		}
	SysErrDat <<  endl;		
	SysErrDat << "BckYields" << endl;	
	SysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
		for(Int_t iPt=iPt0;iPt<NBinsPt;iPt++){
			SysErrDat << iPt << "\t" <<
			BckYields[iPt] << "+-" << BckYieldsError[iPt] << "\t" <<
			BckYieldsWide[iPt] << "+-" << BckYieldsWideError[iPt] << "\t" <<
			BckYieldsNarrow[iPt] << "+-" << BckYieldsNarrowError[iPt] << "\t" <<
			BckYieldsLeft[iPt] << "+-" << BckYieldsLeftError[iPt]<< "\t" <<
			BckYieldsLeftWide[iPt]<< "+-" << BckYieldsLeftWideError[iPt]<< "\t" <<
			BckYieldsLeftNarrow[iPt]<< "+-" << BckYieldsLeftNarrowError[iPt] << endl;
		}
	SysErrDat <<  endl;		
	SysErrDat << "MesonYields" << endl;	
	SysErrDat << "Bin \t Right \t Right Wide \t Right Narr \t Left \t Left Wide \t Left Narr" << endl;
		for(Int_t iPt=iPt0;iPt<NBinsPt;iPt++){
			SysErrDat << iPt << "\t" <<
			MesonYields[iPt] << "+-" << MesonYieldsError[iPt] << "\t" <<
			MesonYieldsWide[iPt] << "+-" << MesonYieldsWideError[iPt] << "\t" <<
			MesonYieldsNarrow[iPt] << "+-" << MesonYieldsNarrowError[iPt] << "\t" <<
			MesonYieldsLeft[iPt]<< "+-" << MesonYieldsLeftError[iPt]<< "\t" <<
			MesonYieldsLeftWide[iPt]<< "+-" << MesonYieldsLeftWideError[iPt]<< "\t" <<
			MesonYieldsLeftNarrow[iPt]<< "+-" << MesonYieldsLeftNarrowError[iPt] << endl;
		}
	if(fIsMC){
		SysErrDat <<  endl;		
		SysErrDat << "TrueYields" << endl;	
		SysErrDat << "Bin \t True \t True Wide \t True Narr " << endl;
		for(Int_t iPt=iPt0;iPt<NBinsPt;iPt++){
			SysErrDat << iPt << "\t" <<
			MesonTrueYields[iPt] << "\t" <<
			MesonTrueYieldsWide[iPt] << "\t" <<
			MesonTrueYieldsNarrow[iPt] << endl;
		}
	}	
	SysErrDat.close();
	//******************************** OUTPUT END ******************************************************
	
	TString nameMeson = Form("%s/%s_%s_MesonWithBck_%s.%s",suffix.Data(),prefix.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data());
	TString nameCanvas = "MesonWithBckCanvas";
	TString namePad = "MesonWithBckPad";
	PlotInvMassInPtBins( Mapping_GG_InvMass_PtBin, Mapping_BackNorm_InvMass_PtBin, nameMeson, nameCanvas, namePad);	
	
	TString nameMesonSub= Form("%s/%s_%s_MesonSubtracted_%s.%s",suffix.Data(),prefix.Data(),prefix2.Data(), fCutSelection.Data(),suffix.Data());
	TString nameCanvasSub= "MesonCanvasSubtracted";
	TString namePadSub= "MesonPadSubtracted";
	PlotWithFitSubtractedInvMassInPtBins( Mapping_Signal_InvMass_PtBin, Signal_InvMassFit_PtBin, nameMesonSub, nameCanvasSub, namePadSub);

	nameMeson= Form("%s/%s_%s_MesonWithBckLeft_%s.%s",suffix.Data(),prefix.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data());
	nameCanvas = "MesonWithBckCanvasLeft";
	namePad = "MesonWithBckPadLeft";
	PlotInvMassInPtBins( Mapping_GG_InvMass_PtBin, Mapping_BackNorm_InvMass_Left_PtBin, nameMeson, nameCanvas, namePad);

	nameMesonSub= Form("%s/%s_%s_MesonSubtractedLeft_%s.%s",suffix.Data(),prefix.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data());
	nameCanvasSub= "MesonCanvasSubtractedLeft";
	namePadSub= "MesonPadSubtractedLeft";
	PlotWithFitSubtractedInvMassInPtBins( Mapping_Signal_InvMass_Left_PtBin, Signal_InvMassFit_Left_PtBin, nameMesonSub, nameCanvasSub, namePadSub);

	CreatePtHistos();
	FillPtHistos();
	PlotPtHistos(suffix,prefix2);

	if(fIsMC){
	FillHistosArrayMC(fMC_Meson_Pt_Eta_withinAcceptance,fMC_Meson_Pt, fdeltaPt);
	CalculateMesonAcceptance();
	cout << "Calculated MesonAcceptance" << endl;

	histoNameEffi="MesonEffiPt";
	CalculateMesonEfficiency(fHistoYieldMeson,histoNameEffi);
	fMC_MesonEffiPt = fMCMesonEffiPt;
	
	histoNameEffi="MesonWideEffiPt";
	CalculateMesonEfficiency(fHistoYieldMesonWide,histoNameEffi);
	fMC_MesonWideEffiPt = fMCMesonEffiPt;
	
	histoNameEffi="MesonNarrowEffiPt";
	CalculateMesonEfficiency(fHistoYieldMesonNarrow,histoNameEffi);
	fMC_MesonNarrowEffiPt = fMCMesonEffiPt;
	
	histoNameEffi="MesonLeftEffiPt";
	CalculateMesonEfficiency(fHistoYieldMesonLeft,histoNameEffi);
	fMC_MesonLeftEffiPt = fMCMesonEffiPt;
	
	histoNameEffi="MesonLeftNarrowEffiPt";
	CalculateMesonEfficiency(fHistoYieldMesonLeftNarrow,histoNameEffi);
	fMC_MesonLeftNarrowEffiPt = fMCMesonEffiPt;
	
	histoNameEffi="MesonLeftWideEffiPt";
	CalculateMesonEfficiency(fHistoYieldMesonLeftWide,histoNameEffi);
	fMC_MesonLeftWideEffiPt = fMCMesonEffiPt;
	// cout<< "True effi-wide::"<< fMC_MesonLeftWideEffiPt<< " "<<fMCMesonEffiPt << endl;
	
	// True Meson (only once case, because no normalization
	histoNameEffi="TrueMesonEffiPt";
	CalculateMesonEfficiency(fHistoYieldTrueMeson,histoNameEffi);
	fMC_TrueMesonEffiPt = fMCMesonEffiPt;
	// cout<< "True effi::"<< fMC_TrueMesonEffiPt<< " "<<fMCMesonEffiPt << endl;
	histoNameEffi="TrueMesonWideEffiPt";
	CalculateMesonEfficiency(fHistoYieldTrueMesonWide,histoNameEffi);
	fMC_TrueMesonWideEffiPt = fMCMesonEffiPt;
	
	histoNameEffi="TrueMesonNarrowEffiPt";
	CalculateMesonEfficiency(fHistoYieldTrueMesonNarrow,histoNameEffi);
	fMC_TrueMesonNarrowEffiPt = fMCMesonEffiPt;
	
/*		new TCanvas;
		fMC_TrueMesonEffiPt->DrawCopy();
		fMC_TrueMesonWideEffiPt->DrawCopy("same");
		fMC_TrueMesonNarrowEffiPt->DrawCopy("same");
		fMC_MesonEffiPt->DrawCopy("same");
		fMC_MesonWideEffiPt->DrawCopy("same");
		fMC_MesonNarrowEffiPt->DrawCopy("same");
		
		new TCanvas;
		fMC_TrueMesonEffiPt->DrawCopy();
		fMC_TrueMesonWideEffiPt->DrawCopy("same");
		fMC_TrueMesonNarrowEffiPt->DrawCopy("same");
		
		fMC_MesonLeftEffiPt->DrawCopy("same");
		fMC_MesonLeftWideEffiPt->DrawCopy("same");
		fMC_MesonLeftNarrowEffiPt->DrawCopy("same");*/
	
	SaveCorrectionHistos(fIsMC, fCutSelection, prefix2);
	}


	SaveHistos(fIsMC, fCutSelection, prefix2);

	//   const char* Outputname = Form("%s_%s_AnalysisResultsWithoutCorrectionDifferentBins_%s.root",prefix.Data(),prefix2.Data(),fCutSelection.Data());
	//   Output = new TFile(Outputname,"RECREATE");		
	//   cout << "file created" << endl;
	//   for (Int_t BinPt = iPt0; BinPt < NBinsPt; BinPt++){
	//     Mapping_Back_InvMass_PtBin[BinPt]->Write();
	//     cout << "saved Mapping_GG_InvMass_PtBin[" << BinPt<< "]" << endl; 
	//     Mapping_BackNorm_InvMass_PtBin[BinPt]->Write();
	//     cout << "saved Mapping_BackNorm_InvMass_PtBin[" << BinPt<< "]" << endl; 
	//     Mapping_Signal_InvMass_PtBin[BinPt]->Write();
	//     cout << "saved Mapping_Signal_InvMass_PtBin[" << BinPt<< "]" << endl; 
	//     Signal_InvMassFit_PtBin[BinPt]->Write();
	//     cout << "saved Signal_InvMassFit_PtBin[" << BinPt<< "]" << endl; 
	//     SignalPeak_InvMassFit_PtBin[BinPt]->Write();
	//     cout << "saved SignalPeak_InvMassFit_PtBin[" << BinPt<< "]" << endl; 
	//     SignalBck_InvMassFit_PtBin[BinPt]->Write();
	//     cout << "saved SignalBck_InvMassFit_PtBin[" << BinPt<< "]" << endl; 
	//     Mapping_BackNorm_InvMass_Left_PtBin[BinPt]->Write();
	//     cout << "saved Mapping_BackNorm_InvMass_Left_PtBin[" << BinPt<< "]" << endl; 
	//     Mapping_Signal_InvMass_Left_PtBin[BinPt]->Write();
	//     cout << "saved Mapping_Signal_InvMass_Left_PtBin[" << BinPt<< "]" << endl; 
	//     Signal_InvMassFit_Left_PtBin[BinPt]->Write();
	//     cout << "saved Signal_InvMassFit_Left_PtBin[" << BinPt<< "]" << endl; 
	//     SignalPeak_InvMassFit_Left_PtBin[BinPt]->Write();
	//     cout << "saved SignalPeak_InvMassFit_Left_PtBin[" << BinPt<< "]" << endl; 
	//     SignalBck_InvMassFit_Left_PtBin[BinPt]->Write();
	//     cout << "saved SignalBck_InvMassFit_Left_PtBin[" << BinPt<< "]" << endl; 
	//   }
	//   Output->Write();
	//   Output->Close();
	//   delete Output;

}



void ProcessEM(TH1D* fGammaGamma, TH1D* fBck, Double_t * BGFitRange)
{
  fBckNorm = (TH1D*)fBck->Clone("fBckNorm");
  fGammaGamma->Sumw2();
  fBck->Sumw2();
  fBckNorm->Sumw2();


  Double_t r= fGammaGamma->Integral(fGammaGamma->GetXaxis()->FindBin(BGFitRange[0]),fGammaGamma->GetXaxis()->FindBin(BGFitRange[1]));
  Double_t b= fBck->Integral(fBck->GetXaxis()->FindBin(BGFitRange[0]),fBck->GetXaxis()->FindBin(BGFitRange[1]));
 

  Double_t norm = 1;
  if(b != 0) norm = r/b;
  fBckNorm->Scale(norm);
  cout<<"r="<<r<<" b="<<b<<" r/b="<<r/b<< " " << endl;
  fSignal = (TH1D*)fGammaGamma->Clone("fSignal");
  fSignal->Sumw2();
  fSignal->Add(fBckNorm,-1.);
  //cout<<"fBckNorm= "<<fBckNorm<<" fSignal="<<fSignal<<endl;
}

			 
void FillMassHistosArray(TH2F* fGammaGammaInvMassVSPt,TH2F *fBckInvMassVSPt)
{

  for(Int_t iPt=iPt0;iPt<NBinsPt;iPt++){
    histonameGG = Form("Mapping_GG_InvMass_in_Pt_Bin%02d", iPt);
    histonameBack = Form("Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
    //    sprintf(histonameGG,"Mapping_GG_InvMass_in_Pt_Bin%02d", iPt);			
    // sprintf(histonameBack,"Mapping_Back_InvMass_in_Pt_Bin%02d", iPt);
    if(Mapping_GG_InvMass_PtBin[iPt]!= NULL){
      delete Mapping_GG_InvMass_PtBin[iPt];
      Mapping_GG_InvMass_PtBin[iPt]=NULL;
    }
    Mapping_GG_InvMass_PtBin[iPt]=new TH1D(histonameGG.Data(),histonameGG.Data(),fGammaGammaInvMassVSPt->GetNbinsX(),0.,1.);

    if(Mapping_Back_InvMass_PtBin[iPt]!= NULL){
      delete Mapping_Back_InvMass_PtBin[iPt];
      Mapping_Back_InvMass_PtBin[iPt]=NULL;
    }
    Mapping_Back_InvMass_PtBin[iPt]=new TH1D(histonameBack.Data(),histonameBack.Data(),fBckInvMassVSPt->GetNbinsX(),0.,1.);
  
  
    Int_t startBin = fGammaGammaInvMassVSPt->GetYaxis()->FindBin(BinsPt[iPt]+0.001);
    Int_t endBin = fGammaGammaInvMassVSPt->GetYaxis()->FindBin(BinsPt[iPt+1]-0.001);

    
    cout<< "bins::"<< startBin<< " " << endBin<<" "<< BinsPt[iPt]<< " "<<BinsPt[iPt+1]<<  endl;
    cout<< "bin values::"<< fGammaGammaInvMassVSPt->GetYaxis()->GetBinCenter(startBin)<< " "
	<< fGammaGammaInvMassVSPt->GetYaxis()->GetBinCenter(endBin)<< endl;

    fGammaGammaInvMassVSPt->ProjectionX(histonameGG.Data(),startBin,endBin);
    fBckInvMassVSPt->ProjectionX(histonameBack.Data(),startBin,endBin);
    
    Mapping_GG_InvMass_PtBin[iPt]=(TH1D*)gDirectory->Get(histonameGG.Data());

    if(fNRebin[iPt]>1){
      Mapping_GG_InvMass_PtBin[iPt]->Rebin(fNRebin[iPt]);
    }
    Mapping_Back_InvMass_PtBin[iPt]=(TH1D*)gDirectory->Get(histonameBack.Data());
    if(fNRebin[iPt]>1){
      Mapping_Back_InvMass_PtBin[iPt]->Rebin(fNRebin[iPt]);
    }
  }

}
void FillMassMCTrueMesonHistosArray(TH2F* fTrueMesonInvMassVSPt)
{
  for(Int_t iPt=iPt0;iPt<NBinsPt;iPt++){
    histonameTrue = Form("Mapping_TrueMeson_InvMass_in_Pt_Bin%02d", iPt);
    //    sprintf(histonameTrue,"Mapping_TrueMeson_InvMass_in_Pt_Bin%02d", iPt);			
    if(Mapping_TrueMeson_InvMass_PtBin[iPt]!= NULL){
      delete Mapping_TrueMeson_InvMass_PtBin[iPt];
      Mapping_TrueMeson_InvMass_PtBin[iPt]=NULL;
    }
    Mapping_TrueMeson_InvMass_PtBin[iPt] = new TH1D(histonameTrue.Data(),histonameTrue.Data(),fTrueMesonInvMassVSPt->GetNbinsX(),0.,1.);

    Int_t startBin = fTrueMesonInvMassVSPt->GetYaxis()->FindBin(BinsPt[iPt]+0.001);
    Int_t endBin = fTrueMesonInvMassVSPt->GetYaxis()->FindBin(BinsPt[iPt+1]-0.001);

    cout<< "bins::"<< startBin<< " " << endBin<<" "<< BinsPt[iPt]<< " "<<BinsPt[iPt+1]<<  endl;
    cout<< "bin values::"<< fTrueMesonInvMassVSPt->GetYaxis()->GetBinCenter(startBin)<< " "
	<< fTrueMesonInvMassVSPt->GetYaxis()->GetBinCenter(endBin)<< endl;

    fTrueMesonInvMassVSPt->ProjectionX(histonameTrue.Data(),startBin,endBin);
    Mapping_TrueMeson_InvMass_PtBin[iPt]=(TH1D*)gDirectory->Get(histonameTrue.Data());
    if(fNRebin[iPt]>1){
      Mapping_TrueMeson_InvMass_PtBin[iPt]->Rebin(fNRebin[iPt]);
    }
    Mapping_TrueMeson_InvMass_PtBin[iPt]->SetLineColor(2);
  }

}


void CreatePtHistos(){

  for(Int_t iPt=iPt0;iPt<NBinsPt+1;iPt++){
    BinsMt[iPt]=pow((BinsPt[iPt]*BinsPt[iPt]+MesonMassExp*MesonMassExp),0.5);
    cout<<iPt<<" "<< BinsPt[iPt]<< " " << BinsMt[iPt]<< endl;
  }

  fdeltaPt = new TH1F("deltaPt","",NBinsPt,BinsPt);
  fdeltaMt = new TH1F("deltaMt","",NBinsPt,BinsMt);

  fHistoYieldMeson = new TH1F("histoYieldMeson","",NBinsPt,BinsPt);
  fHistoYieldTrueMeson = new TH1F("histoYieldTrueMeson","",NBinsPt,BinsPt);
  fHistoYieldMesonPerEvent = new TH1F("histoYieldMesonPerEvent","",NBinsPt,BinsPt);
  fHistoYieldMesonMt = new TH1F("histoYieldMesonMt","",NBinsPt,BinsMt);
  fHistoYieldMesonMtPerEvent = new TH1F("histoYieldMesonMtPerEvent","",NBinsPt,BinsMt);
  fHistoSignMeson = new TH1F("histoSignMeson","",NBinsPt,BinsPt);
  fHistoSBMeson = new TH1F("histoSBMeson","",NBinsPt,BinsPt);
  fHistoMassMeson = new TH1F("histoMassMeson","",NBinsPt,BinsPt);
  fHistoWidthMeson = new TH1F("histoWidthMeson","",NBinsPt,BinsPt);
  fHistoFWHMMeson = new TH1F("histoFWHMMeson","",NBinsPt,BinsPt);
 

  fHistoYieldMesonNarrow = new TH1F("histoYieldMesonNarrow","",NBinsPt,BinsPt);
  fHistoYieldTrueMesonNarrow = new TH1F("histoYieldTrueMesonNarrow","",NBinsPt,BinsPt);
  fHistoYieldMesonPerEventNarrow = new TH1F("histoYieldMesonPerEventNarrow","",NBinsPt,BinsPt);
  fHistoYieldMesonMtNarrow = new TH1F("histoYieldMesonMtNarrow","",NBinsPt,BinsMt);
  fHistoYieldMesonMtPerEventNarrow = new TH1F("histoYieldMesonMtPerEventNarrow","",NBinsPt,BinsMt);
  fHistoSignMesonNarrow = new TH1F("histoSignMesonNarrow","",NBinsPt,BinsPt);
  fHistoSBMesonNarrow = new TH1F("histoSBMesonNarrow","",NBinsPt,BinsPt);


  fHistoYieldMesonWide = new TH1F("histoYieldMesonWide","",NBinsPt,BinsPt);
  fHistoYieldTrueMesonWide = new TH1F("histoYieldTrueMesonWide","",NBinsPt,BinsPt);
  fHistoYieldMesonPerEventWide = new TH1F("histoYieldMesonPerEventWide","",NBinsPt,BinsPt);
  fHistoYieldMesonMtWide = new TH1F("histoYieldMesonMtWide","",NBinsPt,BinsMt);
  fHistoYieldMesonMtPerEventWide = new TH1F("histoYieldMesonMtPerEventWide","",NBinsPt,BinsMt);
  fHistoSignMesonWide = new TH1F("histoSignMesonWide","",NBinsPt,BinsPt);
  fHistoSBMesonWide = new TH1F("histoSBMesonWide","",NBinsPt,BinsPt);

  // Histos for normalization at the left of the peak

  fHistoYieldMesonLeft = new TH1F("histoYieldMesonLeft","",NBinsPt,BinsPt);
  fHistoYieldMesonLeftPerEvent = new TH1F("histoYieldMesonLeftPerEvent","",NBinsPt,BinsPt);
  fHistoYieldMesonLeftMt = new TH1F("histoYieldMesonLeftMt","",NBinsPt,BinsMt);
  fHistoYieldMesonLeftMtPerEvent = new TH1F("histoYieldMesonLeftMtPerEvent","",NBinsPt,BinsMt);
  fHistoSignMesonLeft = new TH1F("histoSignMesonLeft","",NBinsPt,BinsPt);
  fHistoSBMesonLeft = new TH1F("histoSBMesonLeft","",NBinsPt,BinsPt);
  fHistoMassMesonLeft = new TH1F("histoMassMesonLeft","",NBinsPt,BinsPt);
  fHistoWidthMesonLeft = new TH1F("histoWidthMesonLeft","",NBinsPt,BinsPt);
  fHistoFWHMMesonLeft = new TH1F("histoFWHMMesonLeft","",NBinsPt,BinsPt);
 

  fHistoYieldMesonLeftNarrow = new TH1F("histoYieldMesonLeftNarrow","",NBinsPt,BinsPt);
  fHistoYieldMesonLeftPerEventNarrow = new TH1F("histoYieldMesonLeftPerEventNarrow","",NBinsPt,BinsPt);
  fHistoYieldMesonLeftMtNarrow = new TH1F("histoYieldMesonLeftMtNarrow","",NBinsPt,BinsMt);
  fHistoYieldMesonLeftMtPerEventNarrow = new TH1F("histoYieldMesonLeftMtPerEventNarrow","",NBinsPt,BinsMt);
  fHistoSignMesonLeftNarrow = new TH1F("histoSignMesonLeftNarrow","",NBinsPt,BinsPt);
  fHistoSBMesonLeftNarrow = new TH1F("histoSBMesonLeftNarrow","",NBinsPt,BinsPt);


  fHistoYieldMesonLeftWide = new TH1F("histoYieldMesonLeftWide","",NBinsPt,BinsPt);
  fHistoYieldMesonLeftPerEventWide = new TH1F("histoYieldMesonLeftPerEventWide","",NBinsPt,BinsPt);
  fHistoYieldMesonLeftMtWide = new TH1F("histoYieldMesonLeftMtWide","",NBinsPt,BinsMt);
  fHistoYieldMesonLeftMtPerEventWide = new TH1F("histoYieldMesonLeftMtPerEventWide","",NBinsPt,BinsMt);
  fHistoSignMesonLeftWide = new TH1F("histoSignMesonLeftWide","",NBinsPt,BinsPt);
  fHistoSBMesonLeftWide = new TH1F("histoSBMesonLeftWide","",NBinsPt,BinsPt);

}

void FillPtHistos()
{
  for(Int_t iPt=iPt0;iPt<NBinsPt+1;iPt++){

    fdeltaPt->SetBinContent(iPt,BinsPt[iPt]-BinsPt[iPt-1]);
    fdeltaPt->SetBinError(iPt,0);
    fdeltaMt->SetBinContent(iPt,BinsMt[iPt]-BinsMt[iPt-1]);

    fHistoMassMeson->SetBinContent(iPt,MesonMass[iPt-1]);
    fHistoMassMeson->SetBinError(iPt,MesonMassError[iPt-1]);
    // fHistoWidthMeson->SetBinContent(iPt);
    fHistoFWHMMeson->SetBinContent(iPt,MesonFWHM[iPt-1]);
    fHistoFWHMMeson->SetBinError(iPt,MesonFWHMError[iPt-1]);
    

    fHistoSignMeson->SetBinContent(iPt,MesonSign[iPt-1]);
    fHistoSignMeson->SetBinError(iPt,MesonSignError[iPt-1]);
    fHistoSBMeson->SetBinContent(iPt,MesonSB[iPt-1]);
    fHistoSBMeson->SetBinError(iPt,MesonSBError[iPt-1]);

    fHistoYieldMeson->SetBinContent(iPt,MesonYieldsCorResidualBckFunc[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMeson->SetBinError(iPt,MesonYieldsCorResidualBckFuncError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldTrueMeson->SetBinContent(iPt,MesonTrueYields[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldTrueMeson->SetBinError(iPt,0.);


    fHistoYieldMesonPerEvent->SetBinContent(iPt,(1./nEvt)*MesonYieldsCorResidualBckFunc[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonPerEvent->SetBinError(iPt,(1./nEvt)*MesonYieldsCorResidualBckFuncError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));

    // Pay attention dividing by delta pt * m_T/p_T Check that it is correct !!!
    fHistoYieldMesonMt->SetBinContent(iPt,(fHistoYieldMesonMt->GetBinCenter(iPt)/fHistoYieldMeson->GetBinCenter(iPt))*
				      MesonYieldsCorResidualBckFunc[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonMt->SetBinError(iPt,(fHistoYieldMesonMt->GetBinCenter(iPt)/fHistoYieldMeson->GetBinCenter(iPt))*
				    MesonYieldsCorResidualBckFuncError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));

    // Narrow integration window
    fHistoSignMesonNarrow->SetBinContent(iPt,MesonSignNarrow[iPt-1]);
    fHistoSignMesonNarrow->SetBinError(iPt,MesonSignNarrowError[iPt-1]);
    fHistoSBMesonNarrow->SetBinContent(iPt,MesonSBNarrow[iPt-1]);
    fHistoSBMesonNarrow->SetBinError(iPt,MesonSBNarrowError[iPt-1]);

    fHistoYieldMesonNarrow->SetBinContent(iPt,MesonYieldsCorResidualBckFuncNarrow[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonNarrow->SetBinError(iPt,MesonYieldsCorResidualBckFuncNarrowError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
 
    fHistoYieldTrueMesonNarrow->SetBinContent(iPt,MesonTrueYieldsNarrow[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldTrueMesonNarrow->SetBinError(iPt,0.);

    fHistoYieldMesonPerEventNarrow->SetBinContent(iPt,(1./nEvt)*MesonYieldsCorResidualBckFuncNarrow[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonPerEventNarrow->SetBinError(iPt,(1./nEvt)*MesonYieldsCorResidualBckFuncNarrowError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));


    // Pay attention dividing by delta pt * m_T/p_T Check that it is correct !!!
    fHistoYieldMesonMtNarrow->SetBinContent(iPt,(fHistoYieldMesonMt->GetBinCenter(iPt)/fHistoYieldMeson->GetBinCenter(iPt))*
					    MesonYieldsCorResidualBckFuncNarrow[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonMtNarrow->SetBinError(iPt,(fHistoYieldMesonMt->GetBinCenter(iPt)/fHistoYieldMeson->GetBinCenter(iPt))*
					  MesonYieldsCorResidualBckFuncNarrowError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));





    // Wide integration window
    fHistoSignMesonWide->SetBinContent(iPt,MesonSignWide[iPt-1]);
    fHistoSignMesonWide->SetBinError(iPt,MesonSignWideError[iPt-1]);
    fHistoSBMesonWide->SetBinContent(iPt,MesonSBWide[iPt-1]);
    fHistoSBMesonWide->SetBinError(iPt,MesonSBWideError[iPt-1]);

    fHistoYieldMesonWide->SetBinContent(iPt,MesonYieldsCorResidualBckFuncWide[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonWide->SetBinError(iPt,MesonYieldsCorResidualBckFuncWideError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));

    fHistoYieldTrueMesonWide->SetBinContent(iPt,MesonTrueYieldsWide[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldTrueMesonWide->SetBinError(iPt,0.);

    fHistoYieldMesonPerEventWide->SetBinContent(iPt,(1./nEvt)*MesonYieldsCorResidualBckFuncWide[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonPerEventWide->SetBinError(iPt,(1./nEvt)*MesonYieldsCorResidualBckFuncWideError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));


    // Pay attention dividing by delta pt * m_T/p_T Check that it is correct !!!
    fHistoYieldMesonMtWide->SetBinContent(iPt,(fHistoYieldMesonMt->GetBinCenter(iPt)/fHistoYieldMeson->GetBinCenter(iPt))*
					  MesonYieldsCorResidualBckFuncWide[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonMtWide->SetBinError(iPt,(fHistoYieldMesonMt->GetBinCenter(iPt)/fHistoYieldMeson->GetBinCenter(iPt))*
					MesonYieldsCorResidualBckFuncWideError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));


    // Histos for integration at the left of the peak

    fHistoMassMesonLeft->SetBinContent(iPt,MesonMassLeft[iPt-1]);
    fHistoMassMesonLeft->SetBinError(iPt,MesonMassLeftError[iPt-1]);
    // fHistoWidthMeson->SetBinContent(iPt);
    // fHistoFWHMMeson->SetBinContent(iPt);
    fHistoFWHMMesonLeft->SetBinContent(iPt,MesonFWHMLeft[iPt-1]);
    fHistoFWHMMesonLeft->SetBinError(iPt,MesonFWHMLeftError[iPt-1]);
    

    fHistoSignMesonLeft->SetBinContent(iPt,MesonSignLeft[iPt-1]);
    fHistoSignMesonLeft->SetBinError(iPt,MesonSignLeftError[iPt-1]);
    fHistoSBMesonLeft->SetBinContent(iPt,MesonSBLeft[iPt-1]);
    fHistoSBMesonLeft->SetBinError(iPt,MesonSBLeftError[iPt-1]);

    fHistoYieldMesonLeft->SetBinContent(iPt,MesonYieldsCorResidualBckFuncLeft[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonLeft->SetBinError(iPt,MesonYieldsCorResidualBckFuncLeftError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonLeftPerEvent->SetBinContent(iPt,(1./nEvt)*MesonYieldsCorResidualBckFuncLeft[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonLeftPerEvent->SetBinError(iPt,(1./nEvt)*MesonYieldsCorResidualBckFuncLeftError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));

    // Pay attention dividing by delta pt * m_T/p_T Check that it is correct !!!
    fHistoYieldMesonLeftMt->SetBinContent(iPt,(fHistoYieldMesonLeftMt->GetBinCenter(iPt)/fHistoYieldMesonLeft->GetBinCenter(iPt))*
					  MesonYieldsCorResidualBckFuncLeft[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonLeftMt->SetBinError(iPt,(fHistoYieldMesonLeftMt->GetBinCenter(iPt)/fHistoYieldMesonLeft->GetBinCenter(iPt))*
					MesonYieldsCorResidualBckFuncLeftError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));

    // Narrow integration window
    fHistoSignMesonLeftNarrow->SetBinContent(iPt,MesonSignLeftNarrow[iPt-1]);
    fHistoSignMesonLeftNarrow->SetBinError(iPt,MesonSignLeftNarrowError[iPt-1]);
    fHistoSBMesonLeftNarrow->SetBinContent(iPt,MesonSBLeftNarrow[iPt-1]);
    fHistoSBMesonLeftNarrow->SetBinError(iPt,MesonSBLeftNarrowError[iPt-1]);

    fHistoYieldMesonLeftNarrow->SetBinContent(iPt,MesonYieldsCorResidualBckFuncLeftNarrow[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonLeftNarrow->SetBinError(iPt,MesonYieldsCorResidualBckFuncLeftNarrowError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonLeftPerEventNarrow->SetBinContent(iPt,(1./nEvt)*MesonYieldsCorResidualBckFuncLeftNarrow[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonLeftPerEventNarrow->SetBinError(iPt,(1./nEvt)*MesonYieldsCorResidualBckFuncLeftNarrowError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));


    // Pay attention dividing by delta pt * m_T/p_T Check that it is correct !!!
    fHistoYieldMesonLeftMtNarrow->SetBinContent(iPt,(fHistoYieldMesonLeftMt->GetBinCenter(iPt)/fHistoYieldMesonLeft->GetBinCenter(iPt))*
						MesonYieldsCorResidualBckFuncLeftNarrow[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonLeftMtNarrow->SetBinError(iPt,(fHistoYieldMesonLeftMt->GetBinCenter(iPt)/fHistoYieldMesonLeft->GetBinCenter(iPt))*
					      MesonYieldsCorResidualBckFuncLeftNarrowError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));


    // Wide integration window
    fHistoSignMesonLeftWide->SetBinContent(iPt,MesonSignLeftWide[iPt-1]);
    fHistoSignMesonLeftWide->SetBinError(iPt,MesonSignLeftWideError[iPt-1]);
    fHistoSBMesonLeftWide->SetBinContent(iPt,MesonSBLeftWide[iPt-1]);
    fHistoSBMesonLeftWide->SetBinError(iPt,MesonSBLeftWideError[iPt-1]);

    fHistoYieldMesonLeftWide->SetBinContent(iPt,MesonYieldsCorResidualBckFuncLeftWide[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonLeftWide->SetBinError(iPt,MesonYieldsCorResidualBckFuncLeftWideError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonLeftPerEventWide->SetBinContent(iPt,(1./nEvt)*MesonYieldsCorResidualBckFuncLeftWide[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonLeftPerEventWide->SetBinError(iPt,(1./nEvt)*MesonYieldsCorResidualBckFuncLeftWideError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));


    // Pay attention dividing by delta pt * m_T/p_T Check that it is correct !!!
    fHistoYieldMesonLeftMtWide->SetBinContent(iPt,(fHistoYieldMesonLeftMt->GetBinCenter(iPt)/fHistoYieldMesonLeft->GetBinCenter(iPt))*
					      MesonYieldsCorResidualBckFuncLeftWide[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));
    fHistoYieldMesonLeftMtWide->SetBinError(iPt,(fHistoYieldMesonLeftMt->GetBinCenter(iPt)/fHistoYieldMesonLeft->GetBinCenter(iPt))*
					    MesonYieldsCorResidualBckFuncLeftWideError[iPt-1]/(BinsPt[iPt]-BinsPt[iPt-1]));

  }
}

void PlotPtHistos(TString suffix, TString prefix2)
{
  TCanvas * c_data_pt = new TCanvas("c_data_pt","",1400,900);  // gives the page size	
  //   c_data_pt->SetTopMargin(0.02);
  //
  //   c_data_pt->SetRightMargin(0.02);
  //   c_data_pt->SetLeftMargin(0.02);
  c_data_pt->SetBottomMargin(2.72);	
  TPad * pad_data_pt = new TPad("pad_data_pt","",0.0,0.0,1.,1.,0);   // gives the size of the histo areas 

  TString MesonChar;
  
  if(prefix.CompareTo("Pi0") == 0){// means we want to plot values for the pi0
	  MesonChar = "#pi^{0}/N_{ev}";
  }
  else if(prefix.CompareTo("Eta") == 0){
	  MesonChar = "#eta/N_{ev}";
  }
  
  pad_data_pt->SetFillColor(0);
  pad_data_pt->GetFrame()->SetFillColor(0);
  pad_data_pt->SetBorderMode(0);
  pad_data_pt->SetLogy(1);
  c_data_pt->SetLogy(1);
  //  pad_data_pt->Divide(1,3);
  pad_data_pt->SetBottomMargin(2.72);
  pad_data_pt->Draw();
  // pad_data_pt->cd(1);
  DrawGammaHisto( fHistoYieldMesonPerEvent, fCutSelection,
		  Form("Pt spectra"),
		  "p_{T} (GeV/c)", MesonChar,
		  "kFALSE",1.2,0.,"kFALSE",0.,fHistoYieldMeson->GetMaximum(),"kTRUE",
		  BinsPt[0], BinsPt[NBinsPt],0.);
  /*
    fHistoYieldMesonPerEventWide->SetMarkerStyle(21);
    fHistoYieldMesonPerEventWide->SetMarkerColor(2);
    fHistoYieldMesonPerEventWide->DrawCopy("same");
  
    fHistoYieldMesonPerEventNarrow->SetMarkerStyle(25);
    fHistoYieldMesonPerEventNarrow->SetMarkerColor(4);
    fHistoYieldMesonPerEventNarrow->DrawCopy("same");
  */
  fHistoYieldMesonLeftPerEvent->SetMarkerStyle(23);
  fHistoYieldMesonLeftPerEvent->SetMarkerColor(3);
  fHistoYieldMesonLeftPerEvent->DrawCopy("same");
  /*
    fHistoYieldMesonLeftPerEventWide->SetMarkerStyle(22);
    fHistoYieldMesonLeftPerEventWide->SetMarkerColor(3);
    fHistoYieldMesonLeftPerEventWide->DrawCopy("same");

    fHistoYieldMesonLeftPerEventNarrow->SetMarkerStyle(26);
    fHistoYieldMesonLeftPerEventNarrow->SetMarkerColor(5);
    fHistoYieldMesonLeftPerEventNarrow->DrawCopy("same");
  */
  /*

    pad_data_pt->cd(2);
    fHistoMassMeson->SetMinimum(MesonPlotRange[0]);
    fHistoMassMeson->SetMaximum(MesonPlotRange[1]);
    DrawGammaHisto( fHistoMassMeson,
    Form("Mass position"),
    "Mass (GeV/c^2)", "Counts",
    "kFALSE",1.2,0.,"kFALSE",MesonMassRange[0],fHistoMassMeson->GetMaximum(),"kTRUE",
    BinsPt[0], BinsPt[NBinsPt],0.);
    fHistoMassMesonLeft->DrawCopy("same");
  */
  /*
    pad_data_pt->cd(3);
    DrawGammaHisto( fHistoSignMeson,
    Form("Significance"),
    "p_{T} (GeV/c)", "Counts",
    "kFALSE",1.2,0.,"kFALSE",0.,fHistoSignMeson->GetMaximum(),"kTRUE",
    BinsPt[0], BinsPt[NBinsPt],0.);
    fHistoSignMesonLeft->DrawCopy("same");
  */
  c_data_pt->Print(Form("%s/%s_%s_pt_%s.%s",suffix.Data(),prefix.Data(),prefix2.Data(),fCutSelection.Data(),suffix.Data()));
}


void PlotInvMassInPtBins(TH1D** Mapping_GG_InvMass_PtBin, TH1D** Mapping_BackNorm_InvMass_PtBin, TString namePlot, TString nameCanvas, TString namePad)
{

  cout<<"Kenneth"<<endl;

  TCanvas * c_data_spectra = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size	
  c_data_spectra->SetTopMargin(0.02);
  c_data_spectra->SetBottomMargin(0.02);
  c_data_spectra->SetRightMargin(0.02);
  c_data_spectra->SetLeftMargin(0.02);
		
  TPad * pad_data_spectra = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas 
  //if(NBinsPt == 4){
  //		pad_data_spectra = new TPad("pad_data_spectra","",0.0,0.1,1.,0.9,0);   
  //}		
  pad_data_spectra->SetFillColor(0);
  pad_data_spectra->GetFrame()->SetFillColor(0);
  pad_data_spectra->SetBorderMode(0);
  pad_data_spectra->Divide(column,row);
  pad_data_spectra->Draw();
 
  cout<<"Kenneth"<<endl;
  cout<<"column: "<<column<<" row: "<<row<<endl;

  Int_t place = 0;
  // for(Int_t iPt=0;iPt<4;iPt++){
  for(Int_t iPt=iPt0;iPt<NBinsPt;iPt++){
    cout<<"Pt: "<<iPt<<" of "<<NBinsPt<<endl;
    Double_t startpt = BinsPt[iPt];	
    Double_t endpt = BinsPt[iPt+1];	
 
    place = place + 1;						//give the right place in the page
    if(place > column*row) {
      Int_t page = place/column/row;
      place = place - (page*column*row);
    }
    if(place == 0) place = column*row;	

    cout<<"Place: "<<place<<endl;
	
    pad_data_spectra->cd(place);
    cout<<"Place: "<<place<<endl;
    pad_data_spectra->cd(place)->SetTopMargin(0.12);
    cout<<"Place: "<<place<<endl;
    pad_data_spectra->cd(place)->SetBottomMargin(0.15);
    cout<<"Place: "<<place<<endl;
    pad_data_spectra->cd(place)->SetRightMargin(0.05);
    cout<<"Place: "<<place<<endl;
    pad_data_spectra->cd(place)->SetLeftMargin(0.15);
    cout<<"Place: "<<place<<endl;

    cout<<"Balle"<<endl;
    DrawGammaHisto( Mapping_GG_InvMass_PtBin[iPt], fCutSelection,
		    Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startpt,endpt),
		    "M_{#gamma#gamma} (GeV/c^{2})", "Counts",
		    "kFALSE",1.2,0.,"kFALSE",0.,Mapping_GG_InvMass_PtBin[iPt]->GetMaximum(),"kTRUE",
		    MesonMassRange[0],MesonMassRange[1],0.);

    cout<<"Balle"<<endl;
    DrawGammaHisto( Mapping_BackNorm_InvMass_PtBin[iPt], fCutSelection,
		    Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startpt,endpt),
		    "M_{#gamma#gamma} (GeV/c^{2})", "Counts",
		    "kFALSE",1.2,-100.,"kFALSE",-1000.,Mapping_BackNorm_InvMass_PtBin[iPt]->GetMaximum(),"kTRUE",
		    MesonMassRange[0],MesonMassRange[1],1.);
  }
  cout<<"Kenneth"<<endl;
  c_data_spectra->Print(namePlot.Data());
  cout<<"Kenneth"<<endl;
}

void PlotWithFitSubtractedInvMassInPtBins(TH1D ** Mapping_Signal_InvMass_PtBin, TF1 ** Signal_InvMassFit_PtBin, TString namePlot, TString nameCanvas, TString namePad)
{
  TCanvas *c_data_fit = new TCanvas(nameCanvas.Data(),"",1400,900);  // gives the page size	
  c_data_fit->SetTopMargin(0.02);
  c_data_fit->SetBottomMargin(0.02);
  c_data_fit->SetRightMargin(0.02);
  c_data_fit->SetLeftMargin(0.02);
  TPad * pad_data_fit = new TPad(namePad.Data(),"",0.0,0.0,1.,1.,0);   // gives the size of the histo areas 

  //if(NBinsPt == 4){
  //		pad_data_fit = new TPad("pad_data_fit","",0.0,0.1,1.,0.9,0);   
  //}
  pad_data_fit->SetFillColor(0);
  pad_data_fit->GetFrame()->SetFillColor(0);
  pad_data_fit->SetBorderMode(0);
  pad_data_fit->Divide(column,row);
  pad_data_fit->SetLeftMargin(0.3);
  pad_data_fit->Draw();

  Int_t place = 0;		
  //  for(Int_t iPt = 0; iPt < 4; iPt++){
  for(Int_t iPt = iPt0; iPt < NBinsPt; iPt++){
    //for(Int_t iPt = 2; iPt < 3; iPt++){
    Double_t startpt = BinsPt[iPt];	
    Double_t endpt = BinsPt[iPt+1];	
 			
    place = place + 1;						//give the right place in the page
    if(place > column*row) {
      Int_t page = place/column/row;
      place = place - (page*column*row);
    }
    if(place == 0) place = column*row;		
    pad_data_fit->cd(place);
    pad_data_fit->cd(place)->SetTopMargin(0.12);
    pad_data_fit->cd(place)->SetBottomMargin(0.15);
    pad_data_fit->cd(place)->SetRightMargin(0.05);
    pad_data_fit->cd(place)->SetLeftMargin(0.15);
    Mapping_Signal_InvMass_PtBin[iPt]->SetAxisRange(MesonMassRange[0],MesonMassRange[1]);
    cout<<"Maximum::"<<Mapping_Signal_InvMass_PtBin[iPt]->GetMaximum()<<endl;
    DrawGammaHisto( Mapping_Signal_InvMass_PtBin[iPt], fCutSelection,
		    Form("%3.2f GeV/c < p_{t} < %3.2f GeV/c",startpt,endpt),
		    "M_{#gamma#gamma} (GeV/c^{2})", "Counts",
		    "kTRUE",2.2,0.,"kTRUE",-500,2*Mapping_Signal_InvMass_PtBin[iPt]->GetMaximum(),"kTRUE",
		    MesonMassRange[0],MesonMassRange[1],0.);
    Signal_InvMassFit_PtBin[iPt]->DrawCopy("same");
    if(fIsMC) Mapping_TrueMeson_InvMass_PtBin[iPt]->DrawCopy("same");
  }
  c_data_fit->Print(namePlot.Data());
}
void FitAllSubtractedInvMassInPtBins(TH1D** Mapping_Signal_InvMass_PtBin)
{
  for(Int_t iPt = iPt0; iPt < NBinsPt; iPt++){
    //for(Int_t iPt = 2; iPt < 3; iPt++){
    Double_t startpt = BinsPt[iPt];	
    Double_t endpt = BinsPt[iPt+1];	
    FitSubtractedInvMassInPtBins(Mapping_Signal_InvMass_PtBin[iPt],MesonIntRange);
    Signal_InvMassFit_PtBin[iPt]=fReco;   
  }

}

void FitSubtractedInvMassInPtBins(TH1D* Mapping_Signal_InvMass_PtBinSingle,Double_t * MesonIntRange)
{

  cout<<"Start Fitting spectra"<<endl;
  Mapping_Signal_InvMass_PtBinSingle->GetXaxis()->SetRangeUser(MesonMassRange[0],MesonMassRange[1]);
  Double_t MesonAmplitude =Mapping_Signal_InvMass_PtBinSingle->GetMaximum();
  Double_t MesonAmplitudeMin = MesonAmplitude*95./100.;
  Double_t MesonAmplitudeMax = MesonAmplitude*105./100.;


  fReco = new TF1("GaussExpLinear","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1-exp(-0.5*((x-[1])/[2])^2)))+[4]+[5]*x)+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2)+[4]+[5]*x)",MesonFitRange[0],MesonFitRange[1]); 


  fGausExp = new TF1("fGaussExp","(x<[1])*([0]*(exp(-0.5*((x-[1])/[2])^2)+exp((x-[1])/[3])*(1-exp(-0.5*((x-[1])/[2])^2))))+(x>=[1])*([0]*exp(-0.5*((x-[1])/[2])^2))",MesonFitRange[0],MesonFitRange[1]); 

  fLinearBck = new TF1("Linear","[0]+[1]*x",MesonFitRange[0],MesonFitRange[1]); 


  fReco->SetParameter(0,MesonAmplitude);
  fReco->SetParameter(1,MesonMassExp);
  fReco->SetParameter(2,MesonWidthExp);
  fReco->SetParameter(3,MesonLambdaTail);
	
  fReco->SetParLimits(0,MesonAmplitudeMin,MesonAmplitudeMax);
  fReco->SetParLimits(1,MesonMassRange[0],MesonMassRange[1]);
  fReco->SetParLimits(2,MesonWidthRange[0],MesonWidthRange[1]);
  fReco->SetParLimits(3,MesonLambdaTailRange[0],MesonLambdaTailRange[1]);

  Mapping_Signal_InvMass_PtBinSingle->Fit(fReco,"RME0");
  Mapping_Signal_InvMass_PtBinSingle->Fit(fReco,"RME0");

  fReco->SetLineColor(3);
  fReco->SetLineWidth(2);
  fReco->SetLineStyle(1);		

  fGausExp->SetParameter(0,fReco->GetParameter(0));
  fGausExp->SetParameter(1,fReco->GetParameter(1));
  fGausExp->SetParameter(2,fReco->GetParameter(2));
  fGausExp->SetParameter(3,fReco->GetParameter(3));
  
  fGausExp->SetParError(0,fReco->GetParError(0));
  fGausExp->SetParError(1,fReco->GetParError(1));
  fGausExp->SetParError(2,fReco->GetParError(2));
  fGausExp->SetParError(3,fReco->GetParError(3));

  fLinearBck->SetParameter(0,fReco->GetParameter(4));
  fLinearBck->SetParameter(1,fReco->GetParameter(5));

  fLinearBck->SetParError(0,fReco->GetParError(4));
  fLinearBck->SetParError(1,fReco->GetParError(5));



  TVirtualFitter * fitter = TVirtualFitter::GetFitter();
  Int_t nfreepar = fReco->GetNumberFreeParameters();
  double * covMatrix = fitter->GetCovarianceMatrix();

  Float_t IntLinearBck = fLinearBck->GetParameter(0)*(MesonIntRange[1]-MesonIntRange[0])+
    0.5*fLinearBck->GetParameter(1)*(MesonIntRange[1]*MesonIntRange[1]-MesonIntRange[0]*MesonIntRange[0]);
												          

  Float_t errorLinearBck = pow(
			       (pow( (MesonIntRange[1]-MesonIntRange[0])*fReco->GetParError(4),2)+
				pow(0.5*(MesonIntRange[1]*MesonIntRange[1]-MesonIntRange[0]*MesonIntRange[0])*fReco->GetParError(5),2)+
				2*covMatrix[nfreepar*nfreepar-2]*(MesonIntRange[1]-MesonIntRange[0])*
				0.5*(MesonIntRange[1]*MesonIntRange[1]-MesonIntRange[0]*MesonIntRange[0]))
			       ,0.5);


  cout<<"Remaining bck +-error::"<< IntLinearBck/Mapping_Signal_InvMass_PtBinSingle->GetBinWidth(10) <<" +/- "<< errorLinearBck/Mapping_Signal_InvMass_PtBinSingle->GetBinWidth(10)<<endl;
 
  fIntLinearBck = IntLinearBck/Mapping_Signal_InvMass_PtBinSingle->GetBinWidth(10);
  fIntLinearBckError = errorLinearBck/Mapping_Signal_InvMass_PtBinSingle->GetBinWidth(10);

  fReco->DrawCopy("same");

  
}


void IntegrateHistoInvMass(TH1D * Mapping_Signal_InvMass_PtBinSingle, Double_t * MesonIntRange)
{
  Int_t ilowmassMeson_bin = Mapping_Signal_InvMass_PtBinSingle->GetXaxis()->FindBin(MesonIntRange[0]);  
  Int_t ihighmassMeson_bin = Mapping_Signal_InvMass_PtBinSingle->GetXaxis()->FindBin(MesonIntRange[1]);	
  fYields = Mapping_Signal_InvMass_PtBinSingle->IntegralAndError(ilowmassMeson_bin,ihighmassMeson_bin,fYieldsError);


}


void IntegrateFitFunc(TF1 * fFunc, TH1D *  Mapping_Signal_InvMass_PtBinSingle,Double_t * MesonIntRange)
{

  fYieldsFunc = fFunc->Integral(MesonIntRange[0],MesonIntRange[1])/Mapping_Signal_InvMass_PtBinSingle->GetBinWidth(10);

}



void FillHistosArrayMC(TH2F* fMC_Meson_Pt_Eta_withinAcceptance, TH1D * fMC_Meson_Pt, TH1F * fdeltaPt)
{

  Int_t startBin = fMC_Meson_Pt_Eta_withinAcceptance->GetYaxis()->FindBin(-fYMaxMeson+0.001);
  Int_t endBin = fMC_Meson_Pt_Eta_withinAcceptance->GetYaxis()->FindBin(fYMaxMeson-0.001);

  cout<< "Y-bins::"<< startBin<< " " << endBin<<" "<<   endl;
  cout<< "Y- bin values::"<<  fMC_Meson_Pt_Eta_withinAcceptance->GetYaxis()->GetBinCenter(startBin)<< " "
      << fMC_Meson_Pt_Eta_withinAcceptance->GetYaxis()->GetBinCenter(endBin)<< endl;

  Char_t histoname[100] = "fMC_Meson_Pt_Eta_withinAcceptance";
  fMC_Meson_Pt_Eta_withinAcceptance->ProjectionX(histoname,startBin,endBin);

  fMC_Meson_Pt_WithinAcceptance = (TH1D*)gDirectory->Get(histoname);
  fMC_MesonWithinAccepPt = (TH1D*)fMC_Meson_Pt_WithinAcceptance->Rebin(NBinsPt,"",BinsPt);; // Proper bins in Pt
  fMC_MesonWithinAccepPt->Divide(fdeltaPt);
  fMC_MesonPt = (TH1D*)fMC_Meson_Pt->Rebin(NBinsPt,"",BinsPt);; // Proper bins in Pt
  fMC_MesonPt->Divide(fdeltaPt);

}


void CalculateMesonAcceptance()
{
  fMC_MesonAccepPt = new TH1D("fMCMesonAccepPt","",NBinsPt,BinsPt);
  fMC_MesonAccepPt->Sumw2();
  fMC_MesonAccepPt->Divide(fMC_MesonWithinAccepPt,fMC_MesonPt,1.,1.,"B");
  fMC_MesonAccepPt->DrawCopy();
}

void CalculateMesonEfficiency(TH1F* fMC_MesonYieldsPt,TString nameEfi )
{
  fMCMesonEffiPt = new TH1D(nameEfi.Data(),"",NBinsPt,BinsPt);
  fMCMesonEffiPt->Sumw2();
  fMCMesonEffiPt->Divide(fMC_MesonYieldsPt,fMC_MesonWithinAccepPt,1.,1.,"B");
  //fMCMesonEffiPt->DrawCopy();

}

void SaveHistos(Int_t fIsMC, TString fCutID, TString prefix2)
 {
  const char* Outputname = Form("%s_%s_AnalysisResultsWithoutCorrection_%s.root",prefix.Data(),prefix2.Data(),fCutID.Data());
  Output1 = new TFile(Outputname,"RECREATE");		

    fHistoYieldMeson->Write();
      cout << "written fHistoYieldMeson" << endl;
    fHistoYieldMesonPerEvent->Write();
    cout << "written fHistoYieldMesonPerEvent" << endl;
    fHistoYieldMesonMt->Write();
    cout << "written fHistoYieldMesonMt" << endl;
    fHistoYieldMesonMtPerEvent->Write();
    cout << "written fHistoYieldMesonMtPerEvent" << endl;
    fHistoSignMeson->Write();
    cout << "written fHistoSignMeson" << endl;
    fHistoSBMeson->Write();
    cout << "written fHistoSBMeson" << endl;

    fHistoYieldMesonNarrow->Write();
    cout << "written fHistoYieldMesonNarrow" << endl;
    fHistoYieldMesonPerEventNarrow->Write();
    cout << "written fHistoYieldMesonPerEventNarrow" << endl;
    fHistoYieldMesonMtNarrow->Write();
    cout << "written fHistoYieldMesonMtNarrow" << endl;
    fHistoYieldMesonMtPerEventNarrow->Write();
    cout << "written fHistoYieldMesonMtPerEventNarrow" << endl;
    fHistoSignMesonNarrow->Write();
	cout << "written fHistoSignMesonNarrow" << endl;
    fHistoSBMesonNarrow->Write();
    cout << "written fHistoSBMesonNarrow" << endl;

    fHistoYieldMesonWide->Write();
    cout << "written fHistoYieldMesonWide" << endl;
    fHistoYieldMesonPerEventWide->Write();
    cout << "written fHistoYieldMesonPerEventWide" << endl;
    fHistoYieldMesonMtWide->Write();
    cout << "written fHistoYieldMesonMtWide" << endl;
    
    fHistoYieldMesonMtPerEventWide->Write();
    cout << "written fHistoYieldMesonMtPerEventWide" << endl;
    fHistoSignMesonWide->Write();
    cout << "written fHistoSignMesonWide" << endl;
    fHistoSBMesonWide->Write();
    cout << "written fHistoSBMesonWide" << endl;

    fHistoMassMeson->Write();
    cout << "written fHistoMassMeson" << endl;
    fHistoWidthMeson->Write();
    cout << "written fHistoWidthMeson" << endl;
    fHistoFWHMMeson->Write();
    cout << "written fHistoFWHMMeson" << endl;
    fdeltaPt->Write();
    cout << "written fdeltaPt" << endl;
    fdeltaMt->Write();
    cout << "written fdeltaMt" << endl;
    
    fHistoYieldMesonLeft->Write();
    cout << "written fHistoYieldMesonLeft" << endl;
    fHistoYieldMesonLeftPerEvent->Write();
    cout << "written fHistoYieldMesonLeftPerEvent" << endl;
    fHistoYieldMesonLeftMt->Write();
    cout << "written fHistoYieldMesonLeftMt" << endl;
    fHistoYieldMesonLeftMtPerEvent->Write();
    cout << "written fHistoYieldMesonLeftMtPerEvent" << endl;
    fHistoSignMesonLeft->Write();
    cout << "written fHistoSignMesonLeft" << endl;
    fHistoSBMesonLeft->Write();
    cout << "written fHistoSBMesonLeft" << endl;
    
    fHistoYieldMesonLeftNarrow->Write();
    cout << "written fHistoYieldMesonLeftNarrow" << endl;
    fHistoYieldMesonLeftPerEventNarrow->Write();
    cout << "written fHistoYieldMesonLeftPerEventNarrow" << endl;
    fHistoYieldMesonLeftMtNarrow->Write();
    cout << "written fHistoYieldMesonLeftMtNarrow" << endl;
    fHistoYieldMesonLeftMtPerEventNarrow->Write();
    cout << "written fHistoYieldMesonLeftMtPerEventNarrow" << endl;
    fHistoSignMesonLeftNarrow->Write();
    cout << "written fHistoSignMesonLeftNarrow" << endl;
    fHistoSBMesonLeftNarrow->Write();
    cout << "written fHistoSBMesonLeftNarrow" << endl;

    fHistoYieldMesonLeftWide->Write();
    cout << "written fHistoYieldMesonLeftWide" << endl;
    fHistoYieldMesonLeftPerEventWide->Write();
    cout << "written fHistoYieldMesonLeftPerEventWide" << endl;
    fHistoYieldMesonLeftMtWide->Write();
    cout << "written fHistoYieldMesonLeftMtWide" << endl;
    fHistoYieldMesonLeftMtPerEventWide->Write();
    cout << "written fHistoYieldMesonLeftMtPerEventWide" << endl;
    fHistoSignMesonLeftWide->Write();
    cout << "written fHistoSignMesonLeftWide" << endl;
    fHistoSBMesonLeftWide->Write();
    cout << "written fHistoSBMesonLeftWide" << endl;

    fHistoMassMesonLeft->Write();
    cout << "written fHistoMassMesonLeft" << endl;
    fHistoWidthMesonLeft->Write();
    cout << "written fHistoWidthMesonLeft" << endl;
    fHistoFWHMMesonLeft->Write();
    cout << "written fHistoFWHMMesonLeft" << endl;
    
    fNumberOfGoodESDTracksVtx->Write();
    
  Output1->Write();
  Output1->Close();
}

void SaveHistosPtBins(TH1D * Mapping_Signal_InvMass_PtBin,TString fCutID, TString prefix2)
{
  const char* Outputname = Form("%s_%s_AnalysisResultsWithoutCorrection_%s.root",prefix.Data(),prefix2.Data(),fCutID.Data());
  Output1 = new TFile(Outputname,"UPDATE");		
  
	 Mapping_Signal_InvMass_PtBin->Write();  
	 
    Output1->Write();
    Output1->Close();  
}  

void SaveCorrectionHistos(Int_t fIsMC, TString fCutID, TString prefix2)
{
  const char* Outputname = Form("%s_%s_AnalysisResultsCorrectionHistos_%s.root",prefix.Data(),prefix2.Data(),fCutID.Data());
  Output2 = new TFile(Outputname,"RECREATE");		
  
  fMC_Meson_Pt_Eta_withinAcceptance->Write();
  //fMC_Meson_Pt->Write();
  
  fMC_Meson_Pt_WithinAcceptance->Write();
  cout << "written fMC_Meson_Pt_WithinAcceptance" << endl;
  fMC_MesonWithinAccepPt->Write(); // Proper bins in Pt
  cout << "written fMC_MesonWithinAccepPt" << endl;
  fMC_MesonPt->Write(); // Proper bins in Pt
  cout << "written fMC_MesonPt" << endl;

  fHistoYieldTrueMeson->Write();
  cout << "written fHistoYieldTrueMeson" << endl;
  fHistoYieldTrueMesonWide->Write();
  cout << "written fHistoYieldTrueMesonWide" << endl;
  fHistoYieldTrueMesonNarrow->Write();
  cout << "written fHistoYieldTrueMesonNarrow" << endl;
  fMC_MesonAccepPt->Write();
  cout << "written fMC_MesonAccepPt" << endl;
  fTrueMesonInvMassVSPt->Write();
  cout << "written fTrueMesonInvMassVSPt" << endl;
  fMC_MesonEffiPt->Write();
  cout << "written fMC_MesonEffiPt" << endl;
  fMC_MesonNarrowEffiPt->Write();
  cout << "written fMC_MesonNarrowEffiPt" << endl;
  fMC_MesonWideEffiPt->Write();
  cout << "written fMC_MesonWideEffiPt" << endl;
  fMC_MesonLeftEffiPt->Write();
  cout << "written fMC_MesonLeftEffiPt" << endl;
  fMC_MesonLeftNarrowEffiPt->Write();
  cout << "written fMC_MesonLeftNarrowEffiPt" << endl;
  fMC_MesonLeftWideEffiPt->Write();
  cout << "written fMC_MesonLeftWideEffiPt" << endl;
  fMC_TrueMesonEffiPt->Write();
  cout << "written fMC_TrueMesonEffiPt" << endl;
  fMC_TrueMesonNarrowEffiPt->Write();
  cout << "written fMC_TrueMesonNarrowEffiPt" << endl;  
  fMC_TrueMesonWideEffiPt->Write();
  cout << "written fMC_TrueMesonWideEffiPt" << endl;    
  
  Output2->Write();
  Output2->Close();
}

void Initialize(TString setPi0)
{
	if (setPi0.CompareTo("Pi0") == 0){ // means we should make plots for the pi0
    iPt0=1;
    column = 6; 
    row = 4;	
    
    NBinsPt =24;	
    BinsPt= new Double_t[NBinsPt+1];
    BinsPt[0] = 0;
    BinsPt[1] = 0.3;
    BinsPt[2] = 0.4;
    BinsPt[3] = 0.5;
    BinsPt[4] = 0.6;
    BinsPt[5] = 0.8;
    BinsPt[6] = 1.;
    BinsPt[7] = 1.2;
    BinsPt[8] = 1.4;
    BinsPt[9] = 1.6;
    BinsPt[10] = 1.8;
    BinsPt[11] = 2.;
    BinsPt[12] = 2.2;
    BinsPt[13] = 2.4;
    BinsPt[14] = 2.6;
    BinsPt[15] = 2.8;
    BinsPt[16] = 3.;
    BinsPt[17] = 3.4;
    BinsPt[18] = 3.8;
    BinsPt[19] = 4.2;
    BinsPt[20] = 5.;
    BinsPt[21] = 6.;
    BinsPt[22] = 7.;
    BinsPt[23] = 8.;
    BinsPt[24] = 10.;
/*    BinsPt[25] = 12.;
    BinsPt[26] = 15.;*/
    
    BinsMt= new Double_t[NBinsPt+1];
    
    fNRebin = new Int_t[NBinsPt];
    fNRebin[0] = 4;
    fNRebin[1] = 4;
    fNRebin[2] = 2;
    fNRebin[3] = 2;
    fNRebin[4] = 2;
    fNRebin[5] = 2;
    fNRebin[6] = 2;
    fNRebin[7] = 2;
    fNRebin[8] = 2;
    fNRebin[9] = 2;
    fNRebin[10] = 2;
    fNRebin[11] = 2;
    fNRebin[12] = 2;
    fNRebin[13] = 2;
    fNRebin[14] = 2;
    fNRebin[15] = 2;
    fNRebin[16] = 2;
    fNRebin[17] = 4;
    fNRebin[18] = 4;
    fNRebin[19] = 4;
    fNRebin[20] = 4;
    fNRebin[21] = 5;
    fNRebin[22] = 5;
    fNRebin[23] = 5;
//     fNRebin[24] = 10;
//     fNRebin[25] = 10;

    
    BGFitRange = new Double_t[2]; BGFitRange[0]=0.17; BGFitRange[1]=0.3; //eta 0.9

    BGFitRangeLeft = new Double_t[2]; BGFitRangeLeft[0]=0.05; BGFitRangeLeft[1]=0.08;  // eta 09 
    
    MesonPlotRange = new Double_t[2]; MesonPlotRange[0]=0.13; MesonPlotRange[1]=0.138;

    MesonIntRange = new Double_t[2]; MesonIntRange[0] = 0.1; MesonIntRange[1]=0.145;
    
    MesonIntRangeWide = new Double_t[2]; MesonIntRangeWide[0]=0.08; MesonIntRangeWide[1]=0.160;
    
    MesonIntRangeNarrow = new Double_t[2]; MesonIntRangeNarrow[0]=0.12; MesonIntRangeNarrow[1]=0.14;
    
    MesonMassRange = new Double_t[2]; MesonMassRange[0]=0.; MesonMassRange[1]=0.3;

    MesonFitRange = new Double_t[2]; MesonFitRange[0]=0.07; MesonFitRange[1]=0.200;
    
    MesonId=111;
    
    MesonWidthExp = 0.003;
    
    MesonLambdaTail = 0.007;
    
    MesonWidthRange = new Double_t[2]; MesonWidthRange[0]=0.001; MesonWidthRange[1]=0.007;
    
    MesonLambdaTailRange = new Double_t[2]; MesonLambdaTailRange[0]=0.001; MesonLambdaTailRange[1]=0.016;
  }
  else if (setPi0.CompareTo("Eta") == 0){ // means we should make plots for the eta
    iPt0=1;
    column = 4; 
    row = 3;	
    
    NBinsPt =10;	
    BinsPt= new Double_t[NBinsPt+1];
    BinsPt[0] = 0.;
    BinsPt[1] = 0.4;
    BinsPt[2] = 0.6;
    BinsPt[3] = 1.2;
    BinsPt[4] = 1.8;
    BinsPt[5] = 2.4;
    BinsPt[6] = 3.;
    BinsPt[7] = 4.;
    BinsPt[8] = 6.;
    BinsPt[9] = 8.;
    BinsPt[10] = 10.;

    BinsMt= new Double_t[NBinsPt+1]; // remember to include....
    
    fNRebin = new Int_t[NBinsPt];
    fNRebin[0] = 4;
    fNRebin[1] = 4;
    fNRebin[2] = 4;
    fNRebin[3] = 4;
    fNRebin[4] = 4;
    fNRebin[5] = 4;
    fNRebin[6] = 4;
    fNRebin[7] = 8;
    fNRebin[8] = 10;
    fNRebin[9] = 10;

    BGFitRange = new Double_t[2]; BGFitRange[0]=0.58; BGFitRange[1]=0.62; //eta 0.9

    BGFitRangeLeft = new Double_t[2]; BGFitRangeLeft[0]=0.45; BGFitRangeLeft[1]=0.5;  // eta 09 
    
    MesonPlotRange = new Double_t[2]; MesonPlotRange[0]=0.53; MesonPlotRange[1]=0.560;

    MesonIntRange = new Double_t[2]; MesonIntRange[0] = 0.51; MesonIntRange[1]=0.57;
    
    MesonIntRangeWide = new Double_t[2]; MesonIntRangeWide[0]=0.5; MesonIntRangeWide[1]=0.58;
    
    MesonIntRangeNarrow = new Double_t[2]; MesonIntRangeNarrow[0]=0.52; MesonIntRangeNarrow[1]=0.56;
    
    MesonMassRange = new Double_t[2]; MesonMassRange[0]=0.45; MesonMassRange[1]=0.65;

    MesonFitRange = new Double_t[2]; MesonFitRange[0]=0.45; MesonFitRange[1]=0.65;
    
    MesonId=221;
    
    MesonWidthExp = 0.005;
    
    MesonLambdaTail = 0.007;
    
    MesonWidthRange = new Double_t[2]; MesonWidthRange[0]=0.001; MesonWidthRange[1]=0.007;
    
    MesonLambdaTailRange = new Double_t[2]; MesonLambdaTailRange[0]=0.0005; MesonLambdaTailRange[1]=0.026;
  }
  else {
	  iPt0=1;
	  column = 4; 
	  row = 3;	
	  
	  NBinsPt =10;	
	  BinsPt= new Double_t[NBinsPt+1];
	  BinsPt[0] = 0.;
	  BinsPt[1] = 0.4;
	  BinsPt[2] = 0.6;
	  BinsPt[3] = 1.2;
	  BinsPt[4] = 1.8;
	  BinsPt[5] = 2.4;
	  BinsPt[6] = 3.;
	  BinsPt[7] = 4.;
	  BinsPt[8] = 6.;
	  BinsPt[9] = 8.;
	  BinsPt[10] = 10.;
	  
	  BinsMt= new Double_t[NBinsPt+1]; // remember to include....
	  
	  fNRebin = new Int_t[NBinsPt];
	  fNRebin[0] = 2;
	  fNRebin[1] = 2;
	  fNRebin[2] = 2;
	  fNRebin[3] = 2;
	  fNRebin[4] = 2;
	  fNRebin[5] = 2;
	  fNRebin[6] = 2;
	  fNRebin[7] = 2;
	  fNRebin[8] = 2;
	  fNRebin[9] = 2;
	  
	  BGFitRange = new Double_t[2]; BGFitRange[0]=0.17; BGFitRange[1]=0.3; //eta 0.9
	  
	  BGFitRangeLeft = new Double_t[2]; BGFitRangeLeft[0]=0.05; BGFitRangeLeft[1]=0.08;  // eta 09 
	  
	  MesonPlotRange = new Double_t[2]; MesonPlotRange[0]=0.13; MesonPlotRange[1]=0.138;
	  
	  MesonIntRange = new Double_t[2]; MesonIntRange[0] = 0.1; MesonIntRange[1]=0.145;
	  
	  MesonIntRangeWide = new Double_t[2]; MesonIntRangeWide[0]=0.08; MesonIntRangeWide[1]=0.160;
	  
	  MesonIntRangeNarrow = new Double_t[2]; MesonIntRangeNarrow[0]=0.12; MesonIntRangeNarrow[1]=0.14;
	  
	  MesonMassRange = new Double_t[2]; MesonMassRange[0]=0.; MesonMassRange[1]=0.3;
	  
	  MesonFitRange = new Double_t[2]; MesonFitRange[0]=0.07; MesonFitRange[1]=0.200;
	  
	  MesonId=111;
	  
	  MesonWidthExp = 0.003;
	  
	  MesonLambdaTail = 0.007;
	  
	  MesonWidthRange = new Double_t[2]; MesonWidthRange[0]=0.001; MesonWidthRange[1]=0.007;
	  
	  MesonLambdaTailRange = new Double_t[2]; MesonLambdaTailRange[0]=0.001; MesonLambdaTailRange[1]=0.016;	  

  }	  
	  
	  
  GGYields = new Double_t[NBinsPt];
  BckYields = new Double_t[NBinsPt];
  MesonYields = new Double_t[NBinsPt];
  MesonTrueYields = new Double_t[NBinsPt];
  MesonYieldsFunc = new Double_t[NBinsPt];
  MesonYieldsResidualBckFunc = new Double_t[NBinsPt];
  MesonYieldsCorResidualBckFunc = new Double_t[NBinsPt];
  MesonYieldsPerEvent = new Double_t[NBinsPt];
  MesonMass = new Double_t[NBinsPt];
  MesonWidth = new Double_t[NBinsPt];
  MesonSB = new Double_t[NBinsPt];
  MesonSign = new Double_t[NBinsPt];
  MesonFWHM = new Double_t[NBinsPt];
  
  // Normalization at the left of the peak
  GGYieldsLeft = new Double_t[NBinsPt];
  BckYieldsLeft = new Double_t[NBinsPt];
  MesonYieldsLeft = new Double_t[NBinsPt];
  MesonYieldsFuncLeft = new Double_t[NBinsPt];
  MesonYieldsResidualBckFuncLeft = new Double_t[NBinsPt];
  MesonYieldsCorResidualBckFuncLeft = new Double_t[NBinsPt];
  MesonYieldsLeftPerEvent = new Double_t[NBinsPt];
  MesonMassLeft = new Double_t[NBinsPt];
  MesonWidthLeft = new Double_t[NBinsPt];
  MesonSBLeft = new Double_t[NBinsPt];
  MesonSignLeft = new Double_t[NBinsPt];
  MesonFWHMLeft = new Double_t[NBinsPt];
  
  // Narrow Integration Window
  GGYieldsNarrow = new Double_t[NBinsPt];
  BckYieldsNarrow = new Double_t[NBinsPt];
  MesonYieldsNarrow = new Double_t[NBinsPt];
  MesonTrueYieldsNarrow = new Double_t[NBinsPt];
  MesonYieldsFuncNarrow = new Double_t[NBinsPt];
  MesonYieldsResidualBckFuncNarrow = new Double_t[NBinsPt];
  MesonYieldsCorResidualBckFuncNarrow = new Double_t[NBinsPt];
  MesonYieldsPerEventNarrow = new Double_t[NBinsPt];
  MesonSBNarrow = new Double_t[NBinsPt];
  MesonSignNarrow = new Double_t[NBinsPt];

  GGYieldsLeftNarrow = new Double_t[NBinsPt];
  BckYieldsLeftNarrow = new Double_t[NBinsPt];
  MesonYieldsLeftNarrow = new Double_t[NBinsPt];
  MesonYieldsFuncLeftNarrow = new Double_t[NBinsPt];
  MesonYieldsResidualBckFuncLeftNarrow = new Double_t[NBinsPt];
  MesonYieldsCorResidualBckFuncLeftNarrow = new Double_t[NBinsPt];
  MesonYieldsLeftPerEventNarrow = new Double_t[NBinsPt];
  MesonSBLeftNarrow = new Double_t[NBinsPt];
  MesonSignLeftNarrow = new Double_t[NBinsPt];

  // Wide Integration Window
  GGYieldsWide = new Double_t[NBinsPt];
  BckYieldsWide = new Double_t[NBinsPt];
  MesonYieldsWide = new Double_t[NBinsPt];
  MesonTrueYieldsWide = new Double_t[NBinsPt];
  MesonYieldsFuncWide = new Double_t[NBinsPt];
  MesonYieldsResidualBckFuncWide = new Double_t[NBinsPt];
  MesonYieldsCorResidualBckFuncWide = new Double_t[NBinsPt];
  MesonYieldsPerEventWide = new Double_t[NBinsPt];
  MesonSBWide = new Double_t[NBinsPt];
  MesonSignWide = new Double_t[NBinsPt];

  GGYieldsLeftWide = new Double_t[NBinsPt];
  BckYieldsLeftWide = new Double_t[NBinsPt];
  MesonYieldsLeftWide = new Double_t[NBinsPt];
  MesonYieldsFuncLeftWide = new Double_t[NBinsPt];
  MesonYieldsResidualBckFuncLeftWide = new Double_t[NBinsPt];
  MesonYieldsCorResidualBckFuncLeftWide = new Double_t[NBinsPt];
  MesonYieldsLeftPerEventWide = new Double_t[NBinsPt];
  MesonSBLeftWide = new Double_t[NBinsPt];
  MesonSignLeftWide = new Double_t[NBinsPt];

  GGYieldsError = new Double_t[NBinsPt];
  BckYieldsError = new Double_t[NBinsPt];
  MesonYieldsError = new Double_t[NBinsPt];
  MesonYieldsFuncError = new Double_t[NBinsPt];
  MesonYieldsResidualBckFuncError = new Double_t[NBinsPt];
  MesonYieldsCorResidualBckFuncError = new Double_t[NBinsPt];
  MesonYieldsPerEventError = new Double_t[NBinsPt];
  MesonMassError = new Double_t[NBinsPt];
  MesonWidthError = new Double_t[NBinsPt];
  MesonSBError = new Double_t[NBinsPt];
  MesonSignError = new Double_t[NBinsPt];
  MesonFWHMError = new Double_t[NBinsPt];

  GGYieldsLeftError = new Double_t[NBinsPt];
  BckYieldsLeftError = new Double_t[NBinsPt];
  MesonYieldsLeftError = new Double_t[NBinsPt];
  MesonYieldsFuncLeftError = new Double_t[NBinsPt];
  MesonYieldsResidualBckFuncLeftError = new Double_t[NBinsPt];
  MesonYieldsCorResidualBckFuncLeftError = new Double_t[NBinsPt];
  MesonYieldsLeftPerEventError = new Double_t[NBinsPt];
  MesonMassLeftError = new Double_t[NBinsPt];
  MesonWidthLeftError = new Double_t[NBinsPt];
  MesonSBLeftError = new Double_t[NBinsPt];
  MesonSignLeftError = new Double_t[NBinsPt];
  MesonFWHMLeftError = new Double_t[NBinsPt];
  
  // Narrow integration Window
  GGYieldsNarrowError = new Double_t[NBinsPt];
  BckYieldsNarrowError = new Double_t[NBinsPt];
  MesonYieldsNarrowError = new Double_t[NBinsPt];
  MesonYieldsFuncNarrowError = new Double_t[NBinsPt];
  MesonYieldsResidualBckFuncNarrowError = new Double_t[NBinsPt];
  MesonYieldsCorResidualBckFuncNarrowError = new Double_t[NBinsPt];
  MesonYieldsPerEventNarrowError = new Double_t[NBinsPt];
  MesonSBNarrowError = new Double_t[NBinsPt];
  MesonSignNarrowError = new Double_t[NBinsPt];

  GGYieldsLeftNarrowError = new Double_t[NBinsPt];
  BckYieldsLeftNarrowError = new Double_t[NBinsPt];
  MesonYieldsLeftNarrowError = new Double_t[NBinsPt];
  MesonYieldsFuncLeftNarrowError = new Double_t[NBinsPt];
  MesonYieldsResidualBckFuncLeftNarrowError = new Double_t[NBinsPt];
  MesonYieldsCorResidualBckFuncLeftNarrowError = new Double_t[NBinsPt];
  MesonYieldsLeftPerEventNarrowError = new Double_t[NBinsPt];
  MesonSBLeftNarrowError = new Double_t[NBinsPt];
  MesonSignLeftNarrowError = new Double_t[NBinsPt];

  // Wide integration Window
  GGYieldsWideError = new Double_t[NBinsPt];
  BckYieldsWideError = new Double_t[NBinsPt];
  MesonYieldsWideError = new Double_t[NBinsPt];
  MesonYieldsFuncWideError = new Double_t[NBinsPt];
  MesonYieldsResidualBckFuncWideError = new Double_t[NBinsPt];
  MesonYieldsCorResidualBckFuncWideError = new Double_t[NBinsPt];
  MesonYieldsPerEventWideError = new Double_t[NBinsPt];
  MesonSBWideError = new Double_t[NBinsPt];
  MesonSignWideError = new Double_t[NBinsPt];

  GGYieldsLeftWideError = new Double_t[NBinsPt];
  BckYieldsLeftWideError = new Double_t[NBinsPt];
  MesonYieldsLeftWideError = new Double_t[NBinsPt];
  MesonYieldsFuncLeftWideError = new Double_t[NBinsPt];
  MesonYieldsResidualBckFuncLeftWideError = new Double_t[NBinsPt];
  MesonYieldsCorResidualBckFuncLeftWideError = new Double_t[NBinsPt];
  MesonYieldsLeftPerEventWideError = new Double_t[NBinsPt];
  MesonSBLeftWideError = new Double_t[NBinsPt];
  MesonSignLeftWideError = new Double_t[NBinsPt];


  Mapping_TrueMeson_InvMass_PtBin = new TH1D*[NBinsPt];    // array of histos for pt slices
  Mapping_GG_InvMass_PtBin = new TH1D*[NBinsPt];    // array of histos for pt slices
  Mapping_Back_InvMass_PtBin = new TH1D*[NBinsPt];
  Mapping_BackNorm_InvMass_PtBin = new TH1D*[NBinsPt];
  Mapping_Signal_InvMass_PtBin = new TH1D*[NBinsPt];   
  Signal_InvMassFit_PtBin = new TF1*[NBinsPt];   
  SignalPeak_InvMassFit_PtBin = new TF1*[NBinsPt];   
  SignalBck_InvMassFit_PtBin = new TF1*[NBinsPt];   

  // Histograms for normalization on the left of the peak
  Mapping_BackNorm_InvMass_Left_PtBin = new TH1D*[NBinsPt];
  Mapping_Signal_InvMass_Left_PtBin = new TH1D*[NBinsPt];   

  Signal_InvMassFit_Left_PtBin = new TF1*[NBinsPt];   
  SignalPeak_InvMassFit_Left_PtBin = new TF1*[NBinsPt];   
  SignalBck_InvMassFit_Left_PtBin = new TF1*[NBinsPt];   
}

void CalculateFWHM(TF1 * fFunc, Double_t * MesonIntRange)
{
	cout << "\n \n start calc FWHM \n \n " << endl;
	//FWHM
	fFWHMFunc = fFunc->GetX(fFunc->GetParameter(0)*0.5,fFunc->GetParameter(1), MesonFitRange[1]) - fFunc->GetX(fFunc->GetParameter(0)*0.5,MesonFitRange[0],fFunc->GetParameter(1));
	
	//FWHM error +
	TF1* fFunc_plus;
	fFunc_plus = fFunc;
	fFunc_plus->SetParameter(0,fFunc->GetParameter(0) + fFunc->GetParError(0));
	fFunc_plus->SetParameter(1,fFunc->GetParameter(1) + fFunc->GetParError(1));
	fFunc_plus->SetParameter(2,fFunc->GetParameter(2) + fFunc->GetParError(2));
	fFunc_plus->SetParameter(3,fFunc->GetParameter(3) + fFunc->GetParError(3));
	Double_t FWHM_plus = fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,fFunc_plus->GetParameter(1), MesonFitRange[1]) - fFunc_plus->GetX(fFunc_plus->GetParameter(0)*0.5,MesonFitRange[0],fFunc_plus->GetParameter(1));
	
	//FWHM error -	
	TF1* fFunc_minus;
	fFunc_minus = fFunc;
	fFunc_minus->SetParameter(0,fFunc->GetParameter(0) - fFunc->GetParError(0));
	fFunc_minus->SetParameter(1,fFunc->GetParameter(1) - fFunc->GetParError(1));
	fFunc_minus->SetParameter(2,fFunc->GetParameter(2) - fFunc->GetParError(2));
	fFunc_minus->SetParameter(3,fFunc->GetParameter(3) - fFunc->GetParError(3));
	Double_t FWHM_minus =  fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,fFunc_minus->GetParameter(1), MesonFitRange[1]) -fFunc_minus->GetX(fFunc_minus->GetParameter(0)*0.5,MesonFitRange[0],fFunc_minus->GetParameter(1));
	
	Double_t Error1 = TMath::Abs(fFWHMFunc-FWHM_plus);
	Double_t Error2 = TMath::Abs(fFWHMFunc-FWHM_minus);
	
	if(Error1>=Error2) fFWHMFuncError = Error1;
	if(Error1<Error2) fFWHMFuncError = Error2;
	
	//delete fFunc_plus;
	//delete fFunc_minus;
	cout << "\n \n end calc FWHM \n \n " << endl;
	
}

