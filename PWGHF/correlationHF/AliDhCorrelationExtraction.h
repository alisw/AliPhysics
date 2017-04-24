#ifndef AliDhCorrelationExtraction_H
#define AliDhCorrelationExtraction_H

/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//
//  Base class to extract D-h correlation distribution from analysis task
//
//-----------------------------------------------------------------------
//  Author F.Colamaria
//  INFN Bari
//  fabio.colamaria@cern.ch
//-----------------------------------------------------------------------

#include <iostream>
#include "TObject.h"
#include "TMath.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "THnSparse.h"
#include "AliHFMassFitter.h"

class AliDhCorrelationExtraction : public TObject
{

public:
    
    enum DMesonSpecies {kD0toKpi, kDplusKpipi, kDStarD0pi, kDxHFE};
    enum SandBextraction {kSandBFromFit, kSfromBinCount, kBfromBinCount}; //default for paper: kBfromBinCount
    enum SBscaling {kFitScaling, kBinCountScaling}; //default for paper: kBfromBinCount
    enum selectAnType {kSE, kME};
    enum selectRegion {kSign, kSideb};
    enum MCOrigin {kOrigAll, kOrigC, kOrigB, kOrigLF};  
    enum MCmode {kKine, kReco};
    enum DxHFECorrBins {kDxCorrD0Mass, kDxCorrD0Pt, kDxCorrElPt, kDxCorrdPhi, kDxCorrdEta, kDxCorrPoolbin};
    enum DxHFECorrMCBins {kDxCorrMCD0Mass, kDxCorrMCD0Pt, kDxCorrMCElPt, kDxCorrMCdPhi, kDxCorrMCdEta, kDxCorrMCOriginD0, kDxCorrMCOriginEl,kDxCorrMCGeneratorEl};
    enum DxHFED0Bins {kDxD0Pt, kDxD0Phi, kDxD0Mass, kDxD0Eta, kDxD0Mother};
    enum DxHFEElBins {kDxElPt, kDxElPhi, kDxElEta, kDxElMother, kDxElGenerator};
    enum DxHFEElSource {kAll, kAllEl, kHFE, kNonHFE, kcEl, kbEl};
    enum DxHFED0Source {kcD0=1, kbD0};

    AliDhCorrelationExtraction(); // default constructor
    AliDhCorrelationExtraction(const AliDhCorrelationExtraction &source);
    virtual ~AliDhCorrelationExtraction();

    Bool_t SetDmesonSpecie(DMesonSpecies k);
    void SetSandBextraction(SandBextraction k) {fSandBextraction=k;}
    void SetSBscaling(SBscaling k) {fSBscaling=k;}
    void SetInputFilenameMass(TString filenameMass) {fFileNameMass=filenameMass;}    
    void SetInputFilenameSE(TString filenameSE) {fFileNameSE=filenameSE;}
    void SetInputFilenameME(TString filenameME) {fFileNameME=filenameME;}
    void SetDirNameMass(TString dirNameMass) {fDirNameMass=dirNameMass;}
    void SetListNameMass(TString listMassName) {fListNameMass=listMassName;}
    void SetDirNameSE(TString dirNameSE) {fDirNameSE=dirNameSE;}
    void SetListNameSE(TString listNameSE) {fListNameSE=listNameSE;}
    void SetDirNameME(TString dirNameME) {fDirNameME=dirNameME;}
    void SetListNameME(TString listNameME) {fListNameME=listNameME;}
    void SetMassHistoName(TString massHistoName) {fMassHistoName=massHistoName;}
    void SetSECorrelHistoName(TString correlName) {fSECorrelHistoName=correlName;}
    void SetSECorrelHistoName_DstarBkg(TString correlName_DstarBkg) {fSECorrelHistoName_DstarBkg=correlName_DstarBkg;}
    void SetMECorrelHistoNameSuffix(TString suffix) {fMEsuffix=suffix;}
    void IntegratePtBins(Bool_t intPt=kFALSE) {fIntegratePtBins=intPt;}
    void ReadTTreeOutputFiles(Bool_t treeSE, Bool_t treeME) {fReadTreeSE=treeSE; fReadTreeME=treeME;}
    void SetSubtractSoftPiInMEdistr(Bool_t subtractSoftPiME) {fSubtractSoftPiME=subtractSoftPiME;}
    void SetUseMassVsCentPlots(Bool_t mass2D) {fUseMassVsCentPlots=mass2D;}
    void SetCentralitySelection(Double_t min, Double_t max) {fMinCent=min; fMaxCent=max;} //activated only if both values are != 0
    
    void SetRebinMassPlots(Int_t rebinMassPlots) {fRebinMassPlots=rebinMassPlots;}
    void SetNpTbins(Int_t npt) {fNpTbins=npt;}
    void SetFirstpTbin(Int_t ptFirst) {fFirstpTbin=ptFirst;}
    void SetLastpTbin(Int_t ptLast) {fLastpTbin=ptLast;}
    void SetNumberOfSigmasFitter(Double_t nsigma) {fNumberOfSigmasFitter=nsigma;}
    void SetCorrectPoolsSeparately(Bool_t usePools) {fCorrectPoolsSeparately=usePools;}
    void SetNpools(Int_t npools) {fNpools=npools;}
    void SetDeltaEtaRange(Double_t etaLow=-1., Double_t etaHigh=1) {fDeltaEtaMin=etaLow; fDeltaEtaMax=etaHigh;}
    void SetUseMC(Bool_t useMC) {fUseMC=useMC;}
    void SetUseElSource(Int_t elSource) {fElSource=elSource;}
    void SetUseD0Source(Int_t D0Source) {fD0Source=D0Source;}

    void SetFitRanges(Double_t left, Double_t right) {fLeftFitRange=left; fRightFitRange=right;}
    void SetBkgFitFunction(Int_t func=0) {fBkgFitFunction=func;}
    void SetSignalSigmas(Double_t nsigma=2) {fSignalSigmas=nsigma;}
    void SetAutoSBRange(Bool_t autoSB=kFALSE, Double_t inSigma=0., Double_t outSigma=0.) {fAutoSBRange=autoSB; fSBOuterSigmas=outSigma; fSBInnerSigmas=inSigma;}
    void SetAutoSignRange(Bool_t autoSign=kFALSE) {fAutoSignRange=autoSign;}
    void SetSBSingleBin(Bool_t singleSBbin=kFALSE) {fSBSingleBin=singleSBbin;}
    void SetSBRanges(Double_t* rangesSB1L=0x0, Double_t* rangesSB1R=0x0, Double_t* rangesSB2L=0x0, Double_t* rangesSB2R=0x0);
    void SetSignRanges(Double_t* rangesSignL=0x0, Double_t* rangesSignR=0x0);
    void PrintRanges();
    void PrintSandBForNormal();
    void GetSignalAndBackgroundForNorm(Int_t i, TH1F* &histo);  //evaluates signal and background in 'fSignalSigmas', for trigger normalization and SB correlation rescaling
    void GetSBScalingFactor(Int_t i, TH1F* &histo); //estract sideband scaling factor
    TH2D* GetCorrelHisto(Int_t SEorME, Int_t SorSB, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax, Int_t softPiME=1);
    TH2D* GetCorrelHistoDzero(Int_t SEorME, Int_t SorSB, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax, Int_t softPiME=1);
    TH2D* GetCorrelHistoDplus(Int_t SEorME, Int_t SorSB, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax);
    TH2D* GetCorrelHistoDstar(Int_t SEorME, Int_t SorSB, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax);
    TH2D* GetCorrelHistoDxHFE(Int_t SEorME, Int_t SorSB, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax);
    TH2D* GetCorrelHistoDzeroTTree(Int_t SEorME, Int_t SorSB, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax);
    TH2D* GetCorrelHistoDplusTTree(Int_t SEorME, Int_t SorSB, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax);
    TH2D* GetCorrelHistoDstarTTree(Int_t SEorME, Int_t SorSB, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax);
    TH2D* GetCorrelHisto_MC(Int_t SEorME, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax, Int_t orig, Int_t softPiME=1);
    TH2D* GetCorrelHistoDzero_MC(Int_t SEorME, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax, Int_t orig, Int_t softPiME=1);
    TH2D* GetCorrelHistoDplus_MC(Int_t SEorME, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax, Int_t orig);
    TH2D* GetCorrelHistoDstar_MC(Int_t SEorME, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax, Int_t orig);
    TH2D* GetCorrelHistoDxHFE_MC(Int_t SEorME, Int_t pool, Int_t pTbin, Double_t thrMin, Double_t thrMax, Int_t orig);
    void NormalizeMEplot(TH2D* &histoME, TH2D* &histoMEsoftPi); //normalize ME plots to the average value of the 4 'central' bins
    void RescaleSidebandsInMassBinEdges(Int_t i); //readjust SB scaling factor if single bin is used & ranges passed from outside & ranges don't match bin edges
    void MergeMassPlotsVsPt(); //merge spectra from mass-pT bins in larger correlation-pT bins (as if you have a single pT bin)
    void MergeMassPlotsVsPt_MC(); //merge spectra from mass-pT bins in larger correlation-pT bins (as if you have a single pT bin)
    void MergeCorrelPlotsVsPt(THnSparse* &hsparse, Int_t SEorME, Int_t SorSB=0, Int_t pool=0); //merge THnSparse from mass-pT bins in correlation-pT bins (as if you have 1 pT bin)
    void MergeCorrelPlotsVsPt_MC(THnSparse* &hsparse, Int_t SEorME, Int_t pool=0, Int_t orig=0); //merge THnSparse from mass-pT bins in correlation-pT bins (as if you have 1 pT bin)
    void MergeCorrelPlotsVsPtTTree(TH3D* &h3D, Int_t SEorME, Int_t SorSB, Int_t pool, Double_t thrMin, Double_t thrMax); //as before, but for TTree case
    
    void ClearObjects();

    void SetDebugLevel(Int_t debug) {fDebug=debug;}
    void SetHistoStyle(TH1F* &histo, Int_t colour);
    void SetHistoCorrStyle(TH1D* &histo);
    void DefinePaveText(TPaveText* &paveText, Double_t xmin, Double_t ymin, Double_t xmax, Double_t ymax, TString name);
    void SetPtRanges(Int_t ptBin, THnSparse* thn, Int_t dimension);
    void SetElSource(Int_t type, THnSparse* thn, Int_t dimension);
    void SetD0Source(Int_t type, THnSparse* thn, Int_t dimension);

    void AddOriginType(TString suffix, MCOrigin orig) {fMCOriginSuffix.push_back(suffix); fMCOriginType.push_back(orig);} //for MC case
    void SetRecoMode(MCmode mode) {fMCmode=mode;}

    Bool_t ReadInputs(); //reads input files and loads lists of plots
    Bool_t FitInvariantMass(); //method to perform invariant mass fit via AliHFMassFitter
    Bool_t ExtractCorrelations(Double_t thrMin, Double_t thrMax); //method to retrieve the bkg subtracted and ME corrected correlation distributions
    Bool_t ExtractCorrelations_MC(Double_t thrMin, Double_t thrMax); //method to retrieve the ME corrected correlation distributions in MC
    Bool_t ExtractCorrelations_MC_Orig(Double_t thrMin, Double_t thrMax, Int_t orig); //method to retrieve the ME corrected correlation distributions in MC
    Bool_t ExtractNumberOfTriggers_MC();
    Bool_t DoSingleCanvasMCPlot(Double_t thrMin, Double_t thrMax);
    void DrawMCClosure(Int_t nOrig, Int_t binMin, Int_t binMax, Double_t thrMin, Double_t thrMax);
    
private:
    
    TFile *fFileMass; //file containing the mass histograms
    TFile *fFileSE; //file containing the analysis output SE
    TFile *fFileME; //file containing the analysis output ME
    
    TDirectoryFile *fDirMass; // TDirectory for mass histos
    TDirectoryFile *fDirSE; // TDirectory for SE info
    TDirectoryFile *fDirME; // TDirectory for ME info
    
    TList *fMassList; // TList with D mass
    TList *fTracksList; // TList with track properties
    TList *fSECorrelationList; // TList with Correlations from SE
    TList *fMECorrelationList; // TList with Correlations from ME

    DMesonSpecies fDmesonSpecies;
    SandBextraction fSandBextraction;
    SBscaling fSBscaling;
    TString fDmesonLabel;
    TString fFileNameMass;
    TString fFileNameSE;
    TString fFileNameME;  
    TString fDirNameMass;
    TString fListNameMass;
    TString fDirNameSE; 
    TString fListNameSE; 
    TString fDirNameME; 
    TString fListNameME;
    TString fMassHistoName;
    TString fSECorrelHistoName;
    TString fSECorrelHistoName_DstarBkg;
    TString fMEsuffix;

    Int_t fRebinMassPlots;    
    Int_t fNpTbins;
    Int_t fFirstpTbin;
    Int_t fLastpTbin;
    Double_t fNumberOfSigmasFitter;
    Bool_t fCorrectPoolsSeparately;
    Int_t fNpools;
    Bool_t fReadTreeSE;
    Bool_t fReadTreeME;
    Bool_t fSubtractSoftPiME;
    Bool_t fUseMassVsCentPlots;		//don't use histMass plots, but project histMass2D plots (for offline)
    Double_t fMinCent;
    Double_t fMaxCent;

    Double_t *fDmesonFitterSignal;
    Double_t *fDmesonFitterSignalError;
    Double_t *fDmesonFitterBackground;
    Double_t *fDmesonFitterBackgroundError;
    Double_t *fDMesonFitterSBCand;
    Double_t *fDMesonFitterSBCandErr;
    Double_t *fDmesonFitterMean;
    Double_t *fDmesonFitterMeanError;
    Double_t *fDmesonFitterSigma;
    Double_t *fDmesonFitterSigmaError;
    Double_t *fDmesonFitterSignificance;
    Double_t *fDmesonFitterSignificanceError;
    Double_t *fDmesonFitterSOverB;

    Double_t fLeftFitRange;
    Double_t fRightFitRange;
    Int_t fBkgFitFunction;
    Double_t fSignalSigmas;
    Bool_t fAutoSBRange;
    Bool_t fAutoSignRange;
    Double_t fSBOuterSigmas;
    Double_t fSBInnerSigmas;
    Bool_t fSBSingleBin;
    Double_t fDeltaEtaMin;
    Double_t fDeltaEtaMax;
    Bool_t fUseMC;
    Int_t fElSource;
    Int_t fD0Source;

    Double_t *fSignalCorrel;
    Double_t *fBackgrCorrel;
    Double_t *fRangesSignL;
    Double_t *fRangesSignR;
    Double_t *fRangesSB1L;
    Double_t *fRangesSB1R;
    Double_t *fRangesSB2L;
    Double_t *fRangesSB2R;
    Double_t *fScaleFactor;
    Double_t *fSignalCorrelMC_c;
    Double_t *fSignalCorrelMC_b;    

    Bool_t fIntegratePtBins;

    Int_t fDebug;

    TF1 **fMassFit;
    TF1 **fBkgFit;

    TH1F **fMassHisto;

    std::vector<TString>  fMCOriginSuffix;    //container of suffixes for THnSparse for different origins
    std::vector<Int_t>    fMCOriginType;      //container of specificators of origins
    MCmode		  fMCmode;	      //kine or reco analysis (changes just the filenames for output, for now)

    ClassDef(AliDhCorrelationExtraction,4); // class for plotting HF correlations

};

#endif

