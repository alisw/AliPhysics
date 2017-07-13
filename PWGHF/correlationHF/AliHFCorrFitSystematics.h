#ifndef ALIHFCORRFITSYSTEMATICS_H
#define ALIHFCORRFITSYSTEMATICS_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
/*____________________________________________________________
 | 
 | Class for the standard calculation and visualization of 
 | systematic uncertainties on the fit of the azimuthal correlations
 |                                                             |
 |  Author:                                                    |
 |     Sandro Bjelogrlic (sandro.bjelogrlic@cern.ch)                                       |
 |_____________________________________________________________*/


#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include <TGraphAsymmErrors.h>
#include "AliHFCorrFitter.h"

class TH1F;
class TCanvas;
class TF1;

class AliHFCorrFitSystematics{
    
public:
    
    
    //enum
  enum SystematicModes{kFree=0, kLowestPoint = 1, k2PointsAtPiHalf = 2, k4PointsAtPiHalf=4, kTransverse = 5, kNLowest=6, k2GausNS=7, kBinCount = -1, kTransverseUppStatUnc = 10, kTransverseLowStatUnc = 20,  kMinVar=100, kMaxVar=200, kv2Mod = 300};
  enum SystCombination{kSumQuadr = 1, kMax = 2, kRMS =3, kEnvelope_RMS_BaselStat=4};
    
    AliHFCorrFitSystematics();
    virtual ~AliHFCorrFitSystematics();
    Bool_t AddHistoToFit(TString path);
    Bool_t AddHistoToFit(TString path, TString cname, TString  havname, TString havmin, TString havmax);
  
    // setters
    void SetMinCorrelatedSystematicsDeltaPhi(Double_t mincorrsyst){// set min correleted systematics manually
        fVecMinCorrelatedSyst.push_back(mincorrsyst);

    }
    void SetMaxCorrelatedSystematicsDeltaPhi(Double_t maxcorrsyst){// set max correleted systematics
        fVecMaxCorrelatedSyst.push_back(maxcorrsyst);
    }
    void SetDptLowEdge(Double_t ptlow){// set D meson pt lower edges
        fVecLowEdgeDpt.push_back(ptlow);
    }
    void SetDptUpEdge(Double_t ptup){// set D meson pt upper edges
        fVecUpEdgeDpt.push_back(ptup);
    }
    
    void AddSystematicMode(SystematicModes mode, Bool_t isReference);
    
    void AddSystematicMode(SystematicModes mode, Bool_t isReference, Double_t min, Double_t max, Int_t minNpoints=0,AliHFCorrFitter::FunctionType functype=AliHFCorrFitter::kTwoGausPeriodicity);
    void SetIsFittingTemplate(Bool_t isfitTempl){fIsFittingTemplate=isfitTempl;}
    
    
    void AddSystematicMode(SystematicModes mode){
        
        AddSystematicMode(mode,kFALSE);
        //fVecSystMode.push_back(mode);
    }
    
    void Setv2ForSystematics(Double_t v2had, Double_t v2D){fv2had = v2had; fv2Dmeson = v2D; }
    void Setv2ForSystematics(Double_t v2had, Int_t nbins, std::vector<Double_t>& v2Dbins) {fv2had = v2had; for(Int_t i;i<nbins;i++) fv2DmesonVsPt.push_back(v2Dbins.at(i));}
    
    
    void CheckBaselineRanges();
    void SetUseCorrelatedSystematics(Bool_t k) {fUseCorrelatedSystematics = k;}
    void SetIspPb(Bool_t k){fIspPb = k;}
    void SetIsv2DvsPt(Bool_t k){fV2DvsPt = k;}
    void SetUseCorrelatedSystematicsForWidths(Bool_t k){fUseCorrelatedSystematicsForWidths = k;}
    
    void SetPlotV2SystSeparately(Bool_t k){
        if(!fIspPb){std::cout << "Warning:: Not p-Pb " << std::endl;}
        
        fPlotV2SystSeparately = k;}
    
    
    
    void SetCombineSystematicsMode(SystCombination mode){fCombineSystematics = mode;}
    
    void SetFitRange(Double_t min, Double_t max){fMinFitRange = min; fMaxFitRange = max;}
    void SetReferenceBaselineEstimationRange(Double_t min, Double_t max){fMinReferenceBaselineEstimationRange = min; fMaxReferenceBaselineEstimationRange = max;}
  
    
    void SetAssociatedTrackPtMin(Double_t pt){fAssocTrackPtMin = pt;}
    void SetAssociatedTrackPtMax(Double_t pt){fAssocTrackPtMax = pt;}
    
    Int_t GetBinIndex(Int_t ptbinindex, Int_t systematicIndex){
       // cout << "Inputs = " << fVecSize << endl;
        Int_t value = systematicIndex*fVecSize + ptbinindex;
        return value;}
    
    
    
    void CreateOutputFile(TString filename){
        
        fOutputFileName = filename;
         fOutputFile = new TFile(fOutputFileName.Data(),"RECREATE");
        
    }
    
    
    void CheckHisto(Int_t i);
    Bool_t SetUpSystematicsFitter();
    void Fitv2Systematics();
    Bool_t RunFits();
    void DrawTotalSystematicsOnCanvas(TString variable, TLegend * legend);
    void DrawBaselineSystematicsOnCanvas(TH1D * histo, Int_t iSystMode, Double_t * array, TLegend * legend);
    void DrawRMSBaselineSystematicOnCanvas(TH1D * histoinput, Double_t * arraymin, Double_t * arraymax , TLegend * legend, Bool_t isSystBaseline);
    void DrawTotalBaselineSystematicOnCanvas(TH1D * histoinput, Double_t * arraymin, Double_t * arraymax , TLegend * legend, Bool_t isSystBaseline);
    void PrintAllSystematicsOnShell();
    Bool_t DrawFinalPlots(){return DrawFinalPlots(kTRUE,kTRUE,kFALSE,kFALSE,kTRUE);}
    Bool_t DrawFinalPlots(Bool_t drawNSy, Bool_t drawNSsigma,Bool_t drawASy, Bool_t drawASsigma,Bool_t drawPed);
    Bool_t DrawFinalCorrelationPlot(){
        return DrawFinalCorrelationPlot(kTRUE,kTRUE,kFALSE,kFALSE,kTRUE);
    };
    Bool_t DrawFinalCorrelationPlot(Bool_t drawNSy, Bool_t drawNSsigma,Bool_t drawASy, Bool_t drawASsigma,Bool_t drawPed);
    void DefinePaveText();
    void SetHisto(TH1D *&outputhist, TH1D * inputhist);
    //Bool_t ComputeSystematics();// main function that will compute all the systematics
    Bool_t ComputeRatios(TH1D *&historef, TH1D *&histosys, TH1D *&outputhisto);
    
    
    void SetOutputDirectory(TString outputdirectory){fOutputDirectory = outputdirectory;}
    void SaveCanvasesDotC();
    void SaveCanvasesDotRoot();
    void SaveCanvasesDotEps();
    void SaveCanvasesDotPng();
    void SaveCanvasesDotPdf();
    
    
private:
    //functions
    Bool_t CheckSize();
    Bool_t CheckDiffSystematicsSize();

    Bool_t fIsFittingTemplate;
    //standard vectors in size of D meson pt
    std::vector<TH1D*> fVecHisto;
    std::vector<TH1D*> fVecHistoMinVar;
    std::vector<TH1D*> fVecHistoMaxVar;
    std::vector<TCanvas*> fCorrelationCanvas;
    std::vector<TGraphAsymmErrors*> fCorrelationGraphAsymmErrors;
    std::vector<TF1*> fFitFunctions;
    
    std::vector<Double_t> fVecMinCorrelatedSyst;
    std::vector<Double_t> fVecMaxCorrelatedSyst;
    
    std::vector<Double_t> fVecLowEdgeDpt;
    std::vector<Double_t> fVecUpEdgeDpt;
    
    std::vector<Double_t> fv2DmesonVsPt;
    
    //standard vectors in size of Different modes
    std::vector<SystematicModes> fVecSystMode;
    std::vector<AliHFCorrFitter::FunctionType> fFitFuncType;
    std::vector<Double_t> fMinBaselineEstimationRange;
    std::vector<Double_t> fMaxBaselineEstimationRange;
    std::vector<Int_t> fNMinPointsBaselineEstimationRange;
    std::vector<Bool_t> fVecIsReference;
    
    //variables
    Int_t fEntries;
    Int_t fVecSize;
    Int_t fVecSystModesSize;
    Int_t fDim;
    Int_t fReferenceIndex;
    
    Bool_t fIsReferenceAlreadySet;
    Bool_t fUseCorrelatedSystematics;
    Bool_t fUseMaximumVariation; // uses the maximum variation of the systematics on the baseline and the min and max variation, otherwise sums them in quadrature
    Bool_t fIspPb;
    Bool_t fUseCorrelatedSystematicsForWidths;
    Bool_t fPlotV2SystSeparately;
    
    Bool_t fSaveDotC;
    Bool_t fSaveRoot;
    Bool_t fSavePng;
    Bool_t fSaveEps;
    Bool_t fSavePdf;

    Bool_t fV2DvsPt;
    
    Double_t fMinFitRange;
    Double_t fMaxFitRange;
    
    Double_t fv2had;
    Double_t fv2Dmeson;
 
    
    Double_t fAssocTrackPtMin;
    Double_t fAssocTrackPtMax;
    
    Double_t fMinReferenceBaselineEstimationRange;
    Double_t fMaxReferenceBaselineEstimationRange;
    Double_t fMinNPointsReferenceBaseline;
    Double_t *fValueNSYield;
    Double_t *fRatioNSYield;
    Double_t *fValueNSSigma;
    Double_t *fRatioNSSigma;
    Double_t *fValueASYield;
    Double_t *fRatioASYield;
    Double_t *fValueASSigma;
    Double_t *fRatioASSigma;
    Double_t *fValuePedestal;
    Double_t *fRatioPedestal;
    
    
    Double_t *fValuev2NSYield;
    Double_t *fValuev2NSSigma;
    Double_t *fValuev2ASYield;
    Double_t *fValuev2ASSigma;
    Double_t *fValuev2Pedestal;
    
    Double_t *fSystValuev2NSYield;
    Double_t *fSystValuev2NSSigma;
    Double_t *fSystValuev2ASYield;
    Double_t *fSystValuev2ASSigma;
    Double_t *fSystValuev2Pedestal;
    
    Double_t *fValueSystematicBaselineNSYield;
    Double_t *fValueSystematicBaselineNSSigma;
    Double_t *fValueSystematicBaselineASYield;
    Double_t *fValueSystematicBaselineASSigma;
    Double_t *fValueSystematicBaselinePedestal;

    Double_t *fValueSystematicBaseline_FromBaselStatUp_NSYield;
    Double_t *fValueSystematicBaseline_FromBaselStatUp_NSSigma;
    Double_t *fValueSystematicBaseline_FromBaselStatUp_ASYield;
    Double_t *fValueSystematicBaseline_FromBaselStatUp_ASSigma;
    Double_t *fValueSystematicBaseline_FromBaselStatUp_Pedestal;

    Double_t *fValueSystematicBaseline_FromBaselStatLo_NSYield;
    Double_t *fValueSystematicBaseline_FromBaselStatLo_NSSigma;
    Double_t *fValueSystematicBaseline_FromBaselStatLo_ASYield;
    Double_t *fValueSystematicBaseline_FromBaselStatLo_ASSigma;
    Double_t *fValueSystematicBaseline_FromBaselStatLo_Pedestal;
    
    Double_t *fRMSRelative_NSYield;
    Double_t *fRMSRelative_NSSigma;
    Double_t *fRMSRelative_ASYield;
    Double_t *fRMSRelative_ASSigma;
    Double_t *fRMSRelative_Pedestal;    

    Double_t *fValueSystematicNSYieldUp;
    Double_t *fValueSystematicNSSigmaUp;
    Double_t *fValueSystematicASYieldUp;
    Double_t *fValueSystematicASSigmaUp;
    Double_t *fValueSystematicPedestalUp;
    
    Double_t *fValueSystematicNSYieldLow;
    Double_t *fValueSystematicNSSigmaLow;
    Double_t *fValueSystematicASYieldLow;
    Double_t *fValueSystematicASSigmaLow;
    Double_t *fValueSystematicPedestalLow;
    
    TString fOutputFileName;
    TString fOutputDirectory;
    
    //class enums
    SystematicModes fSystMode;
    SystCombination fCombineSystematics;
    Int_t fDmeson;
    
    //(ali)root objects
    AliHFCorrFitter *fFitter;
    
    TFile *fOutputFile;
    TH1D *fReferenceHistoNSYield;
    TH1D *fReferenceHistoNSSigma;
    TH1D *fReferenceHistoASYield;
    TH1D *fReferenceHistoASSigma;
    TH1D *fReferenceHistoPedestal;
    
    
    TH1D *fValueHistoNSYield;
    TH1D *fRatioHistoNSYield;
    TH1D *fValueHistoNSSigma;
    TH1D *fRatioHistoNSSigma;
    TH1D *fValueHistoASYield;
    TH1D *fRatioHistoASYield;
    TH1D *fValueHistoASSigma;
    TH1D *fRatioHistoASSigma;
    TH1D *fValueHistoPedestal;
    TH1D *fRatioHistoPedestal;
    
    
    TH1D **fRMSHistoNSYield;
    TH1D **fRMSHistoNSSigma;
    TH1D **fRMSHistoASYield;
    TH1D **fRMSHistoASSigma;
    TH1D **fRMSHistoPedestal;
    
    
    TH1D *fSystematicSourcesNSYield;
    TH1D *fSystematicSourcesNSSigma;
    TH1D *fSystematicSourcesASYield;
    TH1D *fSystematicSourcesASSigma;
    TH1D *fSystematicSourcesPedestal;
    
    TCanvas *fCanvasSystematicSourcesNSYield;
    TCanvas *fCanvasSystematicSourcesNSSigma;
    TCanvas *fCanvasSystematicSourcesASYield;
    TCanvas *fCanvasSystematicSourcesASSigma;
    TCanvas *fCanvasSystematicSourcesPedestal;
    
    TCanvas *fCanvasTotalSystematicSourcesNSYield;
    TCanvas *fCanvasTotalSystematicSourcesNSSigma;
    TCanvas *fCanvasTotalSystematicSourcesASYield;
    TCanvas *fCanvasTotalSystematicSourcesASSigma;
    TCanvas *fCanvasTotalSystematicSourcesPedestal;
    
    
    TCanvas *fCanvasFinalTrendNSYield;
    TCanvas *fCanvasFinalTrendNSSigma;
    TCanvas *fCanvasFinalTrendASYield;
    TCanvas *fCanvasFinalTrendASSigma;
    TCanvas *fCanvasFinalTrendPedestal;
    TCanvas *fCanvasVariationBaselineTrendPedestal;
    
    TH1D *fFinalTrendNSYield;
    TH1D *fFinalTrendNSSigma;
    TH1D *fFinalTrendASYield;
    TH1D *fFinalTrendASSigma;
    TH1D *fFinalTrendPedestal;
    
    TGraphAsymmErrors *fFullSystematicsNSYield;
    TGraphAsymmErrors *fv2SystematicsNSYield;
    TGraphAsymmErrors *fFullSystematicsNSSigma;
    TGraphAsymmErrors *fv2SystematicsNSSigma;
    TGraphAsymmErrors *fFullSystematicsASYield;
    TGraphAsymmErrors *fv2SystematicsASYield;
    TGraphAsymmErrors *fFullSystematicsASSigma;
    TGraphAsymmErrors *fv2SystematicsASSigma;
    TGraphAsymmErrors *fFullSystematicsPedestal;
    TGraphAsymmErrors *fv2SystematicsPedestal;
    TGraphAsymmErrors *fBaselineVariationSystematicsPedestal;
    TCanvas *fCanvasRefernce;
    TCanvas **fCanvasFitting;

    ClassDef(AliHFCorrFitSystematics,4);    
    
};

#endif

    

