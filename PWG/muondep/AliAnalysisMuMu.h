#ifndef ALIANALYSISMUMU_H
#define ALIANALYSISMUMU_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

///
/// AliAnalysisMuMu : Facade class of all subclass called to digest/plot/massage results from
/// AliAnalysisTaskMuMu
///
/// author : Laurent Aphecetche (Subatech), Javier Martin Blanco, Benjamin Audurier

#include "AliAnalysisMuMuBinning.h"
#include "AliAnalysisMuMuFnorm.h"
#include "TNamed.h"
#include <map>
#include <set>
#include <string>
#include <TString.h>
#include <vector>
#include "RQ_OBJECT.h"
#include "TH1.h"
#include "TH2.h"

class AliAnalysisMuMuConfig;
class AliAnalysisMuMuResult;
class AliAnalysisMuMuJpsiResult;
class AliAnalysisMuMuSpectra;
class AliCounterCollection;
class AliMergeableCollection;
class TF1;
class TFile;
class TGraph;
class TGraphErrors;
class TH1;
class TMap;

class AliAnalysisMuMu : public TObject, public TQObject
{

public:

    AliAnalysisMuMu(const char* filename, AliAnalysisMuMuConfig& config);

    AliAnalysisMuMu(
      const char* filename="",
      const char* associatedSimFileName="",
      const char* associatedSimFileName2="",
      const char* configurationFile="");

    virtual ~AliAnalysisMuMu();

    /* Basic checks */
    void BasicCounts(
      Bool_t detailTrigger=kFALSE,
      ULong64_t* totalNmb=0x0,
      ULong64_t* totalNmsl=0x0,
      ULong64_t* totalNmul=0x0);

    void TriggerCountCoverage(
      const char* triggerList,
      Bool_t compact=kTRUE,
      Bool_t orderByTriggerCount=kTRUE);

    void SelectRunByTrigger(const char* triggerList);

    AliAnalysisMuMuSpectra* FitParticle(
      const char* particle,
      const char* trigger,
      const char* eventType,
      const char* pairCut,
      const char* centrality,
      const AliAnalysisMuMuBinning& binning,
      const char* spectraType="minv",
      Bool_t corrected=kFALSE);

    AliAnalysisMuMuSpectra* CorrectSpectra(
      const char* type,
      const char* flavour="",
      const char* accEffSubResultName="");

    void PrintDistribution(
      const char              * binType="Y",
      const char              * what="NofJPsi",
      const char              * sResName="",
      const char              * ColSys="PbPb",
      Bool_t divideByBinWidth =kTRUE,
      Bool_t AccEffCorr       =kFALSE);

    void PrintFitParam(
      const char* particle ="PSI",
      const char* param = "mJPsi",
      const char* binType="PT",
      const char* subresult="CB2VWG_2.4_4.5_SP1.2",
      const char* printDirectoryPath="",
      Bool_t AccEffCorr =kFALSE
      )const;

    void ComputeDimuonRawCount(
      const Double_t rlow   = 2.8,
      const Double_t rhight = 3.4,
      const char            * binType="pt",
      const char            * binRangeExluded="PT_BENJ_00.00_00.30,PT_BENJ_01.00_08.00,PT_BENJ_00.30_01.00",
      const char            * flavour="BENJ",
      Bool_t corrected      =kFALSE );

    void ComputePPCrossSection(
      const char        * binType="PT",
      const char        * particle ="PSI",
      const char        * what ="CorrNofJPsi",
      const char        * externfile = "",
      const char        * externfile2 = "",
      Bool_t print      =kFALSE,
      Bool_t AccEffCorr =kFALSE);

    TH2* ComputeSPDCorrection(
      const char       * type="oneOverAccEff",
      const char       * eventSel="PSALL",
      const char       * triggerSel="ANY",
      Bool_t bkgReject =kTRUE);

    void ComputeFnorm();

    void ComputeNumberOfEvent();

    void ComputeMeanFnorm(
      const char* triggerCluster="MUON",
      const char* eventSelection="PSALL",
      const char* what="psi",
      const char* quantity="ntrcorr",
      const char* flavour="D2H");

    void ComputeFnormWeightedMeanGraphs(
      AliAnalysisMuMuFnorm::ETriggerType refTrigger = AliAnalysisMuMuFnorm::kMUL,
      const char          * patternOrList= "",
      const char          * graphName= "");

    void ComputeFnormScalers(
      AliAnalysisMuMuFnorm::ETriggerType refTrigger = AliAnalysisMuMuFnorm::kMUL,
      Bool_t PileUpCorr   =kFALSE);

    void ComputeIntFnormFromCounters(
      AliAnalysisMuMuFnorm::ETriggerType refTrigger = AliAnalysisMuMuFnorm::kMUL,
      Bool_t PileUpCorr =kFALSE);

    void PlotYiedWSyst(const char* triggerCluster="MUON");

    void ComputeRelativeValueAndSESystematics(
      const char* quantity,const char* flavour,const char* value2Test,
      const char* binListToExclude,
      const char* fNormType="mean",
      const char* evSelInt="PSALL",
      const char* evSelDiff="PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00",
      const char* triggerCluster="MUON");

    void ComputeJpsiYield(
      const char        * binType="INTEGRATED",
      const char        * what="NofJPsi",
      const char        * externfile1="externFile_PT.txt",
      const char        * externfile2="externFile_CENT.txt",
      const char        * sResName="",
      const char        * beamYear="",
      Bool_t AccEffCorr =kFALSE);

    void ComputeJpsiMPt(
      Bool_t relative=kTRUE,
      const char* evSelInt="PSALL",
      const char* evSelDiff="PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00",
      const char* spectra="PSI-NTRCORR-AccEffCorr-MeanPtVsMinvUS",
      const char* sResName="");

    void ComputeMBXSectionFractionInBins(
      const char* filePileUpCorr="",
      const char* eventSelection="PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00",
      const char* what="psi",
      const char* quantity="ntrcorr",
      const char* flavour="D2H");

    Double_t ErrorPropagationAxBoverCxD(Double_t a,Double_t b,Double_t c,Double_t d);

    void TwikiOutputFnorm(const char* series="FnormOffline2PUPS,FnormScalersPUPS,FnormBest2,RelDifFnormScalersPUPSvsFnormOffline2PUPS,FnormScalersPUVDM,RelDifFnormScalersPUPSvsFnormScalersPUVDM") const;

    AliAnalysisMuMuSpectra* ComputeYield(const char* type, const char* flavour="",const char* accEffSubResultName="PSICB2");

    void CleanAllSpectra();

    void CleanFNorm();

    void RAAasGraphic(
      const char        * particle="PSI",
      const char        * binType="PT",
      const char        * externfile1="externFile_PT.txt",
      const char        * externfile2="externFile_CENT.txt",
      const char        * RefCent ="V0M_00.00_90.00",
      Bool_t print      = kFALSE,
      Bool_t AccEffCorr =kFALSE)const;

    void DrawFitResults(
      const char        * particle="PSI",
      const char        * binType="INTEGRATED",
      const char        * subresults="",
      Bool_t AccEffCorr =kFALSE)const;

    void PrintNofParticle(
      const char        * particle="PSI",
      const char        * what="NofJPsi",
      const char        * binType="PT",
      Bool_t AccEffCorr =kFALSE) const;

    AliAnalysisMuMuSpectra* RABy(const char* type="", const char* direction="pPb");

    TGraph* PlotEventSelectionEvolution(
    const char         * trigger1="CINT7-B-NOPF-MUON",
    const char         * event1="PSALL",
    const char         * trigger2="CINT7-B-NOPF-MUON",
    const char         * event2="PSALL",
    Bool_t drawFills   =kFALSE,
    Bool_t asRejection =kTRUE) const;

    Bool_t Upgrade();

    Bool_t Upgrade(const char* filename);

    TObjArray* CompareJpsiPerCMUUWithBackground(
      const char* jpsiresults="results.root",
      const char* backgroundresults="background.lhc11d.root");

    TGraph* CompareJpsiPerCMUUWithSimu(
      const char* realjpsiresults="results.root",
      const char* simjpsiresults="results.sim.root");


    static TFile* FileOpen(const char* file);

    static TString ExpandPathName(const char* file);

    Bool_t GetCollections(
      const char            * rootfile,
      const char            * subdir,
      AliMergeableCollection*& oc,
      AliCounterCollection  *& cc,
      AliAnalysisMuMuBinning*& bin,
      std                   ::set<int>& runnumbers);

    AliAnalysisMuMuSpectra* GetSpectra(
      const char* what,
      const char* flavour="") const;

    AliAnalysisMuMuSpectra* GetMCSpectra(
      const char* what ,
      const char* EventSelection ="ALL" ,
      const char* DimuonTrigger="ANY",
      const char* Centrality="V0A",
      const char* PairSelectionKey="pALLPAIRYPAIRPTIN0.0-12.0RABSMATCHLOWETAPDCA",
      const char* flavour="BENJ") const;

    TH1* PlotAccEfficiency(const char* whatever="PSI-INTEGRATED");

    UInt_t GetSum(
      AliCounterCollection& cc,
      const char* triggerList,
      const char* eventSelection,
      Int_t runNumber=-1);

    ULong64_t GetTriggerScalerCount(const char* triggerList, Int_t runNumber);

    Int_t Jpsi(
      const char           * what="integrated",
      const char           * binningFlavour="",
      Bool_t fitmPt        =kTRUE,
      Bool_t onlyCorrected =kTRUE);

    Bool_t IsSimulation() const;

    AliMergeableCollection* OC() const { return fMergeableCollection; }
    AliCounterCollection  * CC() const { return fCounterCollection; }
    AliAnalysisMuMuBinning* BIN() const { return fBinning; }

    void Print(Option_t* opt="") const;

    const std::set<int>& RunNumbers() const { return fRunNumbers; }

    void DrawMinv(
      const char* type,
      const char* particle,
      const char* trigger,
      const char* eventType,
      const char* pairCut,
      const char* centrality,
      const char* subresultname="",
      const char* flavour="") const;

    void DrawMinv(
      const char* type="PT",
      const char* particle="PSI",
      const char* flavour="",
      const char* subresultname="") const;

    Bool_t SetCorrectionPerRun(const TGraph& corr, const char* formula="");

    void UnsetCorrectionPerRun();

    void ExecuteCanvasEvent(Int_t event, Int_t px, Int_t py, TObject *sel);

    std::vector<Double_t> GetMCCB2Tails(const AliAnalysisMuMuBinning::Range& bin) const;

    AliAnalysisMuMu* SIM() const { return fAssociatedSimulation; }

    AliAnalysisMuMu* SIM2() const { return fAssociatedSimulation2; }

    AliAnalysisMuMuSpectra* SPECTRA(const char* fullpath) const;

    void SetParticleName(const char* particleName) { fParticleName = particleName; }

    const char* GetParticleName() { return fParticleName; }

    void Update();

    // AliAnalysisMuMuConfig* Config();

    AliAnalysisMuMuConfig* Config() const { return fConfig; }

    void SetConfig(const AliAnalysisMuMuConfig& config);

    void SetCentralitySelectionList(const char* centralitySelectionList);

private:
    AliAnalysisMuMu(const AliAnalysisMuMu& rhs); // not implemented on purpose
    AliAnalysisMuMu& operator=(const AliAnalysisMuMu& rhs); // not implemented on purpose

    Bool_t SetParticleNameFromFileName(const char* filename);

    void ShowList(const char* title, const TString& list, const char separator=',') const;

    TFile* ReOpen(const char* filename, const char* mode) const;

    TString First(const TString& list) const;

    void GetParametersFromMC(TString& fitType, const char* pathCentrPairCut, const char* spectraName, AliAnalysisMuMuBinning::Range* bin) const;
    void GetParametersFromResult(TString& fitType, AliAnalysisMuMuJpsiResult* minvResult) const;


    void GetCollectionsFromAnySubdir(TDirectory& dir,
                                    AliMergeableCollection*& oc,
                                     AliCounterCollection*& cc,
                                     AliAnalysisMuMuBinning*& bin);

    void GetFileNameAndDirectory(const char* filename);

    void LoadStyles();

private:

    void SetNofInputParticles(AliAnalysisMuMuJpsiResult& r);


    TString fFilename; // file containing the result collections (of objects and counters) from AliAnalysisTaskMuMu
    TString fDirectory; // directory, within fFilename, containing the actual objects

    AliCounterCollection* fCounterCollection; // collection of counters in file

    AliAnalysisMuMuBinning* fBinning; // binning

    AliMergeableCollection* fMergeableCollection; // collection of objects in file

    std::set<int> fRunNumbers; // run numbers

    TGraph* fCorrectionPerRun; // correction factor per run

    AliAnalysisMuMu* fAssociatedSimulation; // associated simulations (if any)
    AliAnalysisMuMu* fAssociatedSimulation2; // second associated simulations (if any)

    TString fParticleName; // Name of the simulated particle in the associated simulations

    AliAnalysisMuMuConfig* fConfig; // configuration

    ClassDef(AliAnalysisMuMu,12) // class to analysis results from AliAnalysisTaskMuMuXXX tasks
};

#endif
