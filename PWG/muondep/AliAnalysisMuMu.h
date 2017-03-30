#ifndef ALIANALYSISMUMU_H
#define ALIANALYSISMUMU_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/**
 
  @ingroup pwg_muondep_mumu

  @class AliAnalysisMuMu

  @brief Facade class of classes called to  digest/plot/massage results from AliAnalysisTaskMuMu

  @details This class calls other classes from the offline AliAnalysisMuMu framework to read/analysis/play with the output file from the online AliAnalysisMuMu framework.

  Methods are mostly writen as interlocked loops over event/trigger/centrality/cut and delegate (usually) the task to another class.

 In the case of results from simulations, the file.root should contain the name of the particle ( for instance AnalysisResults.JPSI.root )
 Some considerations  :

   - A configuration file it .txt format is used for more flexibility in the loop process. See AliAnalysisMuMuConfig for details.

   - Associated files ( results from simulation ) can be used. Some methods needs to compare Data vs MC ( for instance AccxEff calculation )

   - By convention, bin type ('pt', 'y', 'integrated' ... ) are minuscule when comming directly from the output file from the online part, and MAJUSCULE when computed from the offline part. See default values.

 
  @author Laurent Aphecetche (Subatech)
  @author Javier Martin Blanco (Subatech)
  @author Benjamin Audurier (Subatech)
*/

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

    /* Fit Proceedure related */

    AliAnalysisMuMuSpectra* FitParticle(
      const char* particle,
      const char* trigger,
      const char* eventType,
      const char* pairCut,
      const char* centrality,
      const AliAnalysisMuMuBinning& binning,
      Bool_t corrected        =kFALSE,
      const TString* fitMethod = 0x0,
      const char* flavour ="",
      const char* histoType ="minv");

    Int_t FitJpsi(
      const char* binType      ="integrated",
      const char* flavour      ="BENJ",
      const TString fitMethod  ="",
      const char* histoType    ="minv");

    void NormMixedMinv(
      const char       * binType="integrated",
      const char       * particle ="psi",
      const char       * flavour="",
      Bool_t corrected =kFALSE,
      Double_t Mmin =2.0,
      Double_t Mmax= 5.0);

    void DivideRawMixHisto(
      const char       * binType="integrated",
      const char       * particle ="psi",
      const char       * flavour="",
      Bool_t corrected =kFALSE,
      Double_t Mmin =2.0,
      Double_t Mmax= 5.0);


    /* Print/draw/plot related */

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

    TObjArray* PrintDistribution(
      const char              * spectraName ="PSI-INTEGRATED",
      const char              * what="NofJPsi",
      const char              * subresultname="",
      Bool_t divideByBinWidth =kTRUE,
      Bool_t AccEffCorr       =kFALSE);

    void PrintFitParam(
      TString spectraName="PT-BENJ",
      const char* subresult="CB2VWG_2.4_4.5_SP1.2",
      const char* param = "FitChi2PerNDF,mJPsi,sJPsi"
    )const;

    void PlotJpsiYield(
      const char* spectraName="INTEGRATED",
      const char* beamYear="",
      const char* subresultname="",
      int NofMUL   =0 ,
      const char* externfile1="externFile_PT.txt",
      const char* externfile2="externFile_CENT.txt");

    void RAAasGraphic(
      const char        * particle="PSI",
      const char        * binType="PT",
      const char        * externfile1="externFile_PT.txt",
      const char        * externfile2="externFile_CENT.txt",
      const char        * RefCent ="V0M_00.00_90.00",
      Bool_t print      = kFALSE,
      Bool_t AccEffCorr =kFALSE)const;

    void DrawFitResults(
      const char              * what="NofJPsi",
      const char              * spectraName="PSI-INTEGRATED",
      const char              * subresults="",
      Bool_t              mix =kFALSE,
      Bool_t AccEffCorr       =kFALSE)const;

    void PrintNofWhat(
      const char        * what="NofJPsi",
      const char        * spectraName="PSI-INTEGRATED",
      Bool_t mix        =kFALSE,
      Bool_t AccEffCorr =kFALSE) const;

    TGraph* PlotEventSelectionEvolution(
      const char         * trigger1="CINT7-B-NOPF-MUON",
      const char         * event1="PSALL",
      const char         * trigger2="CINT7-B-NOPF-MUON",
      const char         * event2="PSALL",
      Bool_t drawFills   =kFALSE,
      Bool_t asRejection =kTRUE) const;

    /* Computing results related */

    void ComputeDimuonRawCount(
      const Double_t rlow   = 2.8,
      const Double_t rhight = 3.4,
      const char            * binType="pt",
      const char            * binRangeExluded="PT_BENJ_00.00_00.30,PT_BENJ_01.00_08.00,PT_BENJ_00.30_01.00",
      const char            * flavour="BENJ",
      Bool_t corrected      =kFALSE );

    void ComputePPCrossSection(
      const char        * spectraName="PSI-PT",
      const char        * externfile = "",
      const char        * externfile2 = "",
      Bool_t print      =kFALSE,
      const char        * what ="CorrNofJPsi");

    TH2* ComputeSPDCorrection(
      const char       * type="oneOverAccEff",
      const char       * eventSel="PSALL",
      const char       * triggerSel="ANY",
      Bool_t bkgReject =kTRUE);

    void ComputeNumberOfEvent();

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

    void ComputeMBXSectionFractionInBins(
      const char* filePileUpCorr="",
      const char* eventSelection="PSALLHASSPDSPDZQA_RES0.25_ZDIF0.50SPDABSZLT10.00",
      const char* what="psi",
      const char* quantity="ntrcorr",
      const char* flavour="D2H");



    /* Cleaning related */

    void CleanAllSpectra();

    void CleanFNorm();

    void CleanMix();

    /* Asoociated File 1 related */

    AliAnalysisMuMuSpectra* ComputeYield(const char* binType, const char* flavour="",const char* accEffSubResultName="PSICB2");

    AliAnalysisMuMuSpectra* CorrectSpectra(
      const char* binType,
      const char* flavour="",
      const char* accEffSubResultName="PSICOUNT");

    TH1* PlotAccEfficiency(const char* SpectraName="PSI-INTEGRATED", const char* subresultname = "PSICOUNT");

    /* internal get */

    Bool_t GetCollections(
      const char            * rootfile,
      const char            * subdir,
      AliMergeableCollection*& oc,
      AliCounterCollection  *& cc,
      AliAnalysisMuMuBinning*& bin,
      std                   ::set<int>& runnumbers);

    AliAnalysisMuMuSpectra* GetSpectraFromConfig(
      const char* binType,
      const char* flavour="") const;

    AliAnalysisMuMuSpectra* GetSpectra(
      const char* binType ,
      const char* EventSelection ="ALL" ,
      const char* DimuonTrigger="ANY",
      const char* Centrality="V0A",
      const char* PairSelectionKey="pALLPAIRYPAIRPTIN0.0-12.0RABSMATCHLOWETAPDCA",
      const char* flavour="BENJ") const;

    UInt_t GetSum(
      AliCounterCollection& cc,
      const char* triggerList,
      const char* eventSelection,
      Int_t runNumber=-1);

    const char* GetParticleName() { return fParticleName; }

    ULong64_t GetTriggerScalerCount(const char* triggerList, Int_t runNumber);

    const std::set<int>& RunNumbers() const { return fRunNumbers; }

    AliAnalysisMuMuSpectra* SPECTRA(const char* fullpath) const;

    AliMergeableCollection* OC()     const { return fMergeableCollection; }
    AliCounterCollection  * CC()     const { return fCounterCollection; }
    AliAnalysisMuMuBinning* BIN()    const { return fBinning; }
    AliAnalysisMuMu       * SIM()    const { return fAssociatedSimulation; }
    AliAnalysisMuMu       * SIM2()   const { return fAssociatedSimulation2; }
    AliAnalysisMuMuConfig * Config() const { return fConfig; }

    /* Other */

    Double_t ErrorPropagationAxBoverCxD(Double_t a,Double_t b,Double_t c,Double_t d);
    void TwikiOutputFnorm(const char* series="FnormOffline2PUPS,FnormScalersPUPS,FnormBest2,RelDifFnormScalersPUPSvsFnormOffline2PUPS,FnormScalersPUVDM,RelDifFnormScalersPUPSvsFnormScalersPUVDM") const;

    void Update();
    Bool_t Upgrade();
    Bool_t Upgrade(const char* filename);

    void SetParticleName(const char* particleName) { fParticleName = particleName; }
    void SetConfig(const AliAnalysisMuMuConfig& config);

    static TFile* FileOpen(const char* file);
    static TString ExpandPathName(const char* file);

    Bool_t IsSimulation() const;

    void Print(Option_t* opt="") const;

    Bool_t SetCorrectionPerRun(const TGraph& corr, const char* formula="");
    void UnsetCorrectionPerRun();

    void ExecuteCanvasEvent(Int_t event, Int_t px, Int_t py, TObject *sel);

    AliAnalysisMuMuSpectra* MergeSpectra(const char* spectraPath1,const char* spectraPath2,const char* spectraname) ;

    void MergeHistoInCentalityBins(const char* binToMerge = "V0M_00.00_10.00,V0M_10.00_20.00", const char* newBin = "V0M_00.00_20.00",Bool_t mix=kFALSE, const char* refbin = "V0M_00.00_10.00" ) const;

private:
    AliAnalysisMuMu(const AliAnalysisMuMu& rhs); // not implemented on purpose
    AliAnalysisMuMu& operator=(const AliAnalysisMuMu& rhs); // not implemented on purpose

    Bool_t SetParticleNameFromFileName(const char* filename);

    void ShowList(const char* title, const TString& list, const char separator=',') const;

    TFile* ReOpen(const char* filename, const char* mode) const;

    TString First(const TString& list) const;

    Bool_t GetParametersFromMC(TString& fitType, const char* pathCentrPairCut, const char* spectraName, AliAnalysisMuMuBinning::Range* bin) const;
    void GetParametersFromResult(TString& fitType, AliAnalysisMuMuJpsiResult* minvResult) const;


    void GetCollectionsFromAnySubdir(TDirectory& dir,
                                    AliMergeableCollection*& oc,
                                     AliCounterCollection*& cc,
                                     AliAnalysisMuMuBinning*& bin);

    void GetFileNameAndDirectory(const char* filename);

private:

    void SetNofInputParticles(AliAnalysisMuMuJpsiResult& r,const char* event, const char* trigger, const char* centrality);


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
