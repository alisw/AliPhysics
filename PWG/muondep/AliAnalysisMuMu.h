#ifndef ALIANALYSISMUMU_H
#define ALIANALYSISMUMU_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

///
/// AliAnalysisMuMu : helper class to digest/plot/massage results from
/// AliAnalysisTaskMuMu
///
/// author : Laurent Aphecetche (Subatech)

#include "AliAnalysisMuMuBinning.h"
#include "TNamed.h"
#include <map>
#include <set>
#include <string>
#include <TString.h>
#include <vector>
#include "RQ_OBJECT.h"

class AliAnalysisMuMuResult;
class AliAnalysisMuMuSpectra;
class AliCounterCollection;
class AliMergeableCollection;
class TF1;
class TFile;
class TGraph;
class TH1;
class TMap;

class AliAnalysisMuMu : public TObject, public TQObject
{

public:
  
  enum EColor
  {
    kBlue=1500,
    kOrange=1501,
    kGreen=1502
  };
  
  AliAnalysisMuMu(const char* filename="LHC12c_muons_AOD000_000179687.saf.root",
                  const char* associatedSimFileName="");
  
  virtual ~AliAnalysisMuMu();
  
  /* Basic checks */
  void BasicCounts(Bool_t detailTrigger=kFALSE,
                   ULong64_t* totalNmb=0x0,
                   ULong64_t* totalNmsl=0x0,
                   ULong64_t* totalNmul=0x0);
  
  static void BasicCountsEvolution(const char* filelist, Bool_t detailTrigger=kFALSE);
  
  void TriggerCountCoverage(const char* triggerList, Bool_t compact=kTRUE,
                            Bool_t orderByTriggerCount=kTRUE);
  
  
  /** Background evolution functions */
  static TObjArray* ComputeBackgroundEvolution(const char* filelist, const char* triggerList,
                                               Double_t ptmin=0.0,
                                               const char* outputfile="background-evolution.root", const char* outputMode="UPDATE");
  
  static void PlotBackgroundEvolution(const char* gfile, const char* triggerList, Double_t ymax=100.0, Bool_t fillBoundaries=kFALSE);
  
  
  /** Fitting */
  static TMap* ComputeJpsiEvolution(const char* filelist, const char* triggerList,
                                    const char* outputFile="jpsi-evolution.root");
  
  static void PlotJpsiEvolution(const char* resultFile, const char* triggerList, Bool_t fillBoundaries=kFALSE,
                                const char* efficiencyFile=0x0, Bool_t simulation=kFALSE);
  
  
  AliAnalysisMuMuSpectra* FitParticle(const char* particle,
                                      const char* trigger,
                                      const char* eventType,
                                      const char* pairCut,
                                      const char* centrality,
                                      const AliAnalysisMuMuBinning& binning);

  AliAnalysisMuMuSpectra* CorrectSpectra(const char* type, const char* flavour="");

  AliAnalysisMuMuSpectra* ComputeYield(const char* type, const char* flavour="");
  
  void CleanAllSpectra();
  
  ///------
  
//  static AliAnalysisMuMuSpectra* ComputeYield(const char* realFile="ds.list.saf.root",
//                                              const char* simFile="ds.sim.list.saf.root",
//                                              const  char* type="PSI-Y VS PT");

  static AliAnalysisMuMuSpectra* RABy(const char* realFile="ds.list.saf.root", const char* simFile="ds.sim.list.saf.root", const char* type="", const char* direction="pPb");

  static TGraph* ResultEvolution(const char* runlist, const char* period="LHC13f", const char* what="Y",
                                 Bool_t forceRecomputation=kFALSE);
  
  ///-------
  
  TGraph* PlotEventSelectionEvolution(const char* trigger1="CINT7-B-NOPF-MUON", const char* event1="PSALL",
                                   const char* trigger2="CINT7-B-NOPF-MUON", const char* event2="PSALL",
                                      Bool_t drawFills=kFALSE,
                                      Bool_t asRejection=kTRUE) const;

  Bool_t Upgrade();
  
  static Bool_t Upgrade(const char* filename);
  
  static void CentralityCheck(const char* filelist);
  
  static TObjArray* CompareJpsiPerCMUUWithBackground(const char* jpsiresults="results.root",
                                                     const char* backgroundresults="background.lhc11d.root");
  
  static TGraph* CompareJpsiPerCMUUWithSimu(const char* realjpsiresults="results.root",
                                            const char* simjpsiresults="results.sim.root");
  
  static Bool_t DecodeFileName(const char* filename, TString& period, int& esdpass, int& aodtrain, int& runnumber);
  
  
  static TFile* FileOpen(const char* file);
  
  
  static Bool_t GetCollections(const char* rootfile,
                               AliMergeableCollection*& mc,
                               AliCounterCollection*& cc,
                               AliAnalysisMuMuBinning*& bin,
                               std::set<int>& runnumbers);
  
  AliAnalysisMuMuSpectra* GetSpectra(const char* what, const char* flavour="") const;

  static UInt_t GetSum(AliCounterCollection& cc, const char* triggerList, const char* eventSelection, Int_t runNumber=-1);
  
  static ULong64_t GetTriggerScalerCount(const char* triggerList, Int_t runNumber);
  
  Int_t Jpsi(const char* what="integrated", const char* binningFlavour="");
  
  Bool_t IsSimulation() const;
  
  static TObjArray* ReadFileList(const char* filelist);
  
  static Int_t RunNumberFromFileName(const char* filename);
  
  TString DimuonTriggerList() const { return fDimuonTriggers; }
  void SetDimuonTriggerList(const char* dimuonTriggerList) { fDimuonTriggers = dimuonTriggerList; }

  TString MuonTriggerList() const { return fMuonTriggers; }
  void SetMuonTriggerList(const char* muonTriggerList) { fMuonTriggers = muonTriggerList; }

  TString MinbiasTriggerList() const { return fMinbiasTriggers; }
  void SetMinbiasTriggerList(const char* minbiasTriggerList) { fMinbiasTriggers = minbiasTriggerList; }
  
  TString EventSelectionList() const { return fEventSelectionList; }
  void SetEventSelectionList(const char* eventSelectionList) { fEventSelectionList = eventSelectionList; }
  
  TString PairSelectionList() const { return fPairSelectionList; }
  void SetPairSelectionList(const char* pairSelectionList) { fPairSelectionList = pairSelectionList; }

  TString CentralitySelectionList() const { return fCentralitySelectionList; }
  void SetCentralitySelectionList(const char* centralitySelectionList) { fCentralitySelectionList = centralitySelectionList; }

  TString FitTypeList() const { return fFitTypeList; }
  void SetFitTypeList(const char* fitTypelist) { fFitTypeList = fitTypelist; }
  
  static void SetDefaultDimuonTriggerList(const char* dimuonTriggerList) { fgDefaultDimuonTriggers = dimuonTriggerList; }
  static void SetDefaultMuonTriggerList(const char* muonTriggerList) { fgDefaultMuonTriggers = muonTriggerList; }
  static void SetDefaultMinbiasTriggerList(const char* minbiasTriggerList) { fgDefaultMinbiasTriggers = minbiasTriggerList; }
  static void SetDefaultEventSelectionList(const char* eventSelectionList) { fgDefaultEventSelectionList = eventSelectionList; }
  static void SetDefaultPairSelectionList(const char* pairSelectionList) { fgDefaultPairSelectionList = pairSelectionList; }
  static void SetDefaultCentralitySelectionList(const char* centralitySelectionList) { fgDefaultCentralitySelectionList = centralitySelectionList; }
  static void SetDefaultFitTypes(const char* fitTypes) { fgDefaultFitTypeList = fitTypes; }
  
  static void SetDefaultEventSelectionForSimulations(const char* value) { fgDefaultEventSelectionForSimulations = value; }
  static void SetDefaultDimuonTriggerForSimulations(const char* value) { fgDefaultDimuonTriggerForSimulations = value; }

  AliMergeableCollection* MC() const { return fMergeableCollection; }
  AliCounterCollection* CC() const { return fCounterCollection; }
  AliAnalysisMuMuBinning* BIN() const { return fBinning; }
  
  static void SetOCDBPath(const char* ocdbPath) { fgOCDBPath = ocdbPath; }
  
  static void SetColorScheme();
  
  static void SetCompactGraphs(Bool_t value=kTRUE) { fgIsCompactGraphs = value; }
  
  static Bool_t CompactGraphs() { return fgIsCompactGraphs; }
  
  static void Compact(TGraph& g);
  
  void Print(Option_t* opt="") const;
  
  void GetMBR(Int_t runNumber, const char* eventSelection, Double_t& value, Double_t& error) const;

  TGraph* MBREvolution(const char* eventSelection1="PSALLNOTZEROPILEUP", const char* eventSelection2="PSALL") const;

  const std::set<int>& RunNumbers() const { return fRunNumbers; }
  
  void DrawMinv(const char* type,
                const char* particle,
                const char* trigger,
                const char* eventType,
                const char* pairCut,
                const char* centrality,
                const char* subresultname="",
                const char* flavour="") const;

  void DrawMinv(const char* type="PT", const char* particle="PSI", const char* flavour="", const char* subresultname="") const;
  
  Bool_t SetCorrectionPerRun(const TGraph& corr, const char* formula="");
  
  void UnsetCorrectionPerRun();
  
  void ExecuteCanvasEvent(Int_t event, Int_t px, Int_t py, TObject *sel);

  std::vector<Double_t> GetMCCB2Tails(const AliAnalysisMuMuBinning::Range& bin) const;
  
  AliAnalysisMuMu* SIM() const { return fAssociatedSimulation; }
  
  AliAnalysisMuMuSpectra* SPECTRA(const char* fullpath) const;
  
  void Update();

private:
  AliAnalysisMuMu(const AliAnalysisMuMu& rhs); // not implemented on purpose
  AliAnalysisMuMu& operator=(const AliAnalysisMuMu& rhs); // not implemented on purpose
  
  TFile* ReOpen(const char* filename, const char* mode) const;

  TString First(const TString& list) const;
  
private:

  void SetNofInputParticles(AliAnalysisMuMuResult& r);

  static TString ExpandPathName(const char* file);
  
  
  TString fFilename; // file containing the result collections (of objects and counters) from AliAnalysisTaskMuMu
  AliCounterCollection* fCounterCollection; // collection of counters in file
  TString fDimuonTriggers; // list of dimuon triggers to consider
  TString fMuonTriggers; // list of single muon triggers to consider
  TString fMinbiasTriggers;   // list of minbias triggers to consider
  TString fEventSelectionList; // list of event types to consider
  TString fPairSelectionList; // list of pair cuts to consider
  TString fCentralitySelectionList; // list of centrality cuts to consider
  TString fFitTypeList; // list of fit types to perform

  AliAnalysisMuMuBinning* fBinning; // binning
  
  static TString fgOCDBPath; // OCDB to be used (raw:// by default)
  
  static TString fgDefaultMuonTriggers; // default list of single muon triggers
  static TString fgDefaultMinbiasTriggers; // default list of MB triggers
  static TString fgDefaultDimuonTriggers; // default list of dimuon triggers
  static TString fgDefaultEventSelectionList; // default list of event selections
  static TString fgDefaultPairSelectionList; // default list of pair selections
  static TString fgDefaultCentralitySelectionList; // default list of centrality selections
  static TString fgDefaultFitTypeList; // default list of fit types
  static TString fgDefaultEventSelectionForSimulations; // default event selection (simulations)
  static TString fgDefaultDimuonTriggerForSimulations; // default dimuon trigger (simulations)
  
  static Bool_t fgIsCompactGraphs; // whether the graph produced should be compact
  
  AliMergeableCollection* fMergeableCollection; // collection of objects in file

  std::set<int> fRunNumbers; // run numbers
  
  TGraph* fCorrectionPerRun; // correction factor per run
  
  AliAnalysisMuMu* fAssociatedSimulation; // associated simulations (if any)

  ClassDef(AliAnalysisMuMu,11) // class to analysis results from AliAnalysisTaskMuMuXXX tasks
};

#endif
