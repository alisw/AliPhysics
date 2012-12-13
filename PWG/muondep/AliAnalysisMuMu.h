#ifndef ALIANALYSISMUMU_H
#define ALIANALYSISMUMU_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/// 
/// AliAnalysisMuMu : helper class to digest results from
/// AliAnalysisTaskMuMu(fromAOD or fromESD)
///
/// author : Laurent Aphecetche (Subatech)

#include "TNamed.h"
#include <string>
#include <TString.h>

class TH1;
class AliHistogramCollection;
class AliCounterCollection;
class TF1;
class TGraph;
class TMap;
class TFile;

class AliAnalysisMuMu : public TObject
{
  
public:
  
  enum EColor {
    kBlue=1500,
    kOrange=1501,
    kGreen=1502
  };
  
  AliAnalysisMuMu(const char* filename="LHC12c_muons_AOD000_000179687.saf.root");
  virtual ~AliAnalysisMuMu();
  
  class Result : public TNamed
  {
    
  public:
    
    enum EFitType
    {
      kJpsi=(1<<0),
      kJpsiPrime=(1<<1),
      kUpsilon=(1<<2),
      kPbPb2010=(1<<3),
      kPbPb2011=(1<<4),
      kpp2011=(1<<5),
      kMatchAny=(1<<6)
    };
    
    Result(TRootIOCtor* /*io*/) :
    TNamed("",""),
    fTriggerName(),
    fEventSelection(),
    fPairSelection(),
    fCentralitySelection(),
    fNofRuns(),
    fNofTriggers(),
    fHC(0x0), fRebin(0),
    fFitTotal(0x0),
    fMinv(0x0), fMinvLS(0), 
    fFitType(0), fMap(0x0)
    {
    }
    
    Result(const char* triggerClass,
           const char* eventSelection,
           const char* pairSelection,
           const char* centSelection,
           Int_t ntriggers, AliHistogramCollection* hc, Int_t nrebin=1, UInt_t fitType=(kJpsi | kJpsiPrime)) :
    TNamed("",""), fTriggerName(triggerClass),
    fEventSelection(eventSelection),
    fPairSelection(pairSelection),
    fCentralitySelection(centSelection),
    fNofRuns(1),
    fNofTriggers(ntriggers), fHC(hc), fRebin(nrebin),
    fFitTotal(0x0),
    fMinv(0x0), fMinvLS(0), 
    fFitType(fitType), fMap(0x0)
    {
      if (hc) Fit(nrebin);
    }
  
    virtual ~Result();
  
    const char* GetName() const { return fTriggerName.Data(); }
    
    AliHistogramCollection* HC() const { return fHC; }
    
    const char* TriggerClass() const { return GetName(); }
    
    const char* EventSelection() const { return fEventSelection.Data(); }
    const char* PairSelection() const { return fPairSelection.Data(); }
    const char* CentralitySelection() const { return fCentralitySelection.Data(); }
    
    Int_t NofTriggers() const { return fNofTriggers; }
    
    Double_t NofJpsi() const { return GetValue("NofJpsi"); }
    
    Double_t ErrorOnNofJpsi() const { return GetError("NofJpsi"); }
    
    Double_t MeanJpsi() const { return GetValue("MeanJpsi"); }
    
    Double_t ErrorOnMeanJpsi() const { return GetError("MeanJpsi"); }
    
    Double_t SigmaJpsi() const { return GetValue("SigmaJpsi"); }
    
    Double_t ErrorOnSigmaJpsi() const { return GetError("SigmaJpsi"); }
    
    Double_t NofJpsiPrime() const { return GetValue("NofJpsiPrime"); }
    
    Double_t ErrorOnNofJpsiPrime() const { return GetError("NofJpsiPrime"); }
    
    Double_t MeanJpsiPrime() const { return GetValue("MeanJpsiPrime"); }
    
    Double_t ErrorOnMeanJpsiPrime() const { return GetError("MeanJpsiPrime"); }
    
    Double_t SigmaJpsiPrime() const { return GetValue("SigmaJpsiPrime"); }
    
    Double_t ErrorOnSigmaJpsiPrime() const { return GetError("SigmaJpsiPrime"); }
    
    Double_t NofUpsilon() const { return GetValue("NofUpsilon"); }
    
    Double_t ErrorOnNofUpsilon() const { return GetError("NofUpsilon"); }
    
    Double_t MeanUpsilon() const { return GetValue("MeanUpsilon"); }
    
    Double_t ErrorOnMeanUpsilon() const { return GetError("MeanUpsilon"); }
    
    Double_t SigmaUpsilon() const { return GetValue("SigmaUpsilon"); }
    
    Double_t ErrorOnSigmaUpsilon() const { return GetError("SigmaUpsilon"); }
    
    void Set(const char* name, double value, double error);
    
    void Print(Option_t* opt="") const;
    
    void Fit(Int_t rebin=1);
    
    void FitJpsiJpsiPrimeCB(TH1& h);
    
    void FitJpsiJpsiPrimeCustom(TH1& h);
    
    void FitJpsi(TH1& h);
    
    void FitJpsiCBE(TH1& h);
    
    void FitJpsiECBE(TH1& h);
    
    void FitJpsiPCBE(TH1& h);
    
    void FitUpsilon(TH1& h);
    
    void FitJpsiGCBE(TH1& h);
    
    Bool_t HasValue(const char* name) const;
    Double_t GetValue(const char* name) const;
    Double_t GetError(const char* name) const;
    
    TH1* Minv() const { return fMinv; }
    TH1* MinvLS() const { return fMinvLS; }
    
    UInt_t FitType() const { return fFitType; }
    
    Int_t NofRuns() const { return fNofRuns; }
    
    void SetNofRuns(int n) { fNofRuns=n; }
    
  private:
    Result(const Result&);
    Result& operator=(const Result&);
    TMap* Map() const;
    
  private:
    TString fTriggerName;
    TString fEventSelection;
    TString fPairSelection;
    TString fCentralitySelection;
    Int_t fNofRuns;
    Int_t fNofTriggers;
    AliHistogramCollection* fHC;
    Int_t fRebin;
    TF1* fFitTotal;
    TH1* fMinv;
    TH1* fMinvLS;
    UInt_t fFitType;
    
    mutable TMap* fMap; // map of results (string,TObjArray) each TObjArray = (value,error) = (TParameter<double>, TParameter<double>)
    
    ClassDef(AliAnalysisMuMu::Result,7)
  };

  /* Basic checks */
  void BasicCounts(Bool_t detailTrigger=kFALSE,
                   ULong64_t* totalNmb=0x0,
                   ULong64_t* totalNmsl=0x0,
                   ULong64_t* totalNmul=0x0);

  static void BasicCountsEvolution(const char* filelist, Bool_t detailTrigger=kFALSE);
  
  void TriggerCountCoverage(const char* triggerList, Bool_t compact=kTRUE);


  /** Background evolution functions */
  static TObjArray* ComputeBackgroundEvolution(const char* filelist, const char* triggerList,
                                               const char* outputfile="background-evolution.root", const char* outputMode="UPDATE");

  static void PlotBackgroundEvolution(const char* gfile, const char* triggerList);
  

  /** Fitting */
  static TMap* ComputeJpsiEvolution(const char* filelist, const char* triggerList,
                                         const char* outputFile="jpsi-evolution.root", const char* outputMode="UPDATE",
                                         Bool_t simulation=kFALSE);
  
  static void PlotJpsiEvolution(const char* resultFile, const char* triggerList, Bool_t fillBoundaries=kFALSE,
                                const char* efficiencyFile=0x0);

  ///------
  
  static void CentralityCheck(const char* filelist);

  static TObjArray* CompareJpsiPerCMUUWithBackground(const char* jpsiresults="results.root",
                                                     const char* backgroundresults="background.lhc11d.root");

  static TGraph* CompareJpsiPerCMUUWithSimu(const char* realjpsiresults="results.root",
                                            const char* simjpsiresults="results.sim.root");

  static Bool_t DecodeFileName(const char* filename, TString& period, int& esdpass, int& aodtrain, int& runnumber);
  

  static TFile* FileOpen(const char* file);
  

  static Bool_t GetCollections(const char* rootfile,
                               AliHistogramCollection*& hc,
                               AliCounterCollection*& cc);
    
  static Result* GetResult(const AliHistogramCollection& hc,
                           AliCounterCollection& cc,
                           const char* base,
                           const char* selection,
                           const char* paircut,
                           const char* centrality,
                           UInt_t fitType,
                           Int_t nrebin=1);
  
  static UInt_t GetSum(AliCounterCollection& cc, const char* triggerList, const char* eventSelection, Int_t runNumber=-1);
  
  static ULong64_t GetTriggerScalerCount(const char* triggerList, Int_t runNumber);

  TObjArray* Jpsi(Bool_t simulation=kFALSE);

  static TObjArray* ReadFileList(const char* filelist);
  
  static Int_t RunNumberFromFileName(const char* filename);
  
  static void SinglePtPlot(const char* rootfile);

  void SetDimuonTriggerList(const char* dimuonTriggerList) { fDimuonTriggers = dimuonTriggerList; }
  void SetMuonTriggerList(const char* muonTriggerList) { fMuonTriggers = muonTriggerList; }
  void SetMinbiasTriggerList(const char* minbiasTriggerList) { fMinbiasTriggers = minbiasTriggerList; }

  TString DimuonTriggerList() const { return fDimuonTriggers; }
  TString MuonTriggerList() const { return fMuonTriggers; }
  TString MinbiasTriggerList() const { return fMinbiasTriggers; }
  
  void SetEventSelectionList(const char* eventSelectionList) { fEventSelectionList = eventSelectionList; }
  void SetPairSelectionList(const char* pairSelectionList) { fPairSelectionList = pairSelectionList; }
  
  static void SetDefaultDimuonTriggerList(const char* dimuonTriggerList) { fgDefaultDimuonTriggers = dimuonTriggerList; }
  static void SetDefaultMuonTriggerList(const char* muonTriggerList) { fgDefaultMuonTriggers = muonTriggerList; }
  static void SetDefaultMinbiasTriggerList(const char* minbiasTriggerList) { fgDefaultMinbiasTriggers = minbiasTriggerList; }
  static void SetDefaultEventSelectionList(const char* eventSelectionList) { fgDefaultEventSelectionList = eventSelectionList; }
  static void SetDefaultPairSelectionList(const char* pairSelectionList) { fgDefaultPairSelectionList = pairSelectionList; }

  
  AliHistogramCollection* HC() const { return fHistogramCollection; }
  AliCounterCollection* CC() const { return fCounterCollection; }
  
  static void SetOCDBPath(const char* ocdbPath) { fgOCDBPath = ocdbPath; }
  
  static void SetColorScheme();

  static void DrawFill(Int_t run1, Int_t run2, double ymin, double ymax, const char* label);

private:
  AliAnalysisMuMu(const AliAnalysisMuMu& rhs); // not implemented on purpose
  AliAnalysisMuMu& operator=(const AliAnalysisMuMu& rhs); // not implemented on purpose
  
private:

  static TString ExpandPathName(const char* file);

  TString fFilename; // file containing the result collections (of histograms and counters) from AliAnalysisTaskMuMu
  AliHistogramCollection* fHistogramCollection; // collection of histograms in file
  AliCounterCollection* fCounterCollection; // collection of counters in file
  TString fDimuonTriggers; // list of dimuon triggers to consider
  TString fMuonTriggers; // list of single muon triggers to consider
  TString fMinbiasTriggers;   // list of minbias triggers to consider
  TString fEventSelectionList; // list of event types to consider
  TString fPairSelectionList; // list of pair cuts to consider

  static TString fgOCDBPath; // OCDB to be used (raw:// by default)
  
  static TString fgDefaultMuonTriggers; // default list of single muon triggers
  static TString fgDefaultMinbiasTriggers; // default list of MB triggers
  static TString fgDefaultDimuonTriggers; // default list of dimuon triggers
  static TString fgDefaultEventSelectionList; // default list of event selections
  static TString fgDefaultPairSelectionList; // default list of pair selections

  ClassDef(AliAnalysisMuMu,3) // class to analysis results from AliAnalysisTaskMuMuXXX tasks
};

#endif
