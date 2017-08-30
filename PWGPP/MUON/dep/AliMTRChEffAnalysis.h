#ifndef ALIMTRCHEFFANALYSIS_h
#define ALIMTRCHEFFANALYSIS_h

/// \class AliMTRChEffAnalysis
/// \brief Output for Trig chamber effieincy
///
/// The class manipulates the output of AliAnalysisTaskTrigChEff
/// in order to build the trigger chamber efficiency object
/// to be plugged in the OCDB for simulations
///
/// \author Diego Stocco <dstocco@cern.ch>, Subatech
/// \date Nov 8, 2015

#include "TObject.h"
#include "TString.h"
#include "AliTrigChEffOutput.h"
#include <map>
#include <string>
#include <vector>

class TObjArray;
class TDirectory;
class TH1;
class TGraphAsymmErrors;
class TArrayI;
class TArrayD;
class TList;

class AliMTRChEffAnalysis : public TObject {
 public:
  AliMTRChEffAnalysis();
  AliMTRChEffAnalysis ( const char *localFileList, const char *outputName = "testMTRChamberEff" );

  TArrayI GetHomogeneousRanges ( Double_t chi2Cut = 3, Int_t maxNRanges = 4, Double_t minEffVariation = 0.005, Bool_t perRPC = kTRUE, TArrayI* forcedChanges = 0x0, Double_t minEff = 0.85, Double_t maxEff = 1.01 );
  TArrayI GetHomogeneousRanges ( TGraphAsymmErrors* trendGraph, Double_t chi2Cut = 3, Int_t maxNRanges = 4, Double_t minEffVariation = 0.005, TArrayI* forcedChanges = 0x0, Bool_t returnIndex = kFALSE );

  void DrawEffTrend ( Int_t itype, Int_t irpc, Double_t maxNsigmasOutliers = -1., Double_t minEff = 0.8, Double_t maxEff = 1.01 ) const;
  void DrawStatContribution ( Int_t itype, Int_t irpc, Double_t maxNsigmaOutliers = -1., Double_t minY = 0., Double_t maxY = 0.15 ) const;

  Double_t GetAverageStat ( Int_t firstRun, Int_t lastRun, Int_t itype = AliTrigChEffOutput::kHboardEff, Bool_t excludePeriphericBoards = kTRUE ) const;

  TGraphAsymmErrors* GetOutliers ( TGraphAsymmErrors* graph, Double_t maxNsigmas = 3. ) const;

  TH1* GetTrend ( Int_t itype, Int_t icount, Int_t ichamber, Int_t idetelem ) const;
  TGraphAsymmErrors* GetTrendEff ( Int_t itype, Int_t icount, Int_t ichamber, Int_t idetelem ) const;

  Int_t CompareEfficiencies ( const char* sources, const char* titles, const char* opt, const char* canvasNameSuffix = "" ) const;
  Int_t CompareEfficiencyMethods ( const char* source, const char* opt, const char* canvasNameSuffix = "" ) const;
  void CompareMergedEfficiencies ( const char* opt ) const;
  Int_t ComputeAndCompareEfficiencies ( const char* sources, const char* titles, const char* opt, const char* canvasNameSuffix = "") const;

  Bool_t AddSystematicCondition ( const char* physSel, const char* trigClassName, const char* centrality, Int_t itrackSel, Int_t imatch, Int_t imethod );
  Bool_t SetDefaultEffConditions ();
  Bool_t SetEffConditions ( const char* physSel, const char* trigClassName, const char* centrality, Int_t itrackSel, Int_t imatch, Int_t imethod );

  Bool_t MergeOutput ( TArrayI runRanges, Double_t averageStatError = 0.01, Bool_t isIndex = kFALSE );

  Bool_t InitFromGrid ( const char *runList, const char *path, const char *pattern, const char* localFileList = "localFileList.txt", const char* outDir = "", const char *directory = "MTR_ChamberEffMap", const char *outputName = "testMTRChamberEff" );
  Bool_t InitFromLocal ( const char *localFileList, const char *outputName = "testMTRChamberEff" );
  Bool_t InitFromWeb ( const char *runList, const char *path, const char* localFileList = "localFileList.txt", const char* outDir = "", const char *directory = "MTR_ChamberEffMap", const char *outputName = "testMTRChamberEff" );


  Bool_t WriteMergedToOCDB ( const char* outputCDB = "CDB", Bool_t writeSystematics = kFALSE ) const;
  Bool_t DrawSystematicEnvelope ( Bool_t perRPC = kFALSE ) const;
  Bool_t BuildSystematicMap ();
  Bool_t RecoverEfficiency ( const char* runList, const char* ocdb, const char* systOcdb, Int_t referenceRun = -1 );

  static void ZoomPad();

  virtual ~AliMTRChEffAnalysis();

 private:

  Bool_t AddToList ( const char *filename, const char *outputName );
  TArrayI BoardsInRPC ( Int_t irpc ) const;
  void CopyDir ( TDirectory *source ) const;
  Bool_t CopyLocally ( const char* runList, const char* path, const char* pattern, const char* localFileList, const char* outDir, const char* directory ) const;
  Int_t Check() const;
  Int_t  CompareEfficiencies ( TObjArray* effHistoLists, const char* titles, const char* opt, const char* canvasNameSuffix ) const;
  TList* CloneEffHistoList ( TList* effHistos ) const;
  Bool_t ExecCommand ( TString command, Bool_t prompt ) const;
  Double_t FitRangesFunc ( Double_t* x, Double_t* par );
  Double_t GetError ( Double_t errLow, Double_t errHigh ) const;
  TList* GetEffHistoList ( AliTrigChEffOutput* trigOut, TObjArray* condition ) const;
  TH1* GetHisto ( TList* effHistoList, Int_t itype, Int_t icount, Int_t ichamber ) const;
  Int_t GetIndexFromRun ( Int_t runNumber ) const;
  Int_t GetRunNumber ( Int_t ipt ) const;
  TList* GetRunList ( const char* runList ) const;
  TString GetId ( const char* condition, Int_t minRun, Int_t maxRun = -1 ) const;

  TString GetShortConditionTitle ( const char* conditionName ) const;
  AliTrigChEffOutput* Namer() const;

  TH1* GetSum ( AliTrigChEffOutput* trigOut, TObjArray* condition, Int_t itype, Int_t icount, Int_t ichamber ) const;

  Double_t GetThreeOfFour ( TArrayD eff, TArrayD effErr, Double_t &probErr ) const;


  Bool_t HasMergedResults () const;

  TArrayI MergeRangesForStat ( TArrayI runRanges, Double_t averageStatError, Bool_t excludePeriphericBoards = kTRUE ) const;

  TList* ReadEffHistoList ( const char* src ) const;

  Bool_t SetCondition ( const char* physSel, const char* trigClassName, const char* centrality, Int_t itrackSel, Int_t imatch, Int_t imethod, Bool_t isBasic );

  /// Dummy
  AliMTRChEffAnalysis(const AliMTRChEffAnalysis&);
  /// Dummy
  AliMTRChEffAnalysis& operator=(const AliMTRChEffAnalysis&);

  TObjArray* fConditions; //!<! List of conditions for trigger efficiency

  mutable AliTrigChEffOutput* fNamer; //!<! Namer for histograms

  class AliMTRChEffInnerObj : public TObject {
  public:
    AliMTRChEffInnerObj ( const char* filename, const char* outputname, Int_t minRun, Int_t maxRun = -1 );
    virtual ~AliMTRChEffInnerObj ();
    TString GetFilename () const { return fFilename; }
    TString GetOutputname () const { return fOutputname; }
    Int_t GetMinRun () const { return fMinRun; }
    Int_t GetMaxRun () const { return fMaxRun; }
    std::map<std::string,TList*> GetEffLists() { return fEffLists; }
    std::vector<std::string> GetSortKeys() { return fSortKeys; }
    TList* GetEffHistoList ( const char* condition ) const;
    Bool_t AddEffHistoList ( const char* condition, TList* effHistoList );
    Bool_t RemoveEffHistoList ( const char* condition );

  private:
    TString fFilename;
    TString fOutputname;
    Int_t fMinRun;
    Int_t fMaxRun;
    std::map<std::string,TList*> fEffLists;
    std::vector<std::string> fSortKeys;
  };

  std::vector<AliMTRChEffAnalysis::AliMTRChEffInnerObj*> fRunMap; //!<! Map of internal objects per run
  std::vector<AliMTRChEffAnalysis::AliMTRChEffInnerObj*> fMergedMap; //!<! Map of merged internal objects

  /// \cond CLASSIMP
  ClassDef(AliMTRChEffAnalysis, 0); // Trigger chamber efficiencies
  /// \endcond
};

#endif
