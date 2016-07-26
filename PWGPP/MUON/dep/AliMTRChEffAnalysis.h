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

class TObjArray;
class TDirectory;
class TH1;
class TGraphAsymmErrors;
class TArrayI;
class TArrayD;

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

  Bool_t AddSystematicCondition ( const char* physSel, const char* trigClassName, const char* centrality, Int_t itrackSel, Int_t imatch, Int_t imethod );
  Bool_t SetDefaultEffConditions ();
  Bool_t SetEffConditions ( const char* physSel, const char* trigClassName, const char* centrality, Int_t itrackSel, Int_t imatch, Int_t imethod );

  Bool_t MergeOutput ( TArrayI runRanges, Double_t averageStatError = 0.01, Bool_t isIndex = kFALSE );

  Bool_t SetResultsFromGrid ( const char *runList, const char *path, const char *pattern, const char* localFileList = "localFileList.txt", const char* outDir = "", const char *directory = "MTR_ChamberEffMap", const char *outputName = "testMTRChamberEff" );
  Bool_t SetResultsFromWeb ( const char *runList, const char *path, const char* localFileList = "localFileList.txt", const char* outDir = "", const char *directory = "MTR_ChamberEffMap", const char *outputName = "testMTRChamberEff" );

  Bool_t WriteMergedToOCDB ( const char* outputCDB = "CDB" ) const;
  Bool_t DrawSystematicEnvelope ( Bool_t perRPC = kFALSE ) const;

  virtual ~AliMTRChEffAnalysis();

 private:

  Bool_t AddToList ( const char *filename, const char *outputName );
  TArrayI BoardsInRPC ( Int_t irpc ) const;
  void CopyDir ( TDirectory *source ) const;
  Bool_t CopyLocally ( const char* runList, const char* path, const char* pattern, const char* localFileList, const char* outDir, const char* directory ) const;
  Int_t  CompareEfficiencies ( TObjArray* effMapList, const char* titles, const char* opt, const char* canvasNameSuffix ) const;
  Bool_t ExecCommand ( TString command, Bool_t prompt ) const;
  Double_t FitRangesFunc ( Double_t* x, Double_t* par );
  Double_t GetError ( Double_t errLow, Double_t errHigh ) const;
  TList* GetEffHistoList ( AliTrigChEffOutput* trigOut, TObjArray* condition ) const;
  TH1* GetHisto ( TList* effHistoList, Int_t itype, Int_t icount, Int_t ichamber ) const;
  TString GetIdentifier ( AliTrigChEffOutput* trigOut, TObjArray* condition, Int_t itype, Int_t icount, Int_t ichamber ) const;
  Int_t GetIndexFromRun ( UInt_t runNumber ) const;
  Int_t GetRunNumber ( Int_t ipt ) const;
  TList* GetRunList ( const char* runList ) const;

  Bool_t GetShortConditionTitles ( AliTrigChEffOutput* trigOut, TObjArray& condTitles ) const;

  TH1* GetSum ( AliTrigChEffOutput* trigOut, TObjArray* condition, Int_t itype, Int_t icount, Int_t ichamber ) const;

  Double_t GetThreeOfFour ( TArrayD eff, TArrayD effErr, Double_t &probErr ) const;


  Bool_t HasMergedResults () const;

  TArrayI MergeRangesForStat ( TArrayI runRanges, Double_t averageStatError, Bool_t excludePeriphericBoards = kTRUE ) const;

  Bool_t SetCondition ( const char* physSel, const char* trigClassName, const char* centrality, Int_t itrackSel, Int_t imatch, Int_t imethod, Bool_t isBasic );

  Bool_t SetOutList ( const char *localFileList, const char *outputName );

  /// Dummy
  AliMTRChEffAnalysis(const AliMTRChEffAnalysis&);
  /// Dummy
  AliMTRChEffAnalysis& operator=(const AliMTRChEffAnalysis&);

  TObjArray* fOutputs; ///!<! List of outputs
  TObjArray* fConditions; ///!<! List of conditions for trigger efficiency
  TObjArray* fMergedOutputs; ///!<! List of merged outputs

  /// \cond CLASSIMP
  ClassDef(AliMTRChEffAnalysis, 0); // Trigger chamber efficiencies
  /// \endcond
};

#endif
