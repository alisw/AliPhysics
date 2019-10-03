#ifndef ALIANALYSISTRIGGERSCALERS_H
#define ALIANALYSISTRIGGERSCALERS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

#ifndef ROOT_TString
#  include "TString.h"
#endif
#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#include "Riostream.h"

#include <vector>
#include <set>
#include <map>

class AliTriggerBCMask;
class AliAnalysisTriggerScalerItem;
class TGraph;
class AliTriggerConfiguration;
class AliTriggerRunScalers;
class AliLHCData;
class AliTriggerScalersRecordESD;

/** 

  @ingroup pwg_muondep_misc

  @class AliAnalysisTriggerScalers
  
  @brief utility class to play with OCDB scalers.

 Can use it to e.g. :

 - get the integrated luminosity for a given period
   (see IntegratedLuminosity method)

 - plot the trigger rates (see PlotTriggerEvolution)

 Please note that this class is doing an OCDB "scan" (for as many run
 as you put in the SetRunList method), so it can be really long
 if you're using the raw:// OCDB. If you're planning on using this
 class for anything more than a brief test, you'd better think
 or making a local copy of the OCDB objects required by this class :

 GRP/CTP/Config
 GRP/CTP/Scalers

 and also (to get the fill and period ranges, mainly for drawing)

 GRP/GRP/Data
 GRP/GRP/LHCData


 Please note also that this class evolved from a set a quick macros to
 follow the luminosity and trigger rates during data taking to this stage.

 Now (Feb 2013) would probably be a good time to regroup a bit, think about it and
 make it more general, more robust and just less "had-hoc"...

 Ideas for improvement :

 - make sure the "what" parameter mean the same thing in all methods
   and so can take the same values in all methods

 - get the lumi trigger name and cross-section from e.g. OADB instead of
   hard-coding them ?

 - organize the class so that the CTP Scalers and Config can be fetched
   from elsewhere (e.g. from a run based container available in some
   distant future in the ESDs/AODs ?)

 - robustify the PAC fraction computation


 @author Laurent Aphecetche (Subatech)
 */

class AliAnalysisTriggerScalers : public TObject
{
public:
  
  AliAnalysisTriggerScalers(const std::vector<int>& runs, const char* source="");
  AliAnalysisTriggerScalers(const std::set<int>& runs, const char* source="");
  AliAnalysisTriggerScalers(Int_t runNumber, const char* source="");
  AliAnalysisTriggerScalers(const char* runlist, const char* source="");
  
  virtual ~AliAnalysisTriggerScalers();
  
  TString CrossSectionUnit() const { return fCrossSectionUnit; }
  void SetCrossSectionUnit(const char* unit="nb");
  
  void DrawFills(Double_t ymin, Double_t ymax, Int_t color=5);
  void DrawFill(Int_t run1, Int_t run2, double ymin, double ymax, const char* label,Int_t color=5);

  void DrawPeriods(Double_t ymin, Double_t ymax, Int_t color=5);

  void GetCTPObjects(Int_t runNumber, AliTriggerConfiguration*& tc, AliTriggerRunScalers*& trs, AliLHCData*& lhc) const;

  Int_t GetFillNumberFromRunNumber(Int_t runNumber);
  
  void GetFillBoundaries(std::map<int, std::pair<int,int> >& fills);

  TString GetLHCPeriodFromRunNumber(Int_t runNumber) const;
  
  void GetLHCPeriodBoundaries(std::map<std::string, std::pair<int,int> >& periods);

  void GetLuminosityTriggerAndCrossSection(Int_t runNumber,
                                           TString& lumiTriggerClassName,
                                           Double_t& lumiTriggerCrossSection,
                                           Double_t& lumiTriggerCrossSectionError) const;
  
  TObject* GetOCDBObject(const char* path, Int_t runNumber) const;
  
  Float_t GetPauseAndConfigCorrection(Int_t runNumber, const char* triggerClassName);

  void GetPileUpFactor(Int_t runNumber, const char* triggerClassName, Double_t purity, Double_t& value, Double_t& error);
  
  void ShowPileUpFactors(const char* triggerClassName, Double_t purity=1.0);
  
  const std::vector<int>& GetRunList() const { return fRunList; }

  Int_t GetTriggerInput(Int_t runNumber, const char* inputname);
  
  AliAnalysisTriggerScalerItem* GetTriggerScaler(Int_t runNumber, const char* level, const char* triggerClassName);
  
  TGraph* IntegratedLuminosityGraph(const char* triggerName, const char* triggerClassNameForPACEstimate="");
  
  TGraph* IntegratedLuminosityGraph(Int_t runNumber, const char* triggerClassName, const char* triggerClassNameForPACEstimate="");
  
  void IntegratedLuminosity(const char* triggerList="",
                            const char* csvOutputFile="",
                            const char sep='\t'
                            );

  TGraph* MakeGraph(const std::vector<int>& vx, const std::vector<int>& vex,
                    const std::vector<double>& vy, const std::vector<double>& vey);

  static Double_t Mu(Double_t L0B, Double_t Nb);

  Int_t NumberOfInteractingBunches(const AliLHCData& lhc, Int_t runNumber, Bool_t mainSat=kFALSE) const;

  TGraph* PlotTrigger(const char* triggerClassName, const char* what, bool draw=kTRUE);
  
  TGraph* PlotTriggerEvolution(const char* triggerClassName,
                               const char* what,
                               bool draw=kTRUE,
                               double* mean=0x0,
                               bool removeZero=kFALSE,
                               std::vector<int>* runlist=0x0);
  
  TGraph* PlotTriggerRatio(const char* triggerClassName1,
                           const char* what1,
                           const char* triggerClassName2,
                           const char* what2,
                           bool draw=kTRUE,
                           Double_t factor=1.0);
  
  TGraph* PlotTriggerRatioEvolution(const char* triggerClassName1,
                                    const char* what1,
                                    const char* triggerClassName2,
                                    const char* what2,
                                    bool draw=kTRUE,
                                    Double_t factor=1.0);
  
  virtual void Print(Option_t* opt="") const;

  void SetRunList(const std::vector<int>& runlist);
  void SetRunList(const std::set<int>& runlist);
  void SetRunList(Int_t runNumber);
  void SetRunList(const char* runlist);
  
  void ShouldCorrectForPileUp(Bool_t flag) { fShouldCorrectForPileUp = flag; }
  Bool_t ShouldCorrectForPileUp() const { return fShouldCorrectForPileUp; }

  TString GetOCDBPath() { return fOCDBPath; }
  
  void GetTriggerClassesForRun(Int_t runNumber, TString& classes) const;
  
  void PurgeRunList(std::vector<int>& purgedList, const char* triggerClassName, Bool_t verbose=kTRUE);
  
private:

  Bool_t CheckRecord(const AliTriggerScalersRecordESD& record,
                     UInt_t index,
                     ULong64_t refa,
                     ULong64_t refb,
                     UInt_t timelapse) const;

  Double_t Millibarn2CrossSectionUnit() const;

  void SetOCDBPath(const char* path);
  
  void CorrectRunScalers(AliTriggerRunScalers* trs) const;

  void CommonRunList(const std::vector<int>& runlist1,
                     const std::vector<int>& runlist2,
                     std::vector<int>& runlist);

  TString DefaultTitle(std::vector<int>* runlist) const;

private:
  std::vector<int> fRunList; // input run list
  TString fOCDBPath; // OCDB path (default raw://)
  Bool_t fShouldCorrectForPileUp; // whether or not to correct scalers by pile-up
  TString fCrossSectionUnit; // cross-section unit (UB = microbarns) by default
  Bool_t fPurgeRunList; // make each PlotXXX call use its own runlist, purged from run w/o the requested trigger class
  
  ClassDef(AliAnalysisTriggerScalers,4) // Utility class to play with scalers
};

class AliAnalysisTriggerScalerItem : public TObject
{
public:
  AliAnalysisTriggerScalerItem(Int_t runNumber, const char* level, const char* dipoleCurrent, const char* triggerClassName, ULong64_t value, AliTriggerBCMask* mask=0x0, Int_t downscalingFactor=1, time_t duration=0)
  : fRunNumber(runNumber),fLevel(level),fDipoleCurrent(dipoleCurrent), fTriggerClassName(triggerClassName), fValue(value), fNofRuns(1), fTriggerBCMask(mask), fDS(downscalingFactor), fDuration(duration)
  {
  }
  
  AliAnalysisTriggerScalerItem(const char* triggerClassName, const char* level, ULong64_t value, AliTriggerBCMask* mask=0x0, Int_t downscalingFactor=1, time_t duration=0)
  : fRunNumber(-1),fLevel(level),fDipoleCurrent("N/A"), fTriggerClassName(triggerClassName), fValue(value), fNofRuns(0), fTriggerBCMask(mask), fDS(downscalingFactor), fDuration(duration)
  {
  }
  
  virtual Int_t	Compare(const TObject* obj) const;
  
  Int_t RunNumber() const { return fRunNumber; }
  
  const char* TriggerClassName() const { return fTriggerClassName.Data(); }
  
  const char* BCMaskName() const;
  
  ULong64_t Value() const { return fValue; }
  
  Double_t ValueCorrectedForDownscale() const { return fValue*DownscalingFactor(); }
  
  Int_t DownscalingFactor() const { return fDS; }
  
  Double_t Rate() const { return fDuration > 0 ? ValueCorrectedForDownscale() / fDuration : 0.0 ; }
  
  void Increment(ULong64_t val) { fValue += val; ++fNofRuns; }
  
  virtual void Print(Option_t* option = "") const;
  
  const char* DipoleCurrent() const { return fDipoleCurrent; }
  
  Int_t NofRuns() const { return fNofRuns; }
  
  Bool_t IsDipoleON() const { return fDipoleCurrent.Contains("6000"); }
  
  AliTriggerBCMask* BCMask() const { return fTriggerBCMask; }

  const char* Level() const { return fLevel.Data(); }
  
  time_t Duration() const { return fDuration; }
  
private:
  AliAnalysisTriggerScalerItem(const AliAnalysisTriggerScalerItem& rhs);
  AliAnalysisTriggerScalerItem& operator=(const AliAnalysisTriggerScalerItem& rhs);
  
private:
  Int_t fRunNumber; // run number for this scaler
  TString fLevel; // L0, L1 or L2
  TString fDipoleCurrent; // dipole current (A)
  TString fTriggerClassName; // trigger class name for this scaler
  ULong64_t fValue; // counter
  Int_t fNofRuns; // number of runs corresponding to counter
  AliTriggerBCMask* fTriggerBCMask; // pointer to BCMasks
  Int_t fDS; // downscaling factor
  time_t fDuration; // duration
  
  ClassDef(AliAnalysisTriggerScalerItem,5) // class to hold information about one scaler for one trigger class
};

#endif
