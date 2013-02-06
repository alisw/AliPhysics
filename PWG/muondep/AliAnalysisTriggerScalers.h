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

class AliAnalysisTriggerScalers : public TObject
{
public:
  AliAnalysisTriggerScalers(const std::vector<int>& runs, const char* ocdbPath="raw://");
  AliAnalysisTriggerScalers(const std::set<int>& runs, const char* ocdbPath="raw://");
  AliAnalysisTriggerScalers(Int_t runNumber, const char* ocdbPath="raw://");
  AliAnalysisTriggerScalers(const char* runlist, const char* ocdbPath="raw://");
  virtual ~AliAnalysisTriggerScalers();
  
  Int_t GetTriggerInput(Int_t runNumber, const char* inputname);

  Int_t GetFillNumberFromRunNumber(Int_t runNumber);

  AliAnalysisTriggerScalerItem* GetTriggerScaler(Int_t runNumber, const char* level, const char* triggerClassName);

  void IntegratedLuminosity(const char* triggerList="",
                            const char* lumiTrigger="C0TVX-B-NOPF-ALLNOTRD",
                            Double_t lumiCrossSection=0.755*2000,
                            const char* csvOutputFile="",
                            const char sep='\t',
                            const char* csUnit="ub");
    
  TGraph* PlotTriggerRatio(const char* triggerClassName1,
                           const char* what1,
                           const char* triggerClassName2,
                           const char* what2);

  TGraph* PlotTrigger(const char* triggerClassName, const char* what);

  TGraph* PlotTriggerEvolution(const char* triggerClassName,
                               const char* what,
                               bool draw=kTRUE,
                               double* mean=0x0,
                               bool removeZero=kFALSE);

  void SetRunList(const std::vector<int>& runlist);
  void SetRunList(const std::set<int>& runlist);
  void SetRunList(Int_t runNumber);
  void SetRunList(const char* runlist);

  virtual void Print(Option_t* opt="") const;
  
  void SetVerbose(Int_t level=1) { fVerbose=level; }
    
  static void ReadIntegers(const char* filename, std::vector<int>& integers, Bool_t resetVector=kTRUE);

  static void PrintIntegers(const std::vector<int>& integers, const char sep = '\n',
                            std::ostream& out = std::cout);

  TObject* GetOCDBObject(const char* path, Int_t runNumber);
  
  const std::vector<int>& GetRunList() const { return fRunList; }

  void DrawFills(Double_t ymin, Double_t ymax);
  
private:
  
  void DrawFill(Int_t run1, Int_t run2, double ymin, double ymax, const char* label);

  void GetFillBoundaries(std::map<int, std::pair<int,int> >& fills);
  
  Float_t GetPauseAndConfigCorrection(Int_t runNumber, const char* triggerClassName);

  TGraph* MakeGraph(const std::vector<int>& vx, const std::vector<int>& vex,
                    const std::vector<double>& vy, const std::vector<double>& vey);

private:
  std::vector<int> fRunList; // input run list
  Int_t fVerbose; // whether or not to be verbose
  TString fOCDBPath; // OCDB path (default raw://)
  
  ClassDef(AliAnalysisTriggerScalers,2) // Utility class to play with scalers
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
