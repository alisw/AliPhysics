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

class AliTriggerBCMask;
class AliAnalysisTriggerScalerItem;
class TGraph;

class AliAnalysisTriggerScalers : public TObject
{
public:
  AliAnalysisTriggerScalers(Int_t runNumber, const char* ocdbPath="raw://");
  AliAnalysisTriggerScalers(const char* runlist, const char* ocdbPath="raw://");
  virtual ~AliAnalysisTriggerScalers();
  
  Int_t GetTriggerInput(Int_t runNumber, const char* inputname);

  AliAnalysisTriggerScalerItem* GetTriggerScaler(Int_t runNumber, const char* level, const char* triggerClassName);

  void IntegratedLuminosity();
    
  TGraph* PlotTriggerRatio(const char* triggerClassName1,
                           const char* what1,
                           const char* triggerClassName2,
                           const char* what2);

  TGraph* PlotTrigger(const char* triggerClassName, const char* what);

  TGraph* PlotTriggerEvolution(const char* triggerClassName,
                               const char* what,
                               bool draw=kTRUE,
                               double* mean=0x0);

  void SetRunList(Int_t runNumber);
  void SetRunList(const char* runlist);

  void SetVerbose(Int_t level=1) { fVerbose=level; }
    
  static void ReadIntegers(const char* filename, std::vector<int>& integers, Bool_t resetVector=kTRUE);

  static void PrintIntegers(std::vector<int>& integers, const char sep = '\n',
                            std::ostream& out = std::cout);

private:
  TObject* GetOCDBObject(const char* path, Int_t runNumber);
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
  Int_t fRunNumber;
  TString fLevel;
  TString fDipoleCurrent;
  TString fTriggerClassName;
  ULong64_t fValue;  
  Int_t fNofRuns;
  AliTriggerBCMask* fTriggerBCMask;
  Int_t fDS;
  time_t fDuration;
  
  ClassDef(AliAnalysisTriggerScalerItem,5)
};

#endif
