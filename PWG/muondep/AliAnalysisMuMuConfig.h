#ifndef ALIANALYSISMUMUPARAMETERS_H
#define ALIANALYSISMUMUPARAMETERS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

///
/// AliAnalysisMuMuConfig : helper class to store steering options
/// for the AliAnalysisMuMu and AliAnalysisMuMuEvolution classes
///
/// author : Laurent Aphecetche (Subatech)

#include "TObject.h"
#include "TString.h"

class TObjArray;

class AliAnalysisMuMuConfig : public TObject
{

public:
  
  enum EColor
  {
    kBlue=1500,
    kOrange=1501,
    kGreen=1502
  };
  
  enum ETypeList
  {
    kDimuonTriggerList=0, // list of dimuon triggers to consider
    kMuonTriggerList=1, // list of single muon triggers to consider
    kMinbiasTriggerList=2,   // list of minbias triggers to consider
    kEventSelectionList=3, // list of event types to consider
    kPairSelectionList=4, // list of pair cuts to consider
    kCentralitySelectionList=5, // list of centrality cuts to consider
    kFitTypeList=7 // list of fit types to perform    
  };
  
  AliAnalysisMuMuConfig(const char* beamYear="pPb2013");
  virtual ~AliAnalysisMuMuConfig();

  TObjArray* GetListElements(ETypeList type, Bool_t simulation) const;

  TString GetList(ETypeList type, Bool_t simulation) const;

  void SetList(ETypeList type, Bool_t simulation, const char* list);

//  TString AppendToList(ETypeList type, Bool_t simulation, const char* list);

  void SetColorScheme();
  
  void SetOCDBPath(const char* ocdbPath) { fOCDBPath = ocdbPath; }

  TString OCDBPath() const { return fOCDBPath; }
  
  void SetCompactGraphs(Bool_t value=kTRUE) { fIsCompactGraphs = value; }
  
  Bool_t CompactGraphs() { return fIsCompactGraphs; }
  
  void Print(Option_t* opt="") const;
  
  void DefineDefaults(const char* beamYear);
  
private:
  
  void ShowLists(const char* title, ETypeList type, const char separator=',', const TString& sopt="ALL") const;
  
  void ShowList(const char* title, const TString& list, const char separator=',') const;

private:

  TObjArray* fLists; // storage for lists
  
  TString fOCDBPath; // OCDB to be used (raw:// by default)
  Bool_t fIsCompactGraphs; // whether the graph produced should be compact

  ClassDef(AliAnalysisMuMuConfig,2) // class to hold configuration of AliAnalysisMuMu(Evolution) class
};

#endif
