#ifndef ALIANALYSISMULTPBCENTRALITYSELECTOR_H
#define ALIANALYSISMULTPBCENTRALITYSELECTOR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                    AliAnalysisMultPbCentralitySelector
// 
// This class selects collision candidates from data runs, applying selection cuts on triggers 
// and background rejection based on the content of the ESD
//
// Author: Michele Floris, CERN
//-------------------------------------------------------------------------

#include <AliAnalysisCuts.h>
#include <AliLog.h>

#define VERBOSE_STAT

class AliESDEvent;
class TH2F;
class TH1F;
class TCollection;
class AliTriggerAnalysis;
class AliAnalysisTaskSE;


class AliAnalysisMultPbCentralitySelector : public AliAnalysisCuts
{
public:

  AliAnalysisMultPbCentralitySelector() : fIsMC (0) {;}
  virtual ~AliAnalysisMultPbCentralitySelector(){}
    
  // AliAnalysisCuts interface
  virtual UInt_t GetSelectionMask(const TObject* obj) { return IsCentralityBinSelected((const AliESDEvent*) obj); }
  virtual Bool_t IsSelected(TList*) { AliFatal("Not implemented"); return kFALSE; }
  virtual Bool_t IsSelected(TObject* obj)  {return IsCentralityBinSelected ( (AliESDEvent*) obj);}
    
  UInt_t IsCentralityBinSelected(const AliESDEvent* aEsd){ return kTRUE;}
    
  void SetAnalyzeMC(Bool_t flag = kTRUE) { fIsMC = flag; }
    
  virtual void Print(Option_t* option = "") const { Printf ("Multiplitity Selector [AliAnalysisMultPbCentralitySelector] [%s]", option);}
  virtual Long64_t Merge(TCollection* list){list->GetEntries();return 0;}
  
protected:
  Bool_t fIsMC;             // flag if MC is analyzed

  ClassDef(AliAnalysisMultPbCentralitySelector, 1)
    
  private:
  AliAnalysisMultPbCentralitySelector(const AliAnalysisMultPbCentralitySelector&); // not implemented
  AliAnalysisMultPbCentralitySelector& operator=(const AliAnalysisMultPbCentralitySelector&); // not implemented
};

#endif
