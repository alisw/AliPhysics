#ifndef ALIANALYSISMUMURESULT_H
#define ALIANALYSISMUMURESULT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

/**

  @ingroup pwg_muondep_mumu

  @class AliAnalysisMuMuResult

  @brief Store results from AliAnalysisTaskMuMu

 Base class to hold a set of results for the same quantity,
 computed using various methods, each with their errors

 @author Laurent Aphecetche (Subatech)
 @author Benjamin Audurier (Subatech)
 */

#include "TNamed.h"
#include <TString.h>

class TH1;
class THashList;
class TMap;

class AliAnalysisMuMuResult : public TNamed
{

public:

  enum EResultMergingMethod { kMean, kSum };

  AliAnalysisMuMuResult(const char* name="", const char* title="", AliAnalysisMuMuResult::EResultMergingMethod mergindMethod=AliAnalysisMuMuResult::kMean);
  AliAnalysisMuMuResult(const AliAnalysisMuMuResult& rhs);
  AliAnalysisMuMuResult& operator=(const AliAnalysisMuMuResult& rhs);

  virtual ~AliAnalysisMuMuResult();

  Bool_t AdoptSubResult(AliAnalysisMuMuResult* r);

  virtual TObject* Clone(const char* newname = "") const;

  void Set(const char* name, Double_t value, Double_t errorStat, Double_t rms=0.0);

  Int_t HasValue(const char* name, const char* subResultName="") const;

  Double_t GetValue(const char* name, const char* subResultName="") const;

  Double_t GetErrorStat(const char* name, const char* subResultName="") const;

  Double_t GetRMS(const char* name, const char* subResultName="") const;

  void Print(Option_t* opt="") const;

  AliAnalysisMuMuResult* SubResult(const char* subResultName) const;

  TObjArray* SubResults() const { return fSubResults; }

  Long64_t Merge(TCollection* list);

  AliAnalysisMuMuResult* Mother() const { return fMother; }

  THashList* Keys() const;

  Double_t Weight() const { return fWeight; }

  void SetWeight(Double_t w) { fWeight=w; }

  static Double_t ErrorAB(Double_t a, Double_t aerr, Double_t b, Double_t berr);

  static Double_t ErrorABC(Double_t a, Double_t aerr, Double_t b, Double_t berr, Double_t c, Double_t cerror);

  static Double_t ErrorABCD(Double_t a, Double_t aerr, Double_t b, Double_t berr, Double_t c, Double_t cerror,
                            Double_t d, Double_t derror);

  static Double_t ErrorABCDE(Double_t a, Double_t aerr, Double_t b, Double_t berr, Double_t c, Double_t cerror,
                             Double_t d, Double_t derror, Double_t e, Double_t eerror);

  void PrintValue(const char* key, const char* opt, Double_t value, Double_t errorStat, Double_t rms=0.0) const;

  void SetAlias(const char* alias) { fAlias = alias; }

  TString Alias() const { if ( fAlias.Length()>0) return fAlias; else return GetName(); }

  void Include(const char* subResultsList);

  void Exclude(const char* subResultsList);

  Bool_t IsIncluded(const TString& alias) const;

  void Scale(Double_t value);

  void SetMergingMethod(AliAnalysisMuMuResult::EResultMergingMethod mergindMethod) { fResultMergingMethod=mergindMethod; }

  Bool_t IsValid() const { return fIsValid; }

  Bool_t Map() const { return fMap; }

  void Invalidate() { fIsValid = kFALSE; }

  void Show(const char* keyPattern);

  void Hide(const char* keyPattern);

  Bool_t IsValidValue(Double_t val) const;

private:

  enum EIndex
  {
    kValue=0,
    kErrorStat=1,
    kRMS=2
  };

  void PrintParticle(const char* particle, const char* opt) const;

  TList* SubResultsToBeIncluded() const;

  TString GetSubResultNameList() const;

  Int_t NofIncludedSubResults(const char* name) const;

private:
  TObjArray            * fSubResults; // TObjArray of AliAnalysisMuMuResult*
  TMap                 * fMap; // internal parameter map
  AliAnalysisMuMuResult* fMother; // mother result
  mutable THashList    * fKeys; //! keys we have in our internal map (or the one of our subresults)
  Double_t fWeight; // weight of this result (default 1.0)
  TString fAlias; // alias name
  mutable TList        * fSubResultsToBeIncluded; // inclusion list
  EResultMergingMethod fResultMergingMethod; // how to merge result (e.g. mean or sum)
  Bool_t fIsValid; // is this result valid ?
  mutable THashList    * fVisibleKeys; // list of keys to show with the Print method

  ClassDef(AliAnalysisMuMuResult,14) // a class to some results (counts, yields, AccxEff, R_AB, etc...)
};

#endif
