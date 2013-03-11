#ifndef ALIANALYSISMUMURESULT_H
#define ALIANALYSISMUMURESULT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// $Id$

///
/// AliAnalysisMuMuResult : helper class to store results from
/// AliAnalysisTaskMuMu
///
/// author : Laurent Aphecetche (Subatech)

#include "TNamed.h"
#include <TString.h>
#include "AliAnalysisMuMuBinning.h"

class TH1;
class THashList;
class TF1;
class TMap;

class AliAnalysisMuMuResult : public TNamed
{
  
public:
  
  AliAnalysisMuMuResult(TRootIOCtor* io);
  
  AliAnalysisMuMuResult(const TH1& hminv);

  AliAnalysisMuMuResult(const TH1& hminv,
                        const char* fitType,
                        Int_t nrebin);

  AliAnalysisMuMuResult(const TH1& hminv,
                        const char* triggerClass,
                        const char* eventSelection,
                        const char* pairSelection,
                        const char* centSelection,
                        const AliAnalysisMuMuBinning::Range& bin);
  
  AliAnalysisMuMuResult(const AliAnalysisMuMuResult& rhs);
  AliAnalysisMuMuResult& operator=(const AliAnalysisMuMuResult& rhs);
  
  virtual ~AliAnalysisMuMuResult();

  virtual TObject* Clone(const char* newname = "") const;
  
  Bool_t Correct(const AliAnalysisMuMuResult& other, const char* particle, const char* subResultName="");
  
  TH1* Minv() const { return fMinv; }
  
  void Set(const char* name, Double_t value, Double_t errorStat);
  
  Bool_t HasValue(const char* name, const char* subResultName="") const;
  
  Double_t GetValue(const char* name, const char* subResultName="") const;
  
  Double_t GetErrorStat(const char* name, const char* subResultName="") const;

  Int_t NofTriggers() const;
  
  void SetNofTriggers(Int_t n);
  
  void Print(Option_t* opt="") const;
  
  Bool_t AddFit(const char* fitType, Int_t npar=0, Double_t* par=0x0);

  AliAnalysisMuMuResult* CountJpsi(TH1& h);

  AliAnalysisMuMuResult*  FitJpsi(TH1& h);

  AliAnalysisMuMuResult* FitJpsiNA48(const TH1& h);
  AliAnalysisMuMuResult* FitJpsiCB2VWG(const TH1& h);
  AliAnalysisMuMuResult* FitJpsi2CB2VWG(const TH1& h, Double_t alphaLow=-1.0, Double_t nLow=-1.0, Double_t alphaUp=-1.0, Double_t nUp=-1.0);
  
  AliAnalysisMuMuResult* FitJpsiGCBE(TH1& h);
  
  Int_t NofRuns() const;
  
  void SetNofRuns(int n);
  
  const AliAnalysisMuMuBinning::Range& Bin() const;

  void SetBin(const AliAnalysisMuMuBinning::Range& bin);
  
  void SetNofInputParticles(const char* particle, int n);

  void SetNofInputParticles(const TH1& hminv);

  void SetMinv(const TH1& hminv);

  AliAnalysisMuMuResult* SubResult(const char* subResultName) const;
  
  TObjArray* SubResults() const { return fSubResults; }
  
  Long64_t Merge(TCollection* list);

  AliAnalysisMuMuResult* Mother() const { return fMother; }
  
  THashList* Keys() const;
  
  Double_t Weight() const { return fWeight > 0  ? fWeight : fNofTriggers; }
  
  void SetWeight(Double_t w) { fWeight=w; }

  static Double_t CountParticle(const TH1& hminv, const char* particle, Double_t sigma=-1.0);
  
  static Double_t ErrorAB(Double_t a, Double_t aerr, Double_t b, Double_t berr);
  
  static Double_t ErrorABC(Double_t a, Double_t aerr, Double_t b, Double_t berr, Double_t c, Double_t cerror);

  static void PrintValue(const char* key, const char* opt, Double_t value, Double_t errorStat);

  void SetAlias(const char* alias) { fAlias = alias; }
  
  TString Alias() const { if ( fAlias.Length()>0) return fAlias; else return GetName(); }
  
private:
  
  enum EIndex
  {
    kValue=0,
    kErrorStat=1
  };
  
  void PrintParticle(const char* particle, const char* opt) const;

private:
  Int_t fNofRuns; // number of runs used to get this result
  Int_t fNofTriggers; // number of trigger analyzed
  TH1* fMinv; // invariant mass spectrum
  AliAnalysisMuMuBinning::Range fBin; // bin range
  TObjArray* fSubResults; // TObjArray of AliAnalysisMuMuResult*
  TMap* fMap; // internal parameter map
  AliAnalysisMuMuResult* fMother; // mother result
  mutable THashList* fKeys; //! keys we have in our internal map (or the one of our subresults)
  Double_t fWeight; // weight of this result (default 1.0)
  Int_t fRebin; // rebin level of minv spectra
  
  TString fTriggerClass; // trigger class for this result
  TString fEventSelection; // event selection for this result
  TString fPairSelection; // pair selection for this result
  TString fCentralitySelection; // centrality selection for this result

  TString fAlias; // alias name
  
  ClassDef(AliAnalysisMuMuResult,8) // a class to hold invariant mass analysis results (counts, yields, AccxEff, R_AB, etc...)
};

#endif
