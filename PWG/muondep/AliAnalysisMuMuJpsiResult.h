#ifndef ALIANALYSISMUMUJPSIRESULT_H
#define ALIANALYSISMUMUJPSIRESULT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

///
/// AliAnalysisMuMuJpsiResult : helper class to store Jpsi results from
/// AliAnalysisTaskMuMu
///
/// author : Laurent Aphecetche (Subatech)

#include "TNamed.h"
#include <TString.h>
#include "AliAnalysisMuMuResult.h"
#include "AliAnalysisMuMuBinning.h"

class TH1;
class THashList;
class TF1;
class TMap;

class AliAnalysisMuMuJpsiResult : public AliAnalysisMuMuResult
{
  
public:
  
  AliAnalysisMuMuJpsiResult(TRootIOCtor* io);
  
  AliAnalysisMuMuJpsiResult(const TH1& hminv);

  AliAnalysisMuMuJpsiResult(const TH1& hminv,
                        const char* fitType,
                        Int_t nrebin);

  AliAnalysisMuMuJpsiResult(const TH1& hminv,
                        const char* triggerClass,
                        const char* eventSelection,
                        const char* pairSelection,
                        const char* centSelection,
                        const AliAnalysisMuMuBinning::Range& bin);
  
  AliAnalysisMuMuJpsiResult(const AliAnalysisMuMuJpsiResult& rhs);
  AliAnalysisMuMuJpsiResult& operator=(const AliAnalysisMuMuJpsiResult& rhs);
  
  virtual ~AliAnalysisMuMuJpsiResult();

  virtual TObject* Clone(const char* newname = "") const;
  
  Bool_t Correct(const AliAnalysisMuMuJpsiResult& other, const char* particle, const char* subResultName="");
  
  TH1* Minv() const { return fMinv; }
  
  Int_t NofTriggers() const;
  
  void SetNofTriggers(Int_t n);
  
  void Print(Option_t* opt="") const;
  
  Bool_t AddFit(const char* fitType, Int_t npar=0, Double_t* par=0x0);

  AliAnalysisMuMuJpsiResult* CountJpsi(TH1& h);

  AliAnalysisMuMuJpsiResult*  FitJpsi(TH1& h);

  AliAnalysisMuMuJpsiResult* FitJpsiNA48(const TH1& h);
  AliAnalysisMuMuJpsiResult* FitJpsiCB2VWG(const TH1& h);
  AliAnalysisMuMuJpsiResult* FitJpsi2CB2VWG(const TH1& h, Double_t alphaLow=-1.0, Double_t nLow=-1.0, Double_t alphaUp=-1.0, Double_t nUp=-1.0);
  
  AliAnalysisMuMuJpsiResult* FitJpsiGCBE(TH1& h);
  
  Int_t NofRuns() const;
  
  void SetNofRuns(int n);
  
  const AliAnalysisMuMuBinning::Range& Bin() const;

  void SetBin(const AliAnalysisMuMuBinning::Range& bin);
  
  void SetNofInputParticles(const char* particle, int n);

  void SetNofInputParticles(const TH1& hminv);

  void SetMinv(const TH1& hminv);

  Long64_t Merge(TCollection* list);

  static Double_t CountParticle(const TH1& hminv, const char* particle, Double_t sigma=-1.0);
  
  virtual AliAnalysisMuMuJpsiResult* Mother() const { return static_cast<AliAnalysisMuMuJpsiResult*>(AliAnalysisMuMuResult::Mother()); }

  void PrintValue(const char* key, const char* opt, Double_t value, Double_t errorStat, Double_t rms=0.0) const;
  
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
  Int_t fRebin; // rebin level of minv spectra
  
  TString fTriggerClass; // trigger class for this result
  TString fEventSelection; // event selection for this result
  TString fPairSelection; // pair selection for this result
  TString fCentralitySelection; // centrality selection for this result

  ClassDef(AliAnalysisMuMuJpsiResult,1) // a class to hold invariant mass analysis results (counts, yields, AccxEff, R_AB, etc...)
};

#endif
