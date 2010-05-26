#ifndef ALIANALYSISTASKDS_H
#define ALIANALYSISTASKDS_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id: $ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Analysis task to produce Ds candidates mass spectra           //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TH2F.h>

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsDstoKKpi.h"

class AliAnalysisTaskSEDs : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEDs();
  AliAnalysisTaskSEDs(const char *name, AliRDHFCutsDstoKKpi* productioncuts, AliRDHFCutsDstoKKpi* analysiscuts);
  virtual ~AliAnalysisTaskSEDs();
  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetMassRange(Double_t rang=0.2){fMassRange=rang;}
  void SetInvMassBinSize(Double_t binsiz=0.002){fMassBinSize=binsiz;}
  void SetPtBins(Int_t n, Float_t* lim);
  void SetProductionCuts(AliRDHFCutsDstoKKpi* cuts){fProdCuts=cuts;}
  void SetAnalysisCuts(AliRDHFCutsDstoKKpi* cuts){fAnalysisCuts=cuts;}
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:
  Int_t GetHistoIndex(Int_t iPtBin) const { return iPtBin*3;}
  Int_t GetSignalHistoIndex(Int_t iPtBin) const { return iPtBin*3+1;}
  Int_t GetBackgroundHistoIndex(Int_t iPtBin) const { return iPtBin*3+2;}

  enum {kMaxPtBins=10};

  AliAnalysisTaskSEDs(const AliAnalysisTaskSEDs &source);
  AliAnalysisTaskSEDs& operator=(const AliAnalysisTaskSEDs& source); 

  TList*  fOutput;                    //! list send on output slot 0
  TH1F*   fHistNEvents;               //! hist. for No. of events  
  TH1F*   fChanHist[3];               //! hist. with KKpi and piKK candidates (sig,bkg,tot)
  TH1F*   fChanHistCuts[3];           //! hist. with KKpi and piKK candidates analysis cuts
  TH1F*   fMassHist[3*kMaxPtBins];    //! hist. of mass spectra (sig,bkg,tot)
  TH1F*   fMassHistCuts[3*kMaxPtBins];//! hist. of mass spectra (sig,bkg,tot) analysis cuts
  TH1F*   fCosPHist[3*kMaxPtBins];    //! hist. of cos pointing angle (sig,bkg,tot)
  TH1F*   fDLenHist[3*kMaxPtBins];    //! hist. of decay length (sig,bkg,tot)
  TH2F*   fDalitz[3*kMaxPtBins];      //! dalitz plot (sig,bkg,tot)
  Bool_t  fReadMC;                    //  flag for access to MC
  UChar_t fNPtBins;                   // number of Pt bins
  TList *fListCuts; //list of cuts
  Float_t fPtLimits[kMaxPtBins+1];    //  limits for pt bins
  Double_t fMassRange;                // range for mass histogram 
  Double_t fMassBinSize;              // bin size for inv. mass histo

  AliRDHFCutsDstoKKpi *fProdCuts;     // Cuts for Analysis
  AliRDHFCutsDstoKKpi *fAnalysisCuts; // Cuts for Analysis
  
  ClassDef(AliAnalysisTaskSEDs,3);    //  AliAnalysisTaskSE for Ds mass spectra
};

#endif

