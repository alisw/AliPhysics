#ifndef ALIANALYSISTASKDPLUS_H
#define ALIANALYSISTASKDPLUS_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskSEDplus
// AliAnalysisTaskSE for the D+ candidates Invariant Mass Histogram and 
//comparison of heavy-flavour decay candidates
// to MC truth (kinematics stored in the AOD)
// Renu Bala, bala@to.infn.it
//*************************************************************************

#include <TROOT.h>
#include <TSystem.h>
#include <TNtuple.h>
#include <TH1F.h>
#include <TArrayD.h>

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"

class AliAnalysisTaskSEDplus : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskSEDplus();
  AliAnalysisTaskSEDplus(const char *name, Bool_t fillNtuple=kFALSE);
  virtual ~AliAnalysisTaskSEDplus();

  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetDoLikeSign(Bool_t dols=kTRUE){fDoLS=dols;}
  void SetMassLimits(Float_t range);
  void SetMassLimits(Float_t lowlimit, Float_t uplimit);
  void SetPtBinLimit(Int_t n, Double_t *limitarray);
  
  Float_t GetUpperMassLimit(){return fUpmasslimit;}
  Float_t GetLowerMassLimit(){return fLowmasslimit;}
  Int_t GetNBinsPt(){return fNPtBins;}
  Double_t GetPtBinLimit(Int_t ibin);
  
  void LSAnalysis(TClonesArray *arrayOppositeSign,TClonesArray *arrayLikeSign,AliAODEvent *aod,AliAODVertex *vtx1, Int_t nDplusOS);

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:

  AliAnalysisTaskSEDplus(const AliAnalysisTaskSEDplus &source);
  AliAnalysisTaskSEDplus& operator=(const AliAnalysisTaskSEDplus& source); 
  Int_t GetHistoIndex(Int_t iPtBin) const { return iPtBin*3;}
  Int_t GetSignalHistoIndex(Int_t iPtBin) const { return iPtBin*3+1;}
  Int_t GetBackgroundHistoIndex(Int_t iPtBin) const { return iPtBin*3+2;}
  Int_t GetLSHistoIndex(Int_t iPtBin)const { return iPtBin*5;}
 
  enum {kMaxPtBins=10};

  TList   *fOutput; //! list send on output slot 0
  TH1F    *fHistNEvents; //!hist. for No. of events
  TH1F *fMassHist[3*kMaxPtBins]; //!hist. for inv mass (LC)
  TH1F *fMassHistTC[3*kMaxPtBins]; //!hist. for inv mass (TC)
  TH1F *fMassHistLS[5*kMaxPtBins];//!hist. for LS inv mass (LC)
  TH1F *fMassHistLSTC[5*kMaxPtBins];//!hist. for LS inv mass (TC)
  TNtuple *fNtupleDplus; //! output ntuple
  Float_t fUpmasslimit;  //upper inv mass limit for histos
  Float_t fLowmasslimit; //lower inv mass limit for histos
  Int_t fNPtBins; //number of bins in Pt for histograms
  Double_t fArrayBinLimits[kMaxPtBins+1]; //limits for the Pt bins
  Bool_t fFillNtuple;   // flag for filling ntuple
  Bool_t fReadMC;    //flag for access to MC
  Bool_t fDoLS;      //flag to do LS analysis
  AliAnalysisVertexingHF *fVHF;  // Vertexer heavy flavour (used to pass the cuts)
  
  ClassDef(AliAnalysisTaskSEDplus,4); // AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
};

#endif

