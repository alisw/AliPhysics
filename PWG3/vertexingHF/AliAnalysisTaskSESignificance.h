#ifndef ALIANALYSISTASKSESIGNIFICANCE_H
#define ALIANALYSISTASKSESIGNIFICANCE_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//*************************************************************************
// AliAnalysisTaskSESignificane to calculate effects on 
// significance of D mesons  cut 
// Authors: G. Ortona, ortona@to.infn.it
// F. Prino, prino@to.infn.it
// Renu Bala, bala@to.infn.it
// Chiara Bianchin, cbianchi@pd.infn.it
//*************************************************************************

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisVertexingHF.h"

class TH1F;
class AliMultiDimVector;
class AliRDHFCuts;

class AliAnalysisTaskSESignificance : public AliAnalysisTaskSE
{
 public:

  enum FeedDownEnum {kBoth,kCharmOnly,kBeautyOnly};

  AliAnalysisTaskSESignificance();
  AliAnalysisTaskSESignificance(const char *name, TList *listMDV,AliRDHFCuts *RDCuts, Int_t decaychannel,Int_t selectionlevel=AliRDHFCuts::kAll);
 
  virtual ~AliAnalysisTaskSESignificance();

  Bool_t CheckConsistency();
  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetBFeedDown(FeedDownEnum flagB);//see enum
  void SetDFromCharmOnly(){SetBFeedDown(kCharmOnly);}
  void SetDFromBeautyOnly(){SetBFeedDown(kBeautyOnly);}
  void SetMassLimits(Float_t range,Int_t pdg);
  void SetMassLimits(Float_t lowlimit, Float_t uplimit);
  void SetNBins(Int_t nbins){fNBins=nbins;}
  void SetFillWithPartAntiPartBoth(Int_t value){fPartOrAndAntiPart=value;}
  //void SetMultiVector(const AliMultiDimVector *MultiDimVec){fMultiDimVec->CopyStructure(MultiDimVec);}
  Float_t GetUpperMassLimit()const {return fUpmasslimit;}
  Float_t GetLowerMassLimit()const {return fLowmasslimit;}
  Int_t GetNBins()const {return fNBins;}
  Int_t GetFillWithPartAntiPartBoth()const {return fPartOrAndAntiPart;}
  Int_t GetBFeedDown()const {return fBFeedDown;}

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void LocalInit();// {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:

  AliAnalysisTaskSESignificance(const AliAnalysisTaskSESignificance &source);
  AliAnalysisTaskSESignificance& operator=(const AliAnalysisTaskSESignificance& source);
  void SetPDGCodes();
  Int_t GetHistoIndex(Int_t iPtBin) const { return iPtBin*3;}
  Int_t GetSignalHistoIndex(Int_t iPtBin) const { return iPtBin*3+1;}
  Int_t GetBackgroundHistoIndex(Int_t iPtBin) const { return iPtBin*3+2;}
  Int_t GetLSHistoIndex(Int_t iPtBin)const { return iPtBin*5;}


  void FillDplus(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index,Int_t isSel);
  void FillD02p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index, Int_t isSel);
  void FillDs(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index,Int_t isSel);
  void FillDstar(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index,Int_t isSel);
  void FillD04p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index,Int_t isSel);
  void FillLambdac(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t index, Int_t isSel);


  enum {kMaxPtBins=5};
  enum {kMaxCutVar=5};
  enum {kMaxSteps=10};
  enum {kMaxNHist=500000};
  enum {kDplustoKpipi,kD0toKpi,kDstartoKpipi,kDstoKKpi,kD0toKpipipi,kLambdactopKpi};

  TList   *fOutput; //! list send on output slot 0
  TList *fCutList; //Multidimvector container
  TH1F *fMassHist[kMaxNHist]; //!hist. for inv mass
  TH1F *fSigHist[kMaxNHist]; //!hist. for inv mass (sig from MC truth)
  TH1F *fBkgHist[kMaxNHist]; //!hist. for inv mass (bkg from MC truth)
  TH1F *fRflHist[kMaxNHist]; //!hist. for inv mass (bkg from MC truth)
  TH1F* fHistNEvents; //! hist of n of aods
  Float_t fUpmasslimit;  //upper inv mass limit for histos
  Float_t fLowmasslimit; //lower inv mass limit for histos
  AliRDHFCuts *fRDCuts;//prong cut values
  Int_t fNPtBins; //number of pt bins
  Bool_t fReadMC;    //flag for access to MC
  FeedDownEnum fBFeedDown; //flag to search for D from B decays
  Int_t fDecChannel; //decay channel identifier
  Int_t fPDGmother;  // PDG code of D meson
  Int_t fNProngs;         // number of prong of the decay channel  
  Int_t fPDGdaughters[4]; // PDG codes of daughters
  TString fBranchName;    // AOD branch name for channel
  Int_t fSelectionlevel;  //selection level: kALL,kTracks,kCandidate
  Int_t   fNVars;         // number of selection variables
  Float_t fVars[kMaxCutVar];       // array with values of cut variables
  Int_t fNBins;  //number of bins in the mass histograms
  Int_t fPartOrAndAntiPart;  //fill histograms with particle only (+1), antiparticle only (-1), both (0)

  ClassDef(AliAnalysisTaskSESignificance,3); // AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
};

#endif


