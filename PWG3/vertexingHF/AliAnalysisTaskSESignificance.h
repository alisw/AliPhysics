#ifndef ALIANALYSISTASKSESIGNIFICANCE_H
#define ALIANALYSISTASKSESIGNIFICANCE_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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

  AliAnalysisTaskSESignificance();
  AliAnalysisTaskSESignificance(const char *name, TList *listMDV,AliRDHFCuts *RDCuts, Int_t decaychannel,Int_t selectionlevel=AliRDHFCuts::kAll);
 
  virtual ~AliAnalysisTaskSESignificance();

  Bool_t CheckConsistency();
  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetMassLimits(Float_t range,Int_t pdg);
  void SetMassLimits(Float_t lowlimit, Float_t uplimit);
  //void SetMultiVector(const AliMultiDimVector *MultiDimVec){fMultiDimVec->CopyStructure(MultiDimVec);}
  Float_t GetUpperMassLimit(){return fUpmasslimit;}
  Float_t GetLowerMassLimit(){return fLowmasslimit;}
  
  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void LocalInit();// {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
    
 private:

  AliAnalysisTaskSESignificance(const AliAnalysisTaskSESignificance &source);
  AliAnalysisTaskSESignificance& operator=(const AliAnalysisTaskSESignificance& source); 
  Int_t GetHistoIndex(Int_t iPtBin) const { return iPtBin*3;}
  Int_t GetSignalHistoIndex(Int_t iPtBin) const { return iPtBin*3+1;}
  Int_t GetBackgroundHistoIndex(Int_t iPtBin) const { return iPtBin*3+2;}
  Int_t GetLSHistoIndex(Int_t iPtBin)const { return iPtBin*5;}
 
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
  Int_t fDecChannel; //decay channel identifier
  Int_t fSelectionlevel;//selection level: kALL,kTracks,kCandidate
   
  ClassDef(AliAnalysisTaskSESignificance,1); // AliAnalysisTaskSE for the MC association of heavy-flavour decay candidates
};

#endif


