/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id: $ */

/////////////////////////////////////////////////
// AliAnalysisTaskCMEV0PID:
// SimpleTask to start Analyzing
// ALICE Data
// Author: Rihan Haque (mhaque@cern.ch)
//
//
/////////////////////////////////////////////////

#ifndef ALIANALYSISCMXQA_H
#define ALIANALYSISCMXQA_H

#include "AliAnalysisTaskSE.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"
#include "TFile.h"
#include "TString.h"
#include "TProfile.h"
#include "TProfile2D.h"

class    AliVEvent;
class    AliVVertex;    
class    AliESDEvent;       
class    AliAODEvent;
class    AliAODTrack;
class    AliPIDResponse;    
class    AliMultSelection;    
class    AliAnalysisUtils;





class AliAnalysisTaskCMXQA : public AliAnalysisTaskSE {

 public:

  //-----> Mandatory Functions:
  AliAnalysisTaskCMXQA();
  AliAnalysisTaskCMXQA(const char *name);
  virtual ~AliAnalysisTaskCMXQA();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * /*option*/);



  
  //-----> User Defined Functions:
  void  SetupEventAndTaskConfigInfo();
  Int_t GetCentralityScaled0to10(Double_t fCent);
  
  /******* Event Cut Ranges ******/
  void SetCentralityPercentileMin(Double_t centMin) {this->fCentralityMin  = centMin;}
  void SetCentralityPercentileMax(Double_t centMax) {this->fCentralityMax  = centMax;}
  void SetCentralityEstimator(TString sEstim)       {this->sCentrEstimator = sEstim;}
  void SetVzRangeMin(Double_t vzMin)                {this->fMinVzCut       =  vzMin;}
  void SetVzRangeMax(Double_t vzMax)                {this->fMaxVzCut       =  vzMax;}

  /******* Track Cut Ranges ******/
  void SetNSigmaCutTPC(Double_t     nSigTPC)     {this->fNSigmaCut  =  nSigTPC;}
  void SetTrackCutdEdxMin(Float_t  ndEdxMin)     {this->fdEdxMin    = ndEdxMin;}
  void SetTrackCutChi2Min(Double_t  chi2Min)     {this->fTrkChi2Min =  chi2Min;}
  void SetFlagUseKinkTracks(Bool_t    bKink)     {this->bUseKinkTracks = bKink;}
  
  void SetTrackCutNclusterMin(Int_t ncl)         {this->fTPCclustMin = ncl;}
  void SetFilterBit(Int_t fb)                    {this->fFilterBit   =  fb;}
  void SetEtaRangeMin(Double_t emn)              {this->fMinEtaCut   = emn;}
  void SetEtaRangeMax(Double_t emx)              {this->fMaxEtaCut   = emx;}
  void SetPtRangeMin(Double_t ptL)               {this->fMinPtCut    = ptL;}
  void SetPtRangeMax(Double_t ptH)               {this->fMaxPtCut    = ptH;}
  //------ End of user defined function -------

 protected:
  // nothing here
  
 private:
  
  AliVEvent             *fVevent;             //! event
  AliESDEvent           *fESD;                //! esd Event
  AliAODEvent           *fAOD;                //! aod Event
  AliPIDResponse        *fPIDResponse;        //! PID Handler
  AliMultSelection      *fMultSelection;      //! For Centrality 
  AliAnalysisUtils      *fAnalysisUtil;       //! Event Selection Options
  TList                 *fListHist;           //! OutputList


  /// Functions for Pile Up Event Removal:
  TF1                   *fV0CutPU;      //!
  TF1                   *fSPDCutPU;     //!
  TF1                   *fMultCutPU;    //!
  TF1                   *fCenCutLowPU;  //!
  TF1                   *fCenCutHighPU; //!




  
  //Set from the AddTask:
  Float_t         fCentralityMin;  //
  Float_t         fCentralityMax;  //
  
  //Track Variables to be used:
  Int_t                fFilterBit;  //
  Int_t              fTPCclustMin;  // 
  Bool_t           bUseKinkTracks;  //
  
  Float_t              fNSigmaCut;  //
  Float_t               fMinPtCut;  //
  Float_t               fMaxPtCut;  //
  Float_t              fMinEtaCut;  //
  Float_t              fMaxEtaCut;  //
  Float_t             fTrkChi2Min;  //
  Float_t                fdEdxMin;  //  
  //Event Variables to be used:
  Float_t               fMinVzCut;  //
  Float_t               fMaxVzCut;  //
  
  TString         sCentrEstimator;  //


  




       
  
  /////////// NOW HISTOGRAMS /////////
  
  //QA histograms:
  TH1F                  *fCentDistBeforCut;         //! Cent before Any Cut
  TH1F                  *fCentDistAfterCut;         //! Cent After All Cuts



  
  TH1F                  *fHistEventCount;           //!   
  /// Arrays:
  TH2F                  *fHistEtaPtwCutChPos[3][4];   //! Eta-Pt  ChPos
  TH2F                  *fHistEtaPhiwCutChPos[3][4];  //! Eta-Phi ChPos
  
  TH2F                  *fHistEtaPtwCutChNeg[3][4];   //! Eta-Pt  ChNeg 
  TH2F                  *fHistEtaPhiwCutChNeg[3][4];  //! Eta-Phi ChNeg

  TH2F                  *fHistEtaVzwCutChPos[3][4];   //! Eta-Vz  ChPos
  TH2F                  *fHistEtaVzwCutChNeg[3][4];   //! Eta-Vz  ChPo

  

  Bool_t  CheckPileUp2018(AliAODEvent *faod);
  Bool_t  CheckPIDofParticle(AliAODTrack* ftrack,Int_t pidToCheck=1);
  void    SetupPileUpRemovalFunctions();
    
  AliAnalysisTaskCMXQA(const AliAnalysisTaskCMXQA &other);
  AliAnalysisTaskCMXQA &operator=(const AliAnalysisTaskCMXQA &other);    
  ClassDef(AliAnalysisTaskCMXQA, 1) 

};

#endif
