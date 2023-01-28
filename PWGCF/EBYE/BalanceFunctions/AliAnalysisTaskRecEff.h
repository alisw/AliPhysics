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

#ifndef ALIANALYSISTASKRECEFF_H
#define ALIANALYSISTASKRECEFF_H

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
class    AliPIDResponse;    
class    AliMultSelection;    
class    AliAnalysisUtils;





class AliAnalysisTaskRecEff : public AliAnalysisTaskSE {

 public:

  //-----> Mandatory Functions:
  AliAnalysisTaskRecEff();
  AliAnalysisTaskRecEff(const char *name);
  virtual ~AliAnalysisTaskRecEff();
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
  void SetDCAXYRangeMax(Double_t dcaxy)             {this->fDCAxyMax    = dcaxy;}
  void SetDCAZRangeMax(Double_t dcaz)              {this->fDCAzMax    =  dcaz;}
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
  Float_t               fDCAxyMax;  //
  Float_t               fDCAzMax;  //
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
  TH1F                  *fCentDistBeforCut;    //! Cent before Any Cut
  TH1F                  *fCentDistAfterCut;    //! Cent After All Cuts
  TH2F                  *fHistEtaPtBeforCut;   //! Eta-Pt before Any Cut
  TH2F                  *fHistEtaPhiBeforCut;  //! Eta-Phi before Any Cut
  TH2F                  *fHistEtaPhiAfterCut;  //! Eta-Phi before Any Cut
    
  //Analysis Result Histograms:
  TH1F                  *fHistPtDistMBwithCut; //! Pt after All Cuts


  //Add Your stuff here: Rihan
 
  TH2F     *fHistPtRecChrgPos;  //!
  TH2F     *fHistPtGenChrgPos;  //!
  TH2F     *fHistPtRecPionPos;  //!
  TH2F     *fHistPtGenPionPos;  //!
  TH2F     *fHistPtRecKaonPos;  //!
  TH2F     *fHistPtGenKaonPos;  //!
  TH2F     *fHistPtRecProtPos;  //!
  TH2F     *fHistPtGenProtPos;  //!

  TH2F     *fHistPtRecChrgNeg;  //!
  TH2F     *fHistPtGenChrgNeg;  //!
  TH2F     *fHistPtRecPionNeg;  //!
  TH2F     *fHistPtGenPionNeg;  //!
  TH2F     *fHistPtRecKaonNeg;  //!
  TH2F     *fHistPtGenKaonNeg;  //!
  TH2F     *fHistPtRecProtNeg;  //!
  TH2F     *fHistPtGenProtNeg;  //!

  
  TH2F     *fHistPtRecPionPosimp;  //!
  TH2F     *fHistPtRecKaonPosimp;  //!
  TH2F     *fHistPtRecProtPosimp;  //!
  TH2F     *fHistPtRecPionNegimp;  //!
  TH2F     *fHistPtRecKaonNegimp;  //!
  TH2F     *fHistPtRecProtNegimp;  //!
  
  


  TH2F     *fHistnSigmaTPCPionCent1;  //!
  TH2F     *fHistnSigmaTPCKaonCent1;  //!
  TH2F     *fHistnSigmaTPCProtCent1;  //!
  TH2F     *fHistnSigmaTOFPionCent1;  //!
  TH2F     *fHistnSigmaTOFKaonCent1;  //!
  TH2F     *fHistnSigmaTOFProtCent1;  //!
  
  // contamination primary and secondary
  TH2F     *fHistPtContaminationPrimariesChrgPos; //!
  TH2F     *fHistPtContaminationPrimariesChrgNeg; //!
  TH2F     *fHistPtContaminationSecondariesChrgPos; //!
  TH2F     *fHistPtContaminationSecondariesChrgNeg; //!
  TH2F     *fHistPtContaminationSecondariesMaterialChrgPos; //!
  TH2F     *fHistPtContaminationSecondariesMaterialChrgNeg; //!
  TH2F     *fHistPtContaminationSecondariesWeakDecChrgPos; //!
  TH2F     *fHistPtContaminationSecondariesWeakDecChrgNeg; //!



  TH2F     *fMCVzDistBeforCut;  //! MC track Vz
  
  Int_t      fParticlePdgCode;  //! particle PDG code
  TH2F     *fHistAChrggenVsCent; //!
  TH2F     *fHistAChrgrecVsCent; //!


  TH2F            *fHistTPConlyVsCL1Before;   //!  
  TH2F            *fHistTPConlyVsV0MBefore;   //!
  TH2F            *fHistCL0VsV0MBefore;   //!       
  TH2F            *fHistTPConlyVsCL1After;   //!    
  TH2F            *fHistTPConlyVsV0MAfter;   //!   
  TH2F            *fHistCL0VsV0MAfter;   //!   
  TH2F            *fHistGlobalVsV0MBefore;   //!   
  TH2F            *fHistGlobalVsV0MAfter;   //
  TH2F            *fTPCvsGlobalTrkBefore; //!  Global vs TPC tracks for QA
  TH2F            *fTPCvsGlobalTrkAfter; //!  Global vs TPC tracks for QA
  TH2F            *fHistTPCVsESDTrkBefore;   //!  
  TH2F            *fHistTPCVsESDTrkAfter;   //
  TH1F            *fHistPileUpCount;   //!


  TH2D *fHistPIDDecayMotherDaughterPhysicalPrimary; //! Mother vs. Daughter
  TH2D *fHistPIDDecayMotherDaughterNotPhysicalPrimary; //! Mother vs. Daughter
  TH1D *fHistPIDDecayDaughterNoMotherPhysicalPrimary; //! Daughter distribution
  TH1D *fHistPIDDecayDaughterNoMotherNotPhysicalPrimary; //! Daughter distribution
  TH1D *fQARecoEventCounter; //!
  TH1D *fQAGeneEventCounter; //!
  
  //Last Hist in the List:
  TH1F *fHistEventCount;      //!   

  TProfile *fHistAch[10]; //!
 
  
    
  Bool_t CheckEventIsPileUp2018(AliAODEvent* faod);

  AliAnalysisTaskRecEff(const AliAnalysisTaskRecEff &other);
  AliAnalysisTaskRecEff &operator=(const AliAnalysisTaskRecEff &other);    
  ClassDef(AliAnalysisTaskRecEff, 1) 

};

#endif
