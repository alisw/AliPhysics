/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id: $ */

////////////////////////////////////////////////
// AliAnalysisTaskCVE:
// Simple CVE AnalysisTask
// PA: Rihan Haque (mhaque@cern.ch, rihanphys@gmail.com)
// Date: Jan 13, 2020.
// Last Mod: Jan 13, 2020.
/////////////////////////////////////////////////



#ifndef ALIANALYSISTASKv2vscent_H
#define ALIANALYSISTASKv2vscent_H

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





class AliAnalysisTaskv2vscent : public AliAnalysisTaskSE {

 public:

  //-----> Mandatory Functions:
  AliAnalysisTaskv2vscent();
  AliAnalysisTaskv2vscent(const char *name);
  virtual ~AliAnalysisTaskv2vscent();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * /*option*/);
  
  //-----> User Defined Functions:
  Int_t GetCentralityScaled0to10(Double_t fCent);
  void  SetupEventAndTaskConfigInfo();
  void  SetListForTrkCorr(TList *flist)      {this->fListTRKCorr = (TList *) flist->Clone(); }
  void  SetListForNUACorr(TList *flist)      {this->fListNUACorr = (TList *) flist->Clone(); }
  void  SetListForV0MCorr(TList *flist)      {this->fListV0MCorr = (TList *) flist->Clone(); }
  void  SetParticle(Int_t part)              {this->fParticle    =  part;}
  void  SetCumulantHarmonic(Int_t harm)      {this->gHarmonic    =  harm;}

  
  /******* Event Cut Ranges ******/
  void SetCentralityPercentileMin(Double_t centMin) {this->fCentralityMin  = centMin;}
  void SetCentralityPercentileMax(Double_t centMax) {this->fCentralityMax  = centMax;}
  void SetCentralityEstimator(TString sEstim)       {this->sCentrEstimator = sEstim;}
  void SetVzRangeMin(Double_t vzMin)                {this->fMinVzCut       =  vzMin;}
  void SetVzRangeMax(Double_t vzMax)                {this->fMaxVzCut       =  vzMax;}

  void SetPileUpCutParam(Float_t m,Float_t c) {this->fPileUpSlopeParm = m;  this->fPileUpConstParm = c;}
  //void SetFlagSkipPileUpCuts(Bool_t b)        {this->bSkipNUA  = b;}
  void SetDataset(Int_t b)        {this->bdataset  = b;}
  

  /******* Track Cut Ranges ******/
  void SetNSigmaCutTPC(Double_t     nSigTPC)     {this->fNSigmaTPCCut  =  nSigTPC;}
  void SetNSigmaCutTOF(Double_t     nSigTOF)     {this->fNSigmaTOFCut  =  nSigTOF;}
  void SetTrackCutdEdxMin(Float_t  ndEdxMin)     {this->fdEdxMin    = ndEdxMin;}
  void SetTrackCutChi2Min(Double_t  chi2Min)     {this->fTrkChi2Min =  chi2Min;}
  void SetFlagUseKinkTracks(Bool_t    bKink)     {this->bUseKinkTracks = bKink;}
  
  void SetTrackCutNclusterMin(Int_t ncl)         {this->fTPCclustMin = ncl;}
  void SetFilterBit(Int_t fb)                    {this->fFilterBit   =  fb;}
  void SetEtaRangeMin(Double_t emn)              {this->fMinEtaCut   = emn;}
  void SetEtaRangeMax(Double_t emx)              {this->fMaxEtaCut   = emx;}
  void SetPtRangeMin(Double_t ptL)               {this->fMinPtCut    = ptL;}
  void SetPtRangeMax(Double_t ptH)               {this->fMaxPtCut    = ptH;}
  void SetDCAXYRangeMax(Double_t dcaxy)          {this->fDCAxyMax    = dcaxy;}
  void SetDCAZRangeMax(Double_t dcaz)            {this->fDCAzMax    =  dcaz;}
  void SetChi2Range(Double_t chi2)               {this->fChi2    =  chi2;}
  void SetEtaNeg(Double_t etaL)                  {this->fEtaGapNeg   = etaL;}
  void SetEtaPos(Double_t etaH)                  {this->fEtaGapPos   = etaH;}

  //------ End of user defined function -------







  
 protected:

 private:


  AliVEvent             *fVevent;             //! event
  AliESDEvent           *fESD;                //! esd Event
  AliAODEvent           *fAOD;                //! aod Event
  AliPIDResponse        *fPIDResponse;        //! PID Handler
  AliMultSelection      *fMultSelection;      //! For Centrality 
  AliAnalysisUtils      *fAnalysisUtil;       //! Event Selection Options
  TList                 *fListHist;           //! OutputList

  TList                 *fListTRKCorr;        //  Supplied from Task
  TList                 *fListNUACorr;        //  Supplied from Task
  TList                 *fListV0MCorr;        //  Supplied from Task

  //Set from the AddTask:
  Float_t         fCentralityMin;  //
  Float_t         fCentralityMax;  //
  
  //Track Variables to be used:
  Int_t                 gHarmonic;  //
  Int_t                 fParticle;  //
  Int_t                fFilterBit;  //
  Int_t              fTPCclustMin;  //
  Bool_t           bUseKinkTracks;  //
  
  Float_t           fNSigmaTPCCut;  //
  Float_t           fNSigmaTOFCut;  //
  Float_t               fMinPtCut;  //
  Float_t               fMaxPtCut;  //
  Float_t               fDCAxyMax;  //                                                                                                        
  Float_t               fDCAzMax;  // 
  Float_t                 fChi2; 
  Double_t         fPileUpSlopeParm;  //
  Double_t         fPileUpConstParm;  //
  //Bool_t             bSkipNUA;  //
  Int_t             bdataset;  //
  Double_t              fEtaGapNeg;  //
  Double_t              fEtaGapPos;  //
  Float_t              fMinEtaCut;  //
  Float_t              fMaxEtaCut;  //
  Float_t             fTrkChi2Min;  //
  Float_t                fdEdxMin;  // 
  //Event Variables to be used:
  Float_t               fMinVzCut;  //
  Float_t               fMaxVzCut;  //
  
  TString         sCentrEstimator;  //


  

  /////////// Default HISTOGRAMS /////////
  
  //QA histograms:
  TH1F                  *fCentDistBeforCut;    //! Cent before Any Cut
  TH1F                  *fCentDistAfterCut;    //! Cent After All Cuts
 


  
  //User Defined Histograms: 


  ///Used For Tracking Efficiency:
  TH1D          *fHCorrectMCposChrg;
  TH1D          *fHCorrectMCposPion;
  TH1D          *fHCorrectMCposKaon;
  TH1D          *fHCorrectMCposProt;
  TH1D          *fHCorrectMCnegChrg;
  TH1D          *fHCorrectMCnegPion;
  TH1D          *fHCorrectMCnegKaon;
  TH1D          *fHCorrectMCnegProt;

  TH3F          *fHCorrectNUAposChrg;   //!  = centrality bins
  TH3F          *fHCorrectNUAnegChrg;   //! 
  TH3F          *fHCorrectNUAposPion;   //! 
  TH3F          *fHCorrectNUAnegPion;   //! 
  TH3F          *fHCorrectNUAposKaon;   //! 
  TH3F          *fHCorrectNUAnegKaon;   //! 
  TH3F          *fHCorrectNUAposProt;   //! 
  TH3F          *fHCorrectNUAnegProt;   //! 

  TProfile      *fHistv2cumCentChrgAll;
  TProfile      *fHistv2cumCentChrgPos;
  TProfile      *fHistv2cumCentChrgNeg;
  TProfile      *fHistv2cumCentChrgPosChrgNeg;
  TProfile      *fHistv2cumCentChrgNegChrgPos;
  TProfile      *fHistv2cumCentPionPos;
  TProfile      *fHistv2cumCentPionNeg;
  TProfile      *fHistv2cumCentPionPosPionNeg;
  TProfile      *fHistv2cumCentPionNegPionPos;
 


  TF1           *fSPDCutPU;     //!
  TF1           *fV0CutPU;      //!
  TF1           *fMultCutPU;    //!
  TF1           *fCenCutLowPU;  //!
  TF1           *fCenCutHighPU; //!
  



  
  //QA and Stepcount 

  
  TH1F            *fHistEventCount;   //!
  TH1F  	  *fHCorrectEVNTWGTChrg;   //!   //eventwgt for charge



  

  ///Custom Functions:
  void  GetNUACorrectionHist(Int_t run=0,Int_t kParticleID=0);
  void  GetEVNTWGTCorrectionHist(Int_t run=0,Int_t kParticleID=0);
  //void  GetV0MCorrectionHist(Int_t run=0);
  void  GetV0MCorrectionHist(Int_t run=0,Int_t kParticleID=0);
  //void  GetMCCorrectionHist(Int_t run=0);
  void  GetMCCorrectionHist(Int_t run=0,Float_t centr=0);
  
  Bool_t CheckEventIsPileUp(AliAODEvent* faod);
  Bool_t CheckEventIsPileUp2018(AliAODEvent* faod);
  Bool_t PileUpMultiVertex(const AliAODEvent* faod);
  double GetWDist(const AliVVertex* v0, const AliVVertex* v1);


  AliAnalysisTaskv2vscent(const AliAnalysisTaskv2vscent &other);
  AliAnalysisTaskv2vscent &operator=(const AliAnalysisTaskv2vscent &other);    
  ClassDef(AliAnalysisTaskv2vscent, 1) 



};

#endif
