/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskZDCGainEq_H
#define AliAnalysisTaskZDCGainEq_H

/////////////////////////////////////////////////
// AliAnalysisTaskZDCGainEq:
// analysis task for ZDC gain Equalization
// Author: Rihan Haque (mhaque@cern.ch)
/////////////////////////////////////////////////

class AliAODEvent;
class AliVVertex;
class AliFlowEventSimple;
class AliMultSelection;
class AliAnalysisUtils;
class TList;


#include "TString.h"
#include "AliAnalysisTaskSE.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TTree.h"
#include "TList.h"


//===============================================================

class AliAnalysisTaskZDCGainEq : public AliAnalysisTaskSE {

public:
  AliAnalysisTaskZDCGainEq();
  AliAnalysisTaskZDCGainEq(const char *name);
  virtual ~AliAnalysisTaskZDCGainEq();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  //void     SetRelDiffMsub(Double_t diff) { this->fRelDiffMsub = diff; }
  //Double_t GetRelDiffMsub() const        { return this->fRelDiffMsub; }
  //void   SetApplyCorrectionForNUA(Bool_t const applyCorrectionForNUA) {this->fApplyCorrectionForNUA = applyCorrectionForNUA;}
  //Bool_t GetApplyCorrectionForNUA() const {return this->fApplyCorrectionForNUA;}

  //
  void  SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;}
  Int_t GetHarmonic()     const           {return     this->fHarmonic;}

  //rihan:
  void    SetRunFlag(Int_t const runnum)              {this->frunflag           =       runnum;}
  void    SetZDCESEList(TList* const kList1)          {this->fListZDCQxy        =       kList1;}
  void    SetFBEffiList(TList* const kList2)          {this->fListHijing        =       kList2;}
  void    SetZDCChWgtList(TList* const kList3)        {this->fListZDCWgt        =       kList3;}

  void    SetDataSet(TString fdataset)                {this->fDataSet           =     fdataset;}
  void    SetAnalysisSet(TString fanalysisSet)        {this->fAnalysisSet       = fanalysisSet;}
  void    SetCentEstimator(TString centEstim)         {this->sCentEstimator     =    centEstim;}

  void    SetRejectPileUpTight(Bool_t const pileupt8) {this->fRejectPileUpTight =     pileupt8;}
  void    SetRejectPileUp(Bool_t const pileup)        {this->fRejectPileUp      =       pileup;}
  void    SetFillCosSin(Bool_t const fillcossin)      {this->bFillCosSin        =   fillcossin;}
  void    SetFillZDCQA(Bool_t const fillzdcQAon)      {this->bFillZDCQAon       =  fillzdcQAon;}
  void    SetFillQnRunAverage(Bool_t const runAvg)    {this->bRunAveragedQn     =       runAvg;}
  void    SetApplyZDCRecenter(Bool_t const brecent)   {this->bApplyRecent       =      brecent;}



  void    SetVxVyVzBinForAvgQx(Int_t vxbin,Int_t vybin,Int_t vzbin){
           this->vxBin =  vxbin;   
           this->vyBin =  vybin;   
           this->vzBin =  vzbin;
          }


  Bool_t                        plpMV(const AliAODEvent* aod);
  double GetWDist(const AliVVertex* v0, const AliVVertex* v1);



private:

  AliAnalysisTaskZDCGainEq(const AliAnalysisTaskZDCGainEq& aAnalysisTask);
  AliAnalysisTaskZDCGainEq& operator=(const AliAnalysisTaskZDCGainEq& aAnalysisTask);

  AliFlowEventSimple*         fEvent;         //! input event
  AliMultSelection*   fMultSelection;         // MultSelection (RUN2 centrality estimator)
  AliAnalysisUtils*    fAnalysisUtil;         // < Event selection


  TList*                 fListHistos;         //! collection of output
  TList*                 fListZDCQxy;         //
  TList*                 fListZDCWgt;         //
  TList*                 fListDummy1;         //!
  TList*                 fListHijing;         //

  TList*                 fListSubRun;

  Int_t                    fHarmonic;    // harmonic


  TString                   fDataSet;    // Dataset: 2010, 2011, or 2015.
  TString               fAnalysisSet;    // Values: recenter1,recenter2,analysis1
  TString             sCentEstimator;
  Bool_t               fRejectPileUp;    //
  Bool_t          fRejectPileUpTight;    //
  Bool_t                 bFillCosSin;    //
  Bool_t                bFillZDCQAon;    //
  Bool_t              bRunAveragedQn;    //
  Bool_t                bApplyRecent;    //
  Int_t                  runNums[90];    //
  Float_t                   VxCut[2];    //
  Float_t                   VyCut[2];    //
  Float_t                   VzCut[2];    //

  Int_t                     frunflag;  //!
  Int_t                      fievent;  //!
  Int_t                   fcheckOnce;  //!
  Int_t                   fOldRunNum;  //!
  Int_t                        vxBin;  //
  Int_t                        vyBin;  //
  Int_t                        vzBin;  //


  TH1F            *fHist_Event_count;  //!
  TH1F                 *fPileUpCount;  //!
  TH1F          *fPileUpMultSelCount;  //!
  TH1F            *fHist_Task_config;  //!
  TH1F          *fHist_CutParameters;  //!
  TH1F          *fHist_Cent_woZDCcut;  //!
  TH1F          *fHist_Cent_wiZDCcut;  //!

  TH1F           *fHist_ChanWgt_ZDCC;  //!
  TH1F           *fHist_ChanWgt_ZDCA;  //!

  TH2F           *fHist_Recenter_ZDCCx[20]; //!
  TH2F           *fHist_Recenter_ZDCCy[20]; //!
  TH2F           *fHist_Recenter_ZDCAx[20]; //!
  TH2F           *fHist_Recenter_ZDCAy[20]; //!

  TH1F               *fHist_Vx_ArrayFinder; //!
  TH1F               *fHist_Vy_ArrayFinder; //!
  TH1F               *fHist_Vz_ArrayFinder; //!

  TH1F                       *fWeight_Cent; //!
  TH1D            *fFB_Efficiency_Cent[10]; //!

  TH1F            *fHist_Psi1_ZDCC_wGainCorr;  //!
  TH1F            *fHist_Psi1_ZDCA_wGainCorr;  //!
  TH1F            *fHist_Psi1_ZDCC_wRectCorr;  //!
  TH1F            *fHist_Psi1_ZDCA_wRectCorr;  //!

  TProfile     *fHist_Qx_vs_Obs_woCorr[4][5];  //!
  TProfile     *fHist_XX_vs_Obs_woCorr[4][5];  //!
  TProfile     *fHist_Qx_vs_Obs_wiCorr[4][5];  //!
  TProfile     *fHist_XX_vs_Obs_wiCorr[4][5];  //!

  TProfile         *fHist_v2xV1_ZDN_Norm_All;  //!
  TProfile         *fHist_v2xV1_ZDN_Refm_All;  //!
  TProfile         *fHist_v2xV1_ZDN_Cent_All;  //!

  TProfile         *fHist_ZDN_resol_Norm_All;  //!
  TProfile         *fHist_ZDN_resol_Refm_All;  //!
  TProfile         *fHist_ZDN_resol_Cent_All;  //!

  TProfile       *fHist_v2xV1_ZDN_Norm_Sep[4];  //! 
  TProfile       *fHist_v2xV1_ZDN_Cent_Sep[4];  //!
  TProfile       *fHist_ZDN_resol_Norm_Sep[2];  //!
  TProfile       *fHist_ZDN_resol_Cent_Sep[2];  //!

  TProfile   *fHist_v2xV1_ZDN_pTDiff_All[10];  //!
  TProfile    *fHist_v1xV1_ZDN_pTDiff[4][10];  //!
  TProfile   *fHist_v1xV1_ZDN_EtaDiff[4][10];  //!

  TProfile               *fHist_Vx_vs_runnum;  //!
  TProfile               *fHist_Vy_vs_runnum;  //!
  TProfile               *fHist_Vz_vs_runnum;  //!

  TProfile2D     *fHist_znCx_V0_VxVy[90][10];  //!
  TProfile2D     *fHist_znCy_V0_VxVy[90][10];  //!
  TProfile2D     *fHist_znAx_V0_VxVy[90][10];  //!
  TProfile2D     *fHist_znAy_V0_VxVy[90][10];  //!

  TProfile2D  *fHist_ZDCC_En_Run[90];    //!
  TProfile2D  *fHist_ZDCA_En_Run[90];    //!

  TProfile2D  *fHist_ZDCC_En_CommonCh[20];   //!
  TProfile2D  *fHist_ZDCA_En_CommonCh[20];   //!

  ClassDef(AliAnalysisTaskZDCGainEq, 3); // example of analysis
};

//==================================================================

#endif
