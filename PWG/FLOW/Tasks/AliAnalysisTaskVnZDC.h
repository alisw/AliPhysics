/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskVnZDC_H
#define AliAnalysisTaskVnZDC_H

/////////////////////////////////////////////////
// AliAnalysisTaskVnZDC:
// analysis task for Scalar Product method
// Author: Naomi van der Kolk (kolk@nikhef.nl)
/////////////////////////////////////////////////

class AliAODEvent;
class AliVVertex;
class AliFlowEventSimple;
class AliMultSelection;
class AliAnalysisUtils;
class TList;


#include "TString.h"
#include "AliAnalysisTaskSE.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TTree.h"
#include "TList.h"


//===============================================================

class AliAnalysisTaskVnZDC : public AliAnalysisTaskSE {

public:
  AliAnalysisTaskVnZDC();
  AliAnalysisTaskVnZDC(const char *name, Bool_t usePhiWeights);
  virtual ~AliAnalysisTaskVnZDC();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void   SetUsePhiWeights(Bool_t const aPhiW) {this->fUsePhiWeights = aPhiW;}
  Bool_t GetUsePhiWeights() const             {return this->fUsePhiWeights;}

  void     SetRelDiffMsub(Double_t diff) { this->fRelDiffMsub = diff; }
  Double_t GetRelDiffMsub() const        { return this->fRelDiffMsub; }

  void   SetApplyCorrectionForNUA(Bool_t const applyCorrectionForNUA) {this->fApplyCorrectionForNUA = applyCorrectionForNUA;}
  Bool_t GetApplyCorrectionForNUA() const {return this->fApplyCorrectionForNUA;}

  void  SetHarmonic(Int_t const harmonic) {this->fHarmonic = harmonic;}
  Int_t GetHarmonic() const {return this->fHarmonic;};
  void  SetBehaveAsEP() { fNormalizationType = 0; }
  void  SetTotalQvector(const char *tqv) { *this->fTotalQvector = tqv; }
  void  SetBookOnlyBasicCCH(Bool_t const aMB) { this->fMinimalBook = aMB; }

  //rihan:
  //setters:
  void    SetRunFlag(Int_t const runnum)              {this->frunflag = runnum; }
  void    SetZDCESEList(TList* const kList)           {this->fZDCESEList = kList;}
  void    SetDataSet(TString fdataset)                {this->fDataSet = fdataset;}
  void    SetAnalysisSet(TString fanalysisSet)        {this->fAnalysisSet = fanalysisSet;}
  void    SetRejectPileUpTight(Bool_t const pileupt8) {this->fRejectPileUpTight = pileupt8;}
  void    SetRejectPileUp(Bool_t const pileup)        {this->fRejectPileUp = pileup;}
  //getters:
  /*
  TString GetDataSet()                    {return this->fDataSet;}
  TString GetAnalysisSet()                {return this->fAnalysisSet;}
  Bool_t  GetRejectPileUpTight() const    {return this->fRejectPileUpTight;}
  Bool_t  GetRejectPileUp()      const    {return this->fRejectPileUp;}
  */

  Bool_t    plpMV(const AliAODEvent* aod);
  double GetWDist(const AliVVertex* v0, const AliVVertex* v1);



private:

  AliAnalysisTaskVnZDC(const AliAnalysisTaskVnZDC& aAnalysisTask);
  AliAnalysisTaskVnZDC& operator=(const AliAnalysisTaskVnZDC& aAnalysisTask);

  AliFlowEventSimple* fEvent;         //input event

//AliFlowAnalysisIDCSP** fSP;         // analysis object
//AliFlowAnalysisIDCSP*  fSP;         // analysis object only one.
  Int_t              fievent;
  TList         *fListHistos;         // collection of output
  TList        *fListWeights;         // list with weights
  Bool_t        fMinimalBook;         // flag to turn off QA and minimize FlowCommonHist
  Bool_t      fUsePhiWeights;         // use phi weights
  Int_t            fHarmonic;         // harmonic
  Int_t   fNormalizationType;         // 0: EP mode || 1: SP mode (default)
  Int_t           fNCentBins;         // Number of Cent bins
  Double_t      fRelDiffMsub;         // the relative difference the two subevent multiplicities can have
  TString     *fTotalQvector;         // total Q-vector is: "QaQb" (means Qa+Qb), "Qa"  or "Qb"
  Bool_t fApplyCorrectionForNUA;      // apply automatic correction for non-uniform acceptance

  //rihan:
  Int_t             frunflag;  //
  Int_t          runNums[89];  //
  Float_t           VxCut[2];  //
  Float_t           VyCut[2];  //
  Float_t           VzCut[2];  //
  Int_t                vxBin;  //
  Int_t                vyBin;  //
  Int_t                vzBin;  //
  Int_t            checkOnce;  //


  TString                   fDataSet;    // Dataset: 2010, 2011, or 2015.
  TString               fAnalysisSet;    // Values: recenter1,recenter2,analysis1
  Bool_t               fRejectPileUp;    //
  Bool_t          fRejectPileUpTight;    //
  AliMultSelection*   fMultSelection;    // MultSelection (RUN2 centrality estimator)
  AliAnalysisUtils*    fAnalysisUtil;    // < Event selection



  TProfile   *fHist_znCx_V0_Cent[89];    //!
  TProfile   *fHist_znCy_V0_Cent[89];    //!
  TProfile   *fHist_znAx_V0_Cent[89];    //!
  TProfile   *fHist_znAy_V0_Cent[89];    //!
  TList                 *fZDCESEList;    //

  TProfile2D  *fHist_ZDCC_En_Run[89];    //!
  TProfile2D  *fHist_ZDCA_En_Run[89];    //!

  TH2F          *fHist_ZDCn_A_XY[89];    //!
  TH2F          *fHist_ZDCn_C_XY[89];    //!
  TProfile2D        *fHist_QnxRecent;    //!
  TProfile2D        *fHist_QnyRecent;    //!
  TH1F              *fHist_vBincount;    //!
  TH1F              *fHist_Psi1_zdnA;    //!
  TH1F              *fHist_Psi1_zdnC;    //!
  TH1F            *fHist_EventperRun;    //!
  TH1F              *fHist_Vertex_Vz;    //!
  TH1F            *fHist_Event_count;    //!
  TH1F            *fHist_Cent_count1;    //!
  TH1F            *fHist_Cent_count2;    //!
  TH2F              *fHist_Vertex_XY;    //!
  TProfile       *fHist_Vx_vs_runnum;    //!
  TProfile       *fHist_Vy_vs_runnum;    //!
  TProfile       *fHist_Vz_vs_runnum;    //!
  TProfile   *fHist_tracks_vs_runnum;    //!
  TH1F                 *fPileUpCount;    //!
  TH1F          *fPileUpMultSelCount;    //!







//These are non written (non-listed) histograms. must be deleted in the constructor
  TH1F         *fHist_Vx_ArrayFinder;    //!
  TH1F         *fHist_Vy_ArrayFinder;    //!
  TH1F         *fHist_Vz_ArrayFinder;    //!
  TProfile2D *fHist_znCx_V0_VxVy[89];    //!
  TProfile2D *fHist_znCy_V0_VxVy[89];    //!
  TProfile2D *fHist_znAx_V0_VxVy[89];    //!
  TProfile2D *fHist_znAy_V0_VxVy[89];    //!




  ClassDef(AliAnalysisTaskVnZDC, 1); // example of analysis
};

//==================================================================

#endif
