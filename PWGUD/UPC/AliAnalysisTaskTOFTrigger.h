/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKTOFTRIGGER_H
#define ALIANALYSISTASKTOFTRIGGER_H

class TList;
class AliTOFTriggerMask;
class TEfficiency;
class TH2F;
class AliESDtrackCuts;
class TGeoMatrix;

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

class AliAnalysisTaskTOFTrigger : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskTOFTrigger();
  AliAnalysisTaskTOFTrigger(const char *name,Float_t lowpt,Float_t highpt,Int_t highmult,TString trgcls,Int_t nBCs,Bool_t useEVS,Int_t fTrackCutSet,Float_t fMaxTrackError,Float_t mintof,Float_t maxtof);
  virtual ~AliAnalysisTaskTOFTrigger();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  void SetupParameters(Float_t lowpt,Float_t highpt,Int_t highmult,TString trgcls,Int_t nBCs,Bool_t useEVS,Int_t cutSet,Int_t maxErr,Float_t mintof,Float_t maxtof)
  			{fMaxPt = highpt; fMinPt = lowpt; fMaxMulti = highmult; fTriggerClass = trgcls;fMaxBCs = nBCs;fUseEventSelection = useEVS;fTrackCutSet = cutSet; fMaxTrackError = maxErr; fMinTOF = mintof; fMaxTOF = maxtof;}
  void GetLTMIndex(const Int_t * const detind, Int_t *indexLTM);

 private:

  TList *fOutputList;		//<

  AliPIDResponse *fPIDResponse;
  AliESDtrackCuts *fTrackCuts;
  AliTOFTriggerMask *fTOFmask;
  TEfficiency *eff_MaxiPadLTM_All;		//!
  TEfficiency *eff_MaxiPadLTM_Mu;		//!
  TEfficiency *eff_MaxiPadLTM_El;		//!
  TEfficiency *eff_MaxiPadLTM_1Trk_All;		//!
  TEfficiency *eff_MaxiPadLTM_1Trk_Mu;		//!
  TEfficiency *eff_MaxiPadLTM_1Trk_El;		//!
  TEfficiency *eff_AverageTracklets;		//!
  TEfficiency *eff_AverageTrackPt;		//!
  TEfficiency *eff_MaxiPadLTM_Around;		//!
  TEfficiency *eff_MaxiPadLTM_OnlyAround;	//!
  TEfficiency *eff_MaxiPadLTM_Clusters;		//!
  TH2F *hTrackDistributionLTM;			//!
  TH2F *hTrackDistribution_Mu;			//!
  TH2F *hTrackDistribution_El;			//!
  TH2F *hTrackDistribution;			//!
  TH2F *hFiredMaxiPad;				//!
  TH2F *hFiredMaxiPadOnlyAround;		//!
  TH2F *hNotFiredMaxiPadCls;			//!
  TH2F *hExtraFiredMaxiPadCls;			//!
  TH2F *hNotFiredMaxiPadTrk;			//!
  TH2F *hExtraFiredMaxiPadTrk;			//!
  TH2F *hTrackPadCorrPhi;			//!
  TH2F *hTrackPadCorrEta;			//!
  TH2F *hNoiseMaxiPad;				//!
  TH1I *hTriggerCounter;			//!
  TH1I *hTriggerCounterIR1;			//!
  TH1I *hTriggerCounterIR2;			//!
  TH1F *hNFiredMaxiPads;			//!
  TH1F *hNFiredMaxiPadsOnlyAround;		//!
  TH1F *hNTracklets;				//!
  TH1I *hDetIn0;				//!
  TH1I *hDetIn1;				//!
  TH1I *hDetIn2;				//!
  TH1I *hDetIn3;				//!
  TH1I *hDetIn4;				//!
  TH1F *hPadDistance;				//!
  TH1F *hTrackPt;				//!
  TH1I *hNMaxiPadIn;				//!
  TH1I *hNCrossTracks;				//!
  TH2I *hBadMaxiPadMask;			//!
  TH1F *hTOFHitTime;				//!

  Bool_t fGeomLoaded;
  TGeoHMatrix matOrig[18]; 
  TGeoHMatrix matCurr[18];
  Float_t fMaxPt;
  Float_t fMinPt;
  Int_t fMaxMulti;
  TString fTriggerClass;
  Int_t fMaxBCs;
  Bool_t fUseEventSelection;
  Int_t fTrackCutSet;
  Float_t fMaxTrackError;
  Float_t fMinTOF;
  Float_t fMaxTOF;
  
  
  AliEventCuts fEventCuts;	


  AliAnalysisTaskTOFTrigger(const AliAnalysisTaskTOFTrigger&); //not implemented
  AliAnalysisTaskTOFTrigger& operator =(const AliAnalysisTaskTOFTrigger&); //not implemented

  ClassDef(AliAnalysisTaskTOFTrigger, 18);
};

#endif
