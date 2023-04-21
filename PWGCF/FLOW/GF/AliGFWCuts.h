/*
Author: Vytautas Vislavicius
Contains the additional event and track selection used within the <AliGFW> framework.
If used, modified, or distributed, please aknowledge the original author of this code.
*/
#ifndef ALIGFWCUTS__H
#define ALIGFWCUTS__H
#include "AliAODEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "TMath.h"
#include "TString.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "TF1.h"

class AliGFWCuts {
 public:
  AliGFWCuts();
  ~AliGFWCuts();
  Int_t AcceptParticle(AliVParticle*, Int_t BitShift=0, Double_t ptLow=-1, Double_t ptHigh=-1);
  Int_t AcceptVertex(AliAODEvent*, Int_t BitShift=0);
  Int_t AcceptVertex(AliESDEvent*, Int_t BitShift=0);
  Int_t AcceptTrack(AliAODTrack*&, Double_t*, const Int_t &BitShift=0, const Bool_t &lDisableDCAxyCheck=kTRUE);
  Int_t AcceptTrack(AliESDtrack*&, Double_t*, const Int_t &BitShift, UInt_t &PrimFlags);
  void SetPtDepDCAXY(TString newval) { if(fPtDepXYCut) delete fPtDepXYCut; fPtDepXYCut = new TF1("ptDepDCAxy",newval.Data(),0.001,100); fPtDepXYCut->SetParameter(0,fDCAxyCut); };
  void ResetCuts();
  void PrintSetup();
  void SetupCuts(Int_t);
  void SetEta(Double_t newval) { fEta = newval; };
  TString *GetFlagDescription(Int_t flag);
  const char *GetSystPF() { return Form("%s",fSystFlag?Form("_SystFlag%i_",fSystFlag):""); };
  Int_t GetSystFlagIndex() { return fSystFlag; };
  Double_t GetChi2TPCConstrained(const AliESDtrack *);
  Int_t fSystFlag;
  Int_t fFilterBit;//=96;
  Double_t fDCAxyCut;//=-1;
  Double_t fDCAzCut;//=2;
  Int_t fTPCNcls;//=70;
  Double_t fTPCChi2PerCluster;//= 4 or 2.5?
  Double_t fVtxZ;// = 10
  Double_t fEta;
  static const Int_t fNTrackFlags;
  static const Int_t fNEventFlags;
  static AliESDtrackCuts *fTCFB32;
  static AliESDtrackCuts *fTCFB64;
  static AliESDtrackCuts *fTCFB256;
  static AliESDtrackCuts *fTCFB512;
  static AliESDtrackCuts *fTCFB16;
  Bool_t NeedsExtraWeight() { return fRequiresExtraWeight; };
 private:
  TF1 *fPtDepXYCut;
  Bool_t fRequiresExtraWeight;
  void SetupTrackCuts(Int_t);
  void SetupEventCuts(Int_t);
  TString *GetTrackFlagDescriptor(Int_t sysflag);
  TString *GetEventFlagDescriptor(Int_t sysflag);
};
#endif
