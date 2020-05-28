/*
Author: Vytautas Vislavicius
Extention of Generic Flow (https://arxiv.org/abs/1312.3572)
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
class AliGFWCuts {
 public:
  AliGFWCuts();
  ~AliGFWCuts();
  Int_t AcceptParticle(AliVParticle*, Int_t BitShift=0, Double_t ptLow=-1, Double_t ptHigh=-1);
  Int_t AcceptVertex(AliAODEvent*, Int_t BitShift=0);
  Int_t AcceptTrack(AliAODTrack*, Double_t*, Int_t BitShift=0);
  void ResetCuts();
  void PrintSetup();
  void SetupCuts(Int_t);
  TString *GetFlagDescription(Int_t flag);
  const char *GetSystPF() { return Form("%s",fSystFlag?Form("_SystFlag%i_",fSystFlag):""); };
  Int_t fSystFlag;
  Int_t fFilterBit;//=96;
  Double_t fDCAxyCut;//=-1;
  Double_t fDCAzCut;//=2;
  Int_t fTPCNcls;//=70;
  Double_t fVtxZ;// = 10
  Double_t fEta;
  static const Int_t fNTrackFlags;
  static const Int_t fNEventFlags;
  Bool_t NeedsExtraWeight() { return fRequiresExtraWeight; };
 private:
  Bool_t fRequiresExtraWeight;
  void SetupTrackCuts(Int_t);
  void SetupEventCuts(Int_t);
  TString *GetTrackFlagDescriptor(Int_t sysflag);
  TString *GetEventFlagDescriptor(Int_t sysflag);
};
#endif
