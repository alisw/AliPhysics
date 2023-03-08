#ifndef ALIGFWFILTER__H
#define ALIGFWFILTER__H
#include <vector>
#include <complex>
#include "TChain.h"
#include "TH1D.h"
#include "TList.h"
#include "TMath.h"
#include "TF1.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliAODInputHandler.h"
#include "AliAODVertex.h"
#include "AliEventCuts.h"
#include "GFWFlags.h"
using namespace std;
using namespace GFWFlags;
class AliGFWFlags: public TNamed
{
 public:
  AliGFWFlags():TNamed("GFWFlags","GFWFlags"),fEventFlag(0),fTrFlags({}) {};
  ~AliGFWFlags() {fTrFlags.clear(); };
  void AddTrackFlags(const Int_t &trInd, const UInt_t &lFlag) {fTrFlags.push_back(make_pair(trInd,lFlag)); };
  void SetEventFlags(const UInt_t &lFlag) { fEventFlag=lFlag; };
  void CleanUp() { fEventFlag=0; fTrFlags.clear(); };
  const UInt_t GetEventFlags() { return fEventFlag; };
  const Int_t GetTrackIndex(const Int_t &ind) { return fTrFlags[ind].first; };
  const UInt_t GetTrackFlag(const Int_t &ind) { return fTrFlags[ind].second; };
  const Int_t GetNFiltered() { return (Int_t)fTrFlags.size(); };
  //Keeping track of ev & tr flags. Not useful right now, so disabling
  // void SetNFlags(Int_t lNEv, Int_t lNTr) {fFlagCount = (lNEv<<8)|lNTr;};
  // Int_t GetNEventFlags() {return (fFlagCount&0xff00)>>8;};
  // Int_t GetNTrackFlags() {return (fFlagCount&0xff); };
  Bool_t CheckEventFlag(const UInt_t &lFlag) { return (fEventFlag&lFlag)==lFlag; };
  Bool_t CheckTrackFlag(const Int_t &ind, const UInt_t &lFlag) { return (fTrFlags[ind].second & lFlag)==lFlag; };
  ClassDef(AliGFWFlags,1);
 private:
   UInt_t fEventFlag;
   // char16_t fFlagCount;
   vector< pair <Int_t,UInt_t> > fTrFlags;
};
class AliGFWFilter
{
  public:
    AliGFWFilter();
    ~AliGFWFilter();
    void CleanUp();
    void SetEventCuts(AliEventCuts *fEvCuts) {fEventCuts = fEvCuts; };
    void CheckEvent(AliVEvent *inEv);
    Bool_t AcceptVertex(AliAODEvent *inEv, Double_t *lvtxXYZ);
    void CreateStandardCutMasks(kLocalTrackFlags lStandardChi2Cut = klTPCchi2PC25);
    void AddCustomCuts(Bool_t cleanFirst=kTRUE, UInt_t lEv=klVtxZ10+klEventCuts, UInt_t lTr=klFB96+klDCAz20+klDCAxy2011+klTPCchi2PC25+klNTPCcls70);
    AliGFWFlags *GetFlags() { return fRetFlags; };
    void SetPt(Double_t ptMin, Double_t ptMax) {fPtMin = ptMin; fPtMax = ptMax; };
    void SetEta(Double_t etaMin, Double_t etaMax) { fEtaMin = etaMin; fEtaMax = etaMax; };
    // static const TString GetSystPF(UInt_t ind) { Int_t lev=ind<gNEventFlags?ind:0; Int_t ltr=ind<gNEventFlags?0:(ind-gNEventFlags); return GetSystPF(lev,ltr); };
  private:
    inline Double_t f_DCAxy2010(Double_t &pt) { return (0.0026+0.0050/TMath::Power(pt,1.01)); };
    inline Double_t f_DCAxy2011(Double_t &pt) { return (0.0015+0.0050/TMath::Power(pt,1.1)); };
    void AddEv(const UInt_t &lFlag) {fEventFlag|=lFlag;};
    void AddTr(const UInt_t &lFlag) {flTrFlag|=lFlag;};
    Bool_t TSB(const UInt_t &lFlag) { return (flTrFlag&lFlag); }; //Test signle bit -- track
    Bool_t TB(const UInt_t &lFlag) {return (flTrFlag&lFlag)==lFlag; }; //Track bit
    Bool_t EB(const UInt_t &lFlag) {return (fEventFlag&lFlag)==lFlag; }; //Event bit
    UInt_t calculateEventFlag();
    UInt_t calculateTrackFlag();
    AliGFWFlags *fRetFlags;
    UInt_t flTrFlag;
    UInt_t fEventFlag;
    vector<UInt_t> fTrackMasks; //!
    Int_t fNTrMasks; //Number of flags added to the filter
    vector<UInt_t> fEventMasks; //!
    Int_t fNEvMasks; //Number of flags added to the filter
    AliAODTrack *fAODTrack; //!
    AliEventCuts *fEventCuts; //!
    Double_t fPtMin;
    Double_t fPtMax;
    Double_t fEtaMin;
    Double_t fEtaMax;
};
#endif
