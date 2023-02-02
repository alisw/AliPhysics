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
using namespace std;
namespace GFWFlags {
  enum kLocalEventFlags {
    klEventCuts =      BIT(0), //Event cuts (AliEventCuts)
    klVtxOK =          BIT(1), //Vertex rezolution
    klVtxZ10 =         BIT(2), //Vtx. z < 10 cm
    klVtxZ5 =          BIT(3), //Vtx. z < 5 cm
    klVtxZ7 =          BIT(4), //Vtx. z < 7 cm
    klVtxZ9 =          BIT(5)  //Vtx. z < 9 cm
  };
  enum kLocalTrackFlags {
    klFB32 =           BIT(0), //FB32
    klFB64 =           BIT(1), //FB64
    klFB256 =          BIT(2), //FB256
    klFB512 =          BIT(3), //FB512
    klFB96 =           BIT(4), //FB96 ( 32 || 64)
    klFB96Tuned =      BIT(5), //FB96 with fraction of shared TPC Clusters < 0.4 (for compatibility with FB768Tuned)
    klFB768 =          BIT(6), //FB768 (256 || 512)
    klFB768Tuned =     BIT(7), //FB768 where 512 FB requires a hit on first SDD (for compatibility with FB96Tuned)
    klSharedClusters = BIT(8), //fraction of shared TPC clusters < 0.4
    klHitOnSDD =       BIT(9), //hit on first layer of SDD
    klNoSPD =          BIT(10), //hit on first layer of SDD
    klDCAz20 =         BIT(11),//DCA z < 2 cm
    klDCAz10 =         BIT(12),//DCA z < 1 cm
    klDCAz05 =         BIT(13),//DCA z < 0.5 cm
    klDCAxy2010 =      BIT(14),//DCAxy cut tuned to 2010
    klDCAxy2011 =      BIT(15),//DCAxy cut tuned to 2011
    klDCAxy8Sigma =    BIT(16),//DCAxy cut on 8 sigma, 2011
    klDCAxy4Sigma =    BIT(17),//DCAxy cut on 4 sigma, 2011
    klDCAxy10Sigma =   BIT(18),//DCAxy cut on 10 sigma, 2011
    klTPCchi2PC25 =    BIT(19),//TPC chi2/cluster < 2.5
    klTPCchi2PC20 =    BIT(20),//TPC chi2/cluster < 2.0
    klTPCchi2PC30 =    BIT(21),//TPC chi2/cluster < 3.0
    klTPCchi2PC40 =    BIT(22),//TPC chi2/cluster < 4.0
    klNTPCcls70 =      BIT(23),//Number of TPC clusters 70
    klNTPCcls80 =      BIT(24),//Number of TPC clusters 80
    klNTPCcls90 =      BIT(25),//Number of TPC clusters 90
    klNTPCcls100 =     BIT(26)//Number of TPC clusters 100
  };
  enum EventFlags { kNominal=0, kVtx9, kVtx7, kVtx5, kAllEvFlags}; //Better keep these as uint_t so that we reuse them as array indeces
  enum TrackFlags { kFB96=0, kFB768,
                    kDCAz10, kDCAz05,
                    kDCA4Sigma, kDCA10Sigma,
                    kChiSq2, kChiSq3,
                    kNTPC80, kNTPC90, kNTPC100,
                    kFB96Tuned, kFB768Tuned, //These are for developing purposes and shouldn't be used!
                    kFB768DCAz,
                    kFB768DCAxyLow,
                    kFB768DCAxyHigh,
                    kFB768ChiSq2, kFB768ChiSq3,
                    kFB768nTPC, kFB96MergedDCA,
                    kChiSq25, //For testing purposes on 15_pass1 data
                    kAllTrFlags
                  };
  static const Int_t BitIndex(const UInt_t &lFlag) {
    if(lFlag==0) return -1;
    for(Int_t i=0;i<sizeof(lFlag)*8;i++) if(lFlag&(1<<i)) return i;
    return -1;
  };
  static const TString GetSystPF(UInt_t lEv, UInt_t lTr) { return TString(Form("_Ev%i_Tr%i",lEv,lTr)); };
  const Int_t gNTrackFlags=21;
  const Int_t gNEventFlags=4;
};
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
