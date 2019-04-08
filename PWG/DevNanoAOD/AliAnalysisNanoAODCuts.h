#ifndef _ALIANALYSISNANOAODCUTSANDSETTERS_H_
#define _ALIANALYSISNANOAODCUTSANDSETTERS_H_

#include "AliAnalysisCuts.h"
#include "AliNanoAODCustomSetter.h"
#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"
#include <map>

class AliEventCuts;

class AliAnalysisNanoAODTrackCuts : public AliAnalysisCuts
{
public:
  AliAnalysisNanoAODTrackCuts();
  virtual ~AliAnalysisNanoAODTrackCuts()  {}
  virtual Bool_t IsSelected(TObject* obj); // TObject should be an AliAODTrack
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }
  UInt_t GetBitMask() { return fBitMask; }
  void  SetBitMask (UInt_t var) { fBitMask = var;}
  Float_t GetMinPt() { return fMinPt; }
  void  SetMinPt (Float_t var) { fMinPt = var;}
  Float_t GetMaxEta() { return fMaxEta; }
  void  SetMaxEta (Float_t var) { fMaxEta = var;}

private:
  UInt_t fBitMask; // Only AOD tracks matching this bit mask are accepted
  Float_t fMinPt; // miminum pt of the tracks
  Float_t fMaxEta; // MaxEta

  ClassDef(AliAnalysisNanoAODTrackCuts,1); // track cut object for nano AOD filtering
};

class AliAnalysisNanoAODV0Cuts : public AliAnalysisCuts
{
public:
  AliAnalysisNanoAODV0Cuts();
  virtual ~AliAnalysisNanoAODV0Cuts() {}
  virtual Bool_t IsSelected(TObject* obj);  // TObject should be an AliAODv0
  virtual Bool_t IsSelected(TList* /* list */) {
    return kTRUE;
  }
  void SetOnFlyStatus(bool flyStat) {
    fSelectOnFly = true;
    fOnFlyStatus = flyStat;
  }
  void Setv0pTMin(float pTMin) {
    fSelectv0pT = true;
    fv0pTMin = pTMin;
  }
  void Setv0EtaMax(float EtaMax) {
    fSelectv0Eta = true;
    fv0EtaMax = EtaMax;
  }
  void SetTransverseRadiusMin(float MinRad) {
    fSelectTransverseRadius = true;
    fTransverseRadiusMin = MinRad;
  }
  void SetTransverseRadiusMax(float MaxRad) {
    fSelectTransverseRadius = true;
    fTransverseRadiusMax = MaxRad;
  }
  void SetCPAMin(float CPAMin) {
    fSelectCPA = true;
    fCPAMin = CPAMin;
  }
  void SetMaxDCADaughtersToV0Vtx(float MaxDCA) {
    fSelectDCADaugv0Vtx = true;
    fDCADaugv0VtxMax = MaxDCA;
  }
  void SetMinDCADaughtersToPrimVtx(float MinDCA) {
    fSelectDCADaugPrimVtx = true;
    fDCADaugPrimVtxMin = MinDCA;
  }
  void SetDaughterEtaMax(float EtaMax) {
    fSelectDaugEta = true;
    fDaughEtaMax = EtaMax;
  }
  void SetMinClsTPCDaughters(int minCls) {
    fSelectClsTPC = true;
    fDaugMinClsTPC = minCls;
  }
  void SetDaughterMaxnSigTPC(float maxnSig) {
    fSelectDaugnPID = true;
    fDaugnSigTPCMax = maxnSig;
  }
  void SetRequireTimingDaughters(bool require) {
    fSelectTimingDaug = true;
    fTimingDaughter = require;
  }
private:
  bool fSelectOnFly;
  bool fOnFlyStatus;
  bool fSelectv0pT;
  float fv0pTMin;
  bool fSelectv0Eta;
  float fv0EtaMax;
  bool fSelectTransverseRadius;
  float fTransverseRadiusMin;
  float fTransverseRadiusMax;
  bool fSelectCPA;
  float fCPAMin;
  bool fSelectDCADaugv0Vtx;
  float fDCADaugv0VtxMax;
  bool fSelectDCADaugPrimVtx;
  float fDCADaugPrimVtxMin;
  bool fSelectDaugEta;
  float fDaughEtaMax;
  bool fSelectClsTPC;
  int fDaugMinClsTPC;
  bool fSelectDaugnPID;
  float fDaugnSigTPCMax;
  bool fSelectTimingDaug;
  bool fTimingDaughter;

  ClassDef(AliAnalysisNanoAODV0Cuts, 1); // track cut object for nano AOD filtering
};

class AliAnalysisNanoAODEventCuts : public AliAnalysisCuts
{
public:
  AliAnalysisNanoAODEventCuts();
  virtual ~AliAnalysisNanoAODEventCuts() {}
  virtual Bool_t IsSelected(TObject* obj); // TObject should be an AliAODEvent
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }
  
  void SetMultiplicityRange(AliAnalysisCuts* cutObject, Int_t minMultiplicity, Int_t maxMultiplicity) { fTrackCut = cutObject; fMinMultiplicity = minMultiplicity; fMaxMultiplicity = maxMultiplicity; }
  AliEventCuts& GetAliEventCuts() { return fEventCuts; }
  
private:
  AliAnalysisCuts* fTrackCut; // track cut object for multiplicity cut
  Int_t fMinMultiplicity;   // minimum number of tracks to accept this event
  Int_t fMaxMultiplicity;   // maximal number of tracks to accept this event
  
  AliEventCuts fEventCuts; // AliEventCut object for Run 2
  
  ClassDef(AliAnalysisNanoAODEventCuts, 3); // event cut object for nano AOD filtering
};

class AliNanoAODSimpleSetter : public AliNanoAODCustomSetter
{
public:
  AliNanoAODSimpleSetter() : fInitialized(kFALSE), fMultMap() {;}
  virtual ~AliNanoAODSimpleSetter(){;}

  virtual void SetNanoAODHeader(const AliAODEvent * event   , AliNanoAODHeader * head ,TString varListHeader  );
  virtual void SetNanoAODTrack (const AliAODTrack * /*aodTrack*/, AliNanoAODTrack * /*spTrack*/);
  
protected:
  void Init(AliNanoAODHeader* head, TString varListHeader);
  
  Bool_t fInitialized;
  std::map<TString,int> fMultMap;

  ClassDef(AliNanoAODSimpleSetter, 2)

};

#endif /* _ALIANALYSISNANOAODCUTSANDSETTERS_H_ */
