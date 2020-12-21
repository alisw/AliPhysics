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
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }

  void SaveOnlyOnTheFly()                            { fSelectOnFly = true; fOnFlyStatus = kFALSE; }
  void SaveOnlyOffline()                             { fSelectOnFly = true; fOnFlyStatus = kTRUE; }
  void Setv0pTMin(Float_t pTMin)                     { fv0pTMin = pTMin;  }
  void Setv0EtaMax(Float_t EtaMax)                   { fv0EtaMax = EtaMax; }
  void SetTransverseRadius(Float_t min, Float_t max) { fTransverseRadiusMin = min; fTransverseRadiusMax = max; }
  void SetCPAMin(Float_t CPAMin)                     { fCPAMin = CPAMin; }
  void SetMaxDCADaughtersToV0Vtx(Float_t MaxDCA)     { fDCADaugv0VtxMax = MaxDCA; }
  void SetMinDCADaughtersToPrimVtx(Float_t MinDCA)   { fDCADaugPrimVtxMin = MinDCA; }
  void SetDaughterEtaMax(Float_t EtaMax)             { fDaughEtaMax = EtaMax; }
  void SetMinClsTPCDaughters(int minCls)             { fDaugMinClsTPC = minCls; }
  void SetLambdaDaugnSigTPCMax(Float_t maxnSig)      { fLambdaDaugnSigTPCMax = maxnSig; }
  void SetRequireTimingDaughters(Bool_t require)     { fCheckDaughterPileup = require; }
  void SetRequireTPCRefitDaughters(Bool_t require)     { fCheckDaughterTPCRefit = require; }
  
private:
  Bool_t fSelectOnFly;
  Bool_t fOnFlyStatus;
  Float_t fv0pTMin;
  Float_t fv0EtaMax;
  Float_t fTransverseRadiusMin;
  Float_t fTransverseRadiusMax;
  Float_t fCPAMin;
  Float_t fDCADaugv0VtxMax;
  Float_t fDCADaugPrimVtxMin;
  Float_t fDaughEtaMax;
  Int_t fDaugMinClsTPC;
  Float_t fLambdaDaugnSigTPCMax;
  Bool_t fCheckDaughterPileup;
  Bool_t fCheckDaughterTPCRefit;

  ClassDef(AliAnalysisNanoAODV0Cuts, 1); // track cut object for nano AOD filtering
};

class AliAnalysisNanoAODV0ParametricCuts : public AliAnalysisCuts
{
public:
    AliAnalysisNanoAODV0ParametricCuts();
    virtual ~AliAnalysisNanoAODV0ParametricCuts() {}
    void SetupDefaultPbPb2015Cuts(); //setup loose cuts that are ok for 2015
    virtual Bool_t IsSelected(TObject* obj);  // TObject should be an AliAODv0
    virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }
    
    void SetOnFlyStatus(Bool_t flyStat)                { fSelectOnFly = true; fOnFlyStatus = flyStat; }
    void Setv0pTMin(Float_t pTMin)                     { fv0pTMin = pTMin;  }
    void SetTransverseRadius(Float_t min, Float_t max) { fTransverseRadiusMin = min; fTransverseRadiusMax = max; }
    void SetCPAMin(Float_t CPAMin)                     { fCPAMin = CPAMin; }
    void SetMaxDCAV0Dau(Float_t MaxDCAV0Dau)           { fMaxDCAV0Dau = MaxDCAV0Dau; }
    void SetMinDCADaughtersToPV(Float_t MinDCADauToPV) { fDCADauToPV = MinDCADauToPV; }
    void SetDaughterEtaMax(Float_t EtaMax)             { fDaughEtaMax = EtaMax; }
    void SetMinClsTPCDaughters(int minCls)             { fDaugMinClsTPC = minCls; }
    void SetLambdaDaugnSigTPCMax(Float_t maxnSig)      { fLambdaDaugnSigTPCMax = maxnSig; }
    void SetRequireTimingDaughters(Bool_t require)     { fCheckDaughterPileup = require; }
    void SetParametricCosPA(Float_t *lVars)     {
        fUseParametricCosPA = kTRUE;
        for(Int_t ii=0; ii<5; ii++) fParCosPA[ii] = lVars[ii];
    }
    
private:
    Bool_t fSelectOnFly;
    Bool_t fOnFlyStatus;
    Float_t fv0pTMin;
    Float_t fTransverseRadiusMin;
    Float_t fTransverseRadiusMax;
    Float_t fCPAMin;
    Float_t fMaxDCAV0Dau;
    Float_t fDCADauToPV;
    Float_t fDaughEtaMax;
    Int_t fDaugMinClsTPC;
    Float_t fLambdaDaugnSigTPCMax;
    Bool_t fCheckDaughterPileup;
    
    //Parametric CosPA
    Bool_t fUseParametricCosPA;
    Float_t fParCosPA[5];
    
    ClassDef(AliAnalysisNanoAODV0ParametricCuts, 1); // track cut object for nano AOD filtering
};

class AliAnalysisNanoAODCascadeCuts : public AliAnalysisCuts
{
 public:
  AliAnalysisNanoAODCascadeCuts();
  virtual ~AliAnalysisNanoAODCascadeCuts() {};
  virtual Bool_t IsSelected(TObject* obj);  // TObject should be an AliAODCascade
  virtual Bool_t IsSelected(TList*   /* list */ )   { return kTRUE; }
  void SetCascpTMin(Float_t pTMin)                  { fCascpTMin = pTMin;  }
  void SetMinDCADaughtersToPrimVtx(Float_t MinDCA)  { fDCADaugPrimVtxMin = MinDCA; }
  void SetCPACascMin(Float_t CPAMin)                { fCPACascMin = CPAMin; }
  void SetCascTransverseRadiusMax(Float_t Max)      { fTransverseRadiusCasc = Max; }
  void SetCPAv0Min(Float_t CPAMin)                  { fCPAv0Min = CPAMin; }
  void Setv0TransverseRadiusMax(Float_t Max)        { fTransverseRadiusv0 = Max; }
  void SetMinDCAv0ToPrimVtx(Float_t MinDCA)         { fDCAv0PrimVtxMin = MinDCA; }
  void SetDaughterEtaMax(Float_t EtaMax)            { fDaughEtaMax = EtaMax; }
  void SetLambdaDaugnSigTPCMax(Float_t maxnSig)     { fCascDaugnSigTPCMax = maxnSig; }
  void SetRequireTimingDaughters(Bool_t require)    { fCheckDaughterPileup = require; }
  void SetRequireTPCRefitDaughters(Bool_t require)     { fCheckDaughterTPCRefit = require; }
 private:
  Float_t fCascpTMin;
  Float_t fDCADaugPrimVtxMin;
  Float_t fCPACascMin;
  Float_t fTransverseRadiusCasc;
  Float_t fCPAv0Min;
  Float_t fTransverseRadiusv0;
  Float_t fDCAv0PrimVtxMin;
  Float_t fDaughEtaMax;
  Float_t fCascDaugnSigTPCMax;
  bool fCheckDaughterPileup;
  bool fCheckDaughterTPCRefit;

  ClassDef(AliAnalysisNanoAODCascadeCuts, 1); // track cut object for nano AOD filtering
};

class AliAnalysisNanoAODCascadeParametricCuts : public AliAnalysisCuts
{
public:
    AliAnalysisNanoAODCascadeParametricCuts();
    virtual ~AliAnalysisNanoAODCascadeParametricCuts() {};
    void SetupDefaultPbPb2015Cuts(); 
    virtual Bool_t IsSelected(TObject* obj);  // TObject should be an AliAODCascade
    virtual Bool_t IsSelected(TList*   /* list */  )   { return kTRUE; }
    void SetCascpTMin(Float_t pTMin)                  { fCascpTMin = pTMin;  }
    void SetMinDCAV0DauToPV(Float_t MinDCA)           { fDCAV0DauToPV = MinDCA; }
    void SetMinDCABachToPV(Float_t MinDCA)            { fDCABachToPV = MinDCA; }
    void SetCPACascMin(Float_t CPAMin)                { fCPACascMin = CPAMin; }
    void SetCPAv0Min(Float_t CPAMin)                  { fCPAv0Min = CPAMin; }
    void SetCascTransverseRadiusMax(Float_t Max)      { fTransverseRadiusCasc = Max; }
    void Setv0TransverseRadiusMax(Float_t Max)        { fTransverseRadiusv0 = Max; }
    void SetMinDCAv0ToPrimVtx(Float_t MinDCA)         { fDCAv0PrimVtxMin = MinDCA; }
    void SetDaughterEtaMax(Float_t EtaMax)            { fDaughEtaMax = EtaMax; }
    void SetLambdaDaugnSigTPCMax(Float_t maxnSig)     { fCascDaugnSigTPCMax = maxnSig; }
    void SetRequireTimingDaughters(Bool_t require)    { fCheckDaughterPileup = require; }
    void SetParametricV0CosPA(Float_t *lVars)     {
        fUseParametricV0CosPA = kTRUE;
        for(Int_t ii=0; ii<5; ii++) fParV0CosPA[ii] = lVars[ii];
    }
    void SetParametricCascCosPA(Float_t *lVars)     {
        fUseParametricCascCosPA = kTRUE;
        for(Int_t ii=0; ii<5; ii++) fParCascCosPA[ii] = lVars[ii];
    }
private:
    Float_t fCascpTMin;
    Float_t fDCAV0DauToPV;
    Float_t fDCABachToPV;
    Float_t fCPACascMin;
    Float_t fCPAv0Min;
    Float_t fTransverseRadiusCasc;
    Float_t fTransverseRadiusv0;
    Float_t fDCAv0PrimVtxMin;
    Float_t fDaughEtaMax;
    Float_t fCascDaugnSigTPCMax;
    bool fCheckDaughterPileup;
    
    //Parametric CosPA
    Bool_t fUseParametricV0CosPA;
    Float_t fParV0CosPA[5];
    Bool_t fUseParametricCascCosPA;
    Float_t fParCascCosPA[5];
    
    ClassDef(AliAnalysisNanoAODCascadeParametricCuts, 1); // track cut object for nano AOD filtering
};

class AliAnalysisNanoAODPhotonCuts : public AliAnalysisCuts
{
public:
  AliAnalysisNanoAODPhotonCuts();
  virtual ~AliAnalysisNanoAODPhotonCuts() {}
  virtual Bool_t IsSelected(TObject* obj);  // TObject should be an AliAODv0
  virtual Bool_t IsSelected(TList*   /* list */ ) { return kTRUE; }

  void SetPhotonEtaMax(float EtaMax)               { fPhotonEtaMax = EtaMax; }
  void SetConversionRadius(float min, Float_t max) { fPhotonConvRadiusMin = min; fPhotonConvRadiusMax = max; }
  void SetPsiPairMax(float psiPairMax)             { fPhotonPsiPairMax = psiPairMax; }

private:
  float fPhotonEtaMax;
  float fPhotonConvRadiusMin;
  float fPhotonConvRadiusMax;
  float fPhotonPsiPairMax;

  ClassDef(AliAnalysisNanoAODPhotonCuts, 1); // track cut object for nano AOD filtering
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

class AliAnalysisNanoAODMCParticleCuts : public AliAnalysisCuts
{
public:
  AliAnalysisNanoAODMCParticleCuts();
  virtual ~AliAnalysisNanoAODMCParticleCuts() {}
  virtual bool IsSelected(TObject* obj);  // TObject should be an AliAODMCParticle
  virtual bool IsSelected(TList* /* list */) {
    return kTRUE;
  }
  bool GetSelectPrimaries() const { return fDoSelectPrimaries; }
  void SetSelectPrimaries(bool doIt) { fDoSelectPrimaries = doIt; }
  bool GetSelectChargedParticles() const { return fDoSelectCharged; }
  void SetSelectChargedParticles(bool doIt) { fDoSelectCharged = doIt; }
  float GetMinPt() const { return fMinPt; }
  void SetMinPt(float var) { fMinPt = var; }
  float GetMaxEta() const { return fMaxEta; }
  void SetMaxEta(float var) { fMaxEta = var; }
  std::vector<int> &GetKeepParticles() { return fPDGToKeep; }
  void SetKeepParticles(std::vector<int> pdgCodes) { fPDGToKeep = pdgCodes; }
  void AddKeepParticles(int pdgCode) { fPDGV0.push_back(pdgCode); }
  std::vector<int> &GetKeepV0s() { return fPDGV0; }
  void SetKeepV0s(std::vector<int> pdgCodes) { fPDGV0 = pdgCodes; }
  void AddKeepV0(int pdgCode) { fPDGV0.push_back(pdgCode); }
  std::vector<int> &GetKeepV0sCascades() { return fPDGV0Cascade; }
  void SetKeepCascadeV0s(std::vector<int> pdgCodes) { fPDGV0Cascade = pdgCodes; }
  void AddKeepCascadeV0(int pdgCode) { fPDGV0Cascade.push_back(pdgCode); }

private:
  bool fDoSelectPrimaries; // switch to only select IsPhysicalPrimary() particles
  bool fDoSelectCharged; // switch to only select charged particles
  float fMinPt;  // miminum pt of the tracks
  float fMaxEta; // MaxEta
  std::vector<int> fPDGToKeep;  // vector of PDG codes that we want to keep anyways
  std::vector<int> fPDGV0;  // vector of PDG codes that we want to match the V0s to
  std::vector<int> fPDGV0Cascade;  // vector of PDG codes that we want to match the V0s to

  ClassDef(AliAnalysisNanoAODMCParticleCuts,2); // MC particle cut object for nano AOD filtering
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
