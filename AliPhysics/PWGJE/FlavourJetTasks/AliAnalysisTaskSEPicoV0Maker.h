#ifndef ALIANALYSISTASKSEPICOV0MAKER_H
#define ALIANALYSISTASKSEPICOV0MAKER_H

//*************************************************************************
// Class AliAnalysisTaskSEPicoV0Maker
// AliAnalysisTaskSE for V0s (K0 short, Lambda... ) filtering
// Author: X-M. Zhang, xmzhang@lbl.gov
//*************************************************************************

#include <TMap.h>
#include <TObjString.h>
#include <TParameter.h>

#include "AliAnalysisTaskSE.h"

class TClonesArray;

class AliAODv0;
class AliESDv0;
class AliAODEvent;
class AliESDEvent;
class AliPIDResponse;
class AliPicoV0RD;
class AliPicoV0MC;

class AliAnalysisTaskSEPicoV0Maker : public AliAnalysisTaskSE {

 public :

  AliAnalysisTaskSEPicoV0Maker();
  AliAnalysisTaskSEPicoV0Maker(const char *name, Bool_t bIsMC);
  virtual ~AliAnalysisTaskSEPicoV0Maker();

  virtual void Init();
  virtual void LocalInit() { Init(); }
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *opt);
  virtual void Terminate(Option_t *opt);
  virtual void NotifyRun();
//=============================================================================

  void SetTriggerMask(UInt_t w)       { fTriggerMask   = w; }
  void SetCollitionType(UInt_t w)     { fCollisionType = w; }

  void OldMultEstimator(const Bool_t b=kTRUE) { fUseMultOld = b; }
  void AddMultEstimator(const TString s) {
    fMultEsti.Add(new TObjString(s.Data()),
                  new TParameter<Float_t>(s.Data(), -999.));
    return;
  }

  void SetMultRange(const Double_t dMin,
                    const Double_t dMax,
                    const TString sEst="V0M") {
    fCutMinMult = dMin;
    fCutMaxMult = dMax;
    fMultEstDef = sEst;
    return;
  }

  void SetUseAnaUtils(const Bool_t b) { fUseAnaUtils = b; }

  void SetSkipFastOnly() { fIsSkipFastOnly = kTRUE; }
  void SetRefitV0ESD()   { fIsRefitV0sESD  = kTRUE; }
  void SetDMPjetMC()     { fIsDPMjetMC     = kTRUE; }
//=============================================================================

  void SetV0PtRange(Double_t  dMin, Double_t dMax) { fCutMinV0Pt  = dMin; fCutMaxV0Pt  = dMax; }
  void SetV0RapRange(Double_t dMin, Double_t dMax) { fCutMinV0Rap = dMin; fCutMaxV0Rap = dMax; }

  void SetDauPtRange(Double_t d) { fCutMinDauPt = d; }
  void SetDauEtaRange(Double_t dMin, Double_t dMax) { fCutMinDauEta = dMin; fCutMaxDauEta = dMax; }
//=============================================================================

 private :

  AliAnalysisTaskSEPicoV0Maker(const AliAnalysisTaskSEPicoV0Maker &);
  AliAnalysisTaskSEPicoV0Maker& operator=(const AliAnalysisTaskSEPicoV0Maker &);
//=============================================================================

  void FillPicoV0s();

  AliPicoV0RD *SelectV0CandidateRD(AliAODv0 const *pV0);
  AliPicoV0RD *SelectV0CandidateRD(AliESDv0 const *pV0);

  AliPicoV0MC *SelectV0CandidateMC(AliAODv0 const *pV0);
  AliPicoV0MC *SelectV0CandidateMC(AliESDv0 const *pV0);

  Bool_t IsEventNotAcpt();
  Bool_t IsEventNotINEL();
  Bool_t IsEventNotMBsa();

  void FillHistogramsEH();
  void FillHistogramsMC();

  void CreateHistogramsEH();
  void CreateHistogramsMC();

  void InitAnalysis();
  void InitParamsPP();
  void InitParamsPA();
  void InitParamsAP();
  void InitParamsAA();
//=============================================================================

  UInt_t fTriggerMask;   //
  UInt_t fCollisionType; //

  Bool_t fIsAnaUseMC; //
  Bool_t fIsDPMjetMC; //

  Bool_t fUseMultOld;
  Bool_t fUseAnaUtils; //

  Bool_t fIsSkipFastOnly; //
  Bool_t fIsRefitV0sESD;  //
//=============================================================================

  Double_t fRapidityShift;  //

  TString  fMultEstDef;
  Double_t fCutMinMult;
  Double_t fCutMaxMult;

  Double_t fCutMinV0Pt;   //
  Double_t fCutMaxV0Pt;   //
  Double_t fCutMinV0Rap;  //
  Double_t fCutMaxV0Rap;  //

  Double_t fCutMinDauPt;   //
  Double_t fCutMinDauEta;  //
  Double_t fCutMaxDauEta;  //

  Double_t fCutMaxV0Chi2;    //
  Double_t fCutMinV0Radius;  //
  Double_t fCutMaxV0Radius;  //

  Double_t fCutMaxDausDCA;      //
  Double_t fCutMinDauDCAtoPV;   //
  Float_t  fCutMinDauXrowsTPC;  //
  Double_t fCutMinDauXrowsOverFindableClusTPC;  //

  Float_t  fCutMaxKshortSigmaTPC;  //
  Double_t fCutMinKshortCosPA;     //
  Double_t fCutMaxKshortCtau;      //
  Double_t fCutMaxKshortArmFrac;   //
  Double_t fCutMinKshortDeltaM;    //

  Float_t  fCutMaxLambdaSigmaTPC;  //
  Double_t fCutMinLambdaCosPA;     //
  Double_t fCutMaxLambdaCtau;      //
  Double_t fCutMaxLambdaArmFrac;   //
  Double_t fCutMinLambdaDeletaM;   //
//=============================================================================

  AliAODEvent    *fEventAOD; //!
  AliESDEvent    *fEventESD; //!
  AliPIDResponse *fRespoPID; //!

  Double_t fPrimaryVtx[3]; //!
  UInt_t fEventAcptMask;   //

  TMap fMultEsti; //
//=============================================================================

  TClonesArray *fPicoV0sClArr;  //!

  TList *fOutputListEH;  //!
  TList *fOutputListMC;  //!
//=============================================================================

  ClassDef(AliAnalysisTaskSEPicoV0Maker, 6)
};

#endif
