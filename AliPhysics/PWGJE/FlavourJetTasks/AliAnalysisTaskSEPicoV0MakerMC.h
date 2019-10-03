#ifndef ALIANALYSISTASKSEPICOV0MAKERMC_H
#define ALIANALYSISTASKSEPICOV0MAKERMC_H

//*************************************************************************
// Class AliAnalysisTaskSEPicoV0MakerMC
// AliAnalysisTaskSE for V0s (K0 short, Lambda... ) filtering
// lite version only for MC
// Author: X-M. Zhang, xmzhang@lbl.gov
//*************************************************************************

#include "AliAnalysisTaskSE.h"

class TClonesArray;

class AliAODv0;
class AliESDv0;
class AliAODEvent;
class AliESDEvent;
class AliPIDResponse;

class AliPicoV0MC;

class AliAnalysisTaskSEPicoV0MakerMC : public AliAnalysisTaskSE {

 public :

  AliAnalysisTaskSEPicoV0MakerMC();
  AliAnalysisTaskSEPicoV0MakerMC(const char *name);
  virtual ~AliAnalysisTaskSEPicoV0MakerMC();

  virtual void Init();
  virtual void LocalInit() { Init(); }
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *opt);
  virtual void Terminate(Option_t *opt);
  virtual void NotifyRun();
//=============================================================================

  void SetTriggerMask(UInt_t w)       { fTriggerMask   = w; }
  void SetCollitionType(UInt_t w)     { fCollisionType = w; }

  void SetMultRange(const Double_t dMin,
                    const Double_t dMax,
                    const TString sEsti="V0M",
                    const Bool_t  bOld=kFALSE) {
    fMultMin = dMin;
    fMultMax = dMax;
    fMultEst = sEsti;
    fMultOld = bOld;
    return;
  }

  void SetRefitV0ESD()   { fIsRefitV0sESD  = kTRUE; }
  void SetSkipFastOnly() { fIsSkipFastOnly = kTRUE; }
  void SetDMPjetMC()     { fIsDPMjetMC     = kTRUE; }
//=============================================================================

  void SetV0PtRange(Double_t  dMin, Double_t dMax) { fCutMinV0Pt  = dMin; fCutMaxV0Pt  = dMax; }
  void SetV0RapRange(Double_t dMin, Double_t dMax) { fCutMinV0Rap = dMin; fCutMaxV0Rap = dMax; }

  void SetDauPtRange(Double_t d) { fCutMinDauPt = d; }
  void SetDauEtaRange(Double_t dMin, Double_t dMax) { fCutMinDauEta = dMin; fCutMaxDauEta = dMax; }

  void SetV0Cuts(Double_t d[14]);
  void SetKaCutMaxKaArmFrac(Double_t d) { fCutMaxKshortArmFrac = d; }
  void SetLaCutMaxLaArmFrac(Double_t d) { fCutMaxLambdaArmFrac = d; }
//=============================================================================

 private :

  AliAnalysisTaskSEPicoV0MakerMC(const AliAnalysisTaskSEPicoV0MakerMC &);
  AliAnalysisTaskSEPicoV0MakerMC& operator=(const AliAnalysisTaskSEPicoV0MakerMC &);
//=============================================================================

  void FillPicoV0s();

  AliPicoV0MC *SelectV0Candidate(AliAODv0 const *pV0);
  AliPicoV0MC *SelectV0Candidate(AliESDv0 const *pV0);

  Bool_t IsEventNotAcpt();
  Bool_t IsEventNotINEL();
  Bool_t IsEventNotMBsa();

  void FillHistograms();
  void CreateHistograms();

  void InitAnalysis();
//=============================================================================

  UInt_t fTriggerMask;   //
  UInt_t fCollisionType; //

  Bool_t fUseAnaUtils; //
  Bool_t fIsDPMjetMC; //

  TString  fMultEst;  //
  Double_t fMultMin;  //
  Double_t fMultMax;  //
  Bool_t   fMultOld;  //

  Bool_t fIsSkipFastOnly; //
  Bool_t fIsRefitV0sESD;  //
//=============================================================================

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

  AliAODEvent      *fEventAOD; //!
  AliESDEvent      *fEventESD; //!
  AliPIDResponse   *fRespoPID; //!

  Double_t fPrimaryVtx[3]; //!
  UInt_t fEventAcptMask;   //
//=============================================================================

  TClonesArray *fPicoV0sClArr;  //!
  TList *fListUserOutputs;  //!
//=============================================================================

  ClassDef(AliAnalysisTaskSEPicoV0MakerMC, 2)
};

#endif
