/*
 * AliOtonOmegaAnalysis.h
 *
 *  Created on: 24 Nov 2017
 *      Author: bernhardhohlweger. Adapted for Proton-Omega by Oton Vazquez Doce.
 */

#ifndef ALIOTONOMEGAANALYSIS_H_
#define ALIOTONOMEGAANALYSIS_H_

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODTrack.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliOtonOmegaCascade.h"
#include "AliOtonOmegaCascadeCuts.h"
#include "AliFemtoDreamControlSample.h"
#include "AliPIDResponse.h"
#include "Rtypes.h"
class AliOtonOmegaAnalysis {
 public:
  AliOtonOmegaAnalysis();
  AliOtonOmegaAnalysis(const AliOtonOmegaAnalysis& analysis);
  AliOtonOmegaAnalysis& operator=(const AliOtonOmegaAnalysis& analysis);
  virtual ~AliOtonOmegaAnalysis();
  void SetMVPileUp(bool mvPileUp) {
    fMVPileUp = mvPileUp;
  }
  ;
  void SetEvtCutQA(bool setQA) {
    fEvtCutQA = setQA;
  }
  ;
  void SetEventCuts(AliFemtoDreamEventCuts *cuts) {
    fEvtCuts = cuts;
  }
  ;
  TList *GetEventCutHists() {
    return fEvtCuts->GetHistList();
  }
  ;
  void SetTrackCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCuts = cuts;
  }
  ;
  TList *GetTrackCutHists() {
    return fTrackCuts->GetQAHists();
  }
  ;
  TList *GetTrackCutHistsMC() {
    return fTrackCuts->GetMCQAHists();
  }
  ;
  void SetAntiTrackCuts(AliFemtoDreamTrackCuts *cuts) {
    fAntiTrackCuts = cuts;
  }
  ;
  TList *GetAntitrackCutHists() {
    return fAntiTrackCuts->GetQAHists();
  }
  ;
  TList *GetAntitrackCutHistsMC() {
    return fAntiTrackCuts->GetMCQAHists();
  }
  ;
  void Setv0Cuts(AliFemtoDreamv0Cuts *cuts) {
    fv0Cuts = cuts;
  }
  ;
  TList *Getv0CutHist() {
    return fv0Cuts->GetQAHists();
  }
  ;
  TList *Getv0MCHist() {
    return fv0Cuts->GetMCQAHists();
  }
  ;
  void SetAntiv0Cuts(AliFemtoDreamv0Cuts *cuts) {
    fAntiv0Cuts = cuts;
  }
  ;
  TList *GetAntiv0CutHist() {
    return fAntiv0Cuts->GetQAHists();
  }
  ;
  TList *GetAntiv0MCHist() {
    return fAntiv0Cuts->GetMCQAHists();
  }
  ;
  void SetCascadeCuts(AliOtonOmegaCascadeCuts *cuts) {
    fCascCuts = cuts;
  }
  ;
  TList *GetCascadeCutHist() {
    return fCascCuts->GetQAHists();
  }
  ;
  TList *GetCascadeMCHist() {
    return fCascCuts->GetMCQAHists();
  }
  ;
  void SetAntiCascadeCuts(AliOtonOmegaCascadeCuts *cuts) {
    fAntiCascCuts = cuts;
  }
  ;
  TList *GetAntiCascadeMCHist() {
    return fAntiCascCuts->GetMCQAHists();
  }
  ;
  TList *GetAntiCascadeCutHist() {
    return fAntiCascCuts->GetQAHists();
  }
  ;



  void SetCascadeOmegaCuts(AliOtonOmegaCascadeCuts *cuts) {
    fCascOmegaCuts = cuts;
  }
  ;
  TList *GetCascadeOmegaCutHist() {
    return fCascOmegaCuts->GetQAHists();
  }
  ;
  TList *GetCascadeOmegaMCHist() {
    return fCascOmegaCuts->GetMCQAHists();
  }
  ;
  void SetAntiCascadeOmegaCuts(AliOtonOmegaCascadeCuts *cuts) {
    fAntiCascOmegaCuts = cuts;
  }
  ;
  TList *GetAntiCascadeOmegaMCHist() {
    return fAntiCascOmegaCuts->GetMCQAHists();
  }
  ;
  TList *GetAntiCascadeOmegaCutHist() {
    return fAntiCascOmegaCuts->GetQAHists();
  }
  ;




  TList *GetResultList() {
    return fPartColl ? fPartColl->GetHistList() : nullptr;
  }
  ;
  TList *GetResultQAList() {
    return fPartColl ? fPartColl->GetQAList() : nullptr;
  }
  ;
  TList *GetResultSampleList() {
    return fControlSample ? fControlSample->GetHistList() : nullptr;
  }
  ;
  TList *GetResultSampleQAList() {
    return fControlSample ? fControlSample->GetQAList() : nullptr;
  }
  ;
  TList *GetQAList() {
    return fQA;
  }
  ;
  void SetTrackBufferSize(int size) {
    fTrackBufferSize = size;
  }
  ;
  void SetCollectionConfig(AliFemtoDreamCollConfig *conf) {
    fConfig = conf;
  }
  ;
  void Init(bool isMonteCarlo, UInt_t trigger);


  //tree stuff:
  TTree *GetTreeCascade() {
   return foTTree;
  }
  TTree *GetTreeOmega() {
   return fomegaTTree;
  }
  void InitializeTreeBooking();
  void InitializeTreeValues();
  Bool_t FillTreeCascade(AliESDEvent *evt, AliESDcascade *casc);
  Bool_t FillTreeTrack(Int_t jj, Int_t idtrack, Int_t V0index, AliESDEvent *evt, AliESDcascade *casc);
  Bool_t FillProtonTrack( AliESDEvent *evt, Int_t idtrack);
  Bool_t FillTreeCascadeAOD(AliAODEvent *evt, AliAODcascade *casc);
  Bool_t FillTreeTrackAOD(Int_t jj, AliAODcascade *casc);
  Bool_t FillProtonTrackAOD( AliAODTrack *Proton);


  TString ClassName() {
    return "AliOtonOmegaAnalysis";
  }
  ;
  void Make(AliAODEvent *evt, bool OmegaTreeFlag);
  void Make(AliESDEvent *evt, AliMCEvent *mcEvent, bool CascadeTreeFlag, bool OmegaTreeFlag, Int_t Cut = 0);
 private:
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fMVPileUp;                           //!
  bool fEvtCutQA;                           //!
  TList *fQA;                               //!
  AliFemtoDreamTrack *fFemtoTrack;          //!
  AliFemtoDreamv0 *fFemtov0;                //!
  AliOtonOmegaCascade *fFemtoCasc;         //!
  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamEventCuts *fEvtCuts;         //!
  AliFemtoDreamTrackCuts *fTrackCuts;       //!
  AliFemtoDreamTrackCuts *fAntiTrackCuts;   //!
  AliFemtoDreamv0Cuts *fv0Cuts;             //!
  AliFemtoDreamv0Cuts *fAntiv0Cuts;         //!
  AliOtonOmegaCascadeCuts *fCascCuts;      //!
  AliOtonOmegaCascadeCuts *fAntiCascCuts;  //!
  AliOtonOmegaCascadeCuts *fCascOmegaCuts;      //!
  AliOtonOmegaCascadeCuts *fAntiCascOmegaCuts;  //!
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamControlSample *fControlSample;   //!
  int fTrackBufferSize;
  AliAODTrack **fGTI;             //!
  AliFemtoDreamCollConfig *fConfig;         //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  AliPIDResponse *fPIDResponse; //! PID response object
  bool fIsMC;                               //!

  
  //more tree stuff:
  //Bool_t FillOmegaTree_tracks;
  //Bool_t FillOmegaTree_lambda;
  Bool_t FillOmegaTree_omega;
  Bool_t FillOmegaTree_aomega;
  TTree* fomegaTTree;
  TTree* foTTree;
  Int_t fTRunNumber;
  Float_t fTVx;
  Float_t fTVy;
  Float_t fTVz;
  Int_t fTMult;

  TRandom3 *frndm;

  const Int_t MAXPROTONS = 150;
  Int_t fTnProton;
  //protons:
  Float_t fTProtonP[150];
  Float_t fTProtonEta[150];
  Float_t fTProtonPx[150];
  Float_t fTProtonPy[150];
  Float_t fTProtonPz[150];
  Float_t fTProtonVPx[150];
  Float_t fTProtonVPy[150];
  Float_t fTProtonVPz[150];
  Float_t fTProtonPt[150];
  Float_t fTProtonmT[150];
  Float_t fTProtonTPCmom[150];
  Short_t fTProtonCharge[150];
  Float_t fTProtonTPCsp[150];
  Float_t fTProtonTOFsp[150];
  Float_t fTProtonDCA[150];
  Int_t fTProtonNcl[150];
  Float_t fTProtonCrF[150];
  Int_t fTProtonShared[150];
  Float_t fTProtonTPCchi2[150];
  Bool_t fTProtonITStime[150];
  Bool_t fTProtonTOFtime[150];
  Bool_t fTProtonTPConly[150];
  Bool_t fTProtonITScomplementary[150];
  Bool_t fTProtonITSpure[150];
  Bool_t fTProtonGLOBAL[150];
  UInt_t fTProtonFilterBit[150];


  const Int_t MAXCASCADES = 300;
  Int_t fTnCascade;
  //cascades:
  Bool_t fTCascadeFlag0[300];
  Bool_t fTCascadeFlag1[300];
  Bool_t fTCascadeFlag2[300];
  Float_t fTCascadeP[300];
  Float_t fTCascadePx[300];
  Float_t fTCascadePy[300];
  Float_t fTCascadePz[300];
  Float_t fTCascadePt[300];
  Float_t fTCascademT[300];
  Short_t fTCascadeCharge[300];
  Float_t fTCascadeDCA[300];
  Float_t fTCascadeDaughtersDCA[300];
  Float_t fTCascadeXiMass[300];
  Float_t fTCascadeOmegaMass[300];
  Float_t fTCascadeVx[300];
  Float_t fTCascadeVy[300];
  Float_t fTCascadeVz[300];
  Float_t fTCascadeVr[300];
  Float_t fTCascadePA[300];
  //lambda
  Float_t fTLambdaP[300];
  Float_t fTLambdaPx[300];
  Float_t fTLambdaPy[300];
  Float_t fTLambdaPz[300];
  Float_t fTLambdaPt[300];
  Float_t fTLambdaDCA[300];
  Float_t fTLambdaDaughtersDCA[300];
  Float_t fTLambdaMass[300];
  Float_t fTLambdaK0Mass[300];
  Float_t fTLambdaVx[300];
  Float_t fTLambdaVy[300];
  Float_t fTLambdaVz[300];
  Float_t fTLambdaVr[300];
  Float_t fTLambdaPA[300];
//tracks: 0 proton, 1 pion, 2 bachelor
  Float_t fTTrackP[300][3];
  Float_t fTTrackPx[300][3];
  Float_t fTTrackPy[300][3];
  Float_t fTTrackPz[300][3];
  Float_t fTTrackPt[300][3];
  Float_t fTTrackTPCmom[300][3];
  Float_t fTTrackEta[300][3];
  Short_t fTTrackCharge[300][3];
  Float_t fTTrackDCA[300][3];
  Float_t fTTrackITSspi[300][3];
  Float_t fTTrackITSsp[300][3];
  Float_t fTTrackITSsk[300][3];
  Float_t fTTrackTPCspi[300][3];
  Float_t fTTrackTPCsp[300][3];
  Float_t fTTrackTPCsk[300][3];
  Float_t fTTrackTOFse[300][3];
  Float_t fTTrackTOFspi[300][3];
  Float_t fTTrackTOFsp[300][3];
  Float_t fTTrackTOFsk[300][3];
  Int_t fTTrackNcl[300][3];
  Float_t fTTrackCrF[300][3];
  Int_t fTTrackShared[300][3];
  Float_t fTTrackTPCchi2[300][3];
  Bool_t fTTrackITStime[300][3];
  Bool_t fTTrackTOFtime[300][3];
  Bool_t fTTrackTPConly[300][3];
  Bool_t fTTrackITScomplementary[300][3];
  Bool_t fTTrackITSpure[300][3];
  Bool_t fTTrackGLOBAL[300][3];
  UInt_t fTTrackFilterBit[300][3];


ClassDef(AliOtonOmegaAnalysis,3)
};

#endif /* ALIOTONOMEGAANALYSIS_H_ */
