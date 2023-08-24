/*
 * AliAnalysisTaskOtonXx.h
 *
 *  reCreated on: 20 Oct 2021
 *      Authors: bernhardhohlweger, Bhawani, Oton
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKOTONXX_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKOTONXX_H_
#include "Rtypes.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCascade.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliOtonOmegaCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "TTree.h"


class AliAnalysisTaskOtonXx : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskOtonXx();
  AliAnalysisTaskOtonXx(const char *name, bool doFDpairing, bool isMC, bool isMCtruth, bool isOmega, bool isPi, bool OnlyXi);
  virtual ~AliAnalysisTaskOtonXx();
  void InitHistograms(AliFemtoDreamTrackCuts *trkCuts, TString trkCutsName, TString MCName);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts = evtCuts;};
  void SetTrackCutsKaon(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsKaon = trkCuts;};
  void SetTrackCutsAntiKaon(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiKaon = trkCuts;};
  void SetXiCuts(AliFemtoDreamCascadeCuts* cascCuts) { fCutsXi = cascCuts; }
  void SetAntiXiCuts(AliFemtoDreamCascadeCuts* cascCuts) { fCutsAntiXi = cascCuts;}
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig = config;};
  void SetRunTaskLightWeight(bool light) {
    fisLightWeight = light;
  }
  Bool_t FillKaon(AliFemtoDreamTrack *TheTrack);
  Bool_t FillXi(AliFemtoDreamCascade *TheCasc, bool isomega=false);

  private:
  AliAnalysisTaskOtonXx(const AliAnalysisTaskOtonXx &task);
  AliAnalysisTaskOtonXx &operator=(const AliAnalysisTaskOtonXx &task);
  void ResetGlobalTrackReference();
//AOD
//  void StoreGlobalTrackReference(AliAODTrack *track);
//NANoAOD
  void StoreGlobalTrackReference(AliVTrack *track);

  bool fisLightWeight;                      //
  int fTrackBufferSize;                     //
  bool fIsMC;                               //
  bool fIsMCtruth;                               //
  bool fdoFDpairing;                               //
  bool fisOmega;                               //
  bool fisPi;                               //
  bool fOnlyXi;                               //
  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamEventCuts *fEventCuts;       //
  AliFemtoDreamTrackCuts *fTrackCutsKaon;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiKaon;  //
  AliFemtoDreamCascadeCuts *fCutsXi;  //
  AliFemtoDreamCascade *fCascade;  //
  AliFemtoDreamCascadeCuts *fCutsAntiXi;  //
  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
//AOD
//  AliAODTrack** fGTI;           //!
//NANoAOD
  AliVTrack **fGTI;  //!

  TList *fEvtList;//!
  TList *fKaonList;//!
  TList* fKaonMCList;//!
  TList *fAntiKaonList;//!
  TList* fAntiKaonMCList;//!
  TList *fXiList;//!
  TList *fAntiXiList;//!
  TList *fResults;                          //!
  TList *fResultsQA;                        //!

  TRandom3 r3;

  TTree* fTree;  //!
  Int_t fTRunNumber;
  Float_t fTVz;
  Int_t fTMult;
  Float_t fTSpher;

  const Int_t MAXKaonS = 300;
  Int_t fTnKaon;
  Float_t fTKaonP[300];
  Float_t fTKaonPx[300];
  Float_t fTKaonPy[300];
  Float_t fTKaonPz[300];
  Float_t fTKaonPTPC[300];
  Float_t fTKaonVPx[300];
  Float_t fTKaonVPy[300];
  Float_t fTKaonVPz[300];
  Float_t fTKaonPt[300];
  Float_t fTKaonmT[300];
  Float_t fTKaonTPCmom[300];
  Short_t fTKaonCharge[300];
  Float_t fTKaonITSsigma_e[300];
  Float_t fTKaonTPCsigma_e[300];
  Float_t fTKaonTOFsigma_e[300];
  Float_t fTKaonITSsigma_pi[300];
  Float_t fTKaonTPCsigma_pi[300];
  Float_t fTKaonTOFsigma_pi[300];
  Float_t fTKaonITSsigma_k[300];
  Float_t fTKaonTPCsigma_k[300];
  Float_t fTKaonTOFsigma_k[300];
  Float_t fTKaonITSsigma_p[300];
  Float_t fTKaonTPCsigma_p[300];
  Float_t fTKaonTOFsigma_p[300];
  Float_t fTKaonITSsigma_d[300];
  Float_t fTKaonTPCsigma_d[300];
  Float_t fTKaonTOFsigma_d[300];
  Float_t fTKaonDCA[300];
  Float_t fTKaonDCAz[300];
  Int_t fTKaonNcl[300];
  Int_t fTKaonShared[300];
  Float_t fTKaonTPCchi2[300];
  Bool_t fTKaonSPDtime[300];
  Bool_t fTKaonITStime[300];
  Bool_t fTKaonTOFtime[300];
  Bool_t fTKaonTPConly[300];
  Bool_t fTKaonITScomplementary[300];
  Bool_t fTKaonITSpure[300];
  Bool_t fTKaonGLOBAL[300];
  Float_t fTKaonPhi[300];
  Int_t fTKaonID[300];
  Bool_t fTKaonIs[300];
  Bool_t fTKaonIsFD[300];
  UInt_t fTKaonFilterBit[300];
  Int_t fTKaonPDG[300];
  Int_t fTKaonMotherWeak[300];
  Short_t fTKaonOrigin[300];
  Int_t fTKaonMotherID[300];

  const Int_t MAXXiS = 10;
  Int_t fTnXi;
  Short_t fTXiCharge[10];
  Float_t fTXiDCA[10];
  Float_t fTXiDaughtersDCA[10];
  Float_t fTXiMass[10];
  Float_t fTXiXiMass[10];
  Float_t fTXiOmegaMass[10];
  Float_t fTXiVr[10];
  Float_t fTXiPA[10];
  Float_t fTXiLambdaDCA[10];
  Float_t fTXiLambdaDaughtersDCA[10];
  Float_t fTXiLambdaMass[10];
  Float_t fTXiLambdaK0Mass[10];
  Float_t fTXiLambdaVr[10];
  Float_t fTXiLambdaPA[10];
  Short_t fTXiTrackCharge[10][3];
  Float_t fTXiTrackPx[10][3];
  Float_t fTXiTrackPy[10][3];
  Float_t fTXiTrackPz[10][3];
  Int_t fTXiTrackID[10][3];
  Float_t fTXiTrackTPCmom[10][3];
  Float_t fTXiTrackDCA[10][3];
  Float_t fTXiTrackTPCsigma[10][3];
  Float_t fTXiTrackTOFsigma[10][3];
  Int_t fTXiTrackNcl[10][3];
  Float_t fTXiTrackCrR[10][3];
  Float_t fTXiTrackCrF[10][3];
  Float_t fTXiTrackTPCchi2[10][3];
  Bool_t fTXiTrackSPDtime[10][3];
  Bool_t fTXiTrackITStime[10][3];
  Bool_t fTXiTrackTOFtime[10][3];
  Int_t fTXiMotherID[10];
   Int_t fTXiPDG[10];
   Int_t fTXiMotherPDG[10];
   Int_t fTXiMotherWeak[10];
   Int_t fTXiOrigin[10];

  // ClassDef 6 ????
  ClassDef(AliAnalysisTaskOtonXx, 6)
};
#endif
