/*
 * AliAnalysisTaskOtonkdAOD.h
 *
 *  reCreated on: 20 Oct 2021
 *      Authors: bernhardhohlweger, Bhawani, Oton
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKOTONKDAOD_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKOTONKDAOD_H_
#include "Rtypes.h"
#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "TTree.h"


class AliAnalysisTaskOtonkdAOD : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskOtonkdAOD();
  AliAnalysisTaskOtonkdAOD(const char *name, bool isMC, bool isMCtruth, bool isIncludeSomeProtons, bool isPions, bool doFDpairing, bool DOpd);
  virtual ~AliAnalysisTaskOtonkdAOD();
  Float_t GetMass2sq(AliFemtoDreamTrack *track)const;
  float MeanTOFMassSqdDeuteron(AliFemtoDreamTrack *track) const;
  float SigmaTOFMassSqdDeuteron(AliFemtoDreamTrack *track) const;
  void InitHistograms(AliFemtoDreamTrackCuts *trkCuts, TString trkCutsName, TString MCName);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts = evtCuts;};
  void SetTrackCutsDeuteron(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsDeuteron = trkCuts;};
//  void SetTrackCutsDeuteronMass(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsDeuteronMass = trkCuts;};
  void SetTrackCutsAntiDeuteron(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiDeuteron = trkCuts;};
//  void SetTrackCutsAntiDeuteronMass(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiDeuteronMass = trkCuts;};
  void SetTrackCutsKaon(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsKaon = trkCuts;};
  void SetTrackCutsAntiKaon(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiKaon = trkCuts;};
  void SetTrackCutsProton(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsProton = trkCuts;};
  void SetTrackCutsAntiProton(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiProton = trkCuts;};
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig = config;};
  void SetRunTaskLightWeight(bool light) {
    fisLightWeight = light;
  }
  void SetSideband(float sigmaUp,float sigmalLow) {
    fSigmaUp = sigmaUp;
    fSigmaLow = sigmalLow;
    fdoSideband= true;
  }
  ;


  Bool_t FillKaon(AliFemtoDreamTrack *TheTrack);
  Bool_t FillProton(AliFemtoDreamTrack *TheTrack);
  Bool_t FillDeuteron(AliFemtoDreamTrack *TheTrack, int IsBkgProton = 0);



  private:
  AliAnalysisTaskOtonkdAOD(const AliAnalysisTaskOtonkdAOD &task);
  AliAnalysisTaskOtonkdAOD &operator=(const AliAnalysisTaskOtonkdAOD &task);
  void ResetGlobalTrackReference();
//AOD
  void StoreGlobalTrackReference(AliAODTrack *track);
//NANoAOD
//  void StoreGlobalTrackReference(AliVTrack *track);

  bool fisLightWeight;                      //
  int fTrackBufferSize;                     //
  bool fIsMC;                               //
  bool fIsMCtruth;                               //
  bool fdoFDpairing;                               //
  bool fDOpd;                               //
  bool fisIncludeSomeProtons;                               //
  bool fisPions;                               //
  bool fdoSideband;                         //
  float fSigmaUp;                           //
  float fSigmaLow;                          //
  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamEventCuts *fEventCuts;       //
  AliFemtoDreamTrackCuts *fTrackCutsDeuteron;  //
//  AliFemtoDreamTrackCuts *fTrackCutsDeuteronMass;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiDeuteron;  //
//  AliFemtoDreamTrackCuts *fTrackCutsAntiDeuteronMass;  //
  AliFemtoDreamTrackCuts *fTrackCutsKaon;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiKaon;  //
  AliFemtoDreamTrackCuts *fTrackCutsProton;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiProton;  //
  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
//AOD
  AliAODTrack** fGTI;           //!
//NANoAOD
//  AliVTrack **fGTI;  //!

  TList *fEvtList;//!
  TList *fKaonList;//!
  TList* fKaonMCList;//!
  TList *fAntiKaonList;//!
  TList* fAntiKaonMCList;//!
  TList *fProtonList;//!
  TList* fProtonMCList;//!
  TList *fAntiProtonList;//!
  TList* fAntiProtonMCList;//!
  TList *fDeuteronList;//!
  TList* fDeuteronMCList;//!
  TList *fAntiDeuteronList;//!
  TList* fAntiDeuteronMCList;//!
//  TList *fDeuteronNoTOFList;//!
//  TList* fDeuteronMCNoTOFList;//!
//  TList *fAntiDeuteronNoTOFList;//!
//  TList* fAntiDeuteronMCNoTOFList;//!
  TList *fResults;                          //!
  TList *fResultsQA;                        //!
  TH2F  *fDeuteronRestMass;                 //!
  TH2F  *fAntiDeuteronRestMass;             //!
//  TH2F  *fDeuteronRestMassNoTOF;            //!
//  TH2F  *fAntiDeuteronRestMassNoTOF;        //!

  TRandom3 r3;

  TTree* fTree;  //!
  Int_t fTRunNumber;
  Float_t fTVz;
  Int_t fTMult;
  Float_t fTSpher;

  const Int_t MAXKaonS = 30;
  Int_t fTnKaon;
  Float_t fTKaonP[30];
  Float_t fTKaonEta[30];
  Float_t fTKaonPx[30];
  Float_t fTKaonPy[30];
  Float_t fTKaonPz[30];
  Float_t fTKaonPTPC[30];
  Float_t fTKaonVPx[30];
  Float_t fTKaonVPy[30];
  Float_t fTKaonVPz[30];
  Float_t fTKaonPt[30];
  Float_t fTKaonmT[30];
  Float_t fTKaonTPCmom[30];
  Short_t fTKaonCharge[30];
  Float_t fTKaonITSsigma_e[30];
  Float_t fTKaonTPCsigma_e[30];
  Float_t fTKaonTOFsigma_e[30];
  Float_t fTKaonITSsigma_pi[30];
  Float_t fTKaonTPCsigma_pi[30];
  Float_t fTKaonTOFsigma_pi[30];
  Float_t fTKaonITSsigma_k[30];
  Float_t fTKaonTPCsigma_k[30];
  Float_t fTKaonTOFsigma_k[30];
  Float_t fTKaonITSsigma_p[30];
  Float_t fTKaonTPCsigma_p[30];
  Float_t fTKaonTOFsigma_p[30];
  Float_t fTKaonITSsigma_d[30];
  Float_t fTKaonTPCsigma_d[30];
  Float_t fTKaonTOFsigma_d[30];
  Float_t fTKaonDCA[30];
  Float_t fTKaonDCAz[30];
  Int_t fTKaonNcl[30];
  Int_t fTKaonShared[30];
  Float_t fTKaonTPCchi2[30];
  Bool_t fTKaonSPDtime[30];
  Bool_t fTKaonITStime[30];
  Bool_t fTKaonTOFtime[30];
  Bool_t fTKaonTPConly[30];
  Bool_t fTKaonITScomplementary[30];
  Bool_t fTKaonITSpure[30];
  Bool_t fTKaonGLOBAL[30];
  Float_t fTKaonPhi[30];
  Int_t fTKaonID[30];
  Bool_t fTKaonIs[30];
  Bool_t fTKaonIsFD[30];
  UInt_t fTKaonFilterBit[30];
  Int_t fTKaonPDG[30];
  Int_t fTKaonMotherWeak[30];
  Short_t fTKaonOrigin[30];


  const Int_t MAXDeuteronS = 10;
  Int_t fTnDeuteron;
  Float_t fTDeuteronP[10];
  Float_t fTDeuteronEta[10];
  Float_t fTDeuteronPx[10];
  Float_t fTDeuteronPy[10];
  Float_t fTDeuteronPz[10];
  Float_t fTDeuteronPTPC[10];
  Float_t fTDeuteronVPx[10];
  Float_t fTDeuteronVPy[10];
  Float_t fTDeuteronVPz[10];
  Float_t fTDeuteronPt[10];
  Float_t fTDeuteronmT[10];
  Float_t fTDeuteronTPCmom[10];
  Short_t fTDeuteronCharge[10];
  Float_t fTDeuteronITSsigma_e[10];
  Float_t fTDeuteronTPCsigma_e[10];
  Float_t fTDeuteronTOFsigma_e[10];
  Float_t fTDeuteronITSsigma_pi[10];
  Float_t fTDeuteronTPCsigma_pi[10];
  Float_t fTDeuteronTOFsigma_pi[10];
  Float_t fTDeuteronITSsigma_k[10];
  Float_t fTDeuteronTPCsigma_k[10];
  Float_t fTDeuteronTOFsigma_k[10];
  Float_t fTDeuteronITSsigma_p[10];
  Float_t fTDeuteronTPCsigma_p[10];
  Float_t fTDeuteronTOFsigma_p[10];
  Float_t fTDeuteronITSsigma_d[10];
  Float_t fTDeuteronTPCsigma_d[10];
  Float_t fTDeuteronTOFsigma_d[10];
  Float_t fTDeuteronDCA[10];
  Float_t fTDeuteronDCAz[10];
  Int_t fTDeuteronNcl[10];
  Int_t fTDeuteronShared[10];
  Float_t fTDeuteronTPCchi2[10];
  Bool_t fTDeuteronSPDtime[10];
  Bool_t fTDeuteronITStime[10];
  Bool_t fTDeuteronTOFtime[10];
  Bool_t fTDeuteronTPConly[10];
  Bool_t fTDeuteronITScomplementary[10];
  Bool_t fTDeuteronITSpure[10];
  Bool_t fTDeuteronGLOBAL[10];
  Float_t fTDeuteronPhi[10];
  Int_t fTDeuteronID[10];
  Float_t fTDeuteronTOFbeta[10];
  Int_t fTDeuteronPDG[10];
  Short_t fTDeuteronOrigin[10];


  ClassDef(AliAnalysisTaskOtonkdAOD, 6)
};
#endif
