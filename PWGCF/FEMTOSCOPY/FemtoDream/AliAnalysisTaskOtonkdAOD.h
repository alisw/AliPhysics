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
  AliAnalysisTaskOtonkdAOD(const char *name, bool isMC);
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
  Bool_t FillDeuteron(AliFemtoDreamTrack *TheTrack);


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


  TTree* fTree;  //!
  Int_t fTRunNumber;
  Float_t fTVz;
  Int_t fTMult;
  Float_t fTSpher;

  const Int_t MAXKaonS = 150;
  Int_t fTnKaon;
  Float_t fTKaonP[150];
  Float_t fTKaonEta[150];
  Float_t fTKaonPx[150];
  Float_t fTKaonPy[150];
  Float_t fTKaonPz[150];
  Float_t fTKaonVPx[150];
  Float_t fTKaonVPy[150];
  Float_t fTKaonVPz[150];
  Float_t fTKaonPt[150];
  Float_t fTKaonmT[150];
  Float_t fTKaonTPCmom[150];
  Short_t fTKaonCharge[150];
  Float_t fTKaonITSsigma_e[150];
  Float_t fTKaonTPCsigma_e[150];
  Float_t fTKaonTOFsigma_e[150];
  Float_t fTKaonITSsigma_pi[150];
  Float_t fTKaonTPCsigma_pi[150];
  Float_t fTKaonTOFsigma_pi[150];
  Float_t fTKaonITSsigma_k[150];
  Float_t fTKaonTPCsigma_k[150];
  Float_t fTKaonTOFsigma_k[150];
  Float_t fTKaonITSsigma_p[150];
  Float_t fTKaonTPCsigma_p[150];
  Float_t fTKaonTOFsigma_p[150];
  Float_t fTKaonITSsigma_d[150];
  Float_t fTKaonTPCsigma_d[150];
  Float_t fTKaonTOFsigma_d[150];
  Float_t fTKaonDCA[150];
  Int_t fTKaonNcl[150];
  Int_t fTKaonShared[150];
  Float_t fTKaonTPCchi2[150];
  Bool_t fTKaonITStime[150];
  Bool_t fTKaonTOFtime[150];
  Bool_t fTKaonTPConly[150];
  Bool_t fTKaonITScomplementary[150];
  Bool_t fTKaonITSpure[150];
  Bool_t fTKaonGLOBAL[150];
  Float_t fTKaonPhi[150];
  Int_t fTKaonID[150];
  Bool_t fTKaonIs[150];
  Bool_t fTKaonIsFD[150];
  UInt_t fTKaonFilterBit[150];

  const Int_t MAXPROTONS = 150;
  Int_t fTnProton;
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
  Float_t fTProtonTPCsigma[150];
  Float_t fTProtonTOFsigma[150];
  Float_t fTProtonDCA[150];
  Int_t fTProtonNcl[150];
  Int_t fTProtonShared[150];
  Float_t fTProtonTPCchi2[150];
  Bool_t fTProtonITStime[150];
  Bool_t fTProtonTOFtime[150];
  Bool_t fTProtonTPConly[150];
  Bool_t fTProtonITScomplementary[150];
  Bool_t fTProtonITSpure[150];
  Bool_t fTProtonGLOBAL[150];
  Float_t fTProtonPhi[150];
  Int_t fTProtonID[150];

  const Int_t MAXDeuteronS = 150;
  Int_t fTnDeuteron;
  Float_t fTDeuteronP[150];
  Float_t fTDeuteronEta[150];
  Float_t fTDeuteronPx[150];
  Float_t fTDeuteronPy[150];
  Float_t fTDeuteronPz[150];
  Float_t fTDeuteronVPx[150];
  Float_t fTDeuteronVPy[150];
  Float_t fTDeuteronVPz[150];
  Float_t fTDeuteronPt[150];
  Float_t fTDeuteronmT[150];
  Float_t fTDeuteronTPCmom[150];
  Short_t fTDeuteronCharge[150];
  Float_t fTDeuteronITSsigma_e[150];
  Float_t fTDeuteronTPCsigma_e[150];
  Float_t fTDeuteronTOFsigma_e[150];
  Float_t fTDeuteronITSsigma_pi[150];
  Float_t fTDeuteronTPCsigma_pi[150];
  Float_t fTDeuteronTOFsigma_pi[150];
  Float_t fTDeuteronITSsigma_k[150];
  Float_t fTDeuteronTPCsigma_k[150];
  Float_t fTDeuteronTOFsigma_k[150];
  Float_t fTDeuteronITSsigma_p[150];
  Float_t fTDeuteronTPCsigma_p[150];
  Float_t fTDeuteronTOFsigma_p[150];
  Float_t fTDeuteronITSsigma_d[150];
  Float_t fTDeuteronTPCsigma_d[150];
  Float_t fTDeuteronTOFsigma_d[150];
  Float_t fTDeuteronDCA[150];
  Int_t fTDeuteronNcl[150];
  Int_t fTDeuteronShared[150];
  Float_t fTDeuteronTPCchi2[150];
  Bool_t fTDeuteronITStime[150];
  Bool_t fTDeuteronTOFtime[150];
  Bool_t fTDeuteronTPConly[150];
  Bool_t fTDeuteronITScomplementary[150];
  Bool_t fTDeuteronITSpure[150];
  Bool_t fTDeuteronGLOBAL[150];
  Float_t fTDeuteronPhi[150];
  Int_t fTDeuteronID[150];
  Float_t fTDeuteronTOFbeta[150];
  Bool_t fTDeuteronIs[150];
  Bool_t fTDeuteronIsFD[150];


  ClassDef(AliAnalysisTaskOtonkdAOD, 6)
};
#endif
