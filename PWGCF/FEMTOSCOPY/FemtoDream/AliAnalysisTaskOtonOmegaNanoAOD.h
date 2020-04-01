/*
 * AliAnalysisTaskOtonOmegaNanoAOD.h
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKOTONOMEGANANOAOD_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKOTONOMEGANANOAOD_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamCascade.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliOtonOmegaCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "TTree.h"

class AliAnalysisTaskOtonOmegaNanoAOD : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskOtonOmegaNanoAOD();
  AliAnalysisTaskOtonOmegaNanoAOD(const char* name, bool OmegaTreeFlag = false);
  //AliAnalysisTaskOtonOmegaNanoAOD(const AliAnalysisTaskOtonOmegaNanoAOD& analysis) = default;
  //AliAnalysisTaskOtonOmegaNanoAOD& operator=(const AliAnalysisTaskOtonOmegaNanoAOD& analysis) = default;
  virtual ~AliAnalysisTaskOtonOmegaNanoAOD();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  void SetRunTaskLightWeight(bool light) {
    fisLightWeight = light;
  }
  void SetEventCuts(AliFemtoDreamEventCuts* evtCuts) {
    fEventCuts = evtCuts;
  }
  void SetProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fProton = trkCuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fAntiProton = trkCuts;
  }

  void SetXiCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fXi = cascCuts;
  }
  void SetAntiXiCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fAntiXi = cascCuts;
  }
  void SetOmegaCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fOmega = cascCuts;
  }
  void SetAntiOmegaCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fAntiOmega = cascCuts;
  }
/*
  void SetXiCuts(AliOtonOmegaCascadeCuts* cascCuts) {
    fXi = cascCuts;
  }
  void SetAntiXiCuts(AliOtonOmegaCascadeCuts* cascCuts) {
    fAntiXi = cascCuts;
  }
  void SetOmegaCuts(AliOtonOmegaCascadeCuts* cascCuts) {
    fOmega = cascCuts;
  }
  void SetAntiOmegaCuts(AliOtonOmegaCascadeCuts* cascCuts) {
    fAntiOmega = cascCuts;
  }
*/
  void SetCorrelationConfig(AliFemtoDreamCollConfig* config) {
    fConfig=config;
  }

  Bool_t FillCascade(AliFemtoDreamCascade *TheCasc);
  Bool_t FillProton(AliFemtoDreamTrack *TheTrack);

 private:
  bool fisLightWeight;//
  bool fOmegaTreeFlag;                        //
  AliFemtoDreamEvent* fEvent;//!
  AliFemtoDreamEventCuts* fEventCuts;//
  TList* fEvtList;//!
  AliFemtoDreamTrack* fTrack;//!
  AliFemtoDreamTrackCuts* fProton;//
  TList* fProtonList;//!
  AliFemtoDreamTrackCuts* fAntiProton;//
  TList* fAntiProtonList;//!
  AliFemtoDreamCascade* fCascade;//!

  AliFemtoDreamCascadeCuts* fXi;//
  TList* fXiList;
  AliFemtoDreamCascadeCuts* fAntiXi;//
  TList* fAntiXiList;
  AliFemtoDreamCascadeCuts* fOmega;//
  TList* fOmegaList;
  AliFemtoDreamCascadeCuts* fAntiOmega;//
  TList* fAntiOmegaList;

/*
  AliOtonOmegaCascadeCuts* fXi;//
  TList* fXiList;
  AliOtonOmegaCascadeCuts* fAntiXi;//
  TList* fAntiXiList;
  AliOtonOmegaCascadeCuts* fOmega;//
  TList* fOmegaList;
  AliOtonOmegaCascadeCuts* fAntiOmega;//
  TList* fAntiOmegaList;
*/

  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  TList *fResults;//!
  TList *fResultsQA;//!
  int fTrackBufferSize;//
  TTree* fOmegaTree;  //!
  AliVTrack **fGTI;  //!

  //Tree Variables:
  /////////////////
  //event
  Int_t fTRunNumber;
  Float_t fTVx;
  Float_t fTVy;
  Float_t fTVz;
  Int_t fTMult;
  //protons:
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
  //cascades:
  const Int_t MAXCASCADES = 300;
  Int_t fTnCascade;
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
  //lambda (from cascade)
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
  //tracks (from cascade): 0 proton, 1 pion, 2 bachelor
  Float_t fTTrackP[300][3];
  Float_t fTTrackPx[300][3];
  Float_t fTTrackPy[300][3];
  Float_t fTTrackPz[300][3];
  Float_t fTTrackPt[300][3];
  Float_t fTTrackTPCmom[300][3];
  Float_t fTTrackEta[300][3];
  Short_t fTTrackCharge[300][3];
  Float_t fTTrackDCA[300][3];
  Float_t fTTrackTPCsigma[300][3];
  Float_t fTTrackTOFsigma[300][3];
  Int_t fTTrackNcl[300][3];
  Float_t fTTrackCrR[300][3];
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
  Float_t fTTrackPhi[300][3];
  Int_t fTTrackID[300][3];

  ClassDef(AliAnalysisTaskOtonOmegaNanoAOD,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskOtonOmegaNanoAOD_H_ */


