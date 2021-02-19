/*
 * AliAnalysisTaskGeorgiosNTuple.h
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKGEORGIOSNTUPLE_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKGEORGIOSNTUPLE_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCascade.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliOtonOmegaCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "TTree.h"

#define MONTECARLO

class AliAnalysisTaskGeorgiosNTuple : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskGeorgiosNTuple();
  AliAnalysisTaskGeorgiosNTuple(const char* name, bool isMC);
  //AliAnalysisTaskGeorgiosNTuple(const AliAnalysisTaskGeorgiosNTuple& analysis) = default;
  //AliAnalysisTaskGeorgiosNTuple& operator=(const AliAnalysisTaskGeorgiosNTuple& analysis) = default;
  virtual ~AliAnalysisTaskGeorgiosNTuple();
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

//Lambda
  void SetLambdaCuts(AliFemtoDreamv0Cuts* v0Cuts){
    fLambda = v0Cuts;
  }
  void SetAntiLambdaCuts(AliFemtoDreamv0Cuts* v0Cuts){
    fAntiLambda = v0Cuts;
  }
//Xi
  void SetXiBGRCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fXiBGR = cascCuts;
  }
  void SetAntiXiBGRCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fAntiXiBGR = cascCuts;
  }
  void SetXiCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fXi = cascCuts;
  }
  void SetAntiXiCuts(AliFemtoDreamCascadeCuts* cascCuts) {
    fAntiXi = cascCuts;
  }
  void SetCorrelationConfig(AliFemtoDreamCollConfig* config) {
    fConfig=config;
  }

  Bool_t Fillv0(AliFemtoDreamv0 *Thev0, int Thev0Charge);
  Bool_t FillCascade(AliFemtoDreamCascade *TheCasc);

 private:
  bool fisLightWeight;//
  bool fisMC;                        //
  AliFemtoDreamEvent* fEvent;//!
  AliFemtoDreamEventCuts* fEventCuts;//
  TList* fEvtList;//!


//Lambda
  AliFemtoDreamv0* fv0;
  AliFemtoDreamv0Cuts* fLambda;
  TList* fLambdaList;
  AliFemtoDreamv0Cuts* fAntiLambda;
  TList* fAntiLambdaList;
//Xi
  AliFemtoDreamCascade* fCascade;//!
  AliFemtoDreamCascadeCuts* fXi;//
  TList* fXiList;
  AliFemtoDreamCascadeCuts* fAntiXi;//
  TList* fAntiXiList;
  AliFemtoDreamCascadeCuts* fXiBGR;//
  TList* fXiBGRList;
  AliFemtoDreamCascadeCuts* fAntiXiBGR;//
  TList* fAntiXiBGRList;

#ifdef MONTECARLO
  TList* fLambdaMCList;
  TList* fAntiLambdaMCList;
  TList* fXiMCList;
  TList* fAntiXiMCList;
  TList* fXiBGRMCList;
  TList* fAntiXiBGRMCList;
#endif
  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  TList *fResults;//!
  TList *fResultsQA;//!
  int fTrackBufferSize;//
  TTree* fGeorgiosTree;  //!
  AliVTrack **fGTI;  //!

  //Tree Variables:
  /////////////////
  //event
  Int_t fTRunNumber;
  Float_t fTVx;
  Float_t fTVy;
  Float_t fTVz;
  Int_t fTMult;

 //v0
  const Int_t MAXv0 = 300;
  Int_t fTnv0;
  Bool_t fTv0Flag0[300];
  Bool_t fTv0Flag1[300];
  //Bool_t fTv0Flag2[300];      //only two daughter particles
  Float_t fTv0Px[300];
  Float_t fTv0Py[300];
  Float_t fTv0Pz[300];
  Float_t fTv0mT[300];
  Short_t fTv0Charge[300];
  Float_t fTv0DCA[300];
  Float_t fTv0DaughtersDCA[300];
//  Float_t fTv0LambdaMass[300];
  //Float_t fTv0SigmaMass[300];
  Float_t fTv0LambdaVr[300];
  Float_t fTv0LambdaPA[300];
  Float_t fTv0Vx[300];
  Float_t fTv0Vy[300];
  Float_t fTv0Vz[300];

//tracks from v0: 0 proton, 1 pion
  Float_t fTTrackv0Px[300][2];
  Float_t fTTrackv0Py[300][2];
  Float_t fTTrackv0Pz[300][2];
  Float_t fTTrackv0TPCmom[300][2];
  //Float_t fTTrackv0Eta[300][2];
  Short_t fTTrackv0Charge[300][2];
  Float_t fTTrackv0DCA[300][2];
  Float_t fTTrackv0TPCsigma[300][2];
  //Float_t fTTrackv0TOFsigma[300][2];

  Int_t fTTrackv0Ncl[300][2];
/* 
  Float_t fTTrackv0CrR[300][2];
  Float_t fTTrackv0CrF[300][2];
  Int_t fTTrackv0Shared[300][2];
  Float_t fTTrackv0TPCchi2[300][2];
  Bool_t fTTrackv0ITStime[300][2];
  Bool_t fTTrackv0TOFtime[300][2];
  Bool_t fTTrackv0TPConly[300][2];
  Bool_t fTTrackv0ITScomplementary[300][2];
  Bool_t fTTrackv0ITSpure[300][2];
  Bool_t fTTrackv0GLOBAL[300][2];
  UInt_t fTTrackv0FilterBit[300][2];
*/
//  Float_t fTTrackv0Phi[300][2];
  Int_t fTTrackv0ID[300][2];

#ifdef MONTECARLO 
   Float_t fTv0LambdaPxMC[300];
   Float_t fTv0LambdaPyMC[300];
   Float_t fTv0LambdaPzMC[300];
   Int_t fTv0LambdaPDG[300];
   Int_t fTv0LambdaMotherPDG[300];
   Int_t fTv0LambdaMotherWeakPDG[300];
   Int_t fTv0LambdaOrigin[300];
   Float_t fTTrackv0PxMC[300][2];
   Float_t fTTrackv0PyMC[300][2];
   Float_t fTTrackv0PzMC[300][2];
   Int_t fTTrackv0PDG[300][2];
   Int_t fTTrackv0MotherPDG[300][2];
   Int_t fTTrackv0MotherWeakPDG[300][2];
   Int_t fTTrackv0Origin[300][2];
#endif

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
  Float_t fTCascademT[300];
  Short_t fTCascadeCharge[300];
  Float_t fTCascadeDCA[300];
  Float_t fTCascadeDaughtersDCA[300];
//  Float_t fTCascadeXiMass[300];     //
//  Float_t fTCascadeOmegaMass[300];
  Float_t fTCascadeVx[300];
  Float_t fTCascadeVy[300];
  Float_t fTCascadeVz[300];
  Float_t fTCascadeVr[300];
  Float_t fTCascadePA[300];
  //lambda (from cascade)
  Float_t fTLambdaDCA[300];
  Float_t fTLambdaDaughtersDCA[300];
//  Float_t fTLambdaMass[300];
//  Float_t fTLambdaK0Mass[300];
  Float_t fTLambdaVx[300];
  Float_t fTLambdaVy[300];
  Float_t fTLambdaVz[300];
  Float_t fTLambdaVr[300];
  Float_t fTLambdaPA[300];
  //tracks (from cascade): 0 proton, 1 pion, 2 bachelor
  Float_t fTTrackPx[300][3];
  Float_t fTTrackPy[300][3];
  Float_t fTTrackPz[300][3];
  Float_t fTTrackTPCmom[300][3];
  //Float_t fTTrackEta[300][3];
  Short_t fTTrackCharge[300][3];
  Float_t fTTrackDCA[300][3];
  Float_t fTTrackTPCsigma[300][3];
//  Float_t fTTrackTOFsigma[300][3];
  Int_t fTTrackNcl[300][3];
/*
 *
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
*/
//  Float_t fTTrackPhi[300][3];
  Int_t fTTrackID[300][3];

#ifdef MONTECARLO 
   Float_t fTCascadePxMC[300];
   Float_t fTCascadePyMC[300];
   Float_t fTCascadePzMC[300];
   Int_t fTCascadePDG[300];
   Int_t fTCascadeMotherPDG[300];
   //Int_t fTCascadeMotherWeakPDG[300];
   Int_t fTCascadeOrigin[300];
   Float_t fTTrackPxMC[300][3];
   Float_t fTTrackPyMC[300][3];
   Float_t fTTrackPzMC[300][3];
   Int_t fTTrackPDG[300][3];
   Int_t fTTrackMotherPDG[300][3];
   Int_t fTTrackMotherWeakPDG[300][3];
   Int_t fTTrackOrigin[300][3];
#endif

  ClassDef(AliAnalysisTaskGeorgiosNTuple,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_AliAnalysisTaskGeorgiosNTuple_H_ */
