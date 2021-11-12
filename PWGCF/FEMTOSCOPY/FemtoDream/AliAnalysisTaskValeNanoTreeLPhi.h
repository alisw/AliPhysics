/*
 * AliAnalysisTaskValeNanoTreeLPhi.h
 *
 *  Created on: Nov 05, 2021
 *      Author: Valentina Mantovani Sarti
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKVALENANOTREELPHI_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKVALENANOTREELPHI_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCascade.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamControlSample.h"
#include "AliFemtoDreamBaseDump.h"

class AliAnalysisTaskValeNanoTreeLPhi : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskValeNanoTreeLPhi();
  AliAnalysisTaskValeNanoTreeLPhi(const char *name, bool isMC);
  virtual ~AliAnalysisTaskValeNanoTreeLPhi();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *){};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) { fEventCuts = evtCuts; };
  void SetLambdaCuts(AliFemtoDreamv0Cuts *cuts) { fLambdaCuts = cuts;}
  void SetAntiLambdaCuts(AliFemtoDreamv0Cuts *cuts){fAntiLambdaCuts = cuts;}
  void SetPosKaonCuts(AliFemtoDreamTrackCuts *trkCuts){fPosKaonCuts = trkCuts;}
  void SetNegKaonCuts(AliFemtoDreamTrackCuts *trkCuts){fNegKaonCuts = trkCuts;}
  void SetPhiCuts(AliFemtoDreamv0Cuts *phiCuts) { fPhiCuts = phiCuts; }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config){fConfig = config;}
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }
  void SetRunTaskLightWeight(bool light){fisLightWeight = light;}

  Bool_t FillKaon(AliFemtoDreamTrack *track);
  Bool_t FillProton(AliFemtoDreamTrack *track);
  Bool_t FillLambda(AliFemtoDreamv0 *v0);

private:
  AliAnalysisTaskValeNanoTreeLPhi(const AliAnalysisTaskValeNanoTreeLPhi &);
  AliAnalysisTaskValeNanoTreeLPhi &operator=(const AliAnalysisTaskValeNanoTreeLPhi &);
  void ResetGlobalTrackReference();
  //AOD
  //  void StoreGlobalTrackReference(AliAODTrack *track);
  //NANoAOD
  void StoreGlobalTrackReference(AliVTrack *track);
  bool fIsMC;                             //
  bool fUseOMixing;                       //
  bool fisLightWeight;                    //
  float fInvMassCutSBdown;                //
  float fInvMassCutSBup;                  //
  UInt_t fTrigger;                        //
  AliVEvent *fInputEvent;                 //! current event
  AliFemtoDreamEvent *fEvent;             //!
  AliFemtoDreamv0 *fLambda;               //!
  AliFemtoDreamv0 *fPhiParticle;          //!
  AliFemtoDreamTrack *fTrack;             //!
  AliFemtoDreamEventCuts *fEventCuts;     //
  AliFemtoDreamTrackCuts *fPosKaonCuts;   //
  AliFemtoDreamTrackCuts *fNegKaonCuts;   //
  AliFemtoDreamv0Cuts *fPhiCuts;          //
  AliFemtoDreamv0Cuts *fLambdaCuts;       //
  AliFemtoDreamv0Cuts *fAntiLambdaCuts;   //
  AliFemtoDreamCollConfig *fConfig;       //
  AliFemtoDreamPairCleaner *fPairCleaner; //!
  AliFemtoDreamPartCollection *fPartColl; //!
  AliFemtoDreamControlSample *fSample;    //!
  //AOD
  //AliAODTrack** fGTI;                   //!
  //NANoAOD
  AliVTrack **fGTI;                       //!
  int fTrackBufferSize;                   //
  TList *fResults;                        //!
  TList *fResultsQA;                      //!
  TList *fQA;                             //!
  TList *fEvtList;                        //!
  TList *fLambdaList;                     //!
  TList *fAntiLambdaList;                 //!
  TList *fKaonPlusList;                   //!
  TList *fKaonMinusList;                  //!
  TList *fPhiList;                        //!
  ///Tree
  TTree *fTree; //!
  Int_t fTRunNumber;
  Float_t fTVz;
  Int_t fTMult;
  
  ClassDef(AliAnalysisTaskValeNanoTreeLPhi, 1)
};

#endif PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKVALENANOTREELPHI_H_
