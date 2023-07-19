/*
 * AliAnalysisTaskFemtoProtonKaonPlus.h
 *
 *  Created on: 2 May 2023
 *  Author: Italiano Luca
 */


#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTOPROTONKAONPLUS_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTOPROTONKAONPLUS_H_
#include "AliAnalysisTaskSE.h"
#include "AliConvEventCuts.h"
#include "AliEventCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamControlSample.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TChain.h"
#include "AliFemtoDreamBaseDump.h"

class AliAnalysisTaskFemtoProtonKaonPlus : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFemtoProtonKaonPlus();
  AliAnalysisTaskFemtoProtonKaonPlus(const char *name, bool isMC);
  virtual ~AliAnalysisTaskFemtoProtonKaonPlus();
  void InitHistograms(AliFemtoDreamTrackCuts *trkCuts, TString trkCutsName, TString MCName);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {};

  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts = evtCuts;};
  void SetTrackCutsKaon(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsKaon = trkCuts;};
  void SetTrackCutsAntiKaon(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiKaon = trkCuts;};
  void SetTrackCutsProton(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsProton = trkCuts;};
  void SetTrackCutsAntiProton(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiProton = trkCuts;};
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig = config;};
  void SetRunTaskLightWeight(bool light) {fisLightWeight = light;};
  void SetDoPairCleaning(bool DoPairCleaning){fDoPairCleaning = DoPairCleaning;};

  private:
  AliAnalysisTaskFemtoProtonKaonPlus(const AliAnalysisTaskFemtoProtonKaonPlus &task);
  AliAnalysisTaskFemtoProtonKaonPlus &operator=(const AliAnalysisTaskFemtoProtonKaonPlus &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fisLightWeight;                      //
  int fTrackBufferSize;                     //
  bool fIsMC;                               //
  bool fDoPairCleaning;			   //

  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamEventCuts *fEventCuts;       //
  AliFemtoDreamTrackCuts *fTrackCutsKaon;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiKaon;  //
  AliFemtoDreamTrackCuts *fTrackCutsProton;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiProton;  //
  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  AliAODTrack** fGTI;           //!
  TList *fEvtList;//!
  TList *fProtonList;//!
  TList* fProtonMCList;//!
  TList *fAntiProtonList;//!
  TList* fAntiProtonMCList;//!
  TList *fKaonList;//!
  TList* fKaonMCList;//!
  TList *fAntiKaonList;//!
  TList* fAntiKaonMCList;//!
  TList *fResults;                          //!
  TList *fResultsQA;                        //!

  ClassDef(AliAnalysisTaskFemtoProtonKaonPlus, 2) 
};
#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTOPROTONKAONPLUS_H_ */
