/*
 * AliAnalysisTaskNanoFemtoProtonPion.h
 *
 *  Created on: 11 Mar 2022
 *  Author: Lesch Marcel
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOFEMTOPROTONPION_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOFEMTOPROTONPION_H_
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
class AliVParticle;
class AliVTrack;

class AliAnalysisTaskNanoFemtoProtonPion : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoFemtoProtonPion();
  AliAnalysisTaskNanoFemtoProtonPion(const char *name, bool isMC);
  virtual ~AliAnalysisTaskNanoFemtoProtonPion();
  void InitHistograms(AliFemtoDreamTrackCuts *trkCuts, TString trkCutsName, TString MCName);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts = evtCuts;};
  void SetTrackCutsPion(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsPion = trkCuts;};
  void SetTrackCutsAntiPion(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiPion = trkCuts;};
  void SetTrackCutsProton(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsProton = trkCuts;};
  void SetTrackCutsAntiProton(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiProton = trkCuts;};
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig = config;};
  void SetRunTaskLightWeight(bool light) {
    fisLightWeight = light;
  }
  void SetDoPairCleaning(bool DoPairCleaning){ 
   fDoPairCleaning = DoPairCleaning;
  }

  private:
  AliAnalysisTaskNanoFemtoProtonPion(const AliAnalysisTaskNanoFemtoProtonPion &task);
  AliAnalysisTaskNanoFemtoProtonPion &operator=(const AliAnalysisTaskNanoFemtoProtonPion &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  bool fisLightWeight;                      //
  int fTrackBufferSize;                     //
  bool fIsMC;                               //
  bool fDoPairCleaning;			   //
  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamEventCuts *fEventCuts;       //
  AliFemtoDreamTrackCuts *fTrackCutsPion;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiPion;  //
  AliFemtoDreamTrackCuts *fTrackCutsProton;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiProton;  //
  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  AliVTrack** fGTI;           //!
  TList *fEvtList;//!
  TList *fProtonList;//!
  TList* fProtonMCList;//!
  TList *fAntiProtonList;//!
  TList* fAntiProtonMCList;//!
  TList *fPionList;//!
  TList* fPionMCList;//!
  TList *fAntiPionList;//!
  TList* fAntiPionMCList;//!
  TList *fResults;                          //!
  TList *fResultsQA;                        //!
  ClassDef(AliAnalysisTaskNanoFemtoProtonPion, 1)
};
#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKNANOFEMTOPROTONPION_H_ */
