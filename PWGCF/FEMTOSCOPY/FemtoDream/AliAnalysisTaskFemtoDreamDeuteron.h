/*
 * AliAnalysisTaskFemtoDreamDeuteron.h
 *
 *  Created on: 21 Mar 2018
 *      Author: bernhardhohlweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTODREAMDEUTERON_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTODREAMDEUTERON_H_
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


class AliAnalysisTaskFemtoDreamDeuteron : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFemtoDreamDeuteron();
  AliAnalysisTaskFemtoDreamDeuteron(const char *name, bool isMC);
  virtual ~AliAnalysisTaskFemtoDreamDeuteron();
  Float_t GetMass2sq(AliFemtoDreamTrack *track)const;
  float MeanTOFMassSqdDeuteron(AliFemtoDreamTrack *track) const;
  float SigmaTOFMassSqdDeuteron(AliFemtoDreamTrack *track) const;
  void InitHistograms(AliFemtoDreamTrackCuts *trkCuts, TString trkCutsName, TString MCName);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {};
  void SetEventCuts(AliFemtoDreamEventCuts *evtCuts) {fEventCuts = evtCuts;};
  void SetTrackCutsDeuteronDCA(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsDeuteronDCA = trkCuts;};
  void SetTrackCutsDeuteronMass(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsDeuteronMass = trkCuts;};
  void SetTrackCutsAntiDeuteronDCA(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiDeuteronDCA = trkCuts;};
  void SetTrackCutsAntiDeuteronMass(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiDeuteronMass = trkCuts;};
  void SetTrackCutsProtonDCA(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsProtonDCA = trkCuts;};
  void SetTrackCutsAntiProtonDCA(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiProtonDCA = trkCuts;};
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
  private:
  AliAnalysisTaskFemtoDreamDeuteron(const AliAnalysisTaskFemtoDreamDeuteron &task);
  AliAnalysisTaskFemtoDreamDeuteron &operator=(const AliAnalysisTaskFemtoDreamDeuteron &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fisLightWeight;                      //
  int fTrackBufferSize;                     //
  bool fIsMC;                               //
  bool fdoSideband;                         //
  float fSigmaUp;                           //
  float fSigmaLow;                          //
  AliFemtoDreamEvent *fEvent;               //!
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamEventCuts *fEventCuts;       //
  AliFemtoDreamTrackCuts *fTrackCutsDeuteronDCA;  //
  AliFemtoDreamTrackCuts *fTrackCutsDeuteronMass;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiDeuteronDCA;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiDeuteronMass;  //
  AliFemtoDreamTrackCuts *fTrackCutsProtonDCA;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiProtonDCA;  //
  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  AliAODTrack** fGTI;           //!
  TList *fEvtList;//!
  TList *fProtonList;//!
  TList* fProtonMCList;//!
  TList *fAntiProtonList;//!
  TList* fAntiProtonMCList;//!
  TList *fDeuteronList;//!
  TList* fDeuteronMCList;//!
  TList *fAntiDeuteronList;//!
  TList* fAntiDeuteronMCList;//!
  TList *fDeuteronNoTOFList;//!
  TList* fDeuteronMCNoTOFList;//!
  TList *fAntiDeuteronNoTOFList;//!
  TList* fAntiDeuteronMCNoTOFList;//!
  TList *fResults;                          //!
  TList *fResultsQA;                        //!
  TH2F  *fDeuteronRestMass;                 //!
  TH2F  *fAntiDeuteronRestMass;             //!
  TH2F  *fDeuteronRestMassNoTOF;            //!
  TH2F  *fAntiDeuteronRestMassNoTOF;        //!
  ClassDef(AliAnalysisTaskFemtoDreamDeuteron, 6)
};
#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTOTUTORIAL_H_ */
