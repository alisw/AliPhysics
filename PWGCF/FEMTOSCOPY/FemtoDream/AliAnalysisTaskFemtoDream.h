/*
 * AliAnalysisTaskFemtoPlotration.h
 *
 *  Created on: 24 Nov 2017
 *      Author: bernhardhohlweger
 */

#ifndef ALIANALYSISTASKFEMTODREAM_H_
#define ALIANALYSISTASKFEMTODREAM_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamAnalysis.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCascadeCuts.h"
#include "TList.h"
class AliAnalysisTaskFemtoDream : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFemtoDream();
  AliAnalysisTaskFemtoDream(const char *name, bool isESD, bool isMC);
  AliAnalysisTaskFemtoDream(const AliAnalysisTaskFemtoDream& task);
  AliAnalysisTaskFemtoDream& operator=(const AliAnalysisTaskFemtoDream& task);
  virtual ~AliAnalysisTaskFemtoDream();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {
  }
  ;
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
  void SetTrackCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCuts = cuts;
  }
  ;
  void SetAntiTrackCuts(AliFemtoDreamTrackCuts *cuts) {
    fAntiTrackCuts = cuts;
  }
  ;
  void Setv0Cuts(AliFemtoDreamv0Cuts *cuts) {
    fv0Cuts = cuts;
  }
  ;
  void SetAntiv0Cuts(AliFemtoDreamv0Cuts *cuts) {
    fAntiv0Cuts = cuts;
  }
  ;
  void SetCascadeCuts(AliFemtoDreamCascadeCuts *cuts) {
    fCascCuts = cuts;
  }
  ;
  void SetAntiCascadeCuts(AliFemtoDreamCascadeCuts *cuts) {
    fAntiCascCuts = cuts;
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
 private:
  int fTrackBufferSize;                     //
  bool fESDAnalysis;                        //
  bool fMinBookingME;                       //
  bool fMinBookingSample;                   //
  bool fMVPileUp;                           //
  bool fEvtCutQA;                           //
  bool fIsMC;                               //
  AliFemtoDreamAnalysis *fAnalysis;         //!
  TList *fQA;                               //!
  AliFemtoDreamEventCuts *fEvtCuts;         //
  TList *fEvtHistList;                      //!
  AliFemtoDreamTrackCuts *fTrackCuts;       //
  TList *fTrackCutHistList;                 //!
  TList *fTrackCutHistMCList;               //!
  AliFemtoDreamTrackCuts *fAntiTrackCuts;   //
  TList *fAntiTrackCutHistList;             //!
  TList *fAntiTrackCutHistMCList;           //!
  AliFemtoDreamv0Cuts *fv0Cuts;             //
  TList *fv0CutHistList;                    //!
  TList *fv0CutHistMCList;                  //!
  AliFemtoDreamv0Cuts *fAntiv0Cuts;         //
  TList *fAntiv0CutHistList;                //!
  TList *fAntiv0CutHistMCList;              //!
  AliFemtoDreamCascadeCuts *fCascCuts;      //
  TList *fCascCutList;                      //!
  TList *fCascCutMCList;                    //!
  AliFemtoDreamCascadeCuts *fAntiCascCuts;  //
  TList *fAntiCascCutList;                  //!
  TList *fAntiCascCutMCList;                //!
  AliFemtoDreamCollConfig *fConfig;         //
  TList *fResults;                          //!
  TList *fResultQA;                         //!
  TList *fResultsSample;                    //!
  TList *fResultQASample;					//!
ClassDef(AliAnalysisTaskFemtoDream,4)
};

#endif /* ALIANALYSISTASKFEMTODREAM_H_ */
