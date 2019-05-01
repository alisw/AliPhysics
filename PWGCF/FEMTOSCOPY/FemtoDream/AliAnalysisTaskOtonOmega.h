/*
 * AliAnalysisTaskFemtoPlotration.h
 *
 *  Created on: 24 Nov 2017
 *      Author: bernhardhohlweger
 */

#ifndef ALIANALYSISTASKOTONOMEGA_H_
#define ALIANALYSISTASKOTONOMEGA_H_
#include "AliAnalysisTaskSE.h"
#include "AliOtonOmegaAnalysis.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliOtonOmegaCascadeCuts.h"
#include "TList.h"
#include "TTree.h"
class AliAnalysisTaskOtonOmega : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskOtonOmega();
  AliAnalysisTaskOtonOmega(const char *name, bool isESD, bool isMC, bool CascadeTreeFlag, bool OmegaTreeFlag);
  AliAnalysisTaskOtonOmega(const AliAnalysisTaskOtonOmega& task);
  AliAnalysisTaskOtonOmega& operator=(const AliAnalysisTaskOtonOmega& task);
  virtual ~AliAnalysisTaskOtonOmega();
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
  void SetCascadeCuts(AliOtonOmegaCascadeCuts *cuts) {
    fCascCuts = cuts;
  }
  ;
  void SetAntiCascadeCuts(AliOtonOmegaCascadeCuts *cuts) {
    fAntiCascCuts = cuts;
  }
  ;

  void SetCascadeOmegaCuts(AliOtonOmegaCascadeCuts *cuts) {
    fCascOmegaCuts = cuts;
  }
  ;
  void SetAntiCascadeOmegaCuts(AliOtonOmegaCascadeCuts *cuts) {
    fAntiCascOmegaCuts = cuts;
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
  bool fCascadeTreeFlag;                        //
  bool fOmegaTreeFlag;                        //
  bool fMinBookingME;                       //
  bool fMinBookingSample;                   //
  bool fMVPileUp;                           //
  bool fEvtCutQA;                           //
  bool fIsMC;                               //
  AliOtonOmegaAnalysis *fAnalysis;         //!
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
  AliOtonOmegaCascadeCuts *fCascCuts;      //
  TList *fCascCutList;                      //!
  TList *fCascCutMCList;                    //!
  AliOtonOmegaCascadeCuts *fAntiCascCuts;  //
  TList *fAntiCascCutList;                  //!
  TList *fAntiCascCutMCList;                //!
  AliOtonOmegaCascadeCuts *fCascOmegaCuts;      //
  TList *fCascOmegaCutList;                      //!
  TList *fCascOmegaCutMCList;                    //!
  AliOtonOmegaCascadeCuts *fAntiCascOmegaCuts;  //
  TList *fAntiCascOmegaCutList;                  //!
  TList *fAntiCascOmegaCutMCList;                //!

  AliFemtoDreamCollConfig *fConfig;         //
  TList *fResults;                          //!
  TList *fResultQA;                         //!
  TList *fResultsSample;                    //!
  TList *fResultQASample;					//!
  TTree* fTreeCascade;  //!
  TTree* fTreeOmega;  //!
ClassDef(AliAnalysisTaskOtonOmega,4)
};

#endif /* ALIANALYSISTASKOTONOMEGA_H_ */
