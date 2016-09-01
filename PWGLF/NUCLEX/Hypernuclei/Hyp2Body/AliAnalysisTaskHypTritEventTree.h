/// \class AliAnalysisTaskHypTritEventTree
/// \brief Hypertriton Analysis in two particle decay channel
///
/// Hypertriton candidates are identified using the on-the-fly V0 finder.
/// Events with Hypertriton candidates are filled in a tree
/// using the AliReducedHypTritEvent class.
///
/// \author Lukas Kreis <lukas.kreis@cern.ch>, GSI
/// \date Sep 1, 2016

#ifndef ALIANALYSISTASKHYPTRITEVENTTREE_H
#define ALIANALYSISTASKHYPTRITEVENTTREE_H

class TH1F;
class TH2F;
class THnSparse;
class AliESDEvent;
class AliESDpid;
class AliESDtrackCuts;
class AliESDv0;
class AliESDVertex;
class AliESDInputHandler;
class AliESDtrack;
class AliReducedHypTritEvent;

#include "AliReducedHypTritEvent.h"
#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"


class AliAnalysisTaskHypTritEventTree : public AliAnalysisTaskSE {

 public:
  AliAnalysisTaskHypTritEventTree();
  AliAnalysisTaskHypTritEventTree(const char *name);
  virtual ~AliAnalysisTaskHypTritEventTree();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(const Option_t*);
  void SetUseEtaCorrection(Bool_t eta = kTRUE) {fEtaCorrection = eta;};
  void SetUseMultiplicityCorrection(Bool_t mult = kTRUE) {fMultiplicityCorrection = mult;};

 private:
  AliESDEvent         *fEvent;
  AliESDInputHandler  *fInputHandler;
  AliESDtrackCuts     *fTrackCutsV0;
  AliESDpid           *fPID;
  Bool_t               fMCtrue;
  TH2F                *fHistdEdx;
  TH2F                *fHistdEdxProton;
  TH2F                *fHistdEdxDeuteron;
  TH2F                *fHistdEdxTriton;
  TH2F                *fHistdEdxHelium3;
  TH2F                *fHistdEdxHypTriton;
  TH2F                *fHistdEdxHypTritonAnti;
  TH1F                *fHistInvMassHypTriton;
  TH1F                *fHistCentrality;
  TH1F                *fHistTrigger;
  TH1F                *fHistNumEvents;
  TTree               *fTree;
  TTree               *fTreeMCGen;
  Double_t             fPosVx;
  Double_t             fPosVy;
  Double_t             fPosVz;
  Int_t                fMCGenRec[40];
  TObjArray           *fMCGenRecArray;
  AliReducedHypTritEvent *fReducedEvent;
  AliReducedHypTritEvent *fReducedEventMCGen;
  Bool_t               fEtaCorrection;
  Bool_t               fMultiplicityCorrection;
  TList               *fHistogramList;

  void MCStackLoop(AliStack *stack);
  Bool_t TriggerSelection();


  AliAnalysisTaskHypTritEventTree(const AliAnalysisTaskHypTritEventTree&);
  AliAnalysisTaskHypTritEventTree &operator=(const AliAnalysisTaskHypTritEventTree&);
  ClassDef(AliAnalysisTaskHypTritEventTree, 2);
};

#endif
