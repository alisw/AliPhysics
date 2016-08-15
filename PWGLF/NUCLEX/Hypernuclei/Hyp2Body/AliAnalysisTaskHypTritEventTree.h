#ifndef AliAnalysisTaskHypTritEventTree_H
#define AliAnalysisTaskHypTritEventTree_H

// Task searching for hypertriton with the V0 finder
// Lukas Kreis

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

 private:
  AliESDEvent         *fEvent;
  AliESDInputHandler  *fInputHandler;
  AliESDtrackCuts     *fTrackCutsV0;
  AliESDpid           *fPID;
  Bool_t               fMCtrue;
  TH2F                *fHistdEdx;
  TH2F                *fHistdEdxDeuteron;
  TH2F                *fHistdEdxTriton;
  TH2F                *fHistdEdxHelium3;
  TH2F                *fHistdEdxHypTriton;
  TH2F                *fHistdEdxHypTritonAnti;
  TH1F                *fHistInvMassHypTriton;
  TH1F                *fHistInvMassHypTritonMC;
  TH1F                *fHistInvMassHypTritonMCAssoc;
  TH1F                *fHistPtHypTriton;
  TH1F                *fHistPtHypTritonMC;
  TH1F                *fHistctHypTritonMC;
  TH1F                *fHistCentrality;
  TH1F                *fHistTrigger;
  TH1F                *fHistPtHypTritonMCAssoc;
  TH2F                *fHistdEdxHelium3NSigma;
  TTree               *fTree;
  TTree               *fTreeMCGen;
  Double_t             fPosVx;
  Double_t             fPosVy;
  Double_t             fPosVz;
  Int_t                fMCGenRec[40];
  TObjArray        *fMCGenRecArray;
  AliReducedHypTritEvent *fReducedEvent;
  AliReducedHypTritEvent *fReducedEventMCGen;

  TList *fOutputContainer;

  void MCStackLoop(AliStack *stack);
  Bool_t TriggerSelection();


  AliAnalysisTaskHypTritEventTree(const AliAnalysisTaskHypTritEventTree&);
  AliAnalysisTaskHypTritEventTree &operator=(const AliAnalysisTaskHypTritEventTree&);
  ClassDef(AliAnalysisTaskHypTritEventTree, 1);
};

#endif
