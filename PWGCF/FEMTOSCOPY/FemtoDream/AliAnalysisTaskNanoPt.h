#ifndef AliAnalysisTaskNanoPt_H
#define AliAnalysisTaskNanoPt_H
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

class AliAnalysisTaskNanoPt : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoPt();
  AliAnalysisTaskNanoPt(const char *name, const bool isMC);
  virtual ~AliAnalysisTaskNanoPt();
  Float_t GetMass2sq(AliFemtoDreamTrack *track) const;
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  void SetEventCuts(AliFemtoDreamEventCuts *cuts) {
    fEvtCuts = cuts;
  }
  void SetProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fProtonTrack = cuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fAntiProtonTrack = cuts;
  }
  void SetDeuteronCuts(AliFemtoDreamTrackCuts *cuts) {
    fDeuteronTrack = cuts;
  }
  void SetAntiDeuteronCuts(AliFemtoDreamTrackCuts *cuts) {
    fAntiDeuteronTrack = cuts;
  }
  void SetDeuteronCutsNoTOF(AliFemtoDreamTrackCuts *cuts) {
    fDeuteronTrackNoTOF = cuts;
  }
  void SetAntiDeuteronCutsNoTOF(AliFemtoDreamTrackCuts *cuts) {
    fAntiDeuteronTrackNoTOF = cuts;
  }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {
    fConfig = config;
  }
  void SetUseDumpster(bool use) {
    fUseDumpster = use;
  }
  void SetMCTruth(bool mct) {
    cout<<"passed argument to MC Truth setter is " << mct<<endl;
    fIsMCTruth = mct;
  }

 private:
  AliAnalysisTaskNanoPt(
    const AliAnalysisTaskNanoPt &task);
  AliAnalysisTaskNanoPt &operator=(
    const AliAnalysisTaskNanoPt &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  int fTrackBufferSize;  //
  bool fIsMC;   //
  bool fIsMCTruth; //
  bool fUseDumpster;  //
  AliVEvent *fInputEvent; //! current event
  AliFemtoDreamEvent *fEvent; //!
  AliFemtoDreamEventCuts *fEvtCuts; //
  AliFemtoDreamTrack *fTrack; //
  AliFemtoDreamTrackCuts *fProtonTrack;  //
  AliFemtoDreamTrackCuts *fAntiProtonTrack; //
  AliFemtoDreamTrackCuts *fDeuteronTrack; //
  AliFemtoDreamTrackCuts *fAntiDeuteronTrack; //
  AliFemtoDreamTrackCuts *fDeuteronTrackNoTOF; //
  AliFemtoDreamTrackCuts *fAntiDeuteronTrackNoTOF; //
  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner; //!
  AliFemtoDreamPartCollection *fPartColl; //!
  AliFemtoDreamDump *fProtonDeuteronDump; //!
  AliFemtoDreamDump *fAntiProtonAntiDeuteronDump; //!
  AliVTrack **fGTI; //!
  TList *fEvtList; //!
  TList *fProtonList; //!
  TList *fProtonMCList; //!
  TList *fAntiProtonList; //!
  TList *fAntiProtonMCList; //!
  TList *fDeuteronList; //!
  TList *fDeuteronMCList; //!
  TList *fAntiDeuteronList; //!
  TList *fAntiDeuteronMCList; //!
  TList *fDeuteronNoTOFList; //!
  TList *fDeuteronMCNoTOFList; //!
  TList *fAntiDeuteronNoTOFList; //!
  TList *fAntiDeuteronMCNoTOFList; //!
  TList *fResults; //!
  TList *fResultsQA; //!
  TList *fDumpster; //!
  TH2F  *fDeuteronRestMass; //!
  TH2F  *fAntiDeuteronRestMass; //!
  TH2F  *fDeuteronRestMassNoTOF; //!
  TH2F  *fAntiDeuteronRestMassNoTOF; //!
  ClassDef(AliAnalysisTaskNanoPt, 6)
};
#endif
