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
  void SetProtonCutsNoTOF(AliFemtoDreamTrackCuts *cuts) {
    fProtonTrackNoTOF = cuts;
  }
  void SetAntiProtonCutsNoTOF(AliFemtoDreamTrackCuts *cuts) {
    fAntiProtonTrackNoTOF = cuts;
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
  void SetUseDumpsterRestPairs(bool use) {
    fUseDumpsterRestPairs = use;
  }

 private:
  AliAnalysisTaskNanoPt(
    const AliAnalysisTaskNanoPt &task);
  AliAnalysisTaskNanoPt &operator=(
    const AliAnalysisTaskNanoPt &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  bool fIsMC;                                         //
  bool fUseDumpster;  //
  bool fUseDumpsterRestPairs;  //
  AliVEvent *fInputEvent;                            //! current event
  AliFemtoDreamEvent *fEvent;                        //!
  AliFemtoDreamEventCuts *fEvtCuts;                  //
  AliFemtoDreamTrack *fTrack;                        //
  AliFemtoDreamTrackCuts *fProtonTrack;              //
  AliFemtoDreamTrackCuts *fAntiProtonTrack;          //
  AliFemtoDreamTrackCuts *fDeuteronTrack;              //
  AliFemtoDreamTrackCuts *fAntiDeuteronTrack;          //
  AliFemtoDreamTrackCuts *fProtonTrackNoTOF;              //
  AliFemtoDreamTrackCuts *fAntiProtonTrackNoTOF;          //
  AliFemtoDreamTrackCuts *fDeuteronTrackNoTOF;              //
  AliFemtoDreamTrackCuts *fAntiDeuteronTrackNoTOF;          //

  int fTrackBufferSize; //
  AliVTrack **fGTI;  //!

  TList *fEvtList;//!
  TList *fProtonList;//!
  TList *fProtonMCList;//!
  TList *fAntiProtonList;//!
  TList *fAntiProtonMCList;//!
  TList *fDeuteronList;//!
  TList *fDeuteronMCList;//!
  TList *fAntiDeuteronList;//!
  TList *fAntiDeuteronMCList;//!
  TList *fProtonNoTOFList;//!
  TList *fProtonMCNoTOFList;//!
  TList *fAntiProtonNoTOFList;//!
  TList *fAntiProtonMCNoTOFList;//!
  TList *fDeuteronNoTOFList;//!
  TList *fDeuteronMCNoTOFList;//!
  TList *fAntiDeuteronNoTOFList;//!
  TList *fAntiDeuteronMCNoTOFList;//!

  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!

  TList *fResults;                          //!
  TList *fResultsQA;                        //!

  TH2F  *fProtonRestMass;                   //!
  TH2F  *fAntiProtonRestMass;               //!
  TH2F  *fDeuteronRestMass;                 //!
  TH2F  *fAntiDeuteronRestMass;             //!
  TH2F  *fProtonRestMassNoTOF;              //!
  TH2F  *fAntiProtonRestMassNoTOF;          //!
  TH2F  *fDeuteronRestMassNoTOF;            //!
  TH2F  *fAntiDeuteronRestMassNoTOF;        //!
  TH2F  *fProtonRestMassMC;                 //!
  TH2F  *fAntiProtonRestMassMC;             //!
  TH2F  *fDeuteronRestMassMC;               //!
  TH2F  *fAntiDeuteronRestMassMC;           //!
  TH2F  *fKaonRestMassMC;                   //!
  TH2F  *fAntiKaonRestMassMC;               //!
  TH2F  *fDProtonRestMassMC;                //!
  TH2F  *fDKaonRestMassMC;                  //!
  TH2F  *fAntiDProtonRestMassMC;            //!
  TH2F  *fAntiDKaonRestMassMC;              //!
  TH2F  *fPionRestMassMC;                   //!
  TH2F  *fAntiPionRestMassMC;               //!
  TH2F  *fDPionRestMassMC;                  //!
  TH2F  *fAntiDPionRestMassMC;              //!
  TH2F  *fProtonBackgroungMC;               //!
  TH2F  *fAntiProtonBackgroundMC;           //!
  TH2F  *fDeuteronBackgroundMC;             //!
  TH2F  *fAntiDeuteronBackgroundMC;         //!

  AliFemtoDreamDump *fProtonProtonDump; //!
  AliFemtoDreamDump *fProtonAntiProtonDump; //!
  AliFemtoDreamDump *fProtonDeuteronDump;   //!
  AliFemtoDreamDump *fProtonAntiDeuteronDump; //!
  AliFemtoDreamDump *fAntiProtonAntiProtonDump; //!
  AliFemtoDreamDump *fAntiProtonDeuteronDump; //!
  AliFemtoDreamDump *fAntiProtonAntiDeuteronDump; //!
  AliFemtoDreamDump *fDeuteronAntiDeuteronDump; //!
  TList* fDumpster; //!
  ClassDef(AliAnalysisTaskNanoPt, 5)
};
#endif
