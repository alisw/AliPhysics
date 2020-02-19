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
#include "AliFemtoDreamBaseDump.h"

class AliAnalysisTaskFemtoDreamDeuteron : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFemtoDreamDeuteron();
  AliAnalysisTaskFemtoDreamDeuteron(const char *name, bool isMC);
  virtual ~AliAnalysisTaskFemtoDreamDeuteron();
  Float_t GetMass2sq(AliFemtoDreamTrack *track);
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
  void SetTrackCutsProtonMass(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsProtonMass = trkCuts;};
  void SetTrackCutsAntiProtonDCA(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiProtonDCA = trkCuts;};
  void SetTrackCutsAntiProtonMass(AliFemtoDreamTrackCuts *trkCuts) {fTrackCutsAntiProtonMass = trkCuts;};
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {fConfig = config;};
 private:
  AliAnalysisTaskFemtoDreamDeuteron(
    const AliAnalysisTaskFemtoDreamDeuteron &task);
  AliAnalysisTaskFemtoDreamDeuteron &operator=(
    const AliAnalysisTaskFemtoDreamDeuteron &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  bool fIsMC;                               //
  AliFemtoDreamEvent *fEvent;               //!  on Runtime
  AliFemtoDreamTrack *fTrack;               //!
  AliFemtoDreamEventCuts *fEventCuts;       //   Stream these bad boys
  AliFemtoDreamTrackCuts *fTrackCutsDeuteronDCA;  //
  AliFemtoDreamTrackCuts *fTrackCutsDeuteronMass;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiDeuteronDCA;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiDeuteronMass;  //
  AliFemtoDreamTrackCuts *fTrackCutsProtonDCA;  //
  AliFemtoDreamTrackCuts *fTrackCutsProtonMass;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiProtonDCA;  //
  AliFemtoDreamTrackCuts *fTrackCutsAntiProtonMass;  //
  TList *fEvtList;//!
  TList *fProtonList;//!
  TList* fProtonMCList;//!
  TList *fAntiProtonList;//!
  TList* fAntiProtonMCList;//!
  TList *fDeuteronList;//!
  TList* fDeuteronMCList;//!
  TList *fAntiDeuteronList;//!
  TList* fAntiDeuteronMCList;//!
  TList *fProtonNoTOFList;//!
  TList* fProtonMCNoTOFList;//!
  TList *fAntiProtonNoTOFList;//!
  TList* fAntiProtonMCNoTOFList;//!
  TList *fDeuteronNoTOFList;//!
  TList* fDeuteronMCNoTOFList;//!
  TList *fAntiDeuteronNoTOFList;//!
  TList* fAntiDeuteronMCNoTOFList;//!
  AliFemtoDreamCollConfig *fConfig;         //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  TList *fResults;                          //!
  TList *fResultsQA;                        //!
  AliAODTrack** fGTI;           //!
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
  TH2F  *fProtonBackgroundMC;               //!
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
  int fTrackBufferSize;                     //
  ClassDef(AliAnalysisTaskFemtoDreamDeuteron, 4)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKFEMTOTUTORIAL_H_ */
