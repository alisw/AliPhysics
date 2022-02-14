#ifndef AliAnalysisTaskNanoAODSigma0Femto_H
#define AliAnalysisTaskNanoAODSigma0Femto_H

#include "AliAnalysisTaskSE.h"
#include "AliConvEventCuts.h"
#include "AliEventCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliSigma0PhotonCuts.h"
#include "AliFemtoDreamControlSample.h"
#include "AliMCEvent.h"
#include "AliSigma0AODPhotonMotherCuts.h"
#include "AliStack.h"
#include "AliV0ReaderV1.h"
#include "TChain.h"
#include "AliFemtoDreamBaseDump.h"

class AliVParticle;
class AliVTrack;

class AliAnalysisTaskNanoAODSigma0Femto : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskNanoAODSigma0Femto();
  AliAnalysisTaskNanoAODSigma0Femto(const char *name, const bool isMC);
  virtual ~AliAnalysisTaskNanoAODSigma0Femto();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  void CastToVector(std::vector<AliFemtoDreamBasePart> &particlesOut,
                    std::vector<AliFemtoDreamBasePart> &particlesIn);

  void SetIsMC(bool isMC) {
    fIsMC = isMC;
  }
  void SetLightweight(bool isLightweight) {
    fIsLightweight = isLightweight;
  }
  void SetUseDumpster(bool use) {
    fUseDumpster = use;
  }
  void SetV0Percentile(float v0perc) {
    fV0PercentileMax = v0perc;
  }
  void SetTrigger(UInt_t trigger) {
    fTrigger = trigger;
  }
  void SetEventCuts(AliFemtoDreamEventCuts *cuts) {
    fEvtCuts = cuts;
  }
  void SetProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCutsPartProton = cuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCutsPartAntiProton = cuts;
  }
  void SetV0Cuts(AliFemtoDreamv0Cuts *cuts) {
    fV0Cuts = cuts;
  }
  void SetAntiV0Cuts(AliFemtoDreamv0Cuts *cuts) {
    fAntiV0Cuts = cuts;
  }
  void SetPhotonCuts(AliSigma0PhotonCuts *cuts) {
    fPhotonCuts = cuts;
  }
  void SetSigmaCuts(AliSigma0AODPhotonMotherCuts *cuts) {
    fSigmaCuts = cuts;
  }
  void SetAntiSigmaCuts(AliSigma0AODPhotonMotherCuts *cuts) {
    fAntiSigmaCuts = cuts;
  }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {
    fConfig = config;
  }

 private:
  AliAnalysisTaskNanoAODSigma0Femto(
      const AliAnalysisTaskNanoAODSigma0Femto &task);
  AliAnalysisTaskNanoAODSigma0Femto &operator=(
      const AliAnalysisTaskNanoAODSigma0Femto &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);

  AliVEvent *fInputEvent;                        //! current event
  AliMCEvent *fMCEvent;                          //! corresponding MC event
  AliSigma0AODPhotonMotherCuts *fSigmaCuts;      //
  AliSigma0AODPhotonMotherCuts *fAntiSigmaCuts;  //

  TRandom3 *fRandom;  //!

  AliFemtoDreamEvent *fEvent;                        //!
  AliFemtoDreamEventCuts *fEvtCuts;                  //
  AliFemtoDreamTrack *fProtonTrack;                  //!
  AliFemtoDreamTrackCuts *fTrackCutsPartProton;      //
  AliFemtoDreamTrackCuts *fTrackCutsPartAntiProton;  //
  AliFemtoDreamv0 *fLambda;                          //!
  AliFemtoDreamv0Cuts *fV0Cuts;                      //
  AliFemtoDreamv0Cuts *fAntiV0Cuts;                  //
  AliSigma0PhotonCuts *fPhotonCuts;                  //
  AliFemtoDreamCollConfig *fConfig;                  //
  AliFemtoDreamPairCleaner *fPairCleaner;            //!
  AliFemtoDreamPartCollection *fPartColl;            //!

  AliFemtoDreamDump *fProtonSigmaDump;               //!
  AliFemtoDreamDump *fAntiProtonAntiSigmaDump;       //!
  AliFemtoDreamDump *fProtonSBDump;                  //!
  AliFemtoDreamDump *fAntiProtonAntiSBDump;          //!

  bool fIsMC;              //
  bool fIsLightweight;     //
  bool fUseDumpster;       //
  float fV0PercentileMax;  //
  UInt_t fTrigger;         //

  TClonesArray *fGammaArray;  //!
  int fTrackBufferSize;
  AliVTrack **fGTI;  //!

  TList *fQA;                      //!
  TList *fEvtHistList;             //!
  TList *fTrackCutHistList;        //!
  TList *fTrackCutHistMCList;      //!
  TList *fAntiTrackCutHistList;    //!
  TList *fAntiTrackCutHistMCList;  //!
  TList *fLambdaHistList;          //!
  TList *fLambdaHistMCList;        //!
  TList *fAntiLambdaHistList;      //!
  TList *fAntiLambdaHistMCList;    //!
  TList *fPhotonHistList;          //!
  TList *fSigmaHistList;           //!
  TList *fAntiSigmaHistList;       //!
  TList *fResultList;              //!
  TList *fResultQAList;            //!
  TList *fDumpster;                //!

ClassDef(AliAnalysisTaskNanoAODSigma0Femto, 8)
};

#endif
