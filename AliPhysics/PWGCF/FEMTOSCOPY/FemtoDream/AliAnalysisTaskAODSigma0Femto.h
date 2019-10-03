#ifndef AliAnalysisTaskAODSigma0Femto_H
#define AliAnalysisTaskAODSigma0Femto_H

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

class AliVParticle;
class AliVTrack;

class AliAnalysisTaskAODSigma0Femto : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskAODSigma0Femto();
  AliAnalysisTaskAODSigma0Femto(const char *name, const bool isMC);
  virtual ~AliAnalysisTaskAODSigma0Femto();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  void CastToVector(std::vector<AliFemtoDreamBasePart> &particlesOut,
                    std::vector<AliFemtoDreamBasePart> &particlesIn);

  void SetV0ReaderName(TString name) { fV0ReaderName = name; }
  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetLightweight(bool isLightweight) { fIsLightweight = isLightweight; }
  void SetCheckDaughterCF(bool doIt) { fCheckDaughterCF = doIt; }
  void SetGoDoThisFemtoJanitor(bool godothis) { fFemtoJanitor = godothis; }
  void SetV0Percentile(float v0perc) { fV0PercentileMax = v0perc; }
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }
  void SetEventCuts(AliFemtoDreamEventCuts *cuts) { fEvtCuts = cuts; }
  void SetProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCutsPartProton = cuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCutsPartAntiProton = cuts;
  }
  void SetV0Cuts(AliFemtoDreamv0Cuts *cuts) { fV0Cuts = cuts; }
  void SetAntiV0Cuts(AliFemtoDreamv0Cuts *cuts) { fAntiV0Cuts = cuts; }
  void SetPhotonCuts(AliSigma0PhotonCuts *cuts) {
    fPhotonCuts = cuts;
  }
  void SetSigmaCuts(AliSigma0AODPhotonMotherCuts *cuts) { fSigmaCuts = cuts; }
  void SetAntiSigmaCuts(AliSigma0AODPhotonMotherCuts *cuts) {
    fAntiSigmaCuts = cuts;
  }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {
    fConfig = config;
  }

 private:
  AliAnalysisTaskAODSigma0Femto(
      const AliAnalysisTaskAODSigma0Femto &task);
  AliAnalysisTaskAODSigma0Femto &operator=(
      const AliAnalysisTaskAODSigma0Femto &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);

  AliVEvent *fInputEvent;                        //! current event
  AliMCEvent *fMCEvent;                          //! corresponding MC event
  AliV0ReaderV1 *fV0Reader;                      //! basic photon Selection Task
  TString fV0ReaderName;                         //

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
  AliFemtoDreamControlSample *fSample;               //!

  bool fIsMC;                //
  bool fIsLightweight;       //
  bool fCheckDaughterCF;     //
  bool fFemtoJanitor;        //
  float fV0PercentileMax;    //
  UInt_t fTrigger;           //

  TClonesArray *fGammaArray;  //!
  int fTrackBufferSize;
  AliAODTrack **fGTI;  //!

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
  TList *fResultsSample;           //!
  TList *fResultsSampleQA;         //!

  ClassDef(AliAnalysisTaskAODSigma0Femto, 3)
};
#endif
