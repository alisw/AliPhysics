#ifndef AliAnalysisTaskSigma0Femto_H
#define AliAnalysisTaskSigma0Femto_H

#include "AliAnalysisTaskSE.h"
#include "AliConvEventCuts.h"
#include "AliEventCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliMCEvent.h"
#include "AliSigma0PhotonMotherCuts.h"
#include "AliSigma0V0Cuts.h"
#include "AliStack.h"
#include "AliV0ReaderV1.h"
#include "TChain.h"

class AliVParticle;
class AliVTrack;

class AliAnalysisTaskSigma0Femto : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSigma0Femto();
  AliAnalysisTaskSigma0Femto(const char *name, const bool isMC);
  virtual ~AliAnalysisTaskSigma0Femto();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  void SetV0ReaderName(TString name) { fV0ReaderName = name; }
  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetLightweight(bool isLightweight) { fIsLightweight = isLightweight; }
  void SetCheckDaughterCF(bool doIt) { fCheckDaughterCF = doIt; }
  void SetV0Percentile(float v0perc) { fV0PercentileMax = v0perc; }
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }
  void SetMultiplicityMode(UInt_t trigger) { fMultMode = trigger; }
  void SetEventCuts(AliFemtoDreamEventCuts *cuts) { fEvtCuts = cuts; }
  void SetProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCutsPartProton = cuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCutsPartAntiProton = cuts;
  }
  void SetV0Cuts(AliSigma0V0Cuts *cuts) { fV0Cuts = cuts; }
  void SetAntiV0Cuts(AliSigma0V0Cuts *cuts) { fAntiV0Cuts = cuts; }
  void SetSigmaCuts(AliSigma0PhotonMotherCuts *cuts) { fSigmaCuts = cuts; }
  void SetAntiSigmaCuts(AliSigma0PhotonMotherCuts *cuts) {
    fAntiSigmaCuts = cuts;
  }
  void SetPhotonLegPileUpCut(bool pileup) { fPhotonLegPileUpCut = pileup; }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {
    fConfig = config;
  }

  void CastToVector(std::vector<AliSigma0ParticleV0> &container,
                    std::vector<AliFemtoDreamBasePart> &particles,
                    const AliMCEvent *mcEvent);
  void CastToVector(std::vector<AliSigma0ParticleV0> &container,
                    const AliVEvent *inputEvent);
  void CastToVector(std::vector<AliSigma0ParticlePhotonMother> &sigmaContainer,
                    std::vector<AliFemtoDreamBasePart> &particles,
                    const AliMCEvent *mcEvent);
  void FillTriggerHisto(TH1F *histo);

 private:
  AliAnalysisTaskSigma0Femto(const AliAnalysisTaskSigma0Femto &task);
  AliAnalysisTaskSigma0Femto &operator=(const AliAnalysisTaskSigma0Femto &task);

  AliVEvent *fInputEvent;                     //! current event
  AliMCEvent *fMCEvent;                       //! corresponding MC event
  AliV0ReaderV1 *fV0Reader;                   //! basic photon Selection Task
  TString fV0ReaderName;                      //
  AliSigma0V0Cuts *fV0Cuts;                   //
  AliSigma0V0Cuts *fAntiV0Cuts;               //
  AliSigma0V0Cuts *fPhotonQA;                 //
  AliSigma0PhotonMotherCuts *fSigmaCuts;      //
  AliSigma0PhotonMotherCuts *fAntiSigmaCuts;  //

  TRandom3 *fRandom;  //!

  AliFemtoDreamEvent *fEvent;                        //!
  AliFemtoDreamEventCuts *fEvtCuts;                  //
  AliFemtoDreamTrack *fProtonTrack;                  //!
  AliFemtoDreamTrackCuts *fTrackCutsPartProton;      //
  AliFemtoDreamTrackCuts *fTrackCutsPartAntiProton;  //
  AliFemtoDreamCollConfig *fConfig;                  //
  AliFemtoDreamPairCleaner *fPairCleaner;            //!
  AliFemtoDreamPartCollection *fPartColl;            //!

  bool fIsMC;                //
  bool fIsLightweight;       //
  bool fCheckDaughterCF;     //
  bool fPhotonLegPileUpCut;  //
  float fV0PercentileMax;    //
  UInt_t fTrigger;           //
  UInt_t fMultMode;          //

  TClonesArray *fGammaArray;  //!

  TList *fQA;                      //!
  TList *fEvtHistList;             //!
  TList *fTrackCutHistList;        //!
  TList *fTrackCutHistMCList;      //!
  TList *fAntiTrackCutHistList;    //!
  TList *fAntiTrackCutHistMCList;  //!
  TList *fLambdaHistList;          //!
  TList *fAntiLambdaHistList;      //!
  TList *fPhotonHistList;          //!
  TList *fSigmaHistList;           //!
  TList *fAntiSigmaHistList;       //!
  TList *fResultList;              //!
  TList *fResultQAList;            //!

  ClassDef(AliAnalysisTaskSigma0Femto, 16)
};
#endif
