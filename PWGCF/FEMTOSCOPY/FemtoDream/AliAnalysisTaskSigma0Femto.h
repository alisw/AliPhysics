#ifndef AliAnalysisTaskSigma0Femto_H
#define AliAnalysisTaskSigma0Femto_H

#include "AliAnalysisTaskSE.h"
#include "AliConvEventCuts.h"
#include "AliEventCuts.h"
#include "AliFemtoDreamCollConfig.h"
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
  AliAnalysisTaskSigma0Femto(const char *name);
  virtual ~AliAnalysisTaskSigma0Femto() {}

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  void SetV0ReaderName(TString name) { fV0ReaderName = name; }
  void SetIsHeavyIon(bool isHeavyIon) { fIsHeavyIon = isHeavyIon; }
  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetLightweight(bool isLightweight) { fIsLightweight = isLightweight; }
  void SetV0Percentile(float v0perc) { fV0PercentileMax = v0perc; }
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }
  void SetProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCutsPartProton = cuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCutsPartAntiProton = cuts;
  }
  void SetV0Cuts(AliSigma0V0Cuts *cuts) { fV0Cuts = cuts; }
  void SetAntiV0Cuts(AliSigma0V0Cuts *cuts) { fAntiV0Cuts = cuts; }
  void SetPhotonV0Cuts(AliSigma0V0Cuts *cuts) { fPhotonV0Cuts = cuts; }
  void SetSigmaCuts(AliSigma0PhotonMotherCuts *cuts) { fSigmaCuts = cuts; }
  void SetAntiSigmaCuts(AliSigma0PhotonMotherCuts *cuts) {
    fAntiSigmaCuts = cuts;
  }
  void SetSigmaPhotonCuts(AliSigma0PhotonMotherCuts *cuts) {
    fSigmaPhotonCuts = cuts;
  }
  void SetAntiSigmaPhotonCuts(AliSigma0PhotonMotherCuts *cuts) {
    fAntiSigmaPhotonCuts = cuts;
  }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {
    fConfig = config;
  }

  bool AcceptEvent(AliVEvent *event);
  void CastToVector(std::vector<AliSigma0ParticleV0> &container,
                    const AliVEvent *inputEvent);
  void CastToVector(std::vector<AliSigma0ParticlePhotonMother> &sigmaContainer,
                    std::vector<AliFemtoDreamBasePart> &particles);
  void FillTriggerHisto(TH1F *histo);

  AliEventCuts fAliEventCuts;

 private:
  AliAnalysisTaskSigma0Femto(const AliAnalysisTaskSigma0Femto &task);
  AliAnalysisTaskSigma0Femto &operator=(const AliAnalysisTaskSigma0Femto &task);

  AliVEvent *fInputEvent;                       //! current event
  AliMCEvent *fMCEvent;                         //! corresponding MC event
  AliV0ReaderV1 *fV0Reader;                     //! basic photon Selection Task
  TString fV0ReaderName;                        //
  AliSigma0V0Cuts *fV0Cuts;                     //
  AliSigma0V0Cuts *fAntiV0Cuts;                 //
  AliSigma0V0Cuts *fPhotonV0Cuts;               //
  AliSigma0PhotonMotherCuts *fSigmaCuts;        //
  AliSigma0PhotonMotherCuts *fAntiSigmaCuts;    //
  AliSigma0PhotonMotherCuts *fSigmaPhotonCuts;  //
  AliSigma0PhotonMotherCuts *fAntiSigmaPhotonCuts;  //

  AliFemtoDreamTrack *fProtonTrack;                  //!
  AliFemtoDreamTrackCuts *fTrackCutsPartProton;      //
  AliFemtoDreamTrackCuts *fTrackCutsPartAntiProton;  //
  AliFemtoDreamCollConfig *fConfig;                  //
  AliFemtoDreamPairCleaner *fPairCleaner;            //!
  AliFemtoDreamPartCollection *fPartColl;            //!

  bool fIsMC;              //
  bool fIsHeavyIon;        //
  bool fIsLightweight;     //
  float fV0PercentileMax;  //
  UInt_t fTrigger;         //

  TClonesArray *fGammaArray;  //!

  // Histograms
  // =====================================================================

  TList *fOutputContainer;                  //!
  TList *fQA;                               //!
  TList *fHistoProton;                      //!
  TList *fHistoAntiProton;                  //!
  TH1F *fHistCutQA;                         //!
  TProfile *fHistRunNumber;                 //!
  TProfile *fHistCutBooking;                //!
  TH1F *fHistCentralityProfileBefore;       //!
  TH1F *fHistCentralityProfileAfter;        //!
  TH1F *fHistCentralityProfileCoarseAfter;  //!
  TH1F *fHistTriggerBefore;                 //!
  TH1F *fHistTriggerAfter;                  //!

  TList *fOutputTree;  //!

  ClassDef(AliAnalysisTaskSigma0Femto, 5)
};
#endif
