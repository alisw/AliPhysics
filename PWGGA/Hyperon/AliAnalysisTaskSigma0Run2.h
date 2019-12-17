#ifndef AliAnalysisTaskSigma0Run2_H
#define AliAnalysisTaskSigma0Run2_H

#include "AliAnalysisTaskSE.h"
#include "AliConvEventCuts.h"
#include "AliEventCuts.h"
#include "AliMCEvent.h"
#include "AliSigma0PhotonMotherCuts.h"
#include "AliSigma0V0Cuts.h"
#include "AliStack.h"
#include "AliV0ReaderV1.h"
#include "TChain.h"

class AliVParticle;
class AliVTrack;

class AliAnalysisTaskSigma0Run2 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskSigma0Run2();
  AliAnalysisTaskSigma0Run2(const char *name);
  virtual ~AliAnalysisTaskSigma0Run2();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  void SetV0ReaderName(TString name) { fV0ReaderName = name; }
  void SetIsHeavyIon(bool isHeavyIon) { fIsHeavyIon = isHeavyIon; }
  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetLightweight(bool isLightweight) { fIsLightweight = isLightweight; }
  void SetIsRun1(bool isRun1) { fIsRun1 = isRun1; }
  void SetV0Percentile(float v0perc) { fV0PercentileMax = v0perc; }
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }
  void SetMultiplicityMode(UInt_t trigger) { fMultMode = trigger; }
  void SetV0Cuts(AliSigma0V0Cuts *cuts) { fV0Cuts = cuts; }
  void SetAntiV0Cuts(AliSigma0V0Cuts *cuts) { fAntiV0Cuts = cuts; }
  void SetSigmaCuts(AliSigma0PhotonMotherCuts *cuts) { fSigmaCuts = cuts; }
  void SetAntiSigmaCuts(AliSigma0PhotonMotherCuts *cuts) {
    fAntiSigmaCuts = cuts;
  }
  void SetPhotonLegPileUpCut(bool pileup) { fPhotonLegPileUpCut = pileup; }
  void SetPhotonDCArCut(bool dcar) { fDoPhotonDCArCut = dcar; }

  bool AcceptEvent(AliVEvent *event);
  bool AcceptEventRun1(AliVEvent *event);
  bool AcceptEventRun2(AliVEvent *event);
  void CastToVector(std::vector<AliSigma0ParticleV0> &container,
                    const AliVEvent *inputEvent);
  void FillTriggerHisto(TH1F *histo);

  AliEventCuts fAliEventCuts;

 private:
  AliAnalysisTaskSigma0Run2(const AliAnalysisTaskSigma0Run2 &task);
  AliAnalysisTaskSigma0Run2 &operator=(const AliAnalysisTaskSigma0Run2 &task);

  AliVEvent *fInputEvent;                     //! current event
  AliMCEvent *fMCEvent;                       //! corresponding MC event
  AliV0ReaderV1 *fV0Reader;                   //! basic photon Selection Task
  TString fV0ReaderName;                      //
  AliSigma0V0Cuts *fV0Cuts;                   //
  AliSigma0V0Cuts *fAntiV0Cuts;               //
  AliSigma0V0Cuts *fPhotonQA;                 //
  AliSigma0PhotonMotherCuts *fSigmaCuts;      //
  AliSigma0PhotonMotherCuts *fAntiSigmaCuts;  //

  bool fIsMC;                //
  bool fIsHeavyIon;          //
  bool fIsLightweight;       //
  bool fIsRun1;              //
  bool fPhotonLegPileUpCut;  //
  bool fDoPhotonDCArCut;    //
  float fV0PercentileMax;    //
  UInt_t fTrigger;           //
  UInt_t fMultMode;          //

  TClonesArray *fGammaArray;  //!

  // Histograms
  // =====================================================================

  TList *fOutputContainer;                  //!
  TList *fQA;                               //!
  TH1F *fHistCutQA;                         //!
  TProfile *fHistRunNumber;                 //!
  TProfile *fHistCutBooking;                //!
  TH1F *fHistCentralityProfileBefore;       //!
  TH1F *fHistCentralityProfileAfter;        //!
  TH1F *fHistCentralityProfileCoarseAfter;  //!
  TH1F *fHistMultiplicityRef08;             //!
  TH1F *fHistTriggerBefore;                 //!
  TH1F *fHistTriggerAfter;                  //!
  TH1I *fHistMultiplicity;                  //!

  ClassDef(AliAnalysisTaskSigma0Run2, 9)
};
#endif
