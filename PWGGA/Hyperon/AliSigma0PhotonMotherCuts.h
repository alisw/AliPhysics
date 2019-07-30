#ifndef AliSigma0PhotonMotherCuts_H
#define AliSigma0PhotonMotherCuts_H

#include <deque>
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliSigma0ParticleBase.h"
#include "AliSigma0ParticlePhotonMother.h"
#include "AliSigma0ParticleV0.h"
#include "AliSigma0V0Cuts.h"
#include "AliV0ReaderV1.h"
#include "AliVEvent.h"
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TObject.h"
#include "TProfile.h"

class AliPIDResponse;

class AliSigma0PhotonMotherCuts : public TObject {
 public:
  AliSigma0PhotonMotherCuts();
  AliSigma0PhotonMotherCuts(const AliSigma0PhotonMotherCuts &);
  AliSigma0PhotonMotherCuts &operator=(const AliSigma0PhotonMotherCuts &);
  virtual ~AliSigma0PhotonMotherCuts() {}

  static AliSigma0PhotonMotherCuts *DefaultCuts();

  void SelectPhotonMother(AliVEvent *inputEvent, AliMCEvent *mcEvent,
                          std::vector<AliSigma0ParticleV0> &photonCandidates,
                          std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void CleanUpClones(std::vector<AliSigma0ParticleV0> &photonCandidates,
                     std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void SingleV0QA(const std::vector<AliSigma0ParticleV0> &photonCandidates,
                  const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void SigmaToLambdaGamma(
      const std::vector<AliSigma0ParticleV0> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  float GetMassSigmaPt(float pt) const;
  void SigmaToLambdaGammaMixedEvent(
      const std::vector<AliSigma0ParticleV0> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void SigmaToLambdaGammaMixedEventBinned(
      const std::vector<AliSigma0ParticleV0> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void FillEventBuffer(
      const std::vector<AliSigma0ParticleV0> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  static int GetMultiplicityBin(float percentile, UInt_t trigger);
  void ProcessMC() const;
  bool CheckDaughters(const AliMCParticle *particle) const;
  bool CheckDaughtersInAcceptance(const AliMCParticle *particle) const;

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetDoCleanUp(bool doCleanUp) { fDoCleanUp = doCleanUp; }
  void SetLightweight(bool isLightweight) { fIsLightweight = isLightweight; }
  void SetIsSpectrum(bool isSpectrum) { fIsSpectrumAnalysis = isSpectrum; }

  void SetSigmaMass(float mass) { fMassSigma = mass; }
  void SetMixingDepth(short mixDepth) { fMixingDepth = mixDepth; }

  // For the containers for Femto
  void SetSigmaMassCut(float cut) { fSigmaMassCut = cut; }
  void SetSigmaSideband(float down, float up) {
    fSidebandCutDown = down, fSidebandCutUp = up;
  }
  void SetMultiplicityMode(UInt_t trigger) { fMultMode = trigger; }

  void SetPhotonMinPt(float minpT) { fPhotonPtMin = minpT; }
  void SetPhotonMaxPt(float maxpT) { fPhotonPtMax = maxpT; }
  void SetMinPt(float minpT) { fPtMin = minpT; }
  void SetArmenterosCut(float qtLow, float qtUp, float alphaLow,
                        float alphaUp) {
    fArmenterosCut = true;
    fArmenterosQtLow = qtLow;
    fArmenterosQtUp = qtUp;
    fArmenterosAlphaLow = alphaLow;
    fArmenterosAlphaUp = alphaUp;
  }
  void SetPDG(const int pdgPID, const int pdgDaughter1,
              const int pdgDaughter2) {
    fPDG = pdgPID;
    fPDGDaughter1 = pdgDaughter1;
    fPDGDaughter2 = pdgDaughter2;
  }
  void SetMCMultThreshold(float multThr) { fMCHighMultThreshold = multThr; }

  void SetLambdaCuts(AliSigma0V0Cuts *lamCut) { fLambdaCuts = lamCut; }
  void SetPhotonCuts(AliSigma0V0Cuts *photCut) { fPhotonCuts = photCut; }
  void SetV0ReaderName(TString name) { fV0ReaderName = name; }
  void SetSigmaMassPt(bool doIt) { fMassWindowPt = doIt; }
  void SetSigmaMassParameters(const float p0, const float p1, const float p2) {
    fMassWindowPt = true;
    fMassWindowP0 = p0;
    fMassWindowP1 = p1;
    fMassWindowP2 = p2;
  }

  void InitCutHistograms(TString appendix = TString(""));
  TList *GetCutHistograms() const { return fHistograms; }

  std::vector<AliSigma0ParticlePhotonMother> &GetSigma() { return fSigma; }
  std::vector<AliSigma0ParticlePhotonMother> &GetSidebandUp() {
    return fSidebandUp;
  }
  std::vector<AliSigma0ParticlePhotonMother> &GetSidebandDown() {
    return fSidebandDown;
  }
  void GetLambda(std::vector<AliSigma0ParticleV0> &vecIn) {
    vecIn.clear();
    for (const auto &sigma : fSigma) {
      vecIn.emplace_back(sigma.GetV0());
    }
  }
  void GetPhoton(std::vector<AliSigma0ParticleV0> &vecIn) {
    vecIn.clear();
    for (const auto &sigma : fSigma) {
      vecIn.emplace_back(sigma.GetPhoton());
    }
  }

 protected:
  TList *fHistograms;    //!
  TList *fHistogramsMC;  //!

  bool fIsMC;                //
  bool fDoCleanUp;           //
  bool fIsLightweight;       //
  bool fIsSpectrumAnalysis;  //

  AliVEvent *fInputEvent;     //!
  AliMCEvent *fMCEvent;       //!
  TDatabasePDG fDataBasePDG;  //!

  std::vector<AliSigma0ParticlePhotonMother> fSigma;         //!
  std::vector<AliSigma0ParticlePhotonMother> fSidebandUp;    //!
  std::vector<AliSigma0ParticlePhotonMother> fSidebandDown;  //!

  deque<vector<AliSigma0ParticleV0> > fLambdaMixed;           //!
  deque<vector<AliSigma0ParticleV0> > fPhotonMixed;           //!
  deque<vector<AliSigma0ParticleV0> > fLambdaMixedBinned[5];  //!
  deque<vector<AliSigma0ParticleV0> > fPhotonMixedBinned[5];  //!

  AliSigma0V0Cuts *fLambdaCuts;  //
  AliSigma0V0Cuts *fPhotonCuts;  //
  AliV0ReaderV1 *fV0Reader;      //! basic photon Selection Task
  TString fV0ReaderName;         //
  UInt_t fMultMode;              //

  short fMixingDepth;  //
  int fPDG;            //
  int fPDGDaughter1;   //
  int fPDGDaughter2;   //

  bool fMassWindowPt;   //
  float fMassWindowP0;  //
  float fMassWindowP1;  //
  float fMassWindowP2;  //

  float fMassSigma;        //
  float fSigmaMassCut;     //
  float fSidebandCutUp;    //
  float fSidebandCutDown;  //
  float fPhotonPtMin;      //
  float fPhotonPtMax;      //
  float fPtMin;            //
  float fRapidityMax;      //

  float fArmenterosCut;       //
  float fArmenterosQtLow;     //
  float fArmenterosQtUp;      //
  float fArmenterosAlphaLow;  //
  float fArmenterosAlphaUp;   //

  float fMCHighMultThreshold;  //

  // Histograms
  // =====================================================================
  TProfile *fHistCutBooking;  //!

  TH1F *fHistNSigma;                                //!
  TH1F *fHistNPhotonBefore;                         //!
  TH1F *fHistNPhotonAfter;                          //!
  TH1F *fHistNLambdaBefore;                         //!
  TH1F *fHistNLambdaAfter;                          //!
  TH1F *fHistNPhotonLabel;                          //!
  TH1F *fHistNLambdaLabel;                          //!
  TH1F *fHistNLambdaGammaLabel;                     //!
  TH1F *fHistMassCutPt;                             //!
  TH1F *fHistInvMass;                               //!
  TH1F *fHistInvMassK0Gamma;                        //!
  TH2F *fHistInvMassSelected;                       //!
  TH2F *fHistInvMassRecPhoton;                      //!
  TH2F *fHistInvMassRecLambda;                      //!
  TH2F *fHistInvMassRec;                            //!
  TH2F *fHistInvMassPt;                             //!
  TH2F *fHistInvMassPtRaw;                          //!
  TH2F *fHistEtaPhi;                                //!
  TH2F *fHistPtRapidity;                            //!
  TH2F *fHistPtMult[5];                             //!
  TH2F *fHistArmenterosBefore;                      //!
  TH2F *fHistArmenterosAfter;                       //!
  TH2F *fHistMixedInvMassPt;                        //!
  TH2F *fHistMixedInvMassBinnedMultPt[5];           //!
  TH2F *fHistDeltaEtaDeltaPhiGammaNegBefore;        //!
  TH2F *fHistDeltaEtaDeltaPhiGammaPosBefore;        //!
  TH2F *fHistDeltaEtaDeltaPhiLambdaNegBefore;       //!
  TH2F *fHistDeltaEtaDeltaPhiLambdaPosBefore;       //!
  TH2F *fHistDeltaEtaDeltaPhiLambdaGammaNegBefore;  //!
  TH2F *fHistDeltaEtaDeltaPhiLambdaGammaPosBefore;  //!
  TH2F *fHistDeltaEtaDeltaPhiGammaNegAfter;         //!
  TH2F *fHistDeltaEtaDeltaPhiGammaPosAfter;         //!
  TH2F *fHistDeltaEtaDeltaPhiLambdaNegAfter;        //!
  TH2F *fHistDeltaEtaDeltaPhiLambdaPosAfter;        //!
  TH2F *fHistDeltaEtaDeltaPhiLambdaGammaNegAfter;   //!
  TH2F *fHistDeltaEtaDeltaPhiLambdaGammaPosAfter;   //!

  TH2F *fHistLambdaPtPhi;   //!
  TH2F *fHistLambdaPtEta;   //!
  TH2F *fHistLambdaMassPt;  //!
  TH2F *fHistPhotonPtPhi;   //!
  TH2F *fHistPhotonPtEta;   //!
  TH2F *fHistPhotonMassPt;  //!

  TH2F *fHistSigmaLambdaPtCorr;  //!
  TH2F *fHistSigmaPhotonPtCorr;  //!
  TH2F *fHistSigmaLambdaPCorr;   //!
  TH2F *fHistSigmaPhotonPCorr;   //!

  TH1F *fHistMCTruthPt;                         //!
  TH1F *fHistMCTruthPtMult[5];                  //!
  TH2F *fHistMCTruthPtY;                        //!
  TH2F *fHistMCTruthDaughterPtY;                //!
  TH2F *fHistMCTruthDaughterPtYAccept;          //!
  TH2F *fHistMCTruthPtYHighMult;                //!
  TH2F *fHistMCTruthDaughterPtYHighMult;        //!
  TH2F *fHistMCTruthDaughterPtYAcceptHighMult;  //!

  TH1F *fHistMCV0Pt;           //!
  TH1F *fHistMCV0Mass;         //!
  TH2F *fHistMCV0Mother;       //!
  TH2F *fHistMCV0Check;        //!
  TH2F *fHistMCV0MotherCheck;  //!

 private:
  ClassDef(AliSigma0PhotonMotherCuts, 29)
};

#endif
