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

  void SelectPhotonMother(
      AliVEvent *inputEvent, AliMCEvent *mcEvent,
      const std::vector<AliSigma0ParticleV0> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void SingleV0QA(const std::vector<AliSigma0ParticleV0> &photonCandidates,
                  const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void SigmaToLambdaGamma(
      const std::vector<AliSigma0ParticleV0> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void SigmaToLambdaGammaMixedEvent(
      const std::vector<AliSigma0ParticleV0> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void SigmaToLambdaGammaMixedEventBinned(
      const std::vector<AliSigma0ParticleV0> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void FillEventBuffer(
      const std::vector<AliSigma0ParticleV0> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  float ComputeRapidity(float pt, float pz, float m) const;
  int GetRapidityBin(float rapidity) const;
  int GetMultiplicityBin(float percentile) const;
  int GetZvertexBin(float zVertex) const;
  void ProcessMC() const;
  bool CheckDaughters(const AliMCParticle *particle) const;
  bool CheckDaughtersInAcceptance(const AliMCParticle *particle) const;

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetLightweight(bool isLightweight) { fIsLightweight = isLightweight; }

  void SetSigmaMass(float mass) { fMassSigma = mass; }
  void SetMixingDepth(short mixDepth) { fMixingDepth = mixDepth; }

  // For the containers for Femto
  void SetSigmaMassCut(float cut) { fSigmaMassCut = cut; }
  void SetSigmaSideband(float down, float up) {
    fSidebandCutDown = down, fSidebandCutUp = up;
  }

  void SetPhotonMinPt(float minpT) { fPhotonPtMin = minpT; }
  void SetPhotonMaxPt(float maxpT) { fPhotonPtMax = maxpT; }
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

  bool fIsMC;           //
  bool fIsLightweight;  //

  AliVEvent *fInputEvent;     //!
  AliMCEvent *fMCEvent;       //!
  TDatabasePDG fDataBasePDG;  //!

  std::vector<AliSigma0ParticlePhotonMother> fSigma;         //!
  std::vector<AliSigma0ParticlePhotonMother> fSidebandUp;    //!
  std::vector<AliSigma0ParticlePhotonMother> fSidebandDown;  //!

  deque<vector<AliSigma0ParticleV0> > fLambdaMixed;               //!
  deque<vector<AliSigma0ParticleV0> > fPhotonMixed;               //!
  deque<vector<AliSigma0ParticleV0> > fLambdaMixedBinned[10][6];  //!
  deque<vector<AliSigma0ParticleV0> > fPhotonMixedBinned[10][6];  //!
  float fTreeVariables[4];                                        //!

  AliSigma0V0Cuts *fLambdaCuts;  //
  AliSigma0V0Cuts *fPhotonCuts;  //
  AliV0ReaderV1 *fV0Reader;      //! basic photon Selection Task
  TString fV0ReaderName;         //

  short fMixingDepth;  //
  int fPDG;            //
  int fPDGDaughter1;   //
  int fPDGDaughter2;   //

  float fMassSigma;        //
  float fSigmaMassCut;     //
  float fSidebandCutUp;    //
  float fSidebandCutDown;  //
  float fPhotonPtMin;      //
  float fPhotonPtMax;      //
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

  TH1F *fHistNSigma;                   //!
  TH1F *fHistMassCutPt;                //!
  TH1F *fHistInvMass;                  //!
  TH1F *fHistInvMassBeforeArmenteros;  //!
  TH2F *fHistInvMassRecPhoton;         //!
  TH2F *fHistInvMassRecLambda;         //!
  TH2F *fHistInvMassRec;               //!
  TH2F *fHistInvMassPt;                //!
  TH2F *fHistInvMassEta;               //!
  TH2F *fHistEtaPhi;                   //!
  TH2F *fHistPtY[22];                  //!
  TH2F *fHistArmenterosBefore;         //!
  TH2F *fHistArmenterosAfter;          //!
  TH2F *fHistMixedPtY[22];             //!
  TH2F *fHistMixedInvMassPt;           //!
  TH2F *fHistMixedInvMassBinnedPt;     //!

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

  TH2F *fHistMCTruthPtY;                          //!
  TH2F *fHistMCTruthPtEta;                        //!
  TH2F *fHistMCTruthDaughterPtY;                  //!
  TH2F *fHistMCTruthDaughterPtEta;                //!
  TH2F *fHistMCTruthDaughterPtYAccept;            //!
  TH2F *fHistMCTruthDaughterPtEtaAccept;          //!
  TH2F *fHistMCTruthPtYHighMult;                  //!
  TH2F *fHistMCTruthPtEtaHighMult;                //!
  TH2F *fHistMCTruthDaughterPtYHighMult;          //!
  TH2F *fHistMCTruthDaughterPtEtaHighMult;        //!
  TH2F *fHistMCTruthDaughterPtYAcceptHighMult;    //!
  TH2F *fHistMCTruthDaughterPtEtaAcceptHighMult;  //!

  TH2F *fHistMCTrueSigmaLambdaPtCorr;  //!
  TH2F *fHistMCTrueSigmaPhotonPtCorr;  //!
  TH2F *fHistMCTrueSigmaLambdaPCorr;   //!
  TH2F *fHistMCTrueSigmaPhotonPCorr;   //!
  TH2F *fHistMCBkgSigmaLambdaPtCorr;   //!
  TH2F *fHistMCBkgSigmaPhotonPtCorr;   //!
  TH2F *fHistMCBkgSigmaLambdaPCorr;    //!
  TH2F *fHistMCBkgSigmaPhotonPCorr;    //!

  TH1F *fHistMCV0Pt;           //!
  TH1F *fHistMCV0Mass;         //!
  TH2F *fHistMCV0Mother;       //!
  TH2F *fHistMCV0Check;        //!
  TH2F *fHistMCV0MotherCheck;  //!

 private:
  ClassDef(AliSigma0PhotonMotherCuts, 13)
};

#endif
