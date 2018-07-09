#ifndef AliSigma0PhotonMotherCuts_H
#define AliSigma0PhotonMotherCuts_H

#include <deque>
#include "AliAODConversionPhoton.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliSigma0ParticleBase.h"
#include "AliSigma0ParticlePhotonMother.h"
#include "AliSigma0ParticleV0.h"
#include "AliVEvent.h"
#include "Riostream.h"
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
  void SigmaToLambdaGamma(
      const std::vector<AliSigma0ParticleV0> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void SigmaToLambdaGammaMixedEvent(
      const std::vector<AliSigma0ParticleV0> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void FillEventBuffer(
      const std::vector<AliSigma0ParticleV0> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  float ComputeRapidity(float pt, float pz, float m) const;
  int GetRapidityBin(float rapidity) const;
  void ProcessMC() const;
  bool CheckDaughters(const AliMCParticle *particle) const;

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetLightweight(bool isLightweight) { fIsLightweight = isLightweight; }
  void SetTreeOutput(bool isTreeOutput) { fIsTreeOutput = isTreeOutput; }

  void SetSigmaMass(float mass) { fMassSigma = mass; }
  void SetMixingDepth(short mixDepth) { fMixingDepth = mixDepth; }
  void SetSigmaMassCut(float cut) { fSigmaMassCut = cut; }
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

  void InitCutHistograms(TString appendix = TString(""));
  TList *GetCutHistograms() const { return fHistograms; }
  TTree *GetSigmaTree() const { return fOutputTree; }

 protected:
  TList *fHistograms;    //!
  TList *fHistogramsMC;  //!

  bool fIsMC;
  bool fIsLightweight;
  bool fIsTreeOutput;

  AliVEvent *fInputEvent;  //!
  AliMCEvent *fMCEvent;    //!

  deque<vector<AliSigma0ParticleV0> > fLambdaMixed;  //!
  deque<vector<AliSigma0ParticleV0> > fPhotonMixed;  //!
  float fTreeVariables[4];                           //!

  short fMixingDepth;  //
  int fPDG;            //
  int fPDGDaughter1;   //
  int fPDGDaughter2;   //

  float fMassSigma;         //
  float fSigmaMassCut;      //
  float fPhotonPtMin;  //
  float fPhotonPtMax;   //

  float fArmenterosCut;       //
  float fArmenterosQtLow;     //
  float fArmenterosQtUp;      //
  float fArmenterosAlphaLow;  //
  float fArmenterosAlphaUp;   //

  // Histograms
  // =====================================================================
  TProfile *fHistCutBooking;  //!

  TH1F *fHistNSigma;                   //!
  TH1F *fHistPt;                       //!
  TH1F *fHistMassCutPt;                //!
  TH1F *fHistInvMass;                  //!
  TH1F *fHistInvMassBeforeArmenteros;  //!
  TH1F *fHistInvMassRec;               //!
  TH2F *fHistInvMassPt;                //!
  TH2F *fHistInvMassEta;               //!
  TH2F *fHistEtaPhi;                   //!
  TH1F *fHistRapidity;                 //!
  TH2F *fHistPtY[22];                  //!
  TH2F *fHistArmenterosBefore;         //!
  TH2F *fHistArmenterosAfter;          //!
  TH1F *fHistMixedPt;                  //!
  TH1F *fHistMixedInvMass;             //!
  TH2F *fHistMixedPtY[22];             //!
  TH2F *fHistMixedInvMassPt;           //!
  TH2F *fHistMixedInvMassEta;          //!

  TH1F *fHistMCTruthPt;             //!
  TH2F *fHistMCTruthPtY;            //!
  TH2F *fHistMCTruthPtEta;          //!
  TH1F *fHistMCTruthDaughterPt;     //!
  TH2F *fHistMCTruthDaughterPtY;    //!
  TH2F *fHistMCTruthDaughterPtEta;  //!

  TH1F *fHistMCV0Pt;    //!
  TH1F *fHistMCV0Mass;  //!

  TTree *fOutputTree;  //!

 private:
  ClassDef(AliSigma0PhotonMotherCuts, 5)
};

#endif
