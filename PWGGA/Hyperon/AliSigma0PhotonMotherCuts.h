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
      const std::vector<AliAODConversionPhoton> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void SigmaToLambdaGamma(
      const std::vector<AliAODConversionPhoton> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void SigmaToLambdaGammaMixedEvent(
      const std::vector<AliAODConversionPhoton> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  void FillEventBuffer(
      const std::vector<AliAODConversionPhoton> &photonCandidates,
      const std::vector<AliSigma0ParticleV0> &lambdaCandidates);
  int GetRapidityBin(float rapidity) const;

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetLightweight(bool isLightweight) { fIsLightweight = isLightweight; }

  void SetSigmaMass(float mass) { fMassSigma = mass; }
  void SetMixingDepth(short mixDepth) { fMixingDepth = mixDepth; }
  void SetSigmaMassCut(float cut) { fSigmaMassCut = cut; }
  void SetSigmaSideband(float low, float up) {
    fSigmaSidebandLow = low;
    fSigmaSidebandUp = up;
  }
  void SetArmenterosCut(float qtLow, float qtUp, float alphaLow,
                        float alphaUp) {
    fArmenterosCut = true;
    fArmenterosQtLow = qtLow;
    fArmenterosQtUp = qtUp;
    fArmenterosAlphaLow = alphaLow;
    fArmenterosAlphaUp = alphaUp;
  }

  void InitCutHistograms(TString appendix = TString(""));
  TList *GetCutHistograms() const { return fHistograms; }

 protected:
  TList *fHistograms;
  TList *fHistogramsMC;

  bool fIsMC;
  bool fIsLightweight;

  AliVEvent *fInputEvent;  //!
  AliMCEvent *fMCEvent;    //!

  deque<vector<AliSigma0ParticleV0> > fLambdaMixed;     //!
  deque<vector<AliAODConversionPhoton> > fPhotonMixed;  //!

  short fMixingDepth;

  float fMassSigma;
  float fSigmaMassCut;
  float fSigmaSidebandLow;
  float fSigmaSidebandUp;

  float fArmenterosCut;
  float fArmenterosQtLow;
  float fArmenterosQtUp;
  float fArmenterosAlphaLow;
  float fArmenterosAlphaUp;

  // Histograms
  // =====================================================================
  TProfile *fHistCutBooking;  //

  TH1F *fHistNSigma;                   //
  TH1F *fHistPt;                       //
  TH1F *fHistMassCutPt;                //
  TH1F *fHistInvMass;                  //
  TH1F *fHistInvMassBeforeArmenteros;  //
  TH1F *fHistInvMassRec;               //
  TH2F *fHistInvMassPt;                //
  TH2F *fHistInvMassEta;               //
  TH2F *fHistEtaPhi;                   //
  TH1F *fHistRapidity;                 //
  TH2F *fHistPtY[22];                  //
  TH2F *fHistArmenterosBefore;         //
  TH2F *fHistArmenterosAfter;          //
  TH1F *fHistMixedPt;                  //
  TH1F *fHistMixedInvMass;             //
  TH2F *fHistMixedPtY[22];             //
  TH2F *fHistMixedInvMassPt;           //
  TH2F *fHistMixedInvMassEta;          //

  TH1F *fHistMCSigmaMassCutPt;                   //
  TH1F *fHistMCTruthSigma0PhotonConvPt;          //
  TH1F *fHistMCTruthSigma0PhotonConvP;           //
  TH1F *fHistMCTruthSigma0PhotonConvInvMass;     //
  TH2F *fHistMCTruthSigma0PhotonConvInvMassPt;   //
  TH2F *fHistMCTruthSigma0PhotonConvPtEta;       //
  TH2F *fHistMCTruthSigma0PhotonConvR;           //
  TH1F *fHistMCTruthSigma0PhotonConvConvPointX;  //
  TH1F *fHistMCTruthSigma0PhotonConvConvPointY;  //
  TH1F *fHistMCTruthSigma0PhotonConvConvPointZ;  //
  TH1F *fHistMCTruthSigma0PhotonConvEleP;        //
  TH1F *fHistMCTruthSigma0PhotonConvElePt;       //
  TH2F *fHistMCTruthSigma0PhotonConvPtY;         //
  TH1F *fHistMCTruthSigma0Pt;                    //
  TH2F *fHistMCTruthSigma0PtY;                   //
  TH2F *fHistMCTruthSigma0PtEta;                 //
  TH1F *fHistMCTruthSigma0PhotonPt;              //
  TH2F *fHistMCTruthSigma0PhotonPtY;             //
  TH2F *fHistMCTruthSigma0PhotonPtEta;           //

 private:
  ClassDef(AliSigma0PhotonMotherCuts, 2)
};

#endif
