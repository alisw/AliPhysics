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
      std::vector<AliAODConversionPhoton> &photonCandidates,
      std::vector<AliSigma0ParticleV0> &lambdaCandidates,
      std::vector<AliSigma0ParticleV0> &antiLambdaCandidates);
  void SigmaToLambdaGamma(
      std::vector<AliAODConversionPhoton> &photonCandidates,
      std::vector<AliSigma0ParticleV0> &lambdaCandidates,
      std::vector<AliSigma0ParticleV0> &antiLambdaCandidates);
  void SigmaToLambdaGammaMixedEvent(
      std::vector<AliAODConversionPhoton> &photonCandidates,
      std::vector<AliSigma0ParticleV0> &lambdaCandidates,
      std::vector<AliSigma0ParticleV0> &antiLambdaCandidates);
  void FillEventBuffer(std::vector<AliAODConversionPhoton> &photonCandidates,
                       std::vector<AliSigma0ParticleV0> &lambdaCandidates,
                       std::vector<AliSigma0ParticleV0> &antiLambdaCandidates);
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
  TList *fHistogramsSigma;
  TList *fHistogramsSigmaMC;
  TList *fHistogramsAntiSigma;
  TList *fHistogramsAntiSigmaMC;

  bool fIsMC;
  bool fIsLightweight;

  AliVEvent *fInputEvent;  //!
  AliMCEvent *fMCEvent;    //!

  deque<vector<AliSigma0ParticleV0> > fLambdaMixed;      //!
  deque<vector<AliSigma0ParticleV0> > fAntiLambdaMixed;  //!
  deque<vector<AliAODConversionPhoton> > fPhotonMixed;   //!

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
  TProfile *fHistCuts;  //

  TH1F *fHistNSigma;                 //
  TH1F *fHistSigmaPt;                //
  TH1F *fHistSigmaMassCutPt;         //
  TH1F *fHistSigmaInvMass;           //
  TH1F *fHistSigmaInvMassRec;        //
  TH2F *fHistSigmaInvMassPt;         //
  TH2F *fHistSigmaInvMassEta;        //
  TH2F *fHistSigmaEtaPhi;            //
  TH1F *fHistSigmaRapidity;          //
  TH2F *fHistSigmaPtY[22];           //
  TH2F *fHistSigmaArmenterosBefore;  //
  TH2F *fHistSigmaArmenterosAfter;   //
  TH1F *fHistSigmaMixedPt;           //
  TH1F *fHistSigmaMixedInvMass;      //
  TH2F *fHistSigmaMixedPtY[22];      //
  TH2F *fHistSigmaMixedInvMassPt;    //
  TH2F *fHistSigmaMixedInvMassEta;   //

  TH1F *fHistNAntiSigma;                 //
  TH1F *fHistAntiSigmaPt;                //
  TH1F *fHistAntiSigmaMassCutPt;         //
  TH1F *fHistAntiSigmaInvMass;           //
  TH1F *fHistAntiSigmaInvMassRec;        //
  TH2F *fHistAntiSigmaInvMassPt;         //
  TH2F *fHistAntiSigmaInvMassEta;        //
  TH2F *fHistAntiSigmaEtaPhi;            //
  TH1F *fHistAntiSigmaRapidity;          //
  TH2F *fHistAntiSigmaPtY[22];           //
  TH2F *fHistAntiSigmaArmenterosBefore;  //
  TH2F *fHistAntiSigmaArmenterosAfter;   //
  TH1F *fHistAntiSigmaMixedPt;           //
  TH1F *fHistAntiSigmaMixedInvMass;      //
  TH2F *fHistAntiSigmaMixedPtY[22];      //
  TH2F *fHistAntiSigmaMixedInvMassPt;    //
  TH2F *fHistAntiSigmaMixedInvMassEta;   //

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

  TH1F *fHistMCAntiSigmaMassCutPt;                   //
  TH1F *fHistMCTruthAntiSigma0PhotonConvPt;          //
  TH1F *fHistMCTruthAntiSigma0PhotonConvP;           //
  TH1F *fHistMCTruthAntiSigma0PhotonConvInvMass;     //
  TH2F *fHistMCTruthAntiSigma0PhotonConvInvMassPt;   //
  TH2F *fHistMCTruthAntiSigma0PhotonConvPtEta;       //
  TH2F *fHistMCTruthAntiSigma0PhotonConvR;           //
  TH1F *fHistMCTruthAntiSigma0PhotonConvConvPointX;  //
  TH1F *fHistMCTruthAntiSigma0PhotonConvConvPointY;  //
  TH1F *fHistMCTruthAntiSigma0PhotonConvConvPointZ;  //
  TH1F *fHistMCTruthAntiSigma0PhotonConvEleP;        //
  TH1F *fHistMCTruthAntiSigma0PhotonConvElePt;       //
  TH2F *fHistMCTruthAntiSigma0PhotonConvPtY;         //
  TH1F *fHistMCTruthAntiSigma0Pt;                    //
  TH2F *fHistMCTruthAntiSigma0PtY;                   //
  TH2F *fHistMCTruthAntiSigma0PtEta;                 //
  TH1F *fHistMCTruthAntiSigma0PhotonPt;              //
  TH2F *fHistMCTruthAntiSigma0PhotonPtY;             //
  TH2F *fHistMCTruthAntiSigma0PhotonPtEta;           //

 private:
  ClassDef(AliSigma0PhotonMotherCuts, 2)
};

#endif
