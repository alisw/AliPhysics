#ifndef AliSigma0PhotonMotherCuts_H
#define AliSigma0PhotonMotherCuts_H

#include "AliAODConversionPhoton.h"
#include "AliMCParticle.h"
#include "AliAODTrack.h"
#include "AliMCEvent.h"
#include "AliSigma0ParticleBase.h"
#include "AliSigma0ParticlePhotonMother.h"
#include "AliSigma0ParticleV0.h"
#include "AliVEvent.h"
#include "Riostream.h"
#include "TObject.h"

#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TProfile.h"

#include <deque>

class AliPIDResponse;

class AliSigma0PhotonMotherCuts : public TObject {
 public:
  AliSigma0PhotonMotherCuts();
  AliSigma0PhotonMotherCuts(const AliSigma0PhotonMotherCuts &);
  AliSigma0PhotonMotherCuts &operator=(const AliSigma0PhotonMotherCuts &);
  virtual ~AliSigma0PhotonMotherCuts() {}

  static AliSigma0PhotonMotherCuts *DefaultCuts();

  void SelectPhotonMother(AliVEvent *inputEvent, AliMCEvent *mcEvent,
      std::vector<AliAODConversionPhoton> &photonCandidates,
      std::vector<AliSigma0ParticleV0> &lambdaCandidates,
      std::vector<AliSigma0ParticleV0> &antiLambdaCandidates,
      std::vector<AliSigma0ParticlePhotonMother> &sigmaCandidates,
      std::vector<AliSigma0ParticlePhotonMother> &antiSigmaCandidates);
  void SigmaToLambdaGamma(
      std::vector<AliAODConversionPhoton> &photonCandidates,
      std::vector<AliSigma0ParticleV0> &lambdaCandidates,
      std::vector<AliSigma0ParticleV0> &antiLambdaCandidates,
      std::vector<AliSigma0ParticlePhotonMother> &sigmaCandidates,
      std::vector<AliSigma0ParticlePhotonMother> &antiSigmaCandidates);
  void SigmaToLambdaGammaGamma(
      std::vector<AliAODConversionPhoton> &photonCandidates,
      std::vector<AliSigma0ParticleV0> &lambdaCandidates,
      std::vector<AliSigma0ParticleV0> &antiLambdaCandidates);
  void TrackCleaner(
      std::vector<AliSigma0ParticlePhotonMother> &sigmaCandidates,
      std::vector<AliSigma0ParticlePhotonMother> &antiSigmaCandidates,
      float massLow, float massUp);
  void FillEventBuffer(std::vector<AliAODConversionPhoton> &photonCandidates,
                       std::vector<AliSigma0ParticleV0> &lambdaCandidates,
                       std::vector<AliSigma0ParticleV0> &antiLambdaCandidates);
  void ProcessMC();
  bool IsConvertedPhoton(TParticle *particle) const;
  float ComputeRapidity(float pt, float pz, float m) const;
  int GetRapidityBin(float rapidity) const;

  void SetIsMC(bool isMC) { fIsMC = isMC; }
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

  void InitCutHistograms();
  TList *GetCutHistograms() const { return fHistograms; }

 protected:
  TList *fHistograms;
  TList *fHistogramsSigma;
  TList *fHistogramsSigmaMC;
  TList *fHistogramsAntiSigma;
  TList *fHistogramsAntiSigmaMC;

  bool fIsMC;

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
  TH1F *fHistSigmaInvMassFinal;      //
  TH1F *fHistSigmaInvMassRec;        //
  TH2F *fHistSigmaInvMassPt;         //
  TH2F *fHistSigmaInvMassEta;        //
  TH2F *fHistSigmaEtaPhi;            //
  TH1F *fHistSigmaRapidity;          //
  TH2F *fHistSigmaPtY[22];           //
  TH1F *fHistSigmaNShared;           //
  TH1F *fHistSigmaMassDiff;          //
  TH2F *fHistSigmaSharedInvMass;     //
  TH2F *fHistSigmaArmenterosBefore;  //
  TH2F *fHistSigmaArmenterosAfter;   //
  TH1F *fHistSigmaMixedPt;           //
  TH1F *fHistSigmaMixedInvMass;      //
  TH2F *fHistSigmaMixedPtY[22];      //
  TH2F *fHistSigmaMixedInvMassPt;    //
  TH2F *fHistSigmaMixedInvMassEta;   //

  TH1F *fHistSigmaPtTwoGamma;                //
  TH1F *fHistSigmaInvMassTwoGamma;           //
  TH2F *fHistSigmaInvMassPtTwoGamma;         //
  TH2F *fHistSigmaInvMassEtaTwoGamma;        //
  TH2F *fHistSigmaArmenterosBeforeTwoGamma;  //
  TH2F *fHistSigmaArmenterosAfterTwoGamma;   //
  TH1F *fHistSigmaMixedPtTwoGamma;           //
  TH1F *fHistSigmaMixedInvMassTwoGamma;      //
  TH2F *fHistSigmaMixedInvMassPtTwoGamma;    //
  TH2F *fHistSigmaMixedInvMassEtaTwoGamma;   //

  TH1F *fHistSigmaPtDiElectron;                //
  TH1F *fHistSigmaInvMassDiElectron;           //
  TH2F *fHistSigmaInvMassPtDiElectron;         //
  TH2F *fHistSigmaInvMassEtaDiElectron;        //
  TH2F *fHistSigmaArmenterosBeforeDiElectron;  //
  TH2F *fHistSigmaArmenterosAfterDiElectron;   //
  TH1F *fHistSigmaMixedPtDiElectron;           //
  TH1F *fHistSigmaMixedInvMassDiElectron;      //
  TH2F *fHistSigmaMixedInvMassPtDiElectron;    //
  TH2F *fHistSigmaMixedInvMassEtaDiElectron;   //

  TH1F *fHistNAntiSigma;                 //
  TH1F *fHistAntiSigmaPt;                //
  TH1F *fHistAntiSigmaMassCutPt;         //
  TH1F *fHistAntiSigmaInvMass;           //
  TH1F *fHistAntiSigmaInvMassFinal;      //
  TH1F *fHistAntiSigmaInvMassRec;        //
  TH2F *fHistAntiSigmaInvMassPt;         //
  TH2F *fHistAntiSigmaInvMassEta;        //
  TH2F *fHistAntiSigmaEtaPhi;            //
  TH1F *fHistAntiSigmaRapidity;          //
  TH2F *fHistAntiSigmaPtY[22];           //
  TH1F *fHistAntiSigmaNShared;           //
  TH1F *fHistAntiSigmaMassDiff;          //
  TH2F *fHistAntiSigmaSharedInvMass;     //
  TH2F *fHistAntiSigmaArmenterosBefore;  //
  TH2F *fHistAntiSigmaArmenterosAfter;   //
  TH1F *fHistAntiSigmaMixedPt;           //
  TH1F *fHistAntiSigmaMixedInvMass;      //
  TH2F *fHistAntiSigmaMixedPtY[22];      //
  TH2F *fHistAntiSigmaMixedInvMassPt;    //
  TH2F *fHistAntiSigmaMixedInvMassEta;   //

  TH1F *fHistAntiSigmaPtTwoGamma;                //
  TH1F *fHistAntiSigmaInvMassTwoGamma;           //
  TH2F *fHistAntiSigmaInvMassPtTwoGamma;         //
  TH2F *fHistAntiSigmaInvMassEtaTwoGamma;        //
  TH2F *fHistAntiSigmaArmenterosBeforeTwoGamma;  //
  TH2F *fHistAntiSigmaArmenterosAfterTwoGamma;   //
  TH1F *fHistAntiSigmaMixedPtTwoGamma;           //
  TH1F *fHistAntiSigmaMixedInvMassTwoGamma;      //
  TH2F *fHistAntiSigmaMixedInvMassPtTwoGamma;    //
  TH2F *fHistAntiSigmaMixedInvMassEtaTwoGamma;   //

  TH1F *fHistAntiSigmaPtDiElectron;                //
  TH1F *fHistAntiSigmaInvMassDiElectron;           //
  TH2F *fHistAntiSigmaInvMassPtDiElectron;         //
  TH2F *fHistAntiSigmaInvMassEtaDiElectron;        //
  TH2F *fHistAntiSigmaArmenterosBeforeDiElectron;  //
  TH2F *fHistAntiSigmaArmenterosAfterDiElectron;   //
  TH1F *fHistAntiSigmaMixedPtDiElectron;           //
  TH1F *fHistAntiSigmaMixedInvMassDiElectron;      //
  TH2F *fHistAntiSigmaMixedInvMassPtDiElectron;    //
  TH2F *fHistAntiSigmaMixedInvMassEtaDiElectron;   //

  TH1F *fHistMCSigmaMassCutPt;                  //
  TH1F *fHistMCTruthSigma0PhotonConvPt;             //
  TH1F *fHistMCTruthSigma0PhotonConvP;              //
  TH1F *fHistMCTruthSigma0PhotonConvInvMass;        //
  TH2F *fHistMCTruthSigma0PhotonConvInvMassPt;      //
  TH2F *fHistMCTruthSigma0PhotonConvPtEta;          //
  TH2F *fHistMCTruthSigma0PhotonConvR;              //
  TH1F *fHistMCTruthSigma0PhotonConvConvPointX;     //
  TH1F *fHistMCTruthSigma0PhotonConvConvPointY;     //
  TH1F *fHistMCTruthSigma0PhotonConvConvPointZ;     //
  TH1F *fHistMCTruthSigma0PhotonConvEleP;           //
  TH1F *fHistMCTruthSigma0PhotonConvElePt;          //
  TH2F *fHistMCTruthSigma0PhotonConvPtY;            //
  TH1F *fHistMCTruthSigma0Pt;                   //
  TH2F *fHistMCTruthSigma0PtY;                  //
  TH2F *fHistMCTruthSigma0PtEta;                //
  TH1F *fHistMCTruthSigma0PhotonPt;                   //
  TH2F *fHistMCTruthSigma0PhotonPtY;                  //
  TH2F *fHistMCTruthSigma0PhotonPtEta;                //
  TH1F *fHistMCTruthSigma0TwoElectronPtPos;     //
  TH1F *fHistMCTruthSigma0TwoElectronPPos;      //
  TH2F *fHistMCTruthSigma0TwoElectronPtEtaPos;  //
  TH1F *fHistMCTruthSigma0TwoElectronPtNeg;     //
  TH1F *fHistMCTruthSigma0TwoElectronPNeg;      //
  TH2F *fHistMCTruthSigma0TwoElectronPtEtaNeg;  //
  TH1F *fHistMCTruthSigma0TwoGammaPt;
  TH1F *fHistMCTruthSigma0TwoGammaP;
  TH2F *fHistMCTruthSigma0TwoGammaPtEta;
  TH1F *fHistMCTruthSigma0TwoGammaPtEle;
  TH1F *fHistMCTruthSigma0TwoGammaPEle;
  TH2F *fHistMCTruthSigma0TwoGammaPtEtaEle;

  TH1F *fHistMCAntiSigmaMassCutPt;                  //
  TH1F *fHistMCTruthAntiSigma0PhotonConvPt;             //
  TH1F *fHistMCTruthAntiSigma0PhotonConvP;              //
  TH1F *fHistMCTruthAntiSigma0PhotonConvInvMass;        //
  TH2F *fHistMCTruthAntiSigma0PhotonConvInvMassPt;      //
  TH2F *fHistMCTruthAntiSigma0PhotonConvPtEta;          //
  TH2F *fHistMCTruthAntiSigma0PhotonConvR;              //
  TH1F *fHistMCTruthAntiSigma0PhotonConvConvPointX;     //
  TH1F *fHistMCTruthAntiSigma0PhotonConvConvPointY;     //
  TH1F *fHistMCTruthAntiSigma0PhotonConvConvPointZ;     //
  TH1F *fHistMCTruthAntiSigma0PhotonConvEleP;           //
  TH1F *fHistMCTruthAntiSigma0PhotonConvElePt;          //
  TH2F *fHistMCTruthAntiSigma0PhotonConvPtY;            //
  TH1F *fHistMCTruthAntiSigma0Pt;                   //
  TH2F *fHistMCTruthAntiSigma0PtY;                  //
  TH2F *fHistMCTruthAntiSigma0PtEta;                //
  TH1F *fHistMCTruthAntiSigma0PhotonPt;                   //
  TH2F *fHistMCTruthAntiSigma0PhotonPtY;                  //
  TH2F *fHistMCTruthAntiSigma0PhotonPtEta;                //
  TH1F *fHistMCTruthAntiSigma0TwoElectronPtPos;     //
  TH1F *fHistMCTruthAntiSigma0TwoElectronPPos;      //
  TH2F *fHistMCTruthAntiSigma0TwoElectronPtEtaPos;  //
  TH1F *fHistMCTruthAntiSigma0TwoElectronPtNeg;     //
  TH1F *fHistMCTruthAntiSigma0TwoElectronPNeg;      //
  TH2F *fHistMCTruthAntiSigma0TwoElectronPtEtaNeg;  //
  TH1F *fHistMCTruthAntiSigma0TwoGammaPt;
  TH1F *fHistMCTruthAntiSigma0TwoGammaP;
  TH2F *fHistMCTruthAntiSigma0TwoGammaPtEta;
  TH1F *fHistMCTruthAntiSigma0TwoGammaPtEle;
  TH1F *fHistMCTruthAntiSigma0TwoGammaPEle;
  TH2F *fHistMCTruthAntiSigma0TwoGammaPtEtaEle;

 private:
  ClassDef(AliSigma0PhotonMotherCuts, 1)
};

#endif
