#ifndef AliSigma0AODPhotonMotherCuts_H
#define AliSigma0AODPhotonMotherCuts_H

#include <deque>
#include "AliAnalysisManager.h"
#include "AliFemtoDreamBasePart.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliVEvent.h"
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TObject.h"
#include "TProfile.h"

class AliPIDResponse;

class AliSigma0AODPhotonMotherCuts : public TObject {
 public:
  AliSigma0AODPhotonMotherCuts();
  AliSigma0AODPhotonMotherCuts(const AliSigma0AODPhotonMotherCuts &);
  AliSigma0AODPhotonMotherCuts &operator=(const AliSigma0AODPhotonMotherCuts &);
  virtual ~AliSigma0AODPhotonMotherCuts() {}

  static AliSigma0AODPhotonMotherCuts *DefaultCuts();

  void SelectPhotonMother(AliVEvent *inputEvent, AliMCEvent *mcEvent,
                          const std::vector<AliFemtoDreamBasePart> &photonCandidates,
                          const std::vector<AliFemtoDreamBasePart> &lambdaCandidates);
  void CleanUpClones();
  void SingleV0QA();
  void SigmaToLambdaGamma();
  float GetMassSigmaPt(float pt) const;
  float GetArmenterosAlpha(const AliFemtoDreamBasePart &gamma,
                           const AliFemtoDreamBasePart &lambda,
                           const AliFemtoDreamBasePart &sigma) const;
  float GetArmenterosQt(const AliFemtoDreamBasePart &gamma,
                        const AliFemtoDreamBasePart &lambda,
                        const AliFemtoDreamBasePart &sigma) const;

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetDoCleanUp(bool doCleanUp) { fDoCleanUp = doCleanUp; }
  void SetLightweight(bool isLightweight) { fIsLightweight = isLightweight; }
  void SetSigmaMass(float mass) { fMassSigma = mass; }
  // For the containers for Femto
  void SetSigmaMassCut(float cut) { fSigmaMassCut = cut; }
  void SetSigmaSideband(float down, float up) {
    fSidebandCutDown = down, fSidebandCutUp = up;
  }
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

  void SetSigmaMassPt(bool doIt) { fMassWindowPt = doIt; }
  void SetSigmaMassParameters(const float p0, const float p1, const float p2) {
    fMassWindowPt = true;
    fMassWindowP0 = p0;
    fMassWindowP1 = p1;
    fMassWindowP2 = p2;
  }

  void SetDeltaPhiEtaMax(float maxVal) {
    fDeltaPhiEtaMax = maxVal;
    fDoDeltaPhiEtaCut = true;
  }

  void InitCutHistograms(TString appendix = TString(""));
  TList *GetCutHistograms() const { return fHistograms; }

  std::vector<AliFemtoDreamBasePart> &GetSigma() { return fSigma; }
  std::vector<AliFemtoDreamBasePart> &GetSidebandUp() { return fSidebandUp; }
  std::vector<AliFemtoDreamBasePart> &GetSidebandDown() {
    return fSidebandDown;
  }
  std::vector<AliFemtoDreamBasePart> &GetLambda() { return fLambda; }
  std::vector<AliFemtoDreamBasePart> &GetPhoton() { return fPhoton; }

 protected:
  TList *fHistograms;    //!
  TList *fHistogramsMC;  //!

  bool fIsMC;           //
  bool fDoCleanUp;      //
  bool fIsLightweight;  //

  AliVEvent *fInputEvent;     //!
  AliMCEvent *fMCEvent;       //!
  TDatabasePDG fDataBasePDG;  //!

  std::vector<AliFemtoDreamBasePart> fSigma;         //!
  std::vector<AliFemtoDreamBasePart> fLambda;         //!
  std::vector<AliFemtoDreamBasePart> fPhoton;         //!
  std::vector<AliFemtoDreamBasePart> fSidebandUp;    //!
  std::vector<AliFemtoDreamBasePart> fSidebandDown;  //!

  std::vector<AliFemtoDreamBasePart> fLambdaCandidates;         //!
  std::vector<AliFemtoDreamBasePart> fPhotonCandidates;         //!

  int fPDG;           //
  int fPDGDaughter1;  //
  int fPDGDaughter2;  //

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
  bool fDoDeltaPhiEtaCut;     //
  float fDeltaPhiEtaMax;      //

  // Histograms
  // =====================================================================
  TProfile *fHistCutBooking;  //!

  TH1F *fHistNSigma;                                //!
  TH1F *fHistNCandidates;                           //!
  TH1F *fHistNPhotonBefore;                         //!
  TH1F *fHistNPhotonAfter;                          //!
  TH1F *fHistNLambdaBefore;                         //!
  TH1F *fHistNLambdaAfter;                          //!
  TH1F *fHistNPhotonLabel;                          //!
  TH1F *fHistNLambdaLabel;                          //!
  TH1F *fHistNLambdaGammaLabel;                     //!
  TH1F *fHistNPhotonSplit;                          //!
  TH1F *fHistNLambdaSplit;                          //!
  TH1F *fHistNLambdaGammaSplit;                     //!
  TH1F *fHistMassCutPt;                             //!
  TH1F *fHistInvMass;                               //!
  TH1F* fHistInvMassOtherChilren[16];               //!
  TH2F *fHistInvMassSelected;                       //!
  TH2F *fHistInvMassRecPhoton;                      //!
  TH2F *fHistInvMassRecLambda;                      //!
  TH2F *fHistInvMassRec;                            //!
  TH2F *fHistInvMassPtRaw;                          //!
  TH2F *fHistEtaPhi;                                //!
  TH2F *fHistArmenterosBefore;                      //!
  TH2F *fHistArmenterosAfter;                       //!
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

  TH1F *fHistMCV0Pt;           //!
  TH1F *fHistMCV0Mass;         //!
  TH2F *fHistMCV0Mother;       //!
  TH2F *fHistMCV0Check;        //!
  TH2F *fHistMCV0MotherCheck;  //!

 private:
  ClassDef(AliSigma0AODPhotonMotherCuts, 6)
};

#endif
