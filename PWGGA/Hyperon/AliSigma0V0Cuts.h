#ifndef AliSigma0V0Cuts_H
#define AliSigma0V0Cuts_H

#include "AliMCEvent.h"
#include "AliPIDResponse.h"
#include "AliSigma0ParticleV0.h"
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TObject.h"

#include "TProfile.h"

class AliPIDResponse;

class AliSigma0V0Cuts : public TObject {
 public:
  enum PileUpRejectionMode {
    BothDaughtersCombined =
        0,  ///< Impose timing criteria on both daughters, ITS or TOF
    OneDaughterCombined =
        1,     ///< Impose timing criteria on both daughters, ITS or TOF
    None = 2,  ///< No timing information required
    BothDaughtersITSonly =
        3,  ///< Impose timing criteria on both daughters, ITS only
    BothDaughtersTOFonly =
        4,  ///< Impose timing criteria on both daughters, TOF only
    OneDaughterITSonly =
        5,  ///< Impose timing criteria on both daughters, ITS only
    OneDaughterTOFonly =
        6,  ///< Impose timing criteria on both daughters, ITS only
  };

  AliSigma0V0Cuts();
  AliSigma0V0Cuts(const AliSigma0V0Cuts &);
  AliSigma0V0Cuts &operator=(const AliSigma0V0Cuts &);
  virtual ~AliSigma0V0Cuts();

  static AliSigma0V0Cuts *LambdaCuts();
  static AliSigma0V0Cuts *PhotonCuts();

  void SelectV0(AliVEvent *inputEvent, AliMCEvent *mcEvent);
  bool V0QualityCuts(const AliESDv0 *v0) const;
  bool V0PID(const AliESDv0 *v0, const AliESDtrack *pos,
             const AliESDtrack *neg) const;
  bool SingleParticlePID(const AliVTrack *track, AliPID::EParticleType particle,
                         float &prob) const;
  void PlotSingleParticlePID(const AliVTrack *track,
                             AliPID::EParticleType particle) const;
  bool SingleParticleQualityCuts(AliESDtrack *track) const;
  bool PileUpRejection(AliESDtrack *pos, AliESDtrack *neg) const;
  bool V0TopologicalSelection(const AliESDv0 *v0) const;
  void PlotMasses(AliESDv0 *v0) const;
  bool LambdaSelection(AliESDv0 *v0) const;
  bool PhotonSelection(AliESDv0 *v0) const;
  float ComputePhotonMass(const AliESDv0 *v0) const;
  float ComputePhotonMassRefit(const AliESDv0 *v0) const;
  float ComputePsiPair(const AliESDv0 *v0) const;
  void PhotonQA(AliVEvent *inputEvent, AliMCEvent *mcEvent, const TClonesArray *photons);

  void SetLightweight(bool isLightweight) { fIsLightweight = isLightweight; }
  void SetCheckCutsMC(bool checkCuts) { fCheckCutsMC = checkCuts; }

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetPileUpRejectionMode(PileUpRejectionMode pileUpRej) {
    fPileUpRejectionMode = pileUpRej;
  }
  void SetPID(int pid) { fPID = pid; }
  void SetPosPID(AliPID::EParticleType pid, const int pdg) {
    fPosPID = pid;
    fPosPDG = pdg;
  }
  void SetNegPID(AliPID::EParticleType pid, const int pdg) {
    fNegPID = pid;
    fNegPDG = pdg;
  }
  void SetV0OnFlyStatus(bool onFly) { fV0OnFly = onFly; }
  void SetV0PtMin(float ptMin) { fV0PtMin = ptMin; }
  void SetV0PtMax(float ptMax) { fV0PtMax = ptMax; }
  void SetV0CosPAMin(float cosPAMin) { fV0CosPAMin = cosPAMin; }
  void SetV0RadiusMax(float rMax) { fV0RadiusMax = rMax; }
  void SetV0RadiusMin(float rMin) { fV0RadiusMin = rMin; }
  void SetV0DecayVertexMax(float dvertexMax) { fV0DecayVertexMax = dvertexMax; }
  void SetPIDnSigma(float nSigma) {
    fPIDnSigma = nSigma;
    fUsePID = true;
  }
  void SetTPCclusterMin(short nTPCcluster) { fTPCclusterMin = nTPCcluster; }
  void SetTPCcrossedRowsMin(short nCrossed) { fTPCnCrossedRowsMin = nCrossed; }
  void SetTPCRatioFindable(float ratio) { fTPCratioFindable = ratio; }
  void SetTPCnClsFindable(short nClsFind) { fTPCfindableMin = nClsFind; }
  void SetTPCnSharedMax(short nSharedMax) { fTPCnSharedMax = nSharedMax; }
  void SetEtaMax(float etaMax) { fEtaMax = etaMax; }
  void SetChi2Max(float chi2Max) { fChi2Max = chi2Max; }
  void SetDaughterDCAMax(float ddcaMax) { fDaughterDCAMax = ddcaMax; }
  void SetDaughterDCAtoPV(float dca2pv) { fDaughterDCAPV = dca2pv; }
  void SetK0Rejection(float low, float up) {
    fK0Rejection = true;
    fK0RejectionLow = low, fK0RejectionUp = up;
  }
  void SetArmenterosCut(float qtLow, float qtUp, float alphaLow,
                        float alphaUp) {
    fUseArmenteros = true;
    fArmenterosQtLow = qtLow;
    fArmenterosQtUp = qtUp;
    fArmenterosAlphaLow = alphaLow;
    fArmenterosAlphaUp = alphaUp;
  }
  void SetLambdaSelection(float low, float up) {
    fLambdaSelectionLow = low, fLambdaSelectionUp = up;
  }
  void SetPsiPairMax(float max) { fPsiPairMax = max; }

  void SetMCMultThreshold(float multThr) { fMCHighMultThreshold = multThr; }

  void ProcessMC() const;
  void CheckCutsMC() const;
  bool CheckDaughters(const AliMCParticle *particle) const;
  bool CheckDaughtersInAcceptance(const AliMCParticle *particle) const;

  void InitCutHistograms(TString appendix = TString(""));
  TList *GetCutHistograms() const { return fHistograms; }

  std::vector<AliSigma0ParticleV0> &GetV0s() { return fV0Container; }

 protected:
  TList *fHistograms;        //!
  TList *fHistogramsMC;      //!
  TList *fHistogramsBefore;  //!
  TList *fHistogramsAfter;   //!
  TList *fHistogramsPos;     //!
  TList *fHistogramsNeg;     //!

  AliESDEvent *fInputEvent;   //!
  AliMCEvent *fMCEvent;       //!
  TDatabasePDG fDataBasePDG;  //!

  std::vector<AliSigma0ParticleV0> fV0Container;  //!

  bool fIsLightweight;  //
  bool fCheckCutsMC;    //

  short fV0cut;  //
  int fPID;      //
  int fPosPDG;   //
  int fNegPDG;   //

  bool fIsMC;                                //
  PileUpRejectionMode fPileUpRejectionMode;  //
  AliPID::EParticleType fPosPID;             //
  AliPID::EParticleType fNegPID;             //
  bool fV0OnFly;                             //
  bool fK0Rejection;                         //
  bool fUsePID;                              //
  bool fUseArmenteros;                       //
  float fV0PtMin;                            //
  float fV0PtMax;                            //
  float fV0CosPAMin;                         //
  float fV0RadiusMax;                        //
  float fV0RadiusMin;                        //
  float fV0DecayVertexMax;                   //
  float fPIDnSigma;                          //
  float fEtaMax;                             //
  float fChi2Max;                            //
  float fTPCclusterMin;                      //
  short fTPCnCrossedRowsMin;                 //
  float fTPCratioFindable;                   //
  short fTPCfindableMin;                     //
  short fTPCnSharedMax;                      //
  float fDaughterDCAMax;                     //
  float fDaughterDCAPV;                      //
  float fK0RejectionLow;                     //
  float fK0RejectionUp;                      //
  float fArmenterosQtLow;                    //
  float fArmenterosQtUp;                     //
  float fArmenterosAlphaLow;                 //
  float fArmenterosAlphaUp;                  //
  float fLambdaSelectionLow;                 //
  float fLambdaSelectionUp;                  //
  float fPsiPairMax;                         //

  float fMCHighMultThreshold;  //

  AliPIDResponse *fPIDResponse;  //!  pid response

  // Histograms
  // =====================================================================
  TProfile *fHistCutBooking;  //!

  TH1F *fHistCuts;  //!
  TH1F *fHistNV0;   //!

  TH1F *fHistLambdaMass;       //!
  TH1F *fHistAntiLambdaMass;   //!
  TH1F *fHistPhotonMass;       //!
  TH1F *fHistPhotonMassRefit;  //!
  TH1F *fHistK0Mass;           //!
  TH1F *fHistV0Pt;             //!
  TH1F *fHistV0Mass;           //!
  TH2F *fHistV0MassPt;         //!
  TH1F *fHistK0MassAfter;      //!
  TH2F *fHistCosPA;            //!
  TH2F *fHistEtaPhi;           //!
  TH2F *fHistPsiPair;          //!

  TH2F *fHistDecayVertexXBefore;      //!
  TH2F *fHistDecayVertexYBefore;      //!
  TH2F *fHistDecayVertexZBefore;      //!
  TH2F *fHistDecayVertexXAfter;       //!
  TH2F *fHistDecayVertexYAfter;       //!
  TH2F *fHistDecayVertexZAfter;       //!
  TH2F *fHistTransverseRadiusBefore;  //!
  TH2F *fHistTransverseRadiusAfter;   //!
  TH2F *fHistCosPABefore;             //!
  TH2F *fHistCosPAAfter;              //!
  TH2F *fHistDCADaughtersBefore;      //!
  TH2F *fHistDCADaughtersAfter;       //!
  TH2F *fHistDCA;                     //!
  TH2F *fHistDecayLength;             //!
  TH2F *fHistArmenterosBefore;        //!
  TH2F *fHistArmenterosAfter;         //!

  TH2F *fHistMCTruthV0PtY;                      //!
  TH2F *fHistMCTruthV0DaughterPtY;              //!
  TH2F *fHistMCTruthV0DaughterPtYAccept;        //!
  TH2F *fHistMCTruthPtYHighMult;                //!
  TH2F *fHistMCTruthDaughterPtYHighMult;        //!
  TH2F *fHistMCTruthDaughterPtYAcceptHighMult;  //!
  TH1F *fHistMCV0Pt;                            //!

  TH2F *fHistV0Mother;                                        //!
  TH2F *fHistV0MotherTrue;                                    //!
  TH2F *fHistV0MassPtTrue;                                    //!
  TH2F *fHistDecayVertexXTrue;                                //!
  TH2F *fHistDecayVertexYTrue;                                //!
  TH2F *fHistDecayVertexZTrue;                                //!
  TH2F *fHistTransverseRadiusTrue;                            //!
  TH2F *fHistCosPATrue;                                       //!
  TH2F *fHistDCADaughtersTrue;                                //!
  TH2F *fHistArmenterosTrue;                                  //!
  TH2F *fHistPsiPairTrue;                                     //!
  TH1F *fHistSingleParticlePtTrue[2];                         //!
  TH2F *fHistSingleParticleDCAtoPVTrue[2];                    //!
  TH2F *fHistSingleParticleNclsTPCTrue[2];                    //!
  TH2F *fHistSingleParticleNclsTPCFindableTrue[2];            //!
  TH2F *fHistSingleParticleNclsTPCRatioFindableTrue[2];       //!
  TH2F *fHistSingleParticleNcrossedTPCTrue[2];                //!
  TH2F *fHistSingleParticleNclsTPCSharedTrue[2];              //!
  TH2F *fHistSingleParticleNclsITSSharedTrue[2];              //!
  TH2F *fHistSingleParticleChi2True[2];                       //!
  TH2F *fHistSingleParticlePIDTrue[2];                        //!
  TH2F *fHistV0MassPtTrueSigma;                               //!
  TH2F *fHistDecayVertexXTrueSigma;                           //!
  TH2F *fHistDecayVertexYTrueSigma;                           //!
  TH2F *fHistDecayVertexZTrueSigma;                           //!
  TH2F *fHistTransverseRadiusTrueSigma;                       //!
  TH2F *fHistCosPATrueSigma;                                  //!
  TH2F *fHistDCADaughtersTrueSigma;                           //!
  TH2F *fHistArmenterosTrueSigma;                             //!
  TH2F *fHistPsiPairTrueSigma;                                //!
  TH1F *fHistSingleParticlePtTrueSigma[2];                    //!
  TH2F *fHistSingleParticleDCAtoPVTrueSigma[2];               //!
  TH2F *fHistSingleParticleNclsTPCTrueSigma[2];               //!
  TH2F *fHistSingleParticleNclsTPCFindableTrueSigma[2];       //!
  TH2F *fHistSingleParticleNclsTPCRatioFindableTrueSigma[2];  //!
  TH2F *fHistSingleParticleNcrossedTPCTrueSigma[2];           //!
  TH2F *fHistSingleParticleNclsTPCSharedTrueSigma[2];         //!
  TH2F *fHistSingleParticleNclsITSSharedTrueSigma[2];         //!
  TH2F *fHistSingleParticleChi2TrueSigma[2];                  //!
  TH2F *fHistSingleParticlePIDTrueSigma[2];                   //!
  TH2F *fHistV0MassPtBkg;                                     //!
  TH2F *fHistDecayVertexXBkg;                                 //!
  TH2F *fHistDecayVertexYBkg;                                 //!
  TH2F *fHistDecayVertexZBkg;                                 //!
  TH2F *fHistTransverseRadiusBkg;                             //!
  TH2F *fHistCosPABkg;                                        //!
  TH2F *fHistDCADaughtersBkg;                                 //!
  TH2F *fHistArmenterosBkg;                                   //!
  TH2F *fHistPsiPairBkg;                                      //!
  TH1F *fHistSingleParticlePtBkg[2];                          //!
  TH2F *fHistSingleParticleDCAtoPVBkg[2];                     //!
  TH2F *fHistSingleParticleNclsTPCBkg[2];                     //!
  TH2F *fHistSingleParticleNclsTPCFindableBkg[2];             //!
  TH2F *fHistSingleParticleNclsTPCRatioFindableBkg[2];        //!
  TH2F *fHistSingleParticleNcrossedTPCBkg[2];                 //!
  TH2F *fHistSingleParticleNclsTPCSharedBkg[2];               //!
  TH2F *fHistSingleParticleNclsITSSharedBkg[2];               //!
  TH2F *fHistSingleParticleChi2Bkg[2];                        //! .
  TH2F *fHistSingleParticlePIDBkg[2];                         //!

  TH1F *fHistSingleParticleCuts[2];                        //!
  TH1F *fHistSingleParticlePt[2];                          //!
  TH2F *fHistSingleParticleEtaBefore[2];                   //!
  TH2F *fHistSingleParticleEtaAfter[2];                    //!
  TH2F *fHistSingleParticleChi2Before[2];                  //!
  TH2F *fHistSingleParticleChi2After[2];                   //!
  TH2F *fHistSingleParticleNclsTPCBefore[2];               //!
  TH2F *fHistSingleParticleNclsTPCAfter[2];                //!
  TH2F *fHistSingleParticleNclsTPCFindableBefore[2];       //!
  TH2F *fHistSingleParticleNclsTPCFindableAfter[2];        //!
  TH2F *fHistSingleParticleNclsTPCRatioFindableBefore[2];  //!
  TH2F *fHistSingleParticleNclsTPCRatioFindableAfter[2];   //!
  TH2F *fHistSingleParticleNcrossedTPCBefore[2];           //!
  TH2F *fHistSingleParticleNcrossedTPCAfter[2];            //!
  TH2F *fHistSingleParticleNclsTPCSharedBefore[2];         //!
  TH2F *fHistSingleParticleNclsITSSharedBefore[2];         //!
  TH2F *fHistSingleParticleNclsTPCSharedAfter[2];          //!
  TH2F *fHistSingleParticleNclsITSSharedAfter[2];          //!
  TH2F *fHistSingleParticleDCAtoPVBefore[2];               //!
  TH2F *fHistSingleParticleDCAtoPVAfter[2];                //!
  TH2F *fHistSingleParticlePileUp[2];                      //!
  TH2F *fHistSingleParticlePID[2];                         //!

 private:
  ClassDef(AliSigma0V0Cuts, 10)
};

#endif
