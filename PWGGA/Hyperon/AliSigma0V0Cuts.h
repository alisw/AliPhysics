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

  void SelectLambda(AliVEvent *inputEvent, AliMCEvent *mcEvent,
                    std::vector<AliSigma0ParticleV0> &V0Container);
  void SelectPhoton(AliVEvent *inputEvent, AliMCEvent *mcEvent,
                    std::vector<AliSigma0ParticleV0> &V0Container,
                    AliPID::EParticleType particle1,
                    AliPID::EParticleType particle2);
  bool V0QualityCuts(const AliESDv0 *v0);
  bool V0PID(const AliESDv0 *v0, const AliESDtrack *pos, const AliESDtrack *neg,
             AliPID::EParticleType particle,
             AliPID::EParticleType antiParticle);
  bool SingleParticlePID(const AliVTrack *track, AliPID::EParticleType particle,
                         float &prob) const;
  void PlotSingleParticlePID(const AliVTrack *track,
                             AliPID::EParticleType particle) const;
  bool SingleParticleQualityCuts(AliESDtrack *track);
  bool PileUpRejection(AliESDtrack *pos, AliESDtrack *neg);
  bool V0TopologicalSelection(const AliESDv0 *v0);
  bool LambdaSelection(AliESDv0 *v0);
  float ComputeRapidity(float pt, float pz, float m) const;
  int GetRapidityBin(float rapidity) const;

  void SetLightweight(bool isLightweight) { fIsLightweight = isLightweight; }

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetPileUpRejectionMode(PileUpRejectionMode pileUpRej) {
    fPileUpRejectionMode = pileUpRej;
  }
  void SetPID(int pid) { fPID = pid; }
  void SetPosPID(AliPID::EParticleType pid) { fPosPID = pid; }
  void SetNegPID(AliPID::EParticleType pid) { fNegPID = pid; }
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
  void SetDaughterDCAMax(float ddcaMax) { fDaughterDCAMax = ddcaMax; }
  void SetDaughterDCAtoPV(float dca2pv) { fDaughterDCAPV = dca2pv; }
  void SetK0Rejection(float low, float up) {
    fK0Rejection = true;
    fK0RejectionLow = low, fK0RejectionUp = up;
  }
  void SetArmenterosCut(float qtLow, float qtUp, float alphaLow,
                        float alphaUp) {
    fUsePID = false;
    fArmenterosQtLow = qtLow;
    fArmenterosQtUp = qtUp;
    fArmenterosAlphaLow = alphaLow;
    fArmenterosAlphaUp = alphaUp;
  }
  void SetLambdaSelection(float low, float up) {
    fLambdaSelectionLow = low, fLambdaSelectionUp = up;
  }

  void InitCutHistograms(TString appendix = TString(""));
  TList *GetCutHistograms() const { return fHistograms; }

 protected:
  TList *fHistograms;
  TList *fHistogramsMC;
  TList *fHistogramsBefore;
  TList *fHistogramsAfter;
  TList *fHistogramsPos;
  TList *fHistogramsNeg;

  AliESDEvent *fInputEvent;   //!
  AliMCEvent *fMCEvent;       //!
  TDatabasePDG fPDGDatabase;  //!

  bool fIsLightweight;  //

  short fV0cut;
  short fAntiV0cut;
  short fPID;

  bool fIsMC;
  PileUpRejectionMode fPileUpRejectionMode;
  AliPID::EParticleType fPosPID;
  AliPID::EParticleType fNegPID;
  bool fV0OnFly;
  bool fK0Rejection;
  bool fUsePID;
  float fV0PtMin;
  float fV0PtMax;
  float fV0CosPAMin;
  float fV0RadiusMax;
  float fV0RadiusMin;
  float fV0DecayVertexMax;
  float fPIDnSigma;
  float fEtaMax;
  float fTPCclusterMin;
  short fTPCnCrossedRowsMin;
  float fTPCratioFindable;
  short fTPCfindableMin;
  short fTPCnSharedMax;
  float fDaughterDCAMax;
  float fDaughterDCAPV;
  float fK0RejectionLow;
  float fK0RejectionUp;
  float fArmenterosQtLow;
  float fArmenterosQtUp;
  float fArmenterosAlphaLow;
  float fArmenterosAlphaUp;
  float fLambdaSelectionLow;
  float fLambdaSelectionUp;

  AliPIDResponse *fPIDResponse;  //! pid response

  // Histograms
  // =====================================================================
  TProfile *fHistCutBooking;  //

  TH1F *fHistCuts;  //
  TH1F *fHistNV0;   //

  TH1F *fHistLambdaMass;       //
  TH1F *fHistLambdaPt;         //
  TH2F *fHistLambdaPtY[20];    //
  TH2F *fHistLambdaMassPt;     //
  TH1F *fHistLambdaMassK0Rej;  //
  TH1F *fHistK0Mass;           //
  TH1F *fHistK0MassAfter;      //
  TH2F *fHistCosPA;            //
  TH2F *fHistEtaPhi;           //

  TH1F *fHistDecayVertexXBefore;      //
  TH1F *fHistDecayVertexYBefore;      //
  TH1F *fHistDecayVertexZBefore;      //
  TH1F *fHistDecayVertexXAfter;       //
  TH1F *fHistDecayVertexYAfter;       //
  TH1F *fHistDecayVertexZAfter;       //
  TH1F *fHistTransverseRadiusBefore;  //
  TH1F *fHistTransverseRadiusAfter;   //
  TH1F *fHistCosPABefore;             //
  TH1F *fHistCosPAAfter;              //
  TH1F *fHistDCADaughtersBefore;      //
  TH1F *fHistDCADaughtersAfter;       //
  TH1F *fHistDCA;                     //
  TH1F *fHistDecayLength;             //
  TH2F *fHistArmenterosBefore;        //
  TH2F *fHistArmenterosAfter;         //

  TH1F *fHistMCTruthV0Pt;               //
  TH2F *fHistMCTruthV0PtY;              //
  TH2F *fHistMCTruthV0PtEta;            //
  TH1F *fHistMCTruthV0ProtonPionPt;     //
  TH2F *fHistMCTruthV0ProtonPionPtY;    //
  TH2F *fHistMCTruthV0ProtonPionPtEta;  //

  TH1F *fHistSingleParticleCuts[2];                        //
  TH1F *fHistSingleParticlePt[2];                          //
  TH1F *fHistSingleParticleEtaBefore[2];                   //
  TH1F *fHistSingleParticleEtaAfter[2];                    //
  TH1F *fHistSingleParticleNclsTPCBefore[2];               //
  TH1F *fHistSingleParticleNclsTPCAfter[2];                //
  TH1F *fHistSingleParticleNclsTPCFindableBefore[2];       //
  TH1F *fHistSingleParticleNclsTPCFindableAfter[2];        //
  TH1F *fHistSingleParticleNclsTPCRatioFindableBefore[2];  //
  TH1F *fHistSingleParticleNclsTPCRatioFindableAfter[2];   //
  TH1F *fHistSingleParticleNcrossedTPCBefore[2];           //
  TH1F *fHistSingleParticleNcrossedTPCAfter[2];            //
  TH1F *fHistSingleParticleNclsTPCShared[2];               //
  TH1F *fHistSingleParticleNclsITSShared[2];               //
  TH1F *fHistSingleParticleDCAtoPVBefore[2];               //
  TH1F *fHistSingleParticleDCAtoPVAfter[2];                //
  TH2F *fHistSingleParticlePileUp[2];                      //
  TH2F *fHistSingleParticlePID[2];                         //

 private:
  ClassDef(AliSigma0V0Cuts, 3)
};

#endif
