#ifndef AliSigma0V0Cuts_H
#define AliSigma0V0Cuts_H

#include "AliMCEvent.h"
#include "AliPIDResponse.h"
#include "AliSigma0ParticleV0.h"
#include "Riostream.h"
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

  static AliSigma0V0Cuts *DefaultCuts();

  void SelectV0s(AliVEvent *inputEvent, AliMCEvent *mcEvent,
                 std::vector<AliSigma0ParticleV0> &V0Container,
                 std::vector<AliSigma0ParticleV0> &AntiV0Container,
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
  TList *fHistogramsV0;
  TList *fHistogramsV0MC;
  TList *fHistogramsV0Before;
  TList *fHistogramsV0After;
  TList *fHistogramsV0Pos;
  TList *fHistogramsV0Neg;
  TList *fHistogramsAntiV0;
  TList *fHistogramsAntiV0MC;
  TList *fHistogramsAntiV0Before;
  TList *fHistogramsAntiV0After;
  TList *fHistogramsAntiV0Pos;
  TList *fHistogramsAntiV0Neg;

  AliESDEvent *fInputEvent;  //!
  AliMCEvent *fMCEvent;      //!

  std::vector<AliSigma0ParticleV0> fV0Vector;      //!
  std::vector<AliSigma0ParticleV0> fAntiV0Vector;  //!

  bool fIsLightweight;  //

  short fV0cut;
  short fAntiV0cut;
  short fPID;

  bool fIsMC;
  PileUpRejectionMode fPileUpRejectionMode;
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
  TProfile *fHistCuts;  //

  TH1F *fHistV0Cuts;  //
  TH1F *fHistNV0;     //

  TH1F *fHistV0LambdaMass;       //
  TH1F *fHistV0LambdaPt;         //
  TH2F *fHistV0LambdaPtY[20];    //
  TH2F *fHistV0LambdaMassPt;     //
  TH1F *fHistV0LambdaMassK0Rej;  //
  TH1F *fHistV0K0Mass;           //
  TH1F *fHistV0K0MassAfter;      //
  TH2F *fHistV0CosPA;            //
  TH2F *fHistV0EtaPhi;           //

  TH1F *fHistV0DecayVertexXBefore;      //
  TH1F *fHistV0DecayVertexYBefore;      //
  TH1F *fHistV0DecayVertexZBefore;      //
  TH1F *fHistV0DecayVertexXAfter;       //
  TH1F *fHistV0DecayVertexYAfter;       //
  TH1F *fHistV0DecayVertexZAfter;       //
  TH1F *fHistV0TransverseRadiusBefore;  //
  TH1F *fHistV0TransverseRadiusAfter;   //
  TH1F *fHistV0CosPABefore;             //
  TH1F *fHistV0CosPAAfter;              //
  TH1F *fHistV0DCADaughtersBefore;      //
  TH1F *fHistV0DCADaughtersAfter;       //
  TH1F *fHistV0DCA;                     //
  TH1F *fHistV0DecayLength;             //
  TH2F *fHistV0ArmenterosBefore;        //
  TH2F *fHistV0ArmenterosAfter;         //

  TH1F *fHistMCTruthV0Pt;               //
  TH2F *fHistMCTruthV0PtY;              //
  TH2F *fHistMCTruthV0PtEta;            //
  TH1F *fHistMCTruthV0ProtonPionPt;     //
  TH2F *fHistMCTruthV0ProtonPionPtY;    //
  TH2F *fHistMCTruthV0ProtonPionPtEta;  //

  TH1F *fHistAntiV0Cuts;  //
  TH1F *fHistNAntiV0;     //

  TH1F *fHistAntiV0LambdaMass;       //
  TH1F *fHistAntiV0LambdaPt;         //
  TH2F *fHistAntiV0LambdaPtY[20];    //
  TH2F *fHistAntiV0LambdaMassPt;     //
  TH1F *fHistAntiV0LambdaMassK0Rej;  //
  TH1F *fHistAntiV0K0Mass;           //
  TH1F *fHistAntiV0K0MassAfter;      //
  TH2F *fHistAntiV0CosPA;            //
  TH2F *fHistAntiV0EtaPhi;           //

  TH1F *fHistAntiV0DecayVertexXBefore;      //
  TH1F *fHistAntiV0DecayVertexYBefore;      //
  TH1F *fHistAntiV0DecayVertexZBefore;      //
  TH1F *fHistAntiV0DecayVertexXAfter;       //
  TH1F *fHistAntiV0DecayVertexYAfter;       //
  TH1F *fHistAntiV0DecayVertexZAfter;       //
  TH1F *fHistAntiV0TransverseRadiusBefore;  //
  TH1F *fHistAntiV0TransverseRadiusAfter;   //
  TH1F *fHistAntiV0CosPABefore;             //
  TH1F *fHistAntiV0CosPAAfter;              //
  TH1F *fHistAntiV0DCADaughtersBefore;      //
  TH1F *fHistAntiV0DCADaughtersAfter;       //
  TH1F *fHistAntiV0DCA;                     //
  TH1F *fHistAntiV0DecayLength;             //
  TH2F *fHistAntiV0ArmenterosBefore;        //
  TH2F *fHistAntiV0ArmenterosAfter;         //

  TH1F *fHistMCTruthAntiV0Pt;     //
  TH2F *fHistMCTruthAntiV0PtY;    //
  TH2F *fHistMCTruthAntiV0PtEta;  //

  TH1F *fHistMCTruthAntiV0ProtonPionPt;     //
  TH2F *fHistMCTruthAntiV0ProtonPionPtY;    //
  TH2F *fHistMCTruthAntiV0ProtonPionPtEta;  //

  TH1F *fHistV0SingleParticleCuts[2];                        //
  TH1F *fHistV0SingleParticlePt[2];                          //
  TH1F *fHistV0SingleParticleEtaBefore[2];                   //
  TH1F *fHistV0SingleParticleEtaAfter[2];                    //
  TH1F *fHistV0SingleParticleNclsTPCBefore[2];               //
  TH1F *fHistV0SingleParticleNclsTPCAfter[2];                //
  TH1F *fHistV0SingleParticleNclsTPCFindableBefore[2];       //
  TH1F *fHistV0SingleParticleNclsTPCFindableAfter[2];        //
  TH1F *fHistV0SingleParticleNclsTPCRatioFindableBefore[2];  //
  TH1F *fHistV0SingleParticleNclsTPCRatioFindableAfter[2];   //
  TH1F *fHistV0SingleParticleNcrossedTPCBefore[2];           //
  TH1F *fHistV0SingleParticleNcrossedTPCAfter[2];            //
  TH1F *fHistV0SingleParticleNclsTPCShared[2];               //
  TH1F *fHistV0SingleParticleNclsITSShared[2];               //
  TH1F *fHistV0SingleParticleDCAtoPVBefore[2];               //
  TH1F *fHistV0SingleParticleDCAtoPVAfter[2];                //
  TH2F *fHistV0SingleParticlePileUp[2];                      //
  TH2F *fHistV0SingleParticlePID[2];                         //

  TH1F *fHistAntiV0SingleParticleCuts[2];                        //
  TH1F *fHistAntiV0SingleParticlePt[2];                          //
  TH1F *fHistAntiV0SingleParticleEtaBefore[2];                   //
  TH1F *fHistAntiV0SingleParticleEtaAfter[2];                    //
  TH1F *fHistAntiV0SingleParticleNclsTPCBefore[2];               //
  TH1F *fHistAntiV0SingleParticleNclsTPCAfter[2];                //
  TH1F *fHistAntiV0SingleParticleNclsTPCFindableBefore[2];       //
  TH1F *fHistAntiV0SingleParticleNclsTPCFindableAfter[2];        //
  TH1F *fHistAntiV0SingleParticleNclsTPCRatioFindableBefore[2];  //
  TH1F *fHistAntiV0SingleParticleNclsTPCRatioFindableAfter[2];   //
  TH1F *fHistAntiV0SingleParticleNcrossedTPCBefore[2];           //
  TH1F *fHistAntiV0SingleParticleNcrossedTPCAfter[2];            //
  TH1F *fHistAntiV0SingleParticleNclsTPCShared[2];               //
  TH1F *fHistAntiV0SingleParticleNclsITSShared[2];               //
  TH1F *fHistAntiV0SingleParticleDCAtoPVBefore[2];               //
  TH1F *fHistAntiV0SingleParticleDCAtoPVAfter[2];                //
  TH2F *fHistAntiV0SingleParticlePileUp[2];                      //
  TH2F *fHistAntiV0SingleParticlePID[2];                         //

 private:
  ClassDef(AliSigma0V0Cuts, 3)
};

#endif
