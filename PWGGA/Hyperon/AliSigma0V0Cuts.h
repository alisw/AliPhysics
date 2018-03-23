#ifndef AliSigma0V0Cuts_H
#define AliSigma0V0Cuts_H

#include "AliAODTrack.h"
#include "AliAODv0.h"
#include "AliMCEvent.h"
#include "AliPIDResponse.h"
#include "AliSigma0ParticleV0.h"
#include "AliSigma0SingleParticleCuts.h"
#include "AliVEvent.h"
#include "Riostream.h"
#include "TObject.h"

#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TProfile.h"

class AliPIDResponse;

class AliSigma0V0Cuts : public TObject {
 public:
  AliSigma0V0Cuts();
  AliSigma0V0Cuts(const AliSigma0V0Cuts &);
  AliSigma0V0Cuts &operator=(const AliSigma0V0Cuts &);
  virtual ~AliSigma0V0Cuts();

  static AliSigma0V0Cuts *DefaultCuts();
  static AliSigma0V0Cuts *Sigma0Cuts();

  void SelectV0s(AliVEvent *inputEvent, AliMCEvent *mcEvent,
                 std::vector<AliSigma0ParticleV0> &V0Container,
                 std::vector<AliSigma0ParticleV0> &AntiV0Container,
                 AliPID::EParticleType particle1,
                 AliPID::EParticleType particle2);
  void ProcessESDs(AliPID::EParticleType particle,
                   AliPID::EParticleType antiParticle);
  void ProcessAODs(AliPID::EParticleType particle,
                   AliPID::EParticleType antiParticle);
  bool V0QualityCuts(const AliESDv0 *v0);
  bool V0QualityCuts(const AliAODv0 *v0);
  bool V0PID(const AliVTrack *pos, const AliVTrack *neg,
             AliPID::EParticleType particle,
             AliPID::EParticleType antiParticle);
  bool SingleParticlePID(const AliVTrack *track, AliPID::EParticleType particle,
                         float &prob) const;
  void PlotSingleParticlePID(const AliVTrack *track,
                             AliPID::EParticleType particle) const;
  bool V0TopologicalSelection(const AliAODv0 *v0);
  bool V0TopologicalSelection(const AliESDv0 *v0);
  bool SingleParticleQualityCuts(AliVTrack *track, const float dcaDaughterToPV);
  template <typename T>
  bool LambdaSelection(const T *v0, float massK0, float massLambda);
  bool CheckIfRealV0(const AliAODv0 *v0, int PDGmother, int PDGdaugh1,
                     int PDGdaugh2) const;
  void ProcessMC() const;
  bool IsLambdaProtonPion(AliMCParticle *particle) const;
  float ComputeRapidity(float pt, float pz, float m) const;
  int GetRapidityBin(float rapidity) const;

  void SetSingleParticleCuts(AliSigma0SingleParticleCuts *cuts) {
    fSingleParticleCuts = cuts;
  }

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetPileUpRejection(bool isPileUpRej) { fPileUpRejection = isPileUpRej; }
  void SetExtendedQA(bool isExtended) { fIsExtendedQA = isExtended; }
  void SetV0OnFlyStatus(bool onFly) { fV0OnFly = onFly; }
  void SetV0PtMin(float ptMin) { fV0PtMin = ptMin; }
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
    fArmenterosCut = false;
    fK0RejectionLow = low, fK0RejectionUp = up;
  }
  void SetArmenterosCut(float qtLow, float qtUp, float alphaLow,
                        float alphaUp) {
    fArmenterosCut = true;
    fArmenterosQtLow = qtLow;
    fArmenterosQtUp = qtUp;
    fArmenterosAlphaLow = alphaLow;
    fArmenterosAlphaUp = alphaUp;
  }
  void SetLambdaSelection(float low, float up) {
    fLambdaSelectionLow = low, fLambdaSelectionUp = up;
  }

  void InitCutHistograms(const char *appendix = "default");
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

  AliVEvent *fInputEvent;  //!
  AliMCEvent *fMCEvent;    //!

  AliSigma0SingleParticleCuts *fSingleParticleCuts;

  std::vector<AliSigma0ParticleV0> fV0Vector;      //!
  std::vector<AliSigma0ParticleV0> fAntiV0Vector;  //!

  std::vector<AliAODTrack *> *fGlobalTrackReference;  //!

  short fV0cut;
  short fAntiV0cut;
  short fPID;

  bool fIsMC;
  bool fPileUpRejection;
  bool fIsExtendedQA;
  bool fV0OnFly;
  bool fArmenterosCut;
  bool fUsePID;
  float fV0PtMin;
  float fV0CosPAMin;
  float fV0RadiusMax;
  float fV0RadiusMin;
  float fV0DecayVertexMax;
  float fPIDnSigma;
  float fEtaMax;
  float fTPCclusterMin;
  short fTPCnCrossedRowsMin;
  float	fTPCratioFindable;
  short fTPCfindableMin;
  short	fTPCnSharedMax;
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

  TH1F *fHistMCTruthV0Pt;     //
  TH2F *fHistMCTruthV0PtY;    //
  TH2F *fHistMCTruthV0PtEta;  //
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

  TH1F *fHistV0SingleParticleCuts[2];           //
  TH1F *fHistV0SingleParticlePt[2];             //
  TH1F *fHistV0SingleParticleEtaBefore[2];      //
  TH1F *fHistV0SingleParticleEtaAfter[2];       //
  TH1F *fHistV0SingleParticleNclsTPCBefore[2];  //
  TH1F *fHistV0SingleParticleNclsTPCAfter[2];   //
  TH1F *fHistV0SingleParticleNclsTPCFindableBefore[2];  //
  TH1F *fHistV0SingleParticleNclsTPCFindableAfter[2];   //
  TH1F *fHistV0SingleParticleNclsTPCRatioFindableBefore[2];  //
  TH1F *fHistV0SingleParticleNclsTPCRatioFindableAfter[2];   //
  TH1F *fHistV0SingleParticleNcrossedTPCBefore[2];  //
  TH1F *fHistV0SingleParticleNcrossedTPCAfter[2];   //
  TH1F *fHistV0SingleParticleNclsTPCShared[2];       //
  TH2F *fHistV0SingleParticleNclsTPCSharedTiming[2]; //
  TH1F *fHistV0SingleParticleNclsITSShared[2];       //
  TH2F *fHistV0SingleParticleNclsITSSharedTiming[2]; //
  TH1F *fHistV0SingleParticleDCAtoPVBefore[2];  //
  TH1F *fHistV0SingleParticleDCAtoPVAfter[2];   //
  TH2F *fHistV0SingleParticlePID[2];            //

  TH1F *fHistAntiV0SingleParticleCuts[2];           //
  TH1F *fHistAntiV0SingleParticlePt[2];             //
  TH1F *fHistAntiV0SingleParticleEtaBefore[2];      //
  TH1F *fHistAntiV0SingleParticleEtaAfter[2];       //
  TH1F *fHistAntiV0SingleParticleNclsTPCBefore[2];  //
  TH1F *fHistAntiV0SingleParticleNclsTPCAfter[2];   //
  TH1F *fHistAntiV0SingleParticleNclsTPCFindableBefore[2];  //
  TH1F *fHistAntiV0SingleParticleNclsTPCFindableAfter[2];   //
  TH1F *fHistAntiV0SingleParticleNclsTPCRatioFindableBefore[2];  //
  TH1F *fHistAntiV0SingleParticleNclsTPCRatioFindableAfter[2];   //
  TH1F *fHistAntiV0SingleParticleNcrossedTPCBefore[2];  //
  TH1F *fHistAntiV0SingleParticleNcrossedTPCAfter[2];   //
  TH1F *fHistAntiV0SingleParticleNclsTPCShared[2];       //
  TH2F *fHistAntiV0SingleParticleNclsTPCSharedTiming[2]; //
  TH1F *fHistAntiV0SingleParticleNclsITSShared[2];       //
  TH2F *fHistAntiV0SingleParticleNclsITSSharedTiming[2]; //
  TH1F *fHistAntiV0SingleParticleDCAtoPVBefore[2];  //
  TH1F *fHistAntiV0SingleParticleDCAtoPVAfter[2];   //
  TH2F *fHistAntiV0SingleParticlePID[2];            //

 private:
  ClassDef(AliSigma0V0Cuts, 2)
};

#endif
