#ifndef AliSigma0PhotonCuts_H
#define AliSigma0PhotonCuts_H

#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliPID.h"
#include "Riostream.h"
#include "TDatabasePDG.h"
#include "TObject.h"
#include "TProfile.h"
#include "TH2F.h"
#include "AliAODConversionPhoton.h"
#include "AliPIDResponse.h"
#include "AliFemtoDreamBasePart.h"

class AliSigma0PhotonCuts : public TObject {
 public:

  AliSigma0PhotonCuts();
  AliSigma0PhotonCuts(const AliSigma0PhotonCuts &);
  AliSigma0PhotonCuts &operator=(const AliSigma0PhotonCuts &);
  virtual ~AliSigma0PhotonCuts();

  static AliSigma0PhotonCuts *PhotonCuts();

  void PhotonCuts(AliAODEvent *inputEvent, AliMCEvent *mcEvent,
                  const TClonesArray *photons,
                  std::vector<AliFemtoDreamBasePart> &container);
  bool ProcessPhoton(AliVEvent* event, AliMCEvent *mcEvent, AliAODConversionPhoton *PhotonCandidate,
                     AliAODv0 *v0, AliVTrack *pos, AliVTrack *neg);
  bool PiZeroRejection(AliAODConversionPhoton*photon,
                       const TClonesArray *photons,
                       int iPhoton);
  void RelabelAODPhotonCandidates(AliAODConversionPhoton *PhotonCandidate, AliVEvent* event);
  float ComputeInvMass(AliAODv0 *v0, AliVTrack *pos, AliVTrack *neg, int pdgPos,
                       int pgdNeg) const;
  AliVTrack *GetTrack(AliVEvent * event, int label);
  void SetLightweight(bool isLightweight) {
    fIsLightweight = isLightweight;
  }

  void SetIsMC(bool isMC) {
    fIsMC = isMC;
  }

  void SetV0ReaderName(TString name) {
    fV0ReaderName = name;
  }
  void SetArmenteros2DCut(bool doIt) {
    f2DArmenterosCut = doIt;
  }
  void SetArmenterosQtMax(float qtmax) {
    fQtMax = qtmax;
  }
  void SetPhotonPtMin(float ptmin) {
    fPhotonPtMin = ptmin;
  }
  void SetPhotonEtaMax(float etamax) {
    fPhotonEtaMax = etamax;
  }
  void SetPsiPair2DCut(bool doIt) {
    f2DPsiPairCut = doIt;
  }
  void SetPsiPairMax(float psimax) {
    fPsiPairMax = psimax;
  }
  void SetChi2ForPsiMax(float chi2max) {
    fChi2MaxFor2DPsiPair = chi2max;
  }
  void SetCPAMin(float cpamin) {
    fCPAMin = cpamin;
  }
  void SetDCAzMax(float dcazmax) {
    fDCAzMax = dcazmax;
  }
  void SetDCArMax(float dcarmax) {
    fDCArMax = dcarmax;
  }

  void SetElectronPtMin(float ptmin) {
    fElectronPtMin = ptmin;
  }
  void SetElectronEtaMax(float etamax) {
    fElectronEtaMax = etamax;
  }
  void SetElectronRatioFindable(float ratioMax) {
    fElectronRatioFindable = ratioMax;
  }
  void SetElectronNSigmaTPCMax(float nsigmamax) {
    fElectronNSigmaTPCMax = nsigmamax;
  }
  void SetElectronNSigmaTPCMin(float nsigmamin) {
    fElectronNSigmaTPCMin = nsigmamin;
  }
  void SetTransverseRadiusRejection(float low, float up) {
    fDoTransvRadRejection = true;
    fTransvRadRejectionLow = low;
    fTransvRadRejectionUp = up;
  }
  void SetPhotonQualityCut(int cut) {
    fDoPhotonQualityCut = true;
    fPhotonQuality = cut;
  }
  void SetPhotonPileUpCut(int cut) {
    fDoPhotonPileupCut = true;
    fPhotonPileupCut = cut;
  }
  void DoPiZeroRejection(bool doIt) {
    fDoExtraPiZeroRejection = doIt;
  }

  void InitCutHistograms(TString appendix = TString(""));
  TList *GetCutHistograms() const {
    return fHistograms;
  }

 protected:
  TList *fHistograms;        //!
  TList *fHistogramsBefore;  //!
  TList *fHistogramsAfter;   //!
  TList *fHistogramsPos;     //!
  TList *fHistogramsNeg;     //!

  AliVEvent *fInputEvent;   //!
  AliMCEvent *fMCEvent;       //!
  TDatabasePDG fDataBasePDG;  //!

  bool fIsLightweight;  //
  TString fV0ReaderName;                         //

  bool fIsMC;                                //
  bool f2DArmenterosCut;                     //
  float fQtMax;                              //
  float fPhotonPtMin;                            //
  float fPhotonEtaMax;  //
  bool f2DPsiPairCut;  //
  float fPsiPairMax;  //
  float fChi2MaxFor2DPsiPair;  //
  float fCPAMin;  //
  float fDCAzMax;  //
  float fDCArMax;  //

  float fElectronPtMin;  //
  float fElectronEtaMax;  //
  float fElectronRatioFindable;  //
  float fElectronNSigmaTPCMax;  //
  float fElectronNSigmaTPCMin;  //

  bool fDoTransvRadRejection;   //
  float fTransvRadRejectionLow; //
  float fTransvRadRejectionUp;  //

  bool fDoPhotonQualityCut; //
  int fPhotonQuality; //

  bool fDoPhotonPileupCut; //
  bool fPhotonPileupCut; //

  bool fDoExtraPiZeroRejection; //

  // Histograms
  // =====================================================================
  TProfile *fHistCutBooking;  //!

  TH1F *fHistCuts;  //!
  TH1F *fHistNV0;   //!

  TH1F *fHistPiZeroMass;       //!
  TH1F *fHistLambdaMass;       //!
  TH1F *fHistAntiLambdaMass;   //!
  TH1F *fHistK0Mass;           //!
  TH1F *fHistV0Pt;             //!
  TH1F *fHistV0Mass;           //!
  TH2F *fHistV0MassPt;         //!
  TH2F *fHistCosPA;            //!
  TH2F *fHistEtaPhi;           //!
  TH1F* fHistMCV0Pt;           //!
  TH2F* fHistV0Mother;         //!
  TH2F *fHistPsiPairBefore;          //!
  TH2F *fHistPsiPairAfter;          //!

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
  TH2F *fHistDCArBefore;              //!
  TH2F *fHistDCArAfter;               //!
  TH2F *fHistDCAzBefore;              //!
  TH2F *fHistDCAzAfterOthersBefore;   //!
  TH2F *fHistDCAzAfterOthersBeforeQuality[4];   //!
  TH2F *fHistDCAzAfter;               //!
  TH2F *fHistDCA;                     //!
  TH2F *fHistDecayLength;             //!
  TH2F *fHistArmenterosBefore;        //!
  TH2F *fHistArmenterosAfter;         //!
  TH2F *fHistQualityBefore;           //!
  TH2F *fHistQualityAfter;            //!

  TH2F *fHistTomography;              //!

  TH1F* fHistMCTruthPhotonPt;         //;
  TH1F* fHistMCTruthPhotonSigmaPt;    //;
  TH1F* fHistMCPhotonPt;              //;
  TH1F* fHistMCPhotonSigmaPt;         //;
  TH1F* fHistMCTruthPhotonP;          //;
  TH1F* fHistMCTruthPhotonSigmaP;     //;
  TH1F* fHistMCPhotonP;               //;
  TH1F* fHistMCPhotonSigmaP;          //;
  TH2F* fHistMCPhotonSource;          //;

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
  TH2F *fHistSingleParticleDCAtoPVBefore[2];               //!
  TH2F *fHistSingleParticleDCAtoPVAfter[2];                //!
  TH2F *fHistSingleParticlePileUp[2];                      //!
  TH2F *fHistSingleParticlePID[2];                         //!

  AliPIDResponse *fPIDResponse;  //!  pid response

 private:
ClassDef(AliSigma0PhotonCuts, 9)
};

#endif
