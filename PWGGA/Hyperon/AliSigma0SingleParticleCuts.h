#ifndef AliSigma0SingleParticleCuts_H
#define AliSigma0SingleParticleCuts_H

#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliPIDResponse.h"
#include "AliSigma0ParticleBase.h"
#include "AliVEvent.h"
#include "Riostream.h"
#include "TObject.h"

#include "TH1.h"
#include "TH2.h"
#include "TList.h"
#include "TProfile.h"

class AliPIDResponse;

class AliSigma0SingleParticleCuts : public TObject {
 public:
  AliSigma0SingleParticleCuts();
  AliSigma0SingleParticleCuts(const AliSigma0SingleParticleCuts &);
  AliSigma0SingleParticleCuts &operator=(const AliSigma0SingleParticleCuts &);
  virtual ~AliSigma0SingleParticleCuts();

  static AliSigma0SingleParticleCuts *DefaultCuts();
  static AliSigma0SingleParticleCuts *ElectronCuts();

  void SelectSingleParticles(
      AliVEvent *inputEvent, AliMCEvent *mcEvent,
      std::vector<AliSigma0ParticleBase> &ParticleContainer,
      std::vector<AliSigma0ParticleBase> &AntiParticleContainer);
  void ProcessESDs();
  void ProcessAODs();
  void StoreGlobalTrackReference();
  bool SingleParticlePID(const AliVTrack *track) const;
  bool SingleParticleQualityCuts(AliVTrack *track) const;
  template <typename T>
  bool GoodTPCFitMapSharedMap(const T *track) const;
  bool SingleParticleDCASelection(AliVTrack *track, float &DCAr,
                                  float &DCAz) const;
  void ProcessMC(AliSigma0ParticleBase *baseParticle,
                 AliMCParticle *mcParticle) const;
  float GetBeta(const AliVTrack *track) const;

  void SetIsMC(bool isMC) { fIsMC = isMC; }
  void SetExtendedQA(bool isExtended) { fIsExtendedQA = isExtended; }
  void SetParticle(AliPID::EParticleType particle) { fParticle = particle; }
  void SetParticleName(TString name) { fParticleName = name; }
  void SetParticleRejection(AliPID::EParticleType part1,
                            AliPID::EParticleType part2,
                            AliPID::EParticleType part3) {
    fParticleVeto1 = part1;
    fParticleVeto2 = part2;
    fParticleVeto3 = part3;
  }
  void SetPDGCode(int pdgCode) { fPDGCode = pdgCode; }
  void SetFilterBit(int filterBit) { fFilterbit = filterBit; }
  void SetTPCclusterMin(short nTPCcluster) { fTPCclusterMin = nTPCcluster; }
  void SetTPCnCrossedRowsMin(short nCrossedRows) {
    fTPCnCrossedRowsMin = nCrossedRows;
  }
  void SetTPCFindableMin(float findable) { fTPCfindableMin = findable; }
  void SetUseSharedMap(bool useSharedMap) { fUseSharedMap = useSharedMap; }
  void SetDCArMax(float dcaRmax) { fDCArMax = dcaRmax; }
  void SetDCAzMax(float dcaZmax) { fDCAzMax = dcaZmax; }
  void SetPtMin(float ptMin) { fPtMin = ptMin; }
  void SetPtMax(float ptMax) { fPtMax = ptMax; }
  void SetEtaMax(float etaMax) { fEtaMax = etaMax; }
  void SetPIDnSigmaStrict(float nSigma) { fPIDnSigmaStrict = nSigma; }
  void SetPIDHypothesisRejection(bool isRej) {
    fPIDHypothesisRejection = isRej;
  }
  void SetPIDMomentumSwitch(float momSwitch) { fPIDMomentumSwitch = momSwitch; }
  void SetPIDTOFif(bool tofif) { fPIDTOFif = tofif; }
  void SetAdditionalPIDRejection(AliPID::EParticleType part, float lower,
                                 float upper, float pTlower = 0.f,
                                 float pTupper = 1E30);

  void InitCutHistograms();
  TList *GetCutHistograms() const { return fHistograms; }
  std::vector<AliAODTrack *> *GetGlobalTrackReference() {
    return &fGlobalTrackReference;
  }

 protected:
  TList *fHistograms;
  TList *fHistogramsSingleParticle;
  TList *fHistogramsSingleParticleMC;
  TList *fHistogramsSingleParticleBefore;
  TList *fHistogramsSingleParticleAfter;
  TList *fHistogramsSingleAntiParticle;
  TList *fHistogramsSingleAntiParticleMC;
  TList *fHistogramsSingleAntiParticleBefore;
  TList *fHistogramsSingleAntiParticleAfter;

  bool fIsMC;

  AliVEvent *fInputEvent;  //!
  AliMCEvent *fMCEvent;    //!

  std::vector<AliSigma0ParticleBase> fSingleParticleVector;      //!
  std::vector<AliSigma0ParticleBase> fSingleAntiParticleVector;  //!
  std::vector<AliAODTrack *> fGlobalTrackReference;              //!
  short fTrackReferenceSize;  //! default value of above vector

  TString fParticleName;
  AliPID::EParticleType fParticle;
  AliPID::EParticleType fParticleVeto1;
  AliPID::EParticleType fParticleVeto2;
  AliPID::EParticleType fParticleVeto3;
  int fPDGCode;
  bool fIsExtendedQA;
  int fFilterbit;
  short fTPCclusterMin;
  short fTPCnCrossedRowsMin;
  float fTPCfindableMin;
  bool fUseSharedMap;
  float fDCArMax;
  float fDCAzMax;
  float fPtMin;
  float fPtMax;
  float fEtaMax;
  float fPIDnSigmaStrict;
  float fPIDMomentumSwitch;
  bool fPIDHypothesisRejection;
  bool fPIDAdditionalRejection;
  bool fPIDTOFif;
  std::vector<AliPID::EParticleType> fParticleAdditionalRejection;
  std::vector<float> fPIDAdditionalRejectionUpperBound;
  std::vector<float> fPIDAdditionalRejectionLowerBound;
  std::vector<float> fPIDAdditionalRejectionUpperpTBound;
  std::vector<float> fPIDAdditionalRejectionLowerpTBound;

  AliPIDResponse *fPIDResponse;  //! pid response

  // Histograms
  // =====================================================================
  TProfile *fHistCuts;  //

  TH1F *fHistSingleParticleCuts;                              //
  TH1F *fHistSingleParticlePhi;                               //
  TH1F *fHistSingleParticlePt;                                //
  TH1F *fHistSingleParticleEta;                               //
  TH2F *fHistSingleParticleEtaPhi;                            //
  TH1F *fHistNSingleParticle;                                 //
  TH2F *fHistSingleParticleDCAxy;                             //
  TH1F *fHistSingleParticleEtaBefore;                         //
  TH1F *fHistSingleParticleEtaAfter;                          //
  TH1F *fHistSingleParticlePtBefore;                          //
  TH1F *fHistSingleParticlePtAfter;                           //
  TH1F *fHistSingleParticleNclsTPCBefore;                     //
  TH1F *fHistSingleParticleNclsTPCAfter;                      //
  TH1F *fHistSingleParticleNclsTPCShared;                     //
  TH2F *fHistSingleParticleNclsTPCSharedTiming;               //
  TH1F *fHistSingleParticleNclsITSShared;                     //
  TH2F *fHistSingleParticleNclsITSSharedTiming;               //
  TH1F *fHistSingleParticleNcrossedTPCBefore;                 //
  TH1F *fHistSingleParticleNcrossedTPCAfter;                  //
  TH1F *fHistSingleParticleFindableTPCBefore;                 //
  TH1F *fHistSingleParticleFindableTPCAfter;                  //
  TH1F *fHistSingleParticleDCAxyBefore;                       //
  TH1F *fHistSingleParticleDCAxyAfterDCAz;                    //
  TH1F *fHistSingleParticleDCAxyAfter;                        //
  TH1F *fHistSingleParticleDCAzBefore;                        //
  TH1F *fHistSingleParticleDCAzAfter;                         //
  TH2F *fHistSingleParticleNsigmaBefore;                      //
  TH2F *fHistSingleParticleNsigmaAfter;                       //
  TH2F *fHistSingleParticleNsigmaTPCBefore;                   //
  TH2F *fHistSingleParticleNsigmaTPCAfter;                    //
  TH2F *fHistSingleParticleNsigmaTOFBefore;                   //
  TH2F *fHistSingleParticleNsigmaTOFAfter;                    //
  TH2F *fHistSingleParticleTPCsignalBefore;                   //
  TH2F *fHistSingleParticleTOFsignalBefore;                   //
  TH2F *fHistSingleParticleTPCsignalAfter;                    //
  TH2F *fHistSingleParticleTOFsignalAfter;                    //
  std::vector<TH2F *> fHistSingleParticleNsigmaTPCRejBefore;  //
  std::vector<TH2F *> fHistSingleParticleNsigmaTPCRejAfter;   //

  TH1F *fHistMCRecSingleParticlePtTruth;   //
  TH2F *fHistMCRecSingleParticleMomentum;  //

  TH1F *fHistSingleAntiParticleCuts;                              //
  TH1F *fHistSingleAntiParticlePhi;                               //
  TH1F *fHistSingleAntiParticlePt;                                //
  TH1F *fHistSingleAntiParticleEta;                               //
  TH2F *fHistSingleAntiParticleEtaPhi;                            //
  TH1F *fHistNSingleAntiParticle;                                 //
  TH2F *fHistSingleAntiParticleDCAxy;                             //
  TH1F *fHistSingleAntiParticleEtaBefore;                         //
  TH1F *fHistSingleAntiParticleEtaAfter;                          //
  TH1F *fHistSingleAntiParticlePtBefore;                          //
  TH1F *fHistSingleAntiParticlePtAfter;                           //
  TH1F *fHistSingleAntiParticleNclsTPCBefore;                     //
  TH1F *fHistSingleAntiParticleNclsTPCAfter;                      //
  TH1F *fHistSingleAntiParticleNclsTPCShared;                     //
  TH2F *fHistSingleAntiParticleNclsTPCSharedTiming;               //
  TH1F *fHistSingleAntiParticleNclsITSShared;                     //
  TH2F *fHistSingleAntiParticleNclsITSSharedTiming;               //
  TH1F *fHistSingleAntiParticleNcrossedTPCBefore;                 //
  TH1F *fHistSingleAntiParticleNcrossedTPCAfter;                  //
  TH1F *fHistSingleAntiParticleFindableTPCBefore;                 //
  TH1F *fHistSingleAntiParticleFindableTPCAfter;                  //
  TH1F *fHistSingleAntiParticleDCAxyBefore;                       //
  TH1F *fHistSingleAntiParticleDCAxyAfterDCAz;                    //
  TH1F *fHistSingleAntiParticleDCAxyAfter;                        //
  TH1F *fHistSingleAntiParticleDCAzBefore;                        //
  TH1F *fHistSingleAntiParticleDCAzAfter;                         //
  TH2F *fHistSingleAntiParticleNsigmaBefore;                      //
  TH2F *fHistSingleAntiParticleNsigmaAfter;                       //
  TH2F *fHistSingleAntiParticleNsigmaTPCBefore;                   //
  TH2F *fHistSingleAntiParticleNsigmaTPCAfter;                    //
  TH2F *fHistSingleAntiParticleNsigmaTOFBefore;                   //
  TH2F *fHistSingleAntiParticleNsigmaTOFAfter;                    //
  TH2F *fHistSingleAntiParticleTPCsignalBefore;                   //
  TH2F *fHistSingleAntiParticleTOFsignalBefore;                   //
  TH2F *fHistSingleAntiParticleTPCsignalAfter;                    //
  TH2F *fHistSingleAntiParticleTOFsignalAfter;                    //
  std::vector<TH2F *> fHistSingleAntiParticleNsigmaTPCRejBefore;  //
  std::vector<TH2F *> fHistSingleAntiParticleNsigmaTPCRejAfter;   //

  TH1F *fHistMCRecSingleAntiParticlePtTruth;   //
  TH2F *fHistMCRecSingleAntiParticleMomentum;  //

 private:
  ClassDef(AliSigma0SingleParticleCuts, 1)
};

#endif
