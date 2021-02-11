#ifndef AliAnalysisTaskCharmingFemto_H
#define AliAnalysisTaskCharmingFemto_H

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliRDHFCuts.h"
#include "AliAODVertex.h"
#include "AliHFMLResponse.h"
#include "TChain.h"

class AliVParticle;
class AliVTrack;

class AliAnalysisTaskCharmingFemto : public AliAnalysisTaskSE {
 public:

  enum DecChannel //more HF particles can be added in the future
  {
    kDplustoKpipi
  };

  enum CollSystem
  {
    kpp5TeV,
    kpp13TeV
  };

  enum MassSelection
  {
    kSignal,
    kSidebandRight,
    kSidebandLeft,
    kStrictCut,
  };

  AliAnalysisTaskCharmingFemto();
  AliAnalysisTaskCharmingFemto(const char *name, const bool isMC);
  virtual ~AliAnalysisTaskCharmingFemto();

  virtual void LocalInit();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  void SetIsMC(bool isMC) {
    fIsMC = isMC;
  }
  void SetLightweight(bool isLightweight) {
    fIsLightweight = isLightweight;
  }
  void SetTrigger(UInt_t trigger) {
    fTrigger = trigger;
  }
  void SetEventCuts(AliFemtoDreamEventCuts *cuts) {
    fEvtCuts = cuts;
  }
  void SetProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCutsPartProton = cuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts *cuts) {
    fTrackCutsPartAntiProton = cuts;
  }
  void SetCollectionConfig(AliFemtoDreamCollConfig *config) {
    fConfig = config;
  }
  void SetSystem(int system) {
    fSystem = system;
  }

  // HF related setters
  void SetDecayChannel(int decayChannel=kDplustoKpipi) {
    fDecChannel = decayChannel;
    if (decayChannel == kDplustoKpipi) {
      fDmesonPDGs.push_back(211);
      fDmesonPDGs.push_back(321);
      fDmesonPDGs.push_back(211);
    } else {
      AliFatal("Decay channel not implemented!");
    }
  }
  void SetHFCuts(AliRDHFCuts* cuts) {
    fRDHFCuts = cuts;
  }
  void SetAODMismatchProtection(int opt=0) {
    fAODProtection = opt;
  }
  void SetDoMLApplication(bool flag = true) {
    fApplyML = flag;
  }
  void SetMLConfigFile(TString path = "") {
    fConfigPath = path;
  }
  void SetMassSelection(int type) {
    fMassSelectionType = type;
  }
  void SetNSigmaSelection(double nSigma = 2) {
    fMassSelectionType = kSignal;
    fNSigmaMass = nSigma;
  }
  void SetSidebandChoice(MassSelection type, double offset, double width) {
    fMassSelectionType = type;
    fNSigmaOffsetSideband = offset;
    fSidebandWidth = width;
  }
  void SetDstarWindow(double lower, double upper) {
    fLowerDstarRemoval = lower;
    fUpperDstarRemoval = upper;
  }
  void SetMassWindow(double lower, double upper) {
    fMassSelectionType = kStrictCut;
    fLowerMassSelection = lower;
    fUpperMassSelection = upper;
  }
  void SetIsDependentOnMLSelector(bool flag=true) {
    fDependOnMLSelector = flag;
  }
  void ScaleMCBeautyFraction(double pythiaBeautyFraction,
                             double desiredFraction) {
    AliInfo("Scaling the beauty fraction in MC activated");
    if (desiredFraction > pythiaBeautyFraction) {
      AliFatal(
          "The scaled fraction cannot be larger than the initial fraction");
    }
    fMCBeautyRejection = true;
    fMCBeautyScalingFactor = (1. - pythiaBeautyFraction)
        / (1. - desiredFraction) * desiredFraction / pythiaBeautyFraction;
    AliInfo(Form("Assumed old fraction: %.3f", pythiaBeautyFraction));
    AliInfo(Form("Desired new fraction: %.3f", desiredFraction));
    AliInfo(Form("Scaling factor: %.3f", fMCBeautyScalingFactor));
  }
  void UseTrueDOnly() {
    fUseTrueDOnly = true;
  }

 private:
  AliAnalysisTaskCharmingFemto(const AliAnalysisTaskCharmingFemto &task);
  AliAnalysisTaskCharmingFemto &operator=(
      const AliAnalysisTaskCharmingFemto &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  int IsCandidateSelected(AliAODRecoDecayHF *&dMeson, int absPdgMom, bool &unsetVtx, bool &recVtx, AliAODVertex *&origOwnVtx, std::vector<double> scores);
  bool MassSelection(const double mass, const double pt, const int pdg);

  // Track / event selection objects
  AliAODEvent* fInputEvent;                          //
  AliFemtoDreamEvent *fEvent;                        //!
  AliFemtoDreamEventCuts *fEvtCuts;                  //
  AliFemtoDreamTrack *fProtonTrack;                  //!
  AliFemtoDreamTrackCuts *fTrackCutsPartProton;      //
  AliFemtoDreamTrackCuts *fTrackCutsPartAntiProton;  //

  // Femto classes
  AliFemtoDreamCollConfig *fConfig;                  //
  AliFemtoDreamPairCleaner *fPairCleaner;            //!
  AliFemtoDreamPartCollection *fPartColl;            //!

  bool fIsMC;              //
  bool fIsLightweight;     //
  UInt_t fTrigger;         //
  int fSystem;             //

  int fTrackBufferSize;
  std::vector<unsigned int> fDmesonPDGs;
  AliAODTrack **fGTI;  //!

  TList *fQA;                      //!
  TList *fEvtHistList;             //!
  TList *fTrackCutHistList;        //!
  TList *fTrackCutHistMCList;      //!
  TList *fAntiTrackCutHistList;    //!
  TList *fAntiTrackCutHistMCList;  //!
  TList *fDChargedHistList;		     //!
  TList *fResultList;              //!
  TList *fResultQAList;            //!

  TH2F *fHistDplusInvMassPt;   //!
  TH2F *fHistDplusInvMassPtSel;   //!
  TH1F *fHistDplusEta;         //!
  TH1F *fHistDplusPhi;         //!
  TH1F *fHistDplusChildPt[5];  //!
  TH1F *fHistDplusChildEta[5]; //!
  TH1F *fHistDplusChildPhi[5]; //!
  TH2F *fHistDplusMCPDGPt;     //!
  TH2F *fHistDplusMCPtRes;     //!
  TH2F *fHistDplusMCPhiRes;    //!
  TH2F *fHistDplusMCThetaRes;  //!
  TH2F *fHistDplusMCOrigin;    //!

  TH2F *fHistDminusInvMassPt;   //!
  TH2F *fHistDminusInvMassPtSel;   //!
  TH1F *fHistDminusEta;         //!
  TH1F *fHistDminusPhi;         //!
  TH1F *fHistDminusChildPt[5];  //!
  TH1F *fHistDminusChildEta[5]; //!
  TH1F *fHistDminusChildPhi[5]; //!
  TH2F *fHistDminusMCPDGPt;     //!
  TH2F *fHistDminusMCPtRes;     //!
  TH2F *fHistDminusMCPhiRes;    //!
  TH2F *fHistDminusMCThetaRes;  //!
  TH2F *fHistDminusMCOrigin;    //!
  
  // HF data members
  int fDecChannel;                                         // HF decay channel
  AliRDHFCuts* fRDHFCuts;                                  // HF cut object
  int fAODProtection;                                      // flag to activate protection against AOD-dAOD mismatch.
                                                           // -1: no protection,  0: check AOD/dAOD nEvents only,  1: check AOD/dAOD nEvents + TProcessID names
  int fMassSelectionType;			           // Switch for the D meson inv. mass selection type
  double fNSigmaMass;					                             // Width of the mass window
  double fNSigmaOffsetSideband;                            // Offset of the mass window from the D inv. mass peak
  double fLowerMassSelection;			                         // Lower boundary of the mass selection
  double fUpperMassSelection;			                         // Upper boundary of the mass selection
  double fSidebandWidth;                                   // Width of the sideband
  double fLowerDstarRemoval;                               // Lower boundary to remove the D*
  double fUpperDstarRemoval;                               // Upper boundary to remove the D*

  bool fMCBeautyRejection;                                 // Switch for scaling the beauty feed-down fraction in MC
  double fMCBeautyScalingFactor;                           // Factor for scaling the beauty feed-down
  bool fUseTrueDOnly;

  // variables for ML application
  bool fApplyML;                                           // flag to enable ML application
  TString fConfigPath;                                     // path to ML config file
  AliHFMLResponse* fMLResponse;                            //!<! object to handle ML response

  bool fDependOnMLSelector;                                // flag to read ML scores from a AliAnalysisTaskSECharmHadronMLSelector task
  std::vector<float> fPtLimsML;                            // pT bins in case application of ML model is done in MLSelector task   
  std::vector<std::vector<double> > fMLScoreCuts;          // score cuts used in case application of ML model is done in MLSelector task   
  std::vector<std::vector<std::string> > fMLOptScoreCuts;  // score cut options (lower, upper) used in case application of ML model is done in MLSelector task   

ClassDef(AliAnalysisTaskCharmingFemto, 11)
};

#endif
