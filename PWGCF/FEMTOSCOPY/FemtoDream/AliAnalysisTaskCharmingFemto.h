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
#include "AliVertexingHFUtils.h"

class AliVParticle;
class AliVTrack;

class AliAnalysisTaskCharmingFemto : public AliAnalysisTaskSE {
 public:

  enum DecChannel //more HF particles can be added in the future
  {
    kDplustoKpipi,
    kDstartoKpipi
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
  AliAnalysisTaskCharmingFemto(const char *name, const bool isMC, const bool isMCtruth);
  virtual ~AliAnalysisTaskCharmingFemto();

  virtual void LocalInit();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);

  void SetIsMC(bool isMC) {
    fIsMC = isMC;
  }
  void SetUseMCTruthReco(bool useMCTruthReco) {
    fUseMCTruthReco = useMCTruthReco;
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
  void CheckProtonSPDHit() {
    fCheckProtonSPDHit = true;
  }

  // HF related setters
  void SetDecayChannel(int decayChannel=kDplustoKpipi) {
    fDecChannel = decayChannel;
    if (decayChannel == kDplustoKpipi) {
      fDmesonPDGs.push_back(211);
      fDmesonPDGs.push_back(321);
      fDmesonPDGs.push_back(211);
    } else if (decayChannel == kDstartoKpipi) {
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
  void SetCutWindowMCTRUTH(float low, float high) {
    fInvMassCutLow = low;
    fInvMassCutHigh = high;
  }
  void SetBuddypTLowMCTRUTH(float pT) {
    fBuddypTlow = pT;
  }
  void SetBuddypTHighMCTRUTH(float pT) {
    fBuddypThigh = pT;
  }
  void SetBuddyEtaMCTRUTH(float eta) {
    fBuddyeta = eta;
  }
  void SetBuddyOriginMCTRUTH(int origin) {
    //0:no selection, 1:Physical Primary, 2:Secondary From Weak Decay, 3:Secondary From Material, 4: Primary part
    fBuddyOrigin = origin;
  }
  bool SelectBuddyOrigin(AliAODMCParticle *mcPart) {
    if(fBuddyOrigin==0) {
      return true;
    }
    else if(fBuddyOrigin==1){
      if(mcPart->IsPhysicalPrimary())
        return true;
      else
        return false;
    }
    else if(fBuddyOrigin==2){
      if(mcPart->IsSecondaryFromWeakDecay())
        return true;
      else 
        return false;
    }
    else if(fBuddyOrigin==3){
      if(mcPart->IsSecondaryFromMaterial())
        return true;
      else 
        return false;
    }
    else if(fBuddyOrigin==4){
      if(mcPart->IsPrimary())
        return true;
      else 
        return false;
    }
    return false;
  }
  void SetDmesonOriginMCTRUTH(int origin) {
    //0:no selection, 1:charm, 2:beauty, 3: D*, 4:B
    fDmesonOrigin = origin;
  }
  bool SelectDmesonOrigin(TClonesArray* arrayMC, AliAODMCParticle *mcPart) {
    if(fDmesonOrigin==0) {
      fDoDorigPlots=true;
      return true;
    }
    else {
      int origin = AliVertexingHFUtils::CheckOrigin(arrayMC, mcPart, true);
      if(fDmesonOrigin==1){
        if(origin==4)
          return true;
        else
          return false;
      }
      else if(fDmesonOrigin==2){
        if(origin==5)
          return true;
        else
          return false;
      }
      else {
        int motherID = mcPart->GetMother();
        AliAODMCParticle *mcMother = (AliAODMCParticle *)arrayMC->At(motherID);
        int PDGCodeMoth= mcMother->GetPdgCode();
        if(origin==4 && fDmesonOrigin==3){
          if(TMath::Abs(PDGCodeMoth) == 413)
            return true;
          else 
            return false;
        }
        if(fDmesonOrigin==4){
          if(TMath::Abs(PDGCodeMoth) == 521)
            return true;
          else 
            return false;
        }
      }
      return false;
    }
  }
  void FillMCtruthPDGDmeson(TClonesArray* arrayMC, AliAODMCParticle *mcPart) {
    int motherID = mcPart->GetMother();
    AliAODMCParticle *mcMother = (AliAODMCParticle *)arrayMC->At(motherID);
    int PDGCodeMoth= mcMother->GetPdgCode();
    int PDGCodePart= mcPart->GetPdgCode();
    float pt = mcPart->Pt();
    if(PDGCodePart==411){
      fHistDplusMCtruthmotherPDG->Fill(pt,TMath::Abs(PDGCodeMoth));
    } else if (PDGCodePart==-411){
      fHistDminusMCtruthmotherPDG->Fill(pt,TMath::Abs(PDGCodeMoth));
    }
    return;
  }
  
  void FillMCtruthQuarkOriginDmeson(TClonesArray* arrayMC, AliAODMCParticle *mcPart) {
    int PDGCodePart= mcPart->GetPdgCode();
    float pt = mcPart->Pt();
    if (AliVertexingHFUtils::CheckOrigin(arrayMC, mcPart, true)==4) { //charm
      if(PDGCodePart==411){
        fHistDplusMCtruthQuarkOrigin->Fill(pt,1);
      } else if (PDGCodePart==-411){
        fHistDminusMCtruthQuarkOrigin->Fill(pt,1);
      }
    }
    else if (AliVertexingHFUtils::CheckOrigin(arrayMC, mcPart, true)==5) { //beauty
      if(PDGCodePart==411){
        fHistDplusMCtruthQuarkOrigin->Fill(pt,2);
      } else if (PDGCodePart==-411){
        fHistDminusMCtruthQuarkOrigin->Fill(pt,2);
      }
    }
    else {
      if(PDGCodePart==411){
        fHistDplusMCtruthQuarkOrigin->Fill(pt,3);
      } else if (PDGCodePart==-411){
        fHistDminusMCtruthQuarkOrigin->Fill(pt,3);
      }
    }
    return;
  }

 private:
  AliAnalysisTaskCharmingFemto(const AliAnalysisTaskCharmingFemto &task);
  AliAnalysisTaskCharmingFemto &operator=(
      const AliAnalysisTaskCharmingFemto &task);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  int IsCandidateSelected(AliAODRecoDecayHF *&dMeson, AliAODRecoDecayHF *&dMesonWithVtx, int absPdgMom, bool &unsetVtx, bool &recVtx, AliAODVertex *&origOwnVtx, std::vector<double> scores);
  bool MassSelection(const double mass, const double pt, const int pdg);

  // Track / event selection objects
  AliAODEvent *fInputEvent;                          //
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
  bool fUseMCTruthReco;    //
  bool fIsMCtruth;         //
  bool fIsLightweight;     //
  UInt_t fTrigger;         //
  int fSystem;             //

  bool fCheckProtonSPDHit; //

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

  bool fDoDorigPlots;            //!
  TH2F *fHistDplusMCtruthmotherPDG;  //!
  TH2F *fHistDplusMCtruthQuarkOrigin;     //!
  TH2F *fHistDminusMCtruthmotherPDG;  //!
  TH2F *fHistDminusMCtruthQuarkOrigin;     //!

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

  //MC Truth Stuff
  float fInvMassCutLow;    //
  float fInvMassCutHigh;   //
  float fBuddypTlow;       //
  float fBuddypThigh;      //
  float fBuddyeta;         //
  int fBuddyOrigin;        //
  int fDmesonOrigin;       //

  // variables for ML application
  bool fApplyML;                                           // flag to enable ML application
  TString fConfigPath;                                     // path to ML config file
  AliHFMLResponse* fMLResponse;                            //!<! object to handle ML response

  bool fDependOnMLSelector;                                // flag to read ML scores from a AliAnalysisTaskSECharmHadronMLSelector task
  std::vector<float> fPtLimsML;                            // pT bins in case application of ML model is done in MLSelector task   
  std::vector<std::vector<double> > fMLScoreCuts;          // score cuts used in case application of ML model is done in MLSelector task   
  std::vector<std::vector<std::string> > fMLOptScoreCuts;  // score cut options (lower, upper) used in case application of ML model is done in MLSelector task   

ClassDef(AliAnalysisTaskCharmingFemto, 13)
};

#endif
