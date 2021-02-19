#ifndef ALIANLYSISTASKGAMMACALO_cxx
#define ALIANLYSISTASKGAMMACALO_cxx

#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliV0ReaderV1.h"
#include "AliKFConversionPhoton.h"
#include "AliGammaConversionAODBGHandler.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliCaloPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "AliConversionPhotonCuts.h"
#include "AliConversionMesonCuts.h"
#include "AliAnalysisTaskConvJet.h"
#include "AliAnalysisTaskJetOutlierRemoval.h"
#include "AliAnalysisManager.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"
#include "TGenPhaseSpace.h"
#include <vector>
#include <map>

// simple struct used for merging studies
struct clusterLabel{
  clusterLabel(): mesonID(0), clusID(0), daughterID(0), daughterPDG(0),
                  EClus(0), EFrac(0), ETrue(0), PtMeson(0), EtaMeson(0), OpeningAngle(0),
                  clusVec()
                  {};
  Int_t mesonID,clusID,daughterID,daughterPDG;
  Float_t EClus,EFrac,ETrue,PtMeson,EtaMeson, OpeningAngle;
  TLorentzVector clusVec;
};

class AliAnalysisTaskGammaCalo : public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskGammaCalo();
    AliAnalysisTaskGammaCalo(const char *name);
    virtual ~AliAnalysisTaskGammaCalo();

    virtual void   UserCreateOutputObjects();
    virtual Bool_t Notify();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);
    void InitBack();

    void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
    void SetCaloTriggerHelperName(TString name){fCaloTriggerHelperName=name; return;}
    void SetIsHeavyIon(Int_t flag){
      fIsHeavyIon = flag;
    }

    // base functions for selecting photon and meson candidates in reconstructed data
    void ProcessClusters();
    void ProcessConversionCandidates();
    void ProcessJets();
    void CalculatePi0Candidates();

    // MC functions
    void SetIsMC(Int_t isMC){fIsMC=isMC;}
    void ProcessMCParticles();
    void ProcessAODMCParticles();
    void ProcessTrueClusterCandidates( AliAODConversionPhoton* TruePhotonCandidate, Double_t clusterM02);
    void ProcessTrueClusterCandidatesAOD( AliAODConversionPhoton* TruePhotonCandidate, Double_t clusterM02);
    void ProcessTrueMesonCandidates( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
    void ProcessTrueMesonCandidatesAOD(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);
    void ProcessAODSphericityParticles();

    // switches for additional analysis streams or outputs
    void SetLightOutput(Bool_t flag){fDoLightOutput = flag;}
    void SetECalibOutput( Bool_t flag ){fDoECalibOutput = flag;}
    void SetDoMesonAnalysis(Bool_t flag){fDoMesonAnalysis = flag;}
    void SetDoMesonQA(Int_t flag){fDoMesonQA = flag;}
    void SetDoClusterQA(Int_t flag){fDoClusterQA = flag;}
    void SetDoTHnSparse(Bool_t flag){fDoTHnSparse = flag;}
    void SetPlotHistsExtQA(Bool_t flag){fSetPlotHistsExtQA = flag;}
    void SetAllowOverlapHeaders( Bool_t allowOverlapHeader ) {fAllowOverlapHeaders = allowOverlapHeader;}
    void SetDoPi0Only(Bool_t flag){fDoPi0Only = flag;}

    void SetInOutTimingCluster(Double_t min, Double_t max){
      fDoInOutTimingCluster = kTRUE; fMinTimingCluster = min; fMaxTimingCluster = max;
      return;
    }

      // Setting the cut lists for the conversion photons
    void SetEventCutList(Int_t nCuts, TList *CutArray){
      fnCuts = nCuts;
      fEventCutArray = CutArray;
    }

      // Setting the cut lists for the calo photons
    void SetCaloCutList(Int_t nCuts, TList *CutArray){
      fnCuts = nCuts;
      fClusterCutArray = CutArray;
    }

    // Setting the cut lists for the meson
    void SetMesonCutList(Int_t nCuts, TList *CutArray){
      fnCuts = nCuts;
      fMesonCutArray = CutArray;
    }

    // BG HandlerSettings
    void CalculateBackground();
    void CalculateBackgroundSwapp();
    void CalculateBackgroundRP();
    void RotateParticle(AliAODConversionPhoton *gamma);
    void RotateParticleAccordingToEP(AliAODConversionPhoton *gamma, Double_t previousEventEP, Double_t thisEventEP);
    void FillPhotonBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode);
    void FillPhotonPlusConversionBackgroundHist(AliAODConversionPhoton *TruePhotonCandidate, Int_t pdgCode);
    void FillPhotonBackgroundM02Hist(AliAODConversionPhoton *TruePhotonCandidate, Double_t clusterM02, Int_t pdgCode);
    void FillPhotonPlusConversionBackgroundM02Hist(AliAODConversionPhoton *TruePhotonCandidate, Double_t  clusterM02, Int_t pdgCode);
    void UpdateEventByEventData();

    // Additional functions for convenience
    void SetLogBinningXTH2(TH2* histoRebin);
    Int_t GetSourceClassification(Int_t daughter, Int_t pdgCode);

    Bool_t CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked);
    void FillMultipleCountMap(map<Int_t,Int_t> &ma, Int_t tobechecked);
    void FillMultipleCountHistoAndClear(map<Int_t,Int_t> &ma, TH1F* hist);
    Int_t WhichDDL(Int_t module=0, Int_t cellx=0);

    // set method to enable EOverP tree
    void SetProduceTreeEOverP(Bool_t b){fProduceTreeEOverP = b;}

    // Function to enable MC label sorting
    void SetEnableSortingOfMCClusLabels(Bool_t enableSort)  { fEnableSortForClusMC   = enableSort; }

    // Function to enable local debugging mode
    void SetLocalDebugFlag(Int_t iF) {fLocalDebugFlag = iF;}

    // Function to set correction task setting
    void SetCorrectionTaskSetting(TString setting) {fCorrTaskSetting = setting;}

    void EventDebugMethod();
    void DebugMethod(AliAODConversionMother *pi0cand, AliAODConversionPhoton *gamma0, AliAODConversionPhoton *gamma1);
    void DebugMethodPrint1(AliAODConversionMother *pi0cand, AliAODConversionPhoton *gamma0, AliAODConversionPhoton *gamma1);

    void SetTrackMatcherRunningMode(Int_t mode){fTrackMatcherRunningMode = mode;}

    void SetSoftAnalysis(Bool_t DoSoft)  {fDoSoftAnalysis = DoSoft;}

    Int_t CheckClustersForMCContribution(Int_t mclabel, Bool_t leading = kFALSE);
    Bool_t CheckSpecificClusterForMCContribution(Int_t mclabel, Int_t cluslabel);
    Int_t CountPhotonsInCluster(Int_t cluslabel);

    void DoClusterMergingStudies(AliVCluster* clus, vector<clusterLabel> &labelvect);
    void DoClusterMergingStudiesAOD(AliVCluster* clus, vector<clusterLabel> &labelvect);

  protected:
    AliV0ReaderV1*        fV0Reader;                                            // basic photon Selection Task
    TString               fV0ReaderName;
    TString               fCaloTriggerHelperName;
    TString               fCorrTaskSetting;
    AliGammaConversionAODBGHandler**  fBGHandler;                               // BG handler for Conversion
    AliVEvent*            fInputEvent;                                          // current event
    AliMCEvent*           fMCEvent;                                             // corresponding MC event
    TList**               fCutFolder;                                           // Array of lists for containers belonging to cut
    TList**               fESDList;                                             // Array of lists with histograms with reconstructed properties
    TList**               fBackList;                                            // Array of lists with BG THnSparseF
    TList**               fMotherList;                                          // Array of lists with Signal THnSparseF
    TList**               fTrueList;                                            // Array of lists with histograms with MC validated reconstructed properties
    TList**               fMCList;                                              // Array of lists with histograms with pure MC information
    TList**               fTreeList;                                            // Array of lists with tree for validated MC
    TList**               fClusterTreeList;                                     // Array of lists with tree for EoverP
    TList*                fOutputContainer;                                     // Output container
    TClonesArray*         fReaderGammas;                                        // Array with conversion photons selected by V0Reader Cut
    TList*                fGammaCandidates;                                     // current list of photon candidates
    TList*                fClusterCandidates;                                   //! current list of cluster candidates
    TList*                fEventCutArray;                                       // List with Event Cuts
    AliConvEventCuts*     fEventCuts;                                           // EventCutObject
    TList*                fClusterCutArray;                                     // List with Cluster Cuts
    AliCaloPhotonCuts*    fCaloPhotonCuts;                                      // CaloPhotonCutObject
    TList*                fMesonCutArray;                                       // List with Meson Cuts
    AliConversionMesonCuts*   fMesonCuts;                                       // MesonCutObject
    AliAnalysisTaskConvJet*   fConvJetReader;                                   // JetReader
    AliAnalysisTaskJetOutlierRemoval*   fOutlierJetReader;                      // JetReader
    AliConversionPhotonCuts*  fConversionCuts;                                  // ConversionPhotonCutObject
    Bool_t                fDoJetAnalysis;                                       // Bool to produce Jet Plots
    Bool_t                fDoJetQA;                                             // Bool to produce Jet QA Plots
    Bool_t                fDoTrueSphericity;                                    // Bool to produce Sphericity correlations
    TList**               fJetHistograms;                                       // Jet Histograms
    TList**               fTrueJetHistograms;                                   // True Jet Histograms
    Int_t                 fJetSector;                                           // Sector of the detector with the maximum pt jet
    Int_t                 fMaxPtNearEMCalPlace;                                 // Place in jet vector of highest pt jet that is near the EMCal
    Bool_t                fJetNearEMCal;                                        // If a jet is near the EMCal in the current event
    Int_t*                fDDLRange_HistoClusGamma;                          //! Min and Max Value for Modules in PHOS for fHistoClusGamma E and Pt
    AliCaloTriggerMimicHelper**     fCaloTriggerMimicHelper;                    //!Array wich points to AliCaloTriggerMimicHelper for each Event Cut
    map<TString, Bool_t>  fSetEventCutsOutputlist;                              //! Store, if Output list for Event Cut has already been added

    //histograms for mesons reconstructed quantities
    TH2F**                fHistoMotherInvMassPt;                                //! array of histogram with signal + BG for same event photon pairs, inv Mass, pt
    THnSparseF**          fSparseMotherInvMassPtZM;                             //! array of THnSparseF with signal + BG for same event photon pairs, inv Mass, pt
    TH2F**                fHistoMotherBackInvMassPt;                            //! array of histogram with BG for mixed event photon pairs, inv Mass, pt
    THnSparseF**          fSparseMotherBackInvMassPtZM;                         //! array of THnSparseF with BG for same event photon pairs, inv Mass, pt
    TH2F**                fHistoMotherPi0PtY;                                   //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, Y
    TH2F**                fHistoMotherEtaPtY;                                   //! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65, pt, Y
    TH2F**                fHistoMotherPi0PtAlpha;                               //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, alpha
    TH2F**                fHistoMotherEtaPtAlpha;                               //! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65, pt, alpha
    TH2F**                fHistoMotherPi0PtOpenAngle;                           //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, openAngle
    TH2F**                fHistoMotherEtaPtOpenAngle;                           //! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65, pt, openAngle
    TH2F**                fHistoMotherPtOpenAngle;                              //! array of histograms pt, openAngle
    TH2F**                fHistoMotherPtOpenAngleBck;                           //! array of histograms pt, openAngle - for Background from eventmixing
    TH2F**                fHistoMotherPi0NGoodESDTracksPt;                 //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, ngoodtrakcs, pt
    TH2F**                fHistoMotherEtaNGoodESDTracksPt;                 //! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65, ngoodtrakcs, pt
    TH2F**		            fHistoMotherInvMassECalib;				//! array of histogram with alpha cut of 0.1 for inv mass vs energy of cluster
    TH2F**		            fHistoMotherBackInvMassECalib;			//! array of histogram with BG for mixed event photon pairs with alpha cut of 0.1, inv mass, energy of cluster

    // histograms for rec photon clusters
    TH1F**                fHistoClusGammaPt;                                    //! array of histos with cluster, pt
    TH1F**                fHistoClusGammaE;                                     //! array of histos with cluster, E
    TH1F**                fHistoClusGammaE_BothBM;                              //! array of histos with cluster, E; Only clusters, which are good on triggered bad map and analysis bad map; Only MB
    TH1F**                fHistoClusGammaE_BothBM_highestE;                     //! array of histos with cluster, E; Only highest cluster in event, which is good on triggered bad map and analysis bad map;  MB and tigger but use only triggered clusters for trigger
    TH1F**                fHistoClusGammaE_AnaBM_highestE;                      //! array of histos with cluster, E; Only highest cluster in event, which is good on analysis bad map; Only MB
    TH1F**                fHistoClusGammaE_onlyTriggered;                       //! array of histos with cluster, E
    TH1F***               fHistoClusGammaE_DDL;                                 //! array of histos with cluster, E
    TH1F***               fHistoClusGammaE_DDL_TrBM;                            //! array of histos with cluster, E, only clusters on Triggered Bad map
    TH1F***               fHistoClusGammaE_DDL_wL0_TrigEv;                      //! array of histos with cluster, E, MB histogram, for events with L0 flag with triggered clusters
    TH1F***               fHistoClusGammaE_DDL_woL0_TrigEv;                     //! array of histos with cluster, E, MB histogram, for events with no L0 flag with triggered clusters
    TH1F***               fHistoClusGammaE_DDL_wL0_TrigEv_TrBM;                 //! array of histos with cluster, E, only clusters on Triggered Bad map
    TH1F***               fHistoClusGammaE_DDL_woL0_TrigEv_TrBM;                //! array of histos with cluster, E, only clusters on Triggered Bad map
    TH1I**                fHistoGoodMesonClusters;                                //! Histograms which stores if Pi0 Clusters Trigger
    TH1F**                fHistoClusOverlapHeadersGammaPt;                      //! array of histos with cluster, pt overlapping with other headers
    TH1F**                fHistoClusAllHeadersGammaPt;                          //! array of histos with cluster, pt all headers
    TH1F**                fHistoClusRejectedHeadersGammaPt;                     //! array of histos with cluster, pt rejected with other headers
    TH2F**                fHistoClusGammaPtM02;                                 //! array of histos with cluster M02 vs. pt
    //histograms for pure MC quantities
    TH1I**                fHistoMCHeaders;                                      //! array of histos for header names
    TH1F**                fHistoMCAllGammaPt;                                   //! array of histos with all gamma, pT
    TH2F**                fHistoMCAllSecondaryGammaPt;                          //! array of histos with all secondary gamma, pT
    TH1F**                fHistoMCDecayGammaPi0Pt;                              //! array of histos with decay gamma from pi0, pT
    TH1F**                fHistoMCDecayGammaRhoPt;                              //! array of histos with decay gamma from rho, pT
    TH1F**                fHistoMCDecayGammaEtaPt;                              //! array of histos with decay gamma from eta, pT
    TH1F**                fHistoMCDecayGammaOmegaPt;                            //! array of histos with decay gamma from omega, pT
    TH1F**                fHistoMCDecayGammaEtapPt;                             //! array of histos with decay gamma from eta', pT
    TH1F**                fHistoMCDecayGammaPhiPt;                              //! array of histos with decay gamma from phi, pT
    TH1F**                fHistoMCDecayGammaSigmaPt;                            //! array of histos with decay gamma from Sigma0, pT
    TH1F**                fHistoMCPi0Pt;                                        //! array of histos with weighted pi0, pT
    TH1F**                fHistoMCPi0WOWeightPt;                                //! array of histos with unweighted pi0, pT
    TH1F**                fHistoMCPi0WOEvtWeightPt;                             //! array of histos without event weights pi0, pT
    TH1F**                fHistoMCEtaPt;                                        //! array of histos with weighted eta, pT
    TH1F**                fHistoMCEtaWOWeightPt;                                //! array of histos with unweighted eta, pT
    TH1F**                fHistoMCEtaWOEvtWeightPt;                             //! array of histos without event weights eta, pT
    TH1F**                fHistoMCPi0InAccPt;                                   //! array of histos with weighted pi0 in acceptance, pT
    TH1F**                fHistoMCEtaInAccPt;                                   //! array of histos with weighted eta in acceptance, pT
    TH1F**                fHistoMCPi0WOEvtWeightInAccPt;                        //! array of histos without evt weight pi0 in acceptance, pT
    TH1F**                fHistoMCEtaWOEvtWeightInAccPt;                        //! array of histos without evt weight eta in acceptance, pT
    TH2F**                fHistoMCPi0PtY;                                       //! array of histos with weighted pi0, pT, Y
    TH2F**                fHistoMCEtaPtY;                                       //! array of histos with weighted eta, pT, Y
    TH2F**                fHistoMCPi0PtAlpha;                                   //! array of histos with weighted pi0, pT, alpha
    TH2F**                fHistoMCEtaPtAlpha;                                   //! array of histos with weighted eta, pT, alpha
    TH2F**                fHistoMCPrimaryPtvsSource;                            //! array of histos with weighted primary particles, pT vs source
    TH2F**                fHistoMCSecPi0PtvsSource;                             //! array of histos with secondary pi0, pT, source
    TH1F**                fHistoMCSecPi0Source;                                 //! array of histos with secondary pi0, source
    TH2F**                fHistoMCSecPi0InAccPtvsSource;                        //! array of histos with secondary pi0 in acceptance, pT, source
    TH1F**                fHistoMCSecEtaPt;                                     //! array of histos with secondary eta, pT
    TH1F**                fHistoMCSecEtaSource;                                 //! array of histos with secondary eta, source
    TH2F**                fHistoMCPi0PtJetPt;                                   //! array of histos with weighted pi0, pT, hardest jet pt
    TH2F**                fHistoMCEtaPtJetPt;                                   //! array of histos with weighted eta, pT, hardest jet pt

    // MC validated reconstructed quantities mesons
    TH2F**                fHistoTruePi0InvMassPt;                               //! array of histos with validated mothers, invMass, pt
    TH2F**                fHistoTruePi0InvMassPtAdditional;                               //! array of histos with validated mothers, invMass, pt
    TH2F**                fHistoTrueEtaInvMassPt;                               //! array of histos with validated mothers, invMass, pt
    TH2F**                fHistoTrueEtaInvMassPtAdditional;                               //! array of histos with validated mothers, invMass, pt
    TH2F**                fHistoTruePi0CaloPhotonInvMassPt;                     //! array of histos with validated mothers, photon leading, invMass, pt
    TH2F**                fHistoTrueEtaCaloPhotonInvMassPt;                     //! array of histos with validated mothers, photon leading, invMass, pt
    TH2F**                fHistoTruePi0CaloConvertedPhotonInvMassPt;            //! array of histos with validated pi0, converted photon leading, invMass, pt
    TH2F**                fHistoTrueEtaCaloConvertedPhotonInvMassPt;            //! array of histos with validated eta, converted photon leading, invMass, pt
    TH2F**                fHistoTruePi0CaloMixedPhotonConvPhotonInvMassPt;      //! array of histos with validated mothers, converted photon leading, invMass, pt
    TH2F**                fHistoTrueEtaCaloMixedPhotonConvPhotonInvMassPt;      //! array of histos with validated mothers, converted photon leading, invMass, pt
    TH2F**                fHistoTruePi0CaloElectronInvMassPt;                   //! array of histos with validated mothers, electron leading, invMass, pt
    TH2F**                fHistoTrueEtaCaloElectronInvMassPt;                   //! array of histos with validated mothers, electron leading, invMass, pt
    TH2F**                fHistoTruePi0CaloMergedClusterInvMassPt;              //! array of histos with validated mothers, merged cluster invMass, pt
    TH2F**                fHistoTrueEtaCaloMergedClusterInvMassPt;              //! array of histos with validated mothers, merged cluster invMass, pt
    TH2F**                fHistoTruePi0CaloMergedClusterPartConvInvMassPt;      //! array of histos with validated mothers, merged cluster part conv, invMass, pt
    TH2F**                fHistoTrueEtaCaloMergedClusterPartConvInvMassPt;      //! array of histos with validated mothers, merged cluster part conv, invMass, pt
    TH2F**                fHistoTruePi0NonMergedElectronPhotonInvMassPt;        //! array of histos with validated mothers, merged cluster invMass, pt
    TH2F**                fHistoTruePi0NonMergedElectronMergedPhotonInvMassPt;  //! array of histos with validated mothers, merged cluster invMass, pt
    TH2F**                fHistoTruePi0Category1;                               //! array of histos with validated pi0, pure real photons
    TH2F**                fHistoTrueEtaCategory1;                               //! array of histos with validated eta, pure real photons
    TH2F**                fHistoTruePi0Category2;                               //! array of histos with validated pi0, 1 real photon, 1 merged converted photon
    TH2F**                fHistoTrueEtaCategory2;                               //! array of histos with validated eta, 1 real photon, 1 merged converted photon
    TH2F**                fHistoTruePi0Category3;                               //! array of histos with validated pi0, 1 real photon, 1 electron from conversion unmerged
    TH2F**                fHistoTrueEtaCategory3;                               //! array of histos with validated eta, 1 real photon, 1 electron from conversion unmerged
    TH2F**                fHistoTruePi0Category4_6;                             //! array of histos with validated pi0, 2 electrons from same conversion
    TH2F**                fHistoTrueEtaCategory4_6;                             //! array of histos with validated eta, 2 electrons from same conversion
    TH2F**                fHistoTruePi0Category5;                               //! array of histos with validated pi0, 2 electrons from different conversions, 2 electrons (unseen)
    TH2F**                fHistoTrueEtaCategory5;                               //! array of histos with validated eta, 2 electrons from different conversions, 2 electrons (unseen)
    TH2F**                fHistoTruePi0Category7;                               //! array of histos with validated pi0, 1 electron from conversion, 2 electrons from other conversion merged, 1 electron (unseen)
    TH2F**                fHistoTrueEtaCategory7;                               //! array of histos with validated eta, 1 electron from conversion, 2 electrons from other conversion merged, 1 electron (unseen)
    TH2F**                fHistoTruePi0Category8;                               //! array of histos with validated pi0, 2 electron from conversion merged, 2 electrons from other conversion merged
    TH2F**                fHistoTrueEtaCategory8;                               //! array of histos with validated eta, 2 electron from conversion merged, 2 electrons from other conversion merged
    TH2F**                fHistoTruePrimaryPi0InvMassPt;                        //! array of histos with validated weighted primary mothers, invMass, pt
    TH2F**                fHistoTruePrimaryEtaInvMassPt;                        //! array of histos with validated weighted primary mothers, invMass, pt
    TH2F**                fHistoTruePrimaryPi0W0WeightingInvMassPt;             //! array of histos with validated unweighted primary mothers, invMass, pt
    TH2F**                fHistoTruePrimaryEtaW0WeightingInvMassPt;             //! array of histos with validated unweighted primary mothers, invMass, pt
    TProfile2D**          fProfileTruePrimaryPi0WeightsInvMassPt;               //! array of profiles with weights for validated primary mothers, invMass, pt
    TProfile2D**          fProfileTruePrimaryEtaWeightsInvMassPt;               //! array of profiles with weights for validated primary mothers, invMass, pt
    TH2F**                fHistoTruePrimaryPi0MCPtResolPt;                      //! array of histos with validated weighted primary pi0, MCpt, resol pt
    TH2F**                fHistoTruePrimaryEtaMCPtResolPt;                      //! array of histos with validated weighted primary eta, MCpt, resol pt
    TH2F**                fHistoTrueSecondaryPi0InvMassPt;                      //! array of histos with validated secondary mothers, invMass, pt
    TH2F**                fHistoTrueSecondaryPi0FromK0sInvMassPt;               //! array of histos with validated secondary mothers from K0s, invMass, pt
    TH1F**                fHistoTrueK0sWithPi0DaughterMCPt;                     //! array of histos with K0s with reconstructed pi0 as daughter, pt
    TH2F**                fHistoTrueSecondaryPi0FromK0lInvMassPt;               //! array of histos with validated secondary mothers from K0l, invMass, pt
    TH1F**                fHistoTrueK0lWithPi0DaughterMCPt;                     //! array of histos with K0l with reconstructed pi0 as daughter, pt
    TH2F**                fHistoTrueSecondaryPi0FromEtaInvMassPt;               //! array of histos with validated secondary mothers from eta, invMass, pt
    TH1F**                fHistoTrueEtaWithPi0DaughterMCPt;                     //! array of histos with eta with reconstructed pi0 as daughter, pt
    TH2F**                fHistoTrueSecondaryPi0FromLambdaInvMassPt;            //! array of histos with validated secondary mothers from Lambda, invMass, pt
    TH1F**                fHistoTrueLambdaWithPi0DaughterMCPt;                  //! array of histos with lambda with reconstructed pi0 as daughter, pt
    TH2F**                fHistoTrueBckGGInvMassPt;                             //! array of histos with pure gamma gamma combinatorial BG, invMass, pt
    TH2F**                fHistoTrueBckFullMesonContainedInOneClusterInvMassPt; //! array of histos with pi0 fully contained in one cluster, invMass, pt
    TH2F**                fHistoTrueBckAsymEClustersInvMassPt;                  //! array of histos with asymmetry energy distributions of clusters, invMass, pt
    TH2F**                fHistoTrueBckContInvMassPt;                           //! array of histos with contamination BG, invMass, pt
    TH2F**                fHistoTruePi0PtY;                                     //! array of histos with validated pi0, pt, Y
    TH2F**                fHistoTrueEtaPtY;                                     //! array of histos with validated eta, pt, Y
    TH2F**                fHistoTruePi0PtAlpha;                                 //! array of histos with validated pi0, pt, alpha
    TH2F**                fHistoTrueEtaPtAlpha;                                 //! array of histos with validated eta, pt, alpha
    TH2F**                fHistoTruePi0PtOpenAngle;                             //! array of histos with validated pi0, pt, openAngle
    TH2F**                fHistoTrueEtaPtOpenAngle;                             //! array of histos with validated eta, pt, openAngle
    // MC validated reconstructed quantities photons
    TH2F**                fHistoClusPhotonBGPt;                                 //! array of histos with cluster photon BG, pt, source
    TH2F**                fHistoClusPhotonPlusConvBGPt;                         //! array of histos with cluster photon plus conv BG, pt, source
    TH2F**                fHistoClustPhotonElectronBGPtM02;                     //! array of histos with cluster photon BG electron, M02 vs. pt, source
    TH2F**                fHistoClustPhotonPionBGPtM02;                         //! array of histos with cluster photon BG pion, M02 vs. pt, source
    TH2F**                fHistoClustPhotonKaonBGPtM02;                         //! array of histos with cluster photon BG kaon, M02 vs. pt, source
    TH2F**                fHistoClustPhotonK0lBGPtM02;                          //! array of histos with cluster photon BG k0l, M02 vs. pt, source
    TH2F**                fHistoClustPhotonNeutronBGPtM02;                      //! array of histos with cluster photon BG neutron, M02 vs. pt, source
    TH2F**                fHistoClustPhotonRestBGPtM02;                         //! array of histos with cluster photon BG rest, M02 vs. pt, source
    TH2F**                fHistoClustPhotonPlusConvElectronBGPtM02;             //! array of histos with cluster photon plus conv BG electron, M02 vs. pt, source
    TH2F**                fHistoClustPhotonPlusConvPionBGPtM02;                 //! array of histos with cluster photon plus conv BG pion, M02 vs. pt, source
    TH2F**                fHistoClustPhotonPlusConvKaonBGPtM02;                 //! array of histos with cluster photon plus conv BG kaon, M02 vs. pt, source
    TH2F**                fHistoClustPhotonPlusConvK0lBGPtM02;                  //! array of histos with cluster photon plus conv BG k0l, M02 vs. pt, source
    TH2F**                fHistoClustPhotonPlusConvNeutronBGPtM02;              //! array of histos with cluster photon plus conv BG neutron, M02 vs. pt, source
    TH2F**                fHistoClustPhotonPlusConvRestBGPtM02;                 //! array of histos with cluster photon plus conv BG rest, M02 vs. pt, source
    TH1F**                fHistoTrueClusGammaPt;                                //! array of histos with validated cluster (electron or photon), pt
    TH1F**                fHistoTrueClusUnConvGammaPt;                          //! array of histos with validated unconverted photon, pt
    TH1F**                fHistoTrueClusUnConvGammaMCPt;                        //! array of histos with validated unconverted photon, pt
    TH2F**                fHistoTrueClusGammaPtM02;                             //! array of histos with validated cluster (electron or photon), M02 vs pt
    TH2F**                fHistoTrueClusUnConvGammaPtM02;                       //! array of histos with validated unconverted photon, M02 vs pt
    TH1F**                fHistoTrueClusElectronPt;                             //! array of histos with validated electron, pt
    TH1F**                fHistoTrueClusConvGammaPt;                            //! array of histos with validated converted photon, pt
    TH1F**                fHistoTrueClusConvGammaMCPt;                          //! array of histos with validated converted photon, pt
    TH1F**                fHistoTrueClusConvGammaFullyPt;                       //! array of histos with validated converted photon, fully contained, pt
    TH1F**                fHistoTrueClusMergedGammaPt;                          //! array of histos with validated merged photons, electrons, dalitz, pt
    TH1F**                fHistoTrueClusMergedPartConvGammaPt;                  //! array of histos with validated merged partially converted photons, pt
    TH1F**                fHistoTrueClusDalitzPt;                               //! array of histos with validated Dalitz decay, pt
    TH1F**                fHistoTrueClusDalitzMergedPt;                         //! array of histos with validated Dalitz decay, more than one decay product in cluster, pt
    TH1F**                fHistoTrueClusPhotonFromElecMotherPt;                 //! array of histos with validated photon from electron, pt
    TH1F**                fHistoTrueClusShowerPt;                               //! array of histos with validated shower, pt
    TH1F**                fHistoTrueClusSubLeadingPt;                           //! array of histos with pi0/eta/eta_prime in subleading contribution
    TH1F**                fHistoTrueClusNParticles;                             //! array of histos with number of different particles (pi0/eta/eta_prime) contributing to cluster
    TH1F**                fHistoTrueClusEMNonLeadingPt;                         //! array of histos with cluster with largest energy by hadron
    TH1F**                fHistoTrueNLabelsInClus;                              //! array of histos with number of labels in cluster
    TH1F**                fHistoTruePrimaryClusGammaPt;                         //! array of histos with validated primary photon cluster, pt
    TH2F**                fHistoTruePrimaryClusGammaESDPtMCPt;                  //! array of histos with validated primary photon cluster, rec Pt, MC pt
    TH1F**                fHistoTruePrimaryClusConvGammaPt;                     //! array of histos with validated primary conv photon cluster, pt
    TH2F**                fHistoTruePrimaryClusConvGammaESDPtMCPt;              //! array of histos with validated primary conv photon cluster, rec Pt, MC pt
    TH2F**                fHistoTrueSecondaryClusGammaPt;                       //! array of histos with validated secondary photon cluster, pt
    TH2F**                fHistoTrueSecondaryClusConvGammaPt;                   //! array of histos with validated secondary converted photon cluster, pt
    TH2F**                fHistoTrueSecondaryClusGammaMCPt;                     //! array of histos with validated secondary photon cluster, MC pt
    TH2F**                fHistoTrueSecondaryClusConvGammaMCPt;                 //! array of histos with validated secondary converted photon cluster, MC pt
    TH2F**                fHistoTrueSecondaryClusGammaFromXFromK0sMCPtESDPt;        //! array of histos with validated secondary photon cluster from X from K0s, MC pt, rec pt
    TH2F**                fHistoTrueSecondaryClusConvGammaFromXFromK0sMCPtESDPt;    //! array of histos with validated secondary converted photon cluster from X from K0s, MC pt, rec pt
    TH2F**                fHistoTrueSecondaryClusGammaFromXFromK0lMCPtESDPt;        //! array of histos with validated secondary photon cluster from X from K0l, MC pt, rec pt
    TH2F**                fHistoTrueSecondaryClusConvGammaFromXFromK0lMCPtESDPt;    //! array of histos with validated secondary converted photon cluster from X from K0l, MC pt, rec pt
    TH2F**                fHistoTrueSecondaryClusGammaFromXFromLambdaMCPtESDPt;     //! array of histos with validated secondary photon cluster from X from Lambda, MC pt, rec pt
    TH2F**                fHistoTrueSecondaryClusConvGammaFromXFromLambdaMCPtESDPt; //! array of histos with validated secondary converted photon cluster from X from Lambda, MC pt, rec pt
    TH2F**                fHistoDoubleCountTruePi0InvMassPt;                    //! array of histos with double counted pi0s, invMass, pT
    TH2F**                fHistoDoubleCountTrueEtaInvMassPt;                    //! array of histos with double counted etas, invMass, pT
    TH2F**                fHistoDoubleCountTrueClusterGammaPt;                  //! array of histos with double counted cluster photons
    vector<Int_t>         fVectorDoubleCountTruePi0s;                           //! vector containing labels of validated pi0
    vector<Int_t>         fVectorDoubleCountTrueEtas;                           //! vector containing labels of validated eta
    vector<Int_t>         fVectorDoubleCountTrueClusterGammas;                  //! vector containing labels of validated cluster photons
    TH1F**                fHistoMultipleCountTrueClusterGamma;                  //! array of histos how often TrueClusterGammas are counted
    map<Int_t,Int_t>      fMapMultipleCountTrueClusterGammas;                   //! map containing cluster photon labels that are counted at least twice
    TH2F**                fHistoTruePi0InvMassPtAlpha;                          //! array of histogram with pure pi0 signal inv Mass, energy of cluster
    TH2F**                fHistoTruePi0PureGammaInvMassPtAlpha;                 //! array of histogram with pure pi0 signal (only pure gammas) inv Mass, energy of cluster
    TH2F**                fHistCellIDvsClusterEnergy;                           //! array of histogram with leading cell ID vs cluster Energy
    TH2F**                fHistCellIDvsClusterEnergyMax;                        //! array of histogram with leading cell ID vs maximum cluster energy in event

    // event histograms
    TH1F**                fHistoNEvents;                                        //! array of histos with event information
    TH1F**                fHistoNEventsWOWeight;                                //! array of histos with event information without event weights
    TH1F**                fHistoNGoodESDTracks;                                 //! array of histos with number of good tracks (2010 Standard track cuts)
    TH1F**                fHistoVertexZ;                                        //! array of histos with vertex z distribution for selected events
    TH1F**                fHistoNGammaCandidates;                               //! array of histos with number of gamma candidates per event
    TH1F**                fHistoNGammaCandidatesBasic;                          //! array of histos with number of gamma candidates per event for basic cluster cut
    TH2F**                fHistoNGoodESDTracksVsNGammaCandidates;               //! array of histos with number of good tracks vs gamma candidates
    TH2F**                fHistoSPDClusterTrackletBackground;                   //! array of histos with SPD tracklets vs SPD clusters for background rejection
    TH1F**                fHistoNV0Tracks;                                      //! array of histos with V0 counts
    TProfile**            fProfileEtaShift;                                     //! array of profiles with eta shift
    TProfile**            fProfileJetJetXSection;                               //! array of profiles with xsection for jetjet
    TH1F**                fHistoJetJetNTrials;                                  //! array of histos with ntrials for jetjet
    TH2F**                fHistoPtHardJJWeight;                                 //! array of histos with ntrials for jetjet
    TH1F**                fHistoEventSphericity;                                //! array of histos with event Sphericity
    TH1F**                fHistoEventSphericityAxis;                            //! array of histos with phi of the event Sphericity axis
    TH2F**                fHistoEventSphericityvsNtracks;                       //! array of histos with event Sphericity vs Ntracks
    TH2F**                fHistoEventSphericityvsNJets;                         //! array of histos with event Sphericity vs NJets
    TH2F**                fHistoEventMultiplicityvsNJets;                       //! array of histos with event Multiplicity vs NJets
    TH2F**                fHistoTrueSphericityvsRecSphericity;                  //! array of histos with true sphericity vs rec. sphericity
    TH2F**                fHistoTrueMultiplicityvsRecMultiplicity;              //! array of histos with true multiplicity vs rec. multiplicity
    TH2F**                fHistoEventSphericityvsHighpt;                        //! array of histos with event Sphericity vs highest pt
    TH2F**                fHistoEventSphericityvsTotalpt;                       //! array of histos with event Sphericity vs total pt
    TH2F**                fHistoEventSphericityvsMeanpt;                        //! array of histos with event Sphericity vs mean pt
    TH1F**                fHistoPionSpectrum;                                   //! array of histos with charged pion spectrum
    TH1F**                fHistoProtonSpectrum;                                 //! array of histos with proton spectrum
    TH1F**                fHistoKaonSpectrum;                                   //! array of histos with charged kaon spectrum
    TH1F**                fHistoNPionSpectrum;                                  //! array of histos with Neutral pion spectrum
    TH1F**                fHistoEtaSpectrum;                                    //! array of histos with Eta spectrum
    TH1F**                fHistoDMesonSpectrum;                                 //! array of histos with D0 meson spectrum
    TTree**               tTreeSphericity;                                      //! array of trees with sphericity correlations
    Float_t               fRecSph;                                              //! Reconstructed sphericity
    Float_t               fTrueSph;                                             //! True Sphericity
    Float_t               fPi0Pt;
    Float_t               fPi0InvMass;

    TH1F**                 fHistoPtJet;                                          // Histogram of Jet Pt
    TH1F**                 fHistoJetEta;                                         // Histogram of Jet Eta
    TH1F**                 fHistoJetPhi;                                         // Histogram of Jet Phi
    TH1F**                 fHistoJetArea;                                        // Histogram of Jet Area
    TH1F**                 fHistoNJets;                                          // Histogram of number of jets
    TH1F**                 fHistoEventwJets;                                     // Histogram of number of events with jets > 0
    TH1F**                 fHistoJetPi0PtRatio;                                  // Histogram of PtPi0/PtJet
    TH1F**                 fHistoDoubleCounting;                                 // Histogram if NM candidates are defined within multiple jets
    TH2F**                 fHistoJetMotherInvMassPt;                             // Histogram of NM candidates with a jet in the event
    TH2F**                 fHistoPi0InJetMotherInvMassPt;                        // Histogram of NM candidates that are inside a jet
    TH2F**                 fHistoMotherBackJetInvMassPt;                         // Histogram of Backgrouns candidates that are involved with jets
    TH2F**                 fHistoRJetPi0Cand;                                    // Histogram of RJetPi0Cand vs Pt
    TH2F**                 fHistoEtaPhiJetPi0Cand;                               // Histogram of delta eta and delta phi distr between jet and NM candidates
    TH2F**                 fHistoEtaPhiJetWithPi0Cand;                           // Histogram of delta eta and delta phi distr when pi0 is inside a jet
    TH2F**                 fHistoJetFragmFunc;                                   // Histogram to determine fragmentation function
    TH2F**                 fHistoJetFragmFuncZInvMass;                           // Histogram of Inv Mass distribution with z
    TH2F**                 fHistoTruevsRecJetPt;                                 // Histogram of true jet pt vs reconstructed jet pt
    TH2F**                 fHistoTruePi0JetMotherInvMassPt;                      // Histogram of true pi0s in an event with a jet
    TH2F**                 fHistoTruePi0InJetMotherInvMassPt;                    // Histogram of true pi0s in a jet
    TH2F**                 fHistoTruePrimaryPi0JetInvMassPt;                     // Histogram of true primary pi0s in an event with a jet
    TH2F**                 fHistoTruePrimaryPi0inJetInvMassPt;                   // Histogram of true primary pi0s in a jet
    TH2F**                 fHistoTruePrimaryPi0InJetInvMassTruePt;               // Histogram of true primary pi0s in a jet with their true pt
    TH1F**                 fHistoTrueDoubleCountingPi0Jet;                       // Histogram of when a true pi0 is defined to be in multiple jets
    TH2F**                 fHistoTrueEtaJetMotherInvMassPt;                      // Histogram of true etas in an event with a jet
    TH2F**                 fHistoTrueEtaInJetMotherInvMassPt;                    // Histogram of true etas in a jet
    TH2F**                 fHistoTruePrimaryEtaJetInvMassPt;                     // Histogram of true primary etas in an event with a jet
    TH2F**                 fHistoTruePrimaryEtainJetInvMassPt;                   // Histogram of true primary etas in a jet
    TH1F**                 fHistoTrueDoubleCountingEtaJet;                       // Histogram of when a true eta is defined to be in multiple jets
    TH2F**                 fHistoTruePi0JetFragmFunc;                            // Histogram to determine true pi0 fragmentation function
    TH2F**                 fHistoTruePi0JetFragmFuncZInvMass;                    // Histogram to determine true pi0 Inv Mass distribution with z
    TH2F**                 fHistoTrueEtaJetFragmFunc;                            // Histogram to determine true eta fragmentation function
    TH2F**                 fHistoTrueEtaJetFragmFuncZInvMass;                    // Histogram to determine true eta Inv Mass distribution with z
    TH2F**                 fHistoMCPi0GenVsNClus;                                // pi0 produced on gen level vs Nclus for merging studies
    TH2F**                 fHistoMCPi0GenFoundInOneCluster;                      // pi0 produced on gen level where both decay photons were found in same cluster (merged)
    TH2F**                 fHistoMCPi0GenFoundInTwoCluster;                      // pi0 produced on gen level where both decay photons were found in different clusters
    TH2F**                fHistoMCEtaGenFoundInOneCluster;                      // pi0 produced on gen level where both decay photons were found in same cluster (merged)
    TH2F**                 fHistoMCEtaGenFoundInTwoCluster;                      // pi0 produced on gen level where both decay photons were found in different clusters
    TH2F**                 fHistoMCGammaConvRvsPt;                      // pi0 produced on gen level where both decay photons were found in different clusters
    TH1F**                 fHistoMCPi0JetInAccPt;                                // Histogram with weighted pi0 in a jet event in acceptance, pT
    TH1F**                 fHistoMCPi0inJetInAccPt;                              // Histogram with weighted pi0 in a jet in acceptance, pT
    TH1F**                 fHistoMCEtaJetInAccPt;                                // Histogram with weighted eta in a jet event in acceptance, pT
    TH1F**                 fHistoMCEtainJetInAccPt;                              // Histogram with weighted eta in a jet in acceptance, pT
    TH1F**                 fHistoMCPi0JetEventGenerated;                         // Histogram with mesons in a jet event generated, pT
    TH1F**                 fHistoMCPi0inJetGenerated;                            // Histogram with mesons in a jet generated, pT
    TH1F**                 fHistoMCEtaJetEventGenerated;                         // Histogram with mesons in a jet event generated, pT
    TH1F**                 fHistoMCEtainJetGenerated;                            // Histogram with mesons in a jet generated, pT
    TH2F**                 fHistoTrueSecondaryPi0FromK0sJetInvMassPt;            // Histogram with validated secondary mothers from K0s in an event with a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0FromK0sinJetInvMassPt;          // Histogram with validated secondary mothers from K0s in a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0FromLambdaJetInvMassPt;         // Histogram with validated secondary mothers from lambda in an event with a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0FromLambdainJetInvMassPt;       // Histogram with validated secondary mothers from lambda in a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0FromK0lJetInvMassPt;            // Histogram with validated secondary mothers from K0l in an event with a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0FromK0linJetInvMassPt;          // Histogram with validated secondary mothers from K0l in a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0InvJetMassPt;                   // Histogram with validated secondary mothers in an event with a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0InvinJetMassPt;                 // Histogram with validated secondary mothers in a jet, invMass, pt
    TH2F**                 fHistoMotherPi0inJetPtY;                              // Histogram with the rapidity of the validated pi0s in jets
    TH2F**                 fHistoMotherEtainJetPtY;                              // Histogram with the rapidity of the validated etas in jets
    TH2F**                 fHistoMotherPi0inJetPtPhi;                            // Histogram with the phi of the validated pi0s in jets
    TH2F**                 fHistoMotherEtainJetPtPhi;                            // Histogram with the phi of the validated etas in jets
    TH1F**                 fNumberOfClusters;                                    // Histogram with total number of clusters
    TH1F**                 fNumberOfClustersinJets;                              // Histogram with number of clusters in jets
    TH2F**                 fEnergyRatio;                                         // Histogram with ratio of energy deposited in cluster
    TH2F**                 fEnergyRatioinJets;                                   // Histogram with ratio of energy deposited in cluster in jets
    TH2F**                 fEnergyRatioGamma1;                                   // Histogram with ratio of energy deposited in cluster when first label is a photon
    TH2F**                 fEnergyRatioGamma1inJets;                             // Histogram with ratio of energy deposited in cluster when first label is a photon in jets
    TH2F**                 fEnergyRatioGammaAnywhere;                            // Histogram with ratio of energy deposited in cluster when a label is photon
    TH2F**                 fEnergyRatioGammaAnywhereinJets;                      // Histogram with ratio of energy deposited in cluster when a label is photon in jets
    TH2F**                 fEnergyDeposit;                                       // Histogram with total energy deposited in cluster
    TH2F**                 fEnergyDepositinJets;                                 // Histogram with total energy deposited in cluster in jets
    TH2F**                 fEnergyDepGamma1;                                     // Histogram with total energy deposited when first label is photon
    TH2F**                 fEnergyDepGamma1inJets;                               // Histogram with total energy deposited when first label is photon in jets
    TH2F**                 fEnergyDepGammaAnywhere;                              // Histogram with total energy deposited when a label is photon
    TH2F**                 fEnergyDepGammaAnywhereinJets;                        // Histogram with total energy deposited when a label is photon in jets
    TH1F**                 fEnergyRatioGamma1Helped;                             // Histogram with energy of photon when it is helped by other particles
    TH1F**                 fEnergyRatioGamma1HelpedinJets;                       // Histogram with energy of photon when it is helped by other particles in jets
    TH2F**                 fClusterEtaPhiJets;                                   // Histogram of phi and eta distribution of clusters in jets
    TH2F**                 fHistoUnfoldingAsData;                                // Histogram to use for jet pi0 unfolding
    TH2F**                 fHistoUnfoldingMissed;                                // Histogram to use for jet pi0 unfolding
    TH2F**                 fHistoUnfoldingReject;                                // Histogram to use for jet pi0 unfolding
    TH2F**                 fHistoUnfoldingAsDataInvMassZ;                        // Histogram to use for jet pi0 unfolding for z inmass
    TH2F**                 fHistoUnfoldingMissedInvMassZ;                        // Histogram to use for jet pi0 unfolding for z inmass
    TH2F**                 fHistoUnfoldingRejectInvMassZ;                        // Histogram to use for jet pi0 unfolding for z inmass

    vector<Double_t>      fVectorJetPt;                                         // Vector of JetPt
    vector<Double_t>      fVectorJetPx;                                         // Vector of JetPx
    vector<Double_t>      fVectorJetPy;                                         // Vector of JetPy
    vector<Double_t>      fVectorJetPz;                                         // Vector of JetPz
    vector<Double_t>      fVectorJetEta;                                        // Vector of JetEta
    vector<Double_t>      fVectorJetPhi;                                        // Vector of JetPhi
    vector<Double_t>      fVectorJetArea;                                       // Vector of JetArea
    vector<Double_t>      fTrueVectorJetPt;                                     // Vector of True JetPt
    vector<Double_t>      fTrueVectorJetPx;                                     // Vector of True JetPx
    vector<Double_t>      fTrueVectorJetPy;                                     // Vector of True JetPy
    vector<Double_t>      fTrueVectorJetPz;                                     // Vector of True JetPz
    vector<Double_t>      fTrueVectorJetEta;                                    // Vector of True JetEta
    vector<Double_t>      fTrueVectorJetPhi;                                    // Vector of True JetPhi

    // tree for identified particle properties
    TTree**               tTrueInvMassROpenABPtFlag;                            //! array of trees with
    Float_t               fInvMass;                                             //! InvMass,
    Float_t               fRconv;                                               //! Rconv,
    Float_t               fOpenRPrim;                                           //! opening angle at R = 0,
    Float_t               fInvMassRTOF;                                         //! InvMass at R=375 cm,
    Float_t               fPt;                                                  //! momentum,
    UChar_t               iFlag;                                                //! flag (0 = gamma, 1 = pi0, 2 = eta)

    // tree for alpha/opening angle studies
    TTree**               tSigInvMassPtAlphaTheta;                              //! array of trees
    TTree**               tBckInvMassPtAlphaTheta;                              //! array of trees
    Float_t               fInvMassTreeInvMass;
    Float_t               fInvMassTreePt;
    Float_t               fInvMassTreeAlpha;
    Float_t               fInvMassTreeTheta;
    Int_t                 fInvMassTreeMixPool;
    Float_t               fInvMassTreeZVertex;
    Float_t               fInvMassTreeEta;

    // tree for E/p studies
    TTree**               tClusterEOverP;                                       //! array of trees with tree for E/p studies
    Float_t               fClusterE;                                            //! cluster energy
    Float_t               fClusterM02;                                          //! cluster M02
    Float_t               fClusterM20;                                          //! cluster M20
    Float_t               fClusterEP;                                           //! cluster-track E/p
    Int_t                 fClusterLeadCellID;                                   //! cellID of leading cell in cluster
    Int_t                 fClusterClassification;                               //! classification of cluster in MC
    Float_t               fDeltaEta;                                            //! matching residual track <-> cluster
    Float_t               fDeltaPhi;                                            //! matching residual track <-> cluster
    Float_t               fTrackPt;                                             //! track Pt
    Int_t                 fTrackPID_e;                                          //! track PID e
    Int_t                 fTrackPID_Pi;                                         //! track PID Pi
    Int_t                 fTrackPID_K;                                          //! track PID K
    Int_t                 fTrackPID_P;                                          //! track PID P
    Float_t               fClusterIsoSumClusterEt;                              //! sum of Et of clusters within R<0.2
    Float_t               fClusterIsoSumTrackEt;                                //! sum of Et of tracks within R<0.2

    // tree for timing cut effi studies
    TTree**               tClusterTimingEff;                                    //! array of trees with tree for timing effi studies
    Float_t               fClusterTimeTag;                                      //! cluster time of tagged photon
    Float_t               fClusterTimeProbe;                                    //! cluster time of probe photon
    Float_t               fClusterETag;                                         //! cluster E of tagged photon
    Float_t               fClusterEProbe;                                       //! cluster E of probe photon

    // hists for nonlineartiy calibration
//    TH2F**                fHistoTruePi0NonLinearity;                            //! E_truth/E_rec vs E_rec for TruePi0s
//    TH2F**                fHistoTrueEtaNonLinearity;                            //! E_truth/E_rec vs E_rec for TrueEtas

    // additional variables
    Double_t              fEventPlaneAngle;                                     // EventPlaneAngle
    TRandom3              fRandom;                                              // random
    Int_t                 fnCuts;                                               // number of cuts to be analysed in parallel
    Int_t                 fiCut;                                                // current cut
    Int_t                 fIsHeavyIon;                                          // switch for pp = 0, PbPb = 1, pPb = 2
    Bool_t                fDoLightOutput;                                       // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    Bool_t                fDoECalibOutput;                                      // switch for running with E-Calib Histograms in Light Output, kFALSE -> no E-Calib Histograms, kTRUE -> with E-Calib Histograms
    Bool_t                fDoMesonAnalysis;                                     // flag for meson analysis
    Int_t                 fDoMesonQA;                                           // flag for meson QA
    Int_t                 fDoClusterQA;                                         // flag for cluster QA
    Bool_t                fIsFromDesiredHeader;                                 // flag for MC headers
    Bool_t                fIsOverlappingWithOtherHeader;                        // flag for particles in MC overlapping between headers
    Int_t                 fIsMC;                                                // flag for MC information
    Bool_t                fDoTHnSparse;                                         // flag for using THnSparses for background estimation
    Bool_t                fSetPlotHistsExtQA;                                   // flag for extended QA hists
    Bool_t                fDoSoftAnalysis;                                      // bool for Sphericity analysis without jets
    Double_t              fWeightJetJetMC;                                      // weight for Jet-Jet MC
    Bool_t                fDoInOutTimingCluster;                                // manual timing cut for cluster to combine cluster within timing cut and without
    Double_t              fMinTimingCluster;                                    // corresponding ranges, min
    Double_t              fMaxTimingCluster;                                    // corresponding ranges, max
    Bool_t                fEnableSortForClusMC;                                 // switch on sorting for MC labels in cluster
    Bool_t                fProduceCellIDPlots;                                  // switch to produce CellID plots for fDoClusterQA==2
    Bool_t                fProduceTreeEOverP;                                   // flag for producing tree for E/p studies
    TTree*                tBrokenFiles;                                         // tree for keeping track of broken files
    TObjString*           fFileNameBroken;                                      // string object for broken file name
    TTree*                tTriggerFiles_wL0;                                    // tree for keeping track of triggering files in MB
    TTree*                tTriggerFiles_woL0;                                   // tree for keeping track of triggering files in MB
    TObjString*           fFileNameTrigger;                                     // string object for triggering filename in MB
    TTree*                tClusterQATree;                                       // tree for specific cluster QA
    TObjString*           fCloseHighPtClusters;                                 // file name to indicate clusters with high pT (>15 GeV/c) very close to each other (<17 mrad)
    TGenPhaseSpace        fGenPhaseSpace;                                       // For generation of decays into two gammas
    TClonesArray*         fAODMCTrackArray;                                     // storage of track array

    Int_t                 fLocalDebugFlag;                                      // debug flag for local running, must be '0' for grid running
    Bool_t                fAllowOverlapHeaders;                                 // enable overlapping headers for cluster selection
    Int_t                 fNCurrentClusterBasic;                                // current number of cluster without minE
    Int_t                 fTrackMatcherRunningMode;                             // CaloTrackMatcher running mode
    Bool_t                fDoPi0Only;                                           // switches ranges of histograms and binnings to pi0 specific analysis
  private:
    AliAnalysisTaskGammaCalo(const AliAnalysisTaskGammaCalo&);                  // Prevent copy-construction
    AliAnalysisTaskGammaCalo &operator=(const AliAnalysisTaskGammaCalo&);       // Prevent assignment

    ClassDef(AliAnalysisTaskGammaCalo, 87);
};

#endif
