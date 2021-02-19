#ifndef ALIANLYSISTASKGAMMACONVCALO_cxx
#define ALIANLYSISTASKGAMMACONVCALO_cxx

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
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskConvJet.h"
#include "AliAnalysisTaskJetOutlierRemoval.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"
#include "THnSparse.h"
#include <vector>
#include <map>

class AliAnalysisTaskGammaConvCalo : public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskGammaConvCalo();
    AliAnalysisTaskGammaConvCalo(const char *name);
    virtual ~AliAnalysisTaskGammaConvCalo();

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
    // Function to set correction task setting
    void SetCorrectionTaskSetting(TString setting) {fCorrTaskSetting = setting;}

    // base functions for selecting photon and meson candidates in reconstructed data
    void ProcessClusters();
    void ProcessPhotonCandidates();
    void CalculatePi0Candidates();

    // MC functions
    void SetIsMC                        ( Int_t isMC)                                       { fIsMC = isMC                              ;}
    void ProcessMCParticles             ();
    void ProcessAODMCParticles          ();
    void RelabelAODPhotonCandidates     ( Bool_t mode);
    void ProcessTruePhotonCandidates    ( AliAODConversionPhoton* TruePhotonCandidate);
    void ProcessTrueClusterCandidates   ( AliAODConversionPhoton* TruePhotonCandidate,
                                          Float_t clusM02);
    void ProcessTrueClusterCandidatesAOD( AliAODConversionPhoton* TruePhotonCandidate,
                                          Float_t clusM02);
    void ProcessTruePhotonCandidatesAOD ( AliAODConversionPhoton* TruePhotonCandidate);
    void ProcessTrueMesonCandidates     ( AliAODConversionMother *Pi0Candidate,
                                          AliAODConversionPhoton *TrueGammaCandidate0,
                                          AliAODConversionPhoton *TrueGammaCandidate1,
                                          Bool_t matched);
    void ProcessTrueMesonCandidatesAOD  ( AliAODConversionMother *Pi0Candidate,
                                          AliAODConversionPhoton *TrueGammaCandidate0,
                                          AliAODConversionPhoton *TrueGammaCandidate1,
                                          Bool_t matched);
    void ProcessConversionPhotonsForMissingTags     ();
    void ProcessConversionPhotonsForMissingTagsAOD  ();
    void ProcessJets();

    // switches for additional analysis streams or outputs
    void SetDoPrimaryTrackMatching      ( Bool_t flag )                                     { fDoPrimaryTrackMatching = flag              ;}
    void SetLightOutput                 ( Int_t flag )                                     { fDoLightOutput = flag                       ;}
    void SetPi0EtaSwitch                ( Int_t flag )                                     { fPi0EtaSwitch = flag                       ;}
    void SetECalibOutput                ( Bool_t flag )                                     { fDoECalibOutput = flag                      ;}
    void SetDoMesonAnalysis             ( Bool_t flag )                                     { fDoMesonAnalysis = flag                     ;}
    void SetDoMesonQA                   ( Int_t flag )                                      { fDoMesonQA = flag                           ;}
    void SetDoPhotonQA                  ( Int_t flag )                                      { fDoPhotonQA = flag                          ;}
    void SetDoClusterQA                 ( Int_t flag )                                      { fDoClusterQA = flag                         ;}
    void SetUseTHnSparse                ( Bool_t flag )                                     { fDoTHnSparse = flag                         ;}
    void SetPlotHistsExtQA              ( Bool_t flag )                                     { fSetPlotHistsExtQA = flag                   ;}
    void SetDoTreeConvGammaShowerShape  ( Bool_t flag )                                     { fDoConvGammaShowerShapeTree = flag          ;}
    void SetDoTreeInvMassShowerShape    ( Bool_t flag )                                     { fDoInvMassShowerShapeTree = flag            ;}
    void SetAllowOverlapHeaders         ( Bool_t allowOverlapHeader )                       { fAllowOverlapHeaders = allowOverlapHeader   ;}
    void SetDoMaterialBudgetWeightingOfGammasForTrueMesons(Bool_t flag)                     { fDoMaterialBudgetWeightingOfGammasForTrueMesons = flag;}
    void SetDoHBTHistoOutput            ( Bool_t flag )                                     { fDoHBTHistoOutput = flag                    ;}

    // Setting the cut lists for the conversion photons
    void SetEventCutList                ( Int_t nCuts,
                                          TList *CutArray)                                  {
                                                                                              fnCuts = nCuts                              ;
                                                                                              fEventCutArray = CutArray                   ;
                                                                                            }

    // Setting the cut lists for the conversion photons
    void SetConversionCutList           ( Int_t nCuts,
                                          TList *CutArray)                                  {
                                                                                              fnCuts = nCuts                              ;
                                                                                              fCutArray = CutArray                        ;
                                                                                            }

      // Setting the cut lists for the calo photons
    void SetCaloCutList                 ( Int_t nCuts,
                                          TList *CutArray)                                  {
                                                                                              fnCuts = nCuts                              ;
                                                                                              fClusterCutArray = CutArray                 ;
                                                                                            }

    // Setting the cut lists for the meson
    void SetMesonCutList                ( Int_t nCuts,
                                          TList *CutArray)                                  {
                                                                                              fnCuts = nCuts                              ;
                                                                                              fMesonCutArray = CutArray                   ;
                                                                                            }

    // BG HandlerSettings
    void CalculateBackground            ();
    void CalculateBackgroundSwapp       ();
    void CalculateBackgroundRP          ();
    void RotateParticle                 ( AliAODConversionPhoton *gamma );
    void RotateParticleAccordingToEP    ( AliAODConversionPhoton *gamma,
                                          Double_t previousEventEP,
                                          Double_t thisEventEP );
    void SetMoveParticleAccordingToVertex       ( Bool_t flag )                             { fMoveParticleAccordingToVertex = flag       ;}
    void FillPhotonCombinatorialBackgroundHist  ( AliAODConversionPhoton *TruePhotonCandidate,
                                                  Int_t pdgCode[] );
    void MoveParticleAccordingToVertex          ( AliAODConversionPhoton* particle,
                                                  const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
    void UpdateEventByEventData         ();

    // Additional functions for convenience
    void SetLogBinningXTH2              ( TH2* histoRebin );
    Int_t GetSourceClassification       (Int_t daughter,
                                         Int_t pdgCode );
    Bool_t CheckVectorOnly              ( vector<Int_t> &vec,
                                          Int_t tobechecked );
    Bool_t CheckVectorForDoubleCount    ( vector<Int_t> &vec,
                                          Int_t tobechecked );

    void FillMultipleCountMap           ( map<Int_t,Int_t> &ma,
                                          Int_t tobechecked );
    void FillMultipleCountHistoAndClear ( map<Int_t,
                                          Int_t> &ma,
                                          TH1F* hist );
    void ProcessAODMCParticlesForHBT    ();
    void CalculateHBTBackgroundMC       ();

        // Function to enable MC label sorting
    void SetEnableSortingOfMCClusLabels (Bool_t enableSort) { fEnableSortForClusMC   = enableSort;}

    void SetTrackMatcherRunningMode(Int_t mode){fTrackMatcherRunningMode = mode;}

  protected:
    AliV0ReaderV1*                      fV0Reader;              //! basic photon Selection Task
    TString                             fV0ReaderName;
    TString                             fCaloTriggerHelperName;
    TString                             fCorrTaskSetting;       // Correction Task Special Name
    AliGammaConversionAODBGHandler**    fBGHandler;             //! BG handler for Conversion
    AliConversionAODBGHandlerRP**       fBGHandlerRP;           //! BG handler for Conversion (possibility to mix with respect to RP)
    AliGammaConversionAODBGHandler**    fBGClusHandler;         //! BG handler for Cluster
    AliConversionAODBGHandlerRP**       fBGClusHandlerRP;       //! BG handler for Cluster (possibility to mix with respect to RP)
    AliGammaConversionAODBGHandler**    fBGHBTTrueGammaHandler; //! BG handler for HBTTrueGamma
    AliGammaConversionAODBGHandler**    fBGHBTGenGammaHandler;  //! BG handler for HBTGenGamma
    AliVEvent*                          fInputEvent;            //! current event
    AliMCEvent*                         fMCEvent;               //! corresponding MC event
    TList**                             fCutFolder;             //! Array of lists for containers belonging to cut
    TList**                             fESDList;               //! Array of lists with histograms with reconstructed properties
    TList**                             fBackList;              //! Array of lists with BG THnSparseF
    TList**                             fMotherList;            //! Array of lists with Signal THnSparseF
    TList**                             fGammaERM02;            //! Array of lists with conv photon shower shape trees
    TList**                             fInvMassShowerShape;    //! Array of lists with invMass shower shape trees
    TList**                             fTrueList;              //! Array of lists with histograms with MC validated reconstructed properties
    TList**                             fMCList;                //! Array of lists with histograms with pure MC information
    TList**                             fClusterOutputList;     //! Array of lists of output histograms for cluster photons
    TList*                              fOutputContainer;       //! Output container
    TClonesArray*                       fReaderGammas;          //! Array with conversion photons selected by V0Reader Cut
    TList*                              fGammaCandidates;       //! current list of photon candidates
    TList*                              fClusterCandidates;     //! current list of cluster candidates
    TList*                              fEventCutArray;         // List with Event Cuts
    TList*                              fCutArray;              // List with Conversion Cuts
    TList*                              fClusterCutArray;       // List with Cluster Cuts
    TList*                              fMesonCutArray;         // List with Meson Cuts
    AliAnalysisTaskConvJet*             fConvJetReader;         //! JetReader
    AliAnalysisTaskJetOutlierRemoval*   fOutlierJetReader;      //! JetReader
    Bool_t                              fDoJetAnalysis;         //! Bool to produce Jet Plots
    Bool_t                              fDoJetQA;               //! Bool to produce Jet QA Plots
    TList**                             fJetHistograms;         //! Jet Histograms
    TList**                             fTrueJetHistograms;     //! True Jet Histograms
    Int_t                               fMaxPtNearEMCalPlace;   //! Place in jet vector of highest pt jet that is near the EMCal
    Bool_t                              fJetNearEMCal;          //! If a jet is near the EMCal in the current event
    TList*                              fTrueGammaCandidatesConv;   //! current list of true Conv photon candidates
    TList*                              fTrueGammaCandidatesCalo;   //! current list of true Calo photon candidates
    TList*                              fMCGammaCandidates;     //! current list of MC generated photon candidates
    TClonesArray*                       fAODMCTrackArray;       //! pointer to track array
    AliCaloTriggerMimicHelper**         fCaloTriggerMimicHelper;//!Array wich points to AliCaloTriggerMimicHelper for each Event Cut
    map<TString, Bool_t>                fSetEventCutsOutputlist;//! Store, if Output list for Event Cut has already been added

    //histograms for Conversions reconstructed quantities
    TH1F**                  fHistoConvGammaPt;                  //! histogram conversion photon pT
    TTree**                 fTreeConvGammaPtDcazCat;            //! tree with dca for conversions
    UShort_t                fPtGamma;                           //! pt of conversion for tree
    Short_t                 fDCAzPhoton;                        //! dcaz of conversion for tree
    UChar_t                 fCharCatPhoton;                     //! category of conversion for tree
    UChar_t                 fCharPhotonMCInfo;                  //! MC info of conversion for tree
                      // 0: garbage,
                      // 1: background
                      // 2: secondary photon not from eta or k0s,
                      // 3: secondary photon from eta,
                      // 4: secondary photon from k0s,
                      // 5: dalitz
                      // 6: primary gamma
    TTree**                 tESDGammaERM02;                     //! tree with E_gamma vs R_conv vs M02 vs M20
    Float_t                 tESDClusE;                          //! energy of conv gamma for tree
    Float_t                 tESDGammaConvR;                     //! conv R of conv gamma for tree
    Float_t                 tESDClusterM02;                     //! M02 of matched cluster for tree
    Float_t                 tESDClusterM20;                     //! M20 of matched cluster for tree
    Float_t                 tESDClusterEta;                     //! Eta of matched cluster for tree
    Float_t                 tESDClusterPhi;                     //! Phi of matched cluster for tree
    Float_t                 tESDClusterNCells;                  //! NCells of matched cluster for tree
    Float_t                 tESDClusterMaxECell;                //! MaxEnergyCell of matched cluster for tree
    Float_t                 tESDClusterNLM;                     //! NLM of matched cluster for tree

    TTree**                 tESDInvMassShowerShape;             //! tree with cluster shower shape for meson cand
    Float_t                 tESDIMMesonInvMass;                 //! invMass of meson cand
    Float_t                 tESDIMMesonPt;                      //! pT of meson cand
    Float_t                 tESDIMClusE;                        //! energy of clus gamma for tree
    Float_t                 tESDIMClusterM02;                   //! M02 of cluster for tree
    Float_t                 tESDIMClusterM20;                   //! M20 of cluster for tree
    Int_t                   tESDIMClusterLeadCellID;            //! leading cellID of cluster
    Int_t                   tESDIMClusterClassification;        //! classification of cluster in MC
    Float_t                 tESDIMClusMatchedTrackPt;           //! pT of closest matched track
    Float_t                 tESDIMClusTrackDeltaEta;            //! matching residual closest track <-> cluster
    Float_t                 tESDIMClusTrackDeltaPhi;            //! matching residual closest track <-> cluster
    Float_t                 tESDIMClusterIsoSumClusterEt;       //! sum of Et of clusters within R<0.2
    Float_t                 tESDIMClusterIsoSumTrackEt;         //! sum of Et of tracks within R<0.2
    map<Long_t,Int_t>       tESDmapIsClusterAcceptedWithoutTrackMatch; //! needed for filling of tree -> compute cluster isolation

    //histograms for mesons reconstructed quantities
    TH2F**                  fHistoMotherInvMassPt;              //! array of histogram with signal + BG for same event photon pairs, inv Mass, pt
    TH2F**                  fHistoMotherMatchedInvMassPt;       //! array of histogram with signal + BG for same event photon pairs, inv Mass, pt
    THnSparseF**            fSparseMotherInvMassPtZM;           //! array of THnSparseF with signal + BG for same event photon pairs, inv Mass, pt
    TH2F**                  fHistoMotherBackInvMassPt;          //! array of histogram with BG for mixed event photon pairs, inv Mass, pt
    THnSparseF**            fSparseMotherBackInvMassPtZM;       //! array of THnSparseF with BG for same event photon pairs, inv Mass, pt
    TH2F**                  fHistoMotherInvMassPtAlpha;         //! array of histograms with alpha cut of 0.1 for inv mass vs pt
    TH2F**                  fHistoMotherPi0PtY;                 //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, Y
    TH2F**                  fHistoMotherEtaPtY;                 //! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65, pt, Y
    TH2F**                  fHistoMotherPi0PtAlpha;             //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, alpha
    TH2F**                  fHistoMotherEtaPtAlpha;             //! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65, pt, alpha
    TH2F**                  fHistoMotherPi0PtOpenAngle;         //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17, pt, openAngle
    TH2F**                  fHistoMotherEtaPtOpenAngle;         //! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65, pt, openAngle
    TH2F**                  fHistoMotherPi0ConvPhotonEtaPhi;    //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17 ,eta/phi of conversion photon
    TH2F**                  fHistoMotherEtaConvPhotonEtaPhi;    //! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65 ,eta/phi of conversion photon
    TH2F**                  fHistoMotherInvMassECalib;          //! array of histogram with signal + BG for same event photon pairs, inv Mass, energy of cluster
    TH2F**                  fHistoMotherBackInvMassECalib;      //! array of histogram with BG for mixed event photon pairs, inv Mass, energy of cluster
    TH2F**                  fHistoMotherInvMassECalibPCM;       //! array of histogram with signal + BG for same event photon pairs, inv Mass, energy of PCM photon
    TH2F**                  fHistoMotherBackInvMassECalibPCM;   //! array of histogram with BG for mixed event photon pairs, inv Mass, energy of PCM photon

    // histograms for rec photons tagged by Calo
    TH2F**                  fHistoPhotonPairPtconv;             //! array of histo for pairs vs. pt of converted photon
    TH2F**                  fHistoPhotonPairMixedEventPtconv;   //! array of histo for pairs vs. pt of converted photon
    // histograms for rec photon clusters
    TH1F**                  fHistoClusGammaPt;                  //! array of histos with cluster, pt
    TH1F**                  fHistoClusGammaE;                   //! array of histos with cluster, E
    TH1F**                  fHistoClusGammaE_BothBM;            //! array of histos with cluster, E; Only clusters, which are good on triggered bad map and analysis bad map; Only MB
    TH1F**                  fHistoClusGammaE_BothBM_highestE;   //! array of histos with cluster, E; Only highest cluster in event, which is good on triggered bad map and analysis bad map;  MB and tigger but use only triggered clusters for trigger
    TH1F**                  fHistoClusGammaE_AnaBM_highestE;    //! array of histos with cluster, E; Only highest cluster in event, which is good on analysis bad map; Only MB
    TH1F**                  fHistoClusGammaE_onlyTriggered;     //! array of histos with cluster, E
    TH1I**                  fHistoGoodMesonClusters;              //! Histograms which stores if Pi0 Clusters Trigger
    TH1F**                  fHistoClusOverlapHeadersGammaPt;    //! array of histos with cluster, pt overlapping with other headers
    TH1F**                  fHistoClusAllHeadersGammaPt;        //! array of histos with cluster, pt all headers
    TH1F**                  fHistoClusRejectedHeadersGammaPt;   //! array of histos with cluster, pt rejected headers

    //histograms for pure MC quantities
    TH1I**                  fHistoMCHeaders;                    //! array of histos for header names
    //histograms for pure MC quantities for gamma part
    TH1F**                  fHistoMCAllGammaPt;                 //! array of histos with all gamma, pT
    TH1F**                  fHistoMCAllGammaEMCALAccPt;         //! array of histos with all gamma in EMCAL acceptance, pT
    TH2F**                  fHistoMCAllSecondaryGammaPt;        //! array of histos with all secondary gamma, pT
    TH1F**                  fHistoMCDecayGammaPi0Pt;            //! array of histos with decay gamma from pi0, pT
    TH1F**                  fHistoMCDecayGammaRhoPt;            //! array of histos with decay gamma from rho, pT
    TH1F**                  fHistoMCDecayGammaEtaPt;            //! array of histos with decay gamma from eta, pT
    TH1F**                  fHistoMCDecayGammaOmegaPt;          //! array of histos with decay gamma from omega, pT
    TH1F**                  fHistoMCDecayGammaEtapPt;           //! array of histos with decay gamma from eta', pT
    TH1F**                  fHistoMCDecayGammaPhiPt;            //! array of histos with decay gamma from phi, pT
    TH1F**                  fHistoMCDecayGammaSigmaPt;          //! array of histos with decay gamma from Sigma0, pT
    TH1F**                  fHistoMCConvGammaPt;                //! array of histos with converted gamma, pT
    TH2F**                  fHistoMCSecondaryConvGammaPt;       //! array of histos wiht secondary converted gamma, pt

    //histograms for pure MC quantities for mesaon analysis
    TH1F**                  fHistoMCPi0Pt;                      //! array of histos with weighted pi0, pT
    TH1F**                  fHistoMCPi0WOWeightPt;              //! array of histos with unweighted pi0, pT
    TH1F**                  fHistoMCPi0WOEvtWeightPt;           //! array of histos without event weights pi0, pT
    TH1F**                  fHistoMCEtaPt;                      //! array of histos with weighted eta, pT
    TH1F**                  fHistoMCEtaWOWeightPt;              //! array of histos with unweighted eta, pT
    TH1F**                  fHistoMCEtaWOEvtWeightPt;           //! array of histos without evt weight eta, pT
    TH1F**                  fHistoMCPi0InAccPt;                 //! array of histos with weighted pi0 in acceptance, pT
    TH1F**                  fHistoMCPi0WOWeightInAccPt;         //! array of histos without weight pi0 in acceptance, pT
    TH1F**                  fHistoMCPi0WOEvtWeightInAccPt;      //! array of histos without evt weight pi0 in acceptance, pT
    TH1F**                  fHistoMCEtaInAccPt;                 //! array of histos with weighted eta in acceptance, pT
    TH1F**                  fHistoMCEtaWOWeightInAccPt;         //! array of histos without weight eta in acceptance, pT
    TH1F**                  fHistoMCEtaWOEvtWeightInAccPt;      //! array of histos without evt weight eta in acceptance, pT
    TH2F**                  fHistoMCPi0PtY;                     //! array of histos with weighted pi0, pT, Y
    TH2F**                  fHistoMCEtaPtY;                     //! array of histos with weighted eta, pT, Y
    TH2F**                  fHistoMCPi0PtAlpha;                 //! array of histos with weighted pi0, pT, alpha
    TH2F**                  fHistoMCEtaPtAlpha;                 //! array of histos with weighted eta, pT, alpha
    TH2F**                  fHistoMCPrimaryPtvsSource;          //! array of histos with weighted primary particles, pT vs source
    TH2F**                  fHistoMCSecPi0PtvsSource;           //! array of histos with secondary pi0, pT, source
    TH1F**                  fHistoMCSecPi0Source;               //! array of histos with secondary pi0, source
    TH2F**                  fHistoMCSecPi0InAccPtvsSource;      //! array of histos with secondary pi0, pT, source in acceptance
    TH1F**                  fHistoMCSecEtaPt;                   //! array of histos with secondary eta, pT
    TH1F**                  fHistoMCSecEtaSource;               //! array of histos with secondary eta, source
    TH2F**                  fHistoMCPi0PtJetPt;                 //! array of histos with weighted pi0, pT, hardest jet pt
    TH2F**                  fHistoMCEtaPtJetPt;                 //! array of histos with weighted eta, pT, hardest jet pt

    // histograms for pure MC quantites for dir gamma analysis part
    TH2F**                  fHistoMCPi0PtGammaLeg;              //! array of histos with weighted pi0, pT gamma, different legs
    TH2F**                  fHistoMCPi0WOWeightPtGammaLeg;      //! array of histos with unweighted pi0, pT gamma, different legs
    TH2F**                  fHistoMCPi0InAccPtGammaLeg;         //! array of histos with weighted pi0, pT gamma, different legs
    TH2F**                  fHistoMCPi0WOWeightInAccPtGammaLeg; //! array of histos with unweighted pi0, pT gamma, different legs
    TH2F**                  fHistoMCSecPi0PtGamma1vsSource;     //! array of histos with secondary pi0, pT gamma 1, source
    TH2F**                  fHistoMCSecPi0InAccPtGamma1vsSource;//! array of histos with secondary pi0, pT gamma 1, source in acceptance
    TH2F**                  fHistoMCSecPi0PtGamma2vsSource;     //! array of histos with secondary pi0, pT gamma 2, source
    TH2F**                  fHistoMCSecPi0InAccPtGamma2vsSource;//! array of histos with secondary pi0, pT gamma 2, source in acceptance

    // MC validated reconstructed quantities mesons
    TH2F**                  fHistoTruePi0InvMassPt;                             //! array of histos with validated pi0, invMass, pt
    TH2F**                  fHistoTrueEtaInvMassPt;                             //! array of histos with validated eta, invMass, pt
    TH2F**                  fHistoTruePi0InvMassPtAdditional;                             //! array of histos with validated pi0, invMass, pt
    TH2F**                  fHistoTrueEtaInvMassPtAdditional;                             //! array of histos with validated eta, invMass, pt
    TH2F**                  fHistoTruePi0MatchedInvMassPt;                      //! array of histos with rejected pi0, invMass, pt
    TH2F**                  fHistoTrueEtaMatchedInvMassPt;                      //! array of histos with rejected eta, invMass, pt
    TH2F**                  fHistoTruePi0CaloPhotonInvMassPt;                   //! array of histos with validated pi0, photon leading, invMass, pt
    TH2F**                  fHistoTrueEtaCaloPhotonInvMassPt;                   //! array of histos with validated eta, photon leading, invMass, pt
    TH2F**                  fHistoTruePi0CaloConvertedPhotonInvMassPt;          //! array of histos with validated pi0, converted photon leading, invMass, pt
    TH2F**                  fHistoTruePi0CaloConvertedPhotonMatchedInvMassPt;   //! array of histos with validated pi0 matched with conv photon, converted photon leading, invMass, pt
    TH2F**                  fHistoTruePi0CaloConvPhotonConvRPt;                 //!
    TH2F**                  fHistoTruePi0CaloConvPhotonConvRAlpha;              //!
    TH2F**                  fHistoTruePi0CaloConvPhotonPtAlpha;                 //!
    TH2F**                  fHistoTrueEtaCaloConvertedPhotonInvMassPt;          //! array of histos with validated eta, converted photon leading, invMass, pt
    TH2F**                  fHistoTrueEtaCaloConvertedPhotonMatchedInvMassPt;   //! array of histos with validated eta matched with conv photon, converted photon leading, invMass, pt
    TH2F**                  fHistoTrueEtaCaloConvPhotonConvRPt;                 //!
    TH2F**                  fHistoTrueEtaCaloConvPhotonConvRAlpha;              //!
    TH2F**                  fHistoTrueEtaCaloConvPhotonPtAlpha;                 //!
    TH2F**                  fHistoTruePi0CaloElectronInvMassPt;                 //! array of histos with validated mothers, electron leading, invMass, pt
    TH2F**                  fHistoTrueEtaCaloElectronInvMassPt;                 //! array of histos with validated mothers, electron leading, invMass, pt
    TH2F**                  fHistoTruePi0CaloMergedClusterInvMassPt;            //! array of histos with validated mothers, merged cluster invMass, pt
    TH2F**                  fHistoTrueEtaCaloMergedClusterInvMassPt;            //! array of histos with validated mothers, merged cluster invMass, pt
    TH2F**                  fHistoTrueMotherCaloEMNonLeadingInvMassPt;          //! array of histos with validated mothers, EM non leading, invMass, pt
    TH2F**                  fHistoTruePi0CaloMergedClusterPartConvInvMassPt;    //! array of histos with validated mothers, merged cluster part conv, invMass, pt
    TH2F**                  fHistoTrueEtaCaloMergedClusterPartConvInvMassPt;    //! array of histos with validated mothers, merged cluster part conv, invMass, pt
    TH2F**                  fHistoTruePrimaryPi0InvMassPt;                      //! array of histos with validated weighted primary mothers, invMass, pt
    TH2F**                  fHistoTruePrimaryEtaInvMassPt;                      //! array of histos with validated weighted primary mothers, invMass, pt
    TH2F**                  fHistoTruePrimaryPi0W0WeightingInvMassPt;           //! array of histos with validated unweighted primary mothers, invMass, pt
    TH2F**                  fHistoTruePrimaryEtaW0WeightingInvMassPt;           //! array of histos with validated unweighted primary mothers, invMass, pt
    TProfile2D**            fProfileTruePrimaryPi0WeightsInvMassPt;             //! array of profiles with weights for validated primary mothers, invMass, pt
    TProfile2D**            fProfileTruePrimaryEtaWeightsInvMassPt;             //! array of profiles with weights for validated primary mothers, invMass, pt
    TH2F**                  fHistoTruePrimaryPi0MCPtResolPt;                    //! array of histos with validated weighted primary pi0, MCpt, resol pt
    TH2F**                  fHistoTruePrimaryEtaMCPtResolPt;                    //! array of histos with validated weighted primary eta, MCpt, resol pt
    TH2F**                  fHistoTrueMotherPi0ConvPhotonEtaPhi;                //! array of histograms with invariant mass cut of 0.05 && pi0cand->M() < 0.17 ,eta/phi of conversion photon
    TH2F**                  fHistoTrueMotherEtaConvPhotonEtaPhi;                //! array of histograms with invariant mass cut of 0.45 && pi0cand->M() < 0.65 ,eta/phi of conversion photon
    TH2F**                  fHistoTrueSecondaryPi0InvMassPt;                    //! array of histos with validated secondary mothers, invMass, pt
    TH2F**                  fHistoTrueSecondaryPi0FromK0sInvMassPt;             //! array of histos with validated secondary mothers from K0s, invMass, pt
    TH1F**                  fHistoTrueK0sWithPi0DaughterMCPt;                   //! array of histos with K0s with reconstructed pi0 as daughter, pt
    TH2F**                  fHistoTrueSecondaryPi0FromK0lInvMassPt;             //! array of histos with validated secondary mothers from K0l, invMass, pt
    TH1F**                  fHistoTrueK0lWithPi0DaughterMCPt;                   //! array of histos with K0l with reconstructed pi0 as daughter, pt
    TH2F**                  fHistoTrueSecondaryPi0FromEtaInvMassPt;             //! array of histos with validated secondary mothers from eta, invMass, pt
    TH1F**                  fHistoTrueEtaWithPi0DaughterMCPt;                   //! array of histos with eta with reconstructed pi0 as daughter, pt
    TH2F**                  fHistoTrueSecondaryPi0FromLambdaInvMassPt;          //! array of histos with validated secondary mothers from Lambda, invMass, pt
    TH1F**                  fHistoTrueLambdaWithPi0DaughterMCPt;                //! array of histos with lambda with reconstructed pi0 as daughter, pt
    TH2F**                  fHistoTrueBckGGInvMassPt;                           //! array of histos with pure gamma gamma combinatorial BG, invMass, pt
    TH2F**                  fHistoTrueBckFullMesonContainedInOneClusterInvMassPt; //! array of histos with pi0 fully contained in one cluster, invMass, pt
    TH2F**                  fHistoTrueBckAsymEClustersInvMassPt;                  //! array of histos with asymmetry energy distributions of clusters, invMass, pt
    TH2F**                  fHistoTrueBckContInvMassPt;                         //! array of histos with contamination BG, invMass, pt
    TH2F**                  fHistoTruePi0PtY;                                   //! array of histos with validated pi0, pt, Y
    TH2F**                  fHistoTrueEtaPtY;                                   //! array of histos with validated eta, pt, Y
    TH2F**                  fHistoTruePi0PtAlpha;                               //! array of histos with validated pi0, pt, alpha
    TH2F**                  fHistoTrueEtaPtAlpha;                               //! array of histos with validated eta, pt, alpha
    TH2F**                  fHistoTruePi0PtOpenAngle;                           //! array of histos with validated pi0, pt, openAngle
    TH2F**                  fHistoTrueEtaPtOpenAngle;                           //! array of histos with validated eta, pt, openAngle

    // MC validated reconstructed quantities photons
    TH1F**                  fHistoTrueConvGammaPt;                              //! array of histos with validated conversion photon, pt
    TH2F**                  fHistoCombinatorialPt;                              //! array of histos with combinatorial BG, pt, source
    TH1F**                  fHistoTruePrimaryConvGammaPt;                       //! array of histos with validated primary
    TH2F**                  fHistoTruePrimaryConvGammaESDPtMCPt;                //! array of histos with validated primary conversion photon, rec pt, mc pt conversion photon, pt
    TH2F**                  fHistoTrueSecondaryConvGammaPt;                     //! array of histos with validated secondary conversion photon, rec pt for different sources
    TH2F**                  fHistoTrueSecondaryConvGammaMCPt;                   //! array of histos with validated secondary conversion photon, MC pt  for different sources
    TH2F**                  fHistoTrueSecondaryConvGammaFromXFromK0sMCPtESDPt;  //! array of histos with validated secondary conversion photon from K0s, rec pt, mc pt conversion photon, pt
    TH2F**                  fHistoTrueSecondaryConvGammaFromXFromK0lMCPtESDPt;  //! array of histos with validated secondary conversion photon from K0l, rec pt, mc pt conversion photon, pt
    TH2F**                  fHistoTrueSecondaryConvGammaFromXFromLambdaMCPtESDPt; //! array of histos with validated secondary conversion photon from Lambda, rec pt, mc pt conversion photon, pt

    TH1F**                  fHistoTrueClusGammaPt;                              //! array of histos with validated cluster (electron or photon), pt
    TH1F**                  fHistoTrueClusElectronPt;                           //! array of histos with validated electron, pt
    TH1F**                  fHistoTrueClusConvGammaPt;                          //! array of histos with validated converted photon, pt
    TH1F**                  fHistoTrueClusConvGammaFullyPt;                     //! array of histos with validated converted photon, fully contained, pt
    TH1F**                  fHistoTrueClusMergedGammaPt;                        //! array of histos with validated merged photons, electrons, dalitz, pt
    TH1F**                  fHistoTrueClusMergedPartConvGammaPt;                //! array of histos with validated merged partially converted photons, pt
    TH1F**                  fHistoTrueClusDalitzPt;                             //! array of histos with validated Dalitz decay, pt
    TH1F**                  fHistoTrueClusDalitzMergedPt;                       //! array of histos with validated Dalitz decay, more than one decay product in cluster, pt
    TH1F**                  fHistoTrueClusPhotonFromElecMotherPt;               //! array of histos with validated photon from electron, pt
    TH1F**                  fHistoTrueClusShowerPt;                             //! array of histos with validated shower, pt
    TH1F**                  fHistoTrueClusSubLeadingPt;                         //! array of histos with pi0/eta/eta_prime in subleading contribution
    TH1F**                  fHistoTrueClusNMothers;                             //! array of histos with number of different particles (pi0/eta/eta_prime) contributing to cluster
    TH1F**                  fHistoTrueClusEMNonLeadingPt;                       //! array of histos with cluster with largest energy by hadron
    TH2F**                  fHistoTrueNLabelsInClusPt;                          //! array of histos with number of labels in cluster
    TH1F**                  fHistoTruePrimaryClusGammaPt;                       //! array of histos with validated primary cluster, photons, pt
    TH2F**                  fHistoTruePrimaryClusGammaESDPtMCPt;                //! array of histos with validated primary cluster, photons, rec Pt, MC pt
    TH1F**                  fHistoTruePrimaryClusConvGammaPt;                   //! array of histos with validated primary cluster, converted photons, pt
    TH2F**                  fHistoTruePrimaryClusConvGammaESDPtMCPt;            //! array of histos with validated primary cluster, converted photons, rec Pt, MC pt
    TH1F**                  fHistoTrueSecondaryClusGammaPt;                     //! array of histos with validated secondary cluster, photons, pt
    TH1F**                  fHistoTrueSecondaryClusGammaFromK0sPt;              //! array of histos with validated secondary cluster from K0s, photons, pt
    TH1F**                  fHistoTrueSecondaryClusGammaFromK0lPt;              //! array of histos with validated secondary cluster from K0l, photons, pt
    TH1F**                  fHistoTrueSecondaryClusGammaFromLambdaPt;           //! array of histos with validated secondary cluster from Lambda, photons, pt

    // MC validated reconstructed quantities for direct photon reco
    TH2F**                  fHistoTruePrimaryPi0PhotonPairPtconv;               //! array of histos with validated primary pi0's vs conversion photon pT
    TH2F**                  fHistoTruePrimaryPi0W0WeightsPhotonPairPtconv;      //! array of histos with validated primary pi0's vs conversion photon pT
    TH1F**                  fHistoTruePrimaryPi0DCPtconv;                       //! array of histos with validated primary pi0's vs conversion photon pT, double counting
    TH1F**                  fHistoTruePrimaryPi0MissingPtconv;                  //! array of histos with validated primary pi0's vs conversion photon pT, missing
    TH2F**                  fHistoTrueSecondaryPi0PhotonPairPtconv;             //! array of histos with validated secondary pi0's vs conversion photon pT
    TH2F**                  fHistoTrueSecondaryPi0FromK0sPhotonPairPtconv;      //! array of histos with validated secondary pi0's from K0s vs conversion photon pT
    TH2F**                  fHistoTrueSecondaryPi0FromK0lPhotonPairPtconv;      //! array of histos with validated secondary pi0's from K0l vs conversion photon pT
    TH2F**                  fHistoTrueSecondaryPi0FromLambdaPhotonPairPtconv;   //! array of histos with validated secondary pi0's from Lambda vs conversion photon pT
    TH2F**                  fHistoTrueSecondaryPi0DCPtconvSource;               //! array of histos with validated secondary pi0's vs conversion photon pT, double counting vs source
    TH2F**                  fHistoTrueSecondaryPi0MissingPtconvSource;          //! array of histos with validated secondary pi0's vs conversion photon pT, missing vs source
    TH2F**                  fHistoTrueGammaInvMassPt;                           //! array of histos with Inv mass of validated gammas
    TH2F**                  fHistoMCGammaHBTInvMassPt;                          //! array of histos with Inv mass of generated MC gammas
    TH2F**                  fHistoBckTrueGammaInvMassPt;                        //! array of histos with Bck Inv mass of validated gammas
    TH2F**                  fHistoBckMCGammaHBTInvMassPt;                       //! array of histos with Bck Inv mass of generated MC gammas
    TH2F**                  fHistoHBTOpeningAnglePt;                            //! array of histos with Bck Inv mass of generated MC gammas
    TH2F**                  fHistoBckHBTOpeningAnglePt;                         //! array of histos with Bck Inv mass of generated MC gammas
    TH2F**                  fHistoHBTDeltaEPt;                                  //! array of histos with Bck Inv mass of generated MC gammas
    TH2F**                  fHistoBckHBTDeltaEPt;                               //! array of histos with Bck Inv mass of generated MC gammas
    TH2F**                  fHistoMCGammaHBTOpeningAnglePt;                     //! array of histos with Bck Inv mass of generated MC gammas
    TH2F**                  fHistoBckMCGammaHBTOpeningAnglePt;                  //! array of histos with Bck Inv mass of generated MC gammas
    TH2F**                  fHistoMCGammaHBTDeltaEPt;                           //! array of histos with Bck Inv mass of generated MC gammas
    TH2F**                  fHistoBckMCGammaHBTDeltaEPt;                        //! array of histos with Bck Inv mass of generated MC gammas

    TH1F**                 fHistoPtJet;                                          //! Histogram of Jet Pt
    TH1F**                 fHistoJetEta;                                         //! Histogram of Jet Eta
    TH1F**                 fHistoJetPhi;                                         //! Histogram of Jet Phi
    TH1F**                 fHistoJetArea;                                        //! Histogram of Jet Area
    TH1F**                 fHistoNJets;                                          //! Histogram of number of jets
    TH1F**                 fHistoEventwJets;                                     //! Histogram of number of events with jets > 0
    TH1F**                 fHistoJetPi0PtRatio;                                  //! Histogram of PtPi0/PtJet
    TH1F**                 fHistoDoubleCounting;                                 //! Histogram if NM candidates are defined within multiple jets
    TH2F**                 fHistoJetMotherInvMassPt;                             //! Histogram of NM candidates with a jet in the event
    TH2F**                 fHistoPi0InJetMotherInvMassPt;                        //! Histogram of NM candidates that are inside a jet
    TH2F**                 fHistoMotherBackJetInvMassPt;                         //! Histogram of Backgrouns candidates that are involved with jets
    TH2F**                 fHistoRJetPi0Cand;                                    //! Histogram of RJetPi0Cand vs Pt
    TH2F**                 fHistoEtaPhiJetPi0Cand;                               //! Histogram of delta eta and delta phi distr between jet and NM candidates
    TH2F**                 fHistoEtaPhiJetWithPi0Cand;                           //! Histogram of delta eta and delta phi distr when pi0 is inside a jet
    TH2F**                 fHistoJetFragmFunc;                                   //! Histogram to determine fragmentation function
    TH2F**                 fHistoJetFragmFuncZInvMass;                           //! Histogram of Inv Mass distribution with z
    TH2F**                 fHistoTruevsRecJetPt;                                 //! Histogram of true jet pt vs reconstructed jet pt
    TH2F**                 fHistoTruePi0JetMotherInvMassPt;                      //! Histogram of true pi0s in an event with a jet
    TH2F**                 fHistoTruePi0InJetMotherInvMassPt;                    //! Histogram of true pi0s in a jet
    TH2F**                 fHistoTruePrimaryPi0JetInvMassPt;                     //! Histogram of true primary pi0s in an event with a jet
    TH2F**                 fHistoTruePrimaryPi0inJetInvMassPt;                   //! Histogram of true primary pi0s in a jet
    TH2F**                 fHistoTruePrimaryPi0InJetInvMassTruePt;               //! Histogram of true primary pi0s in a jet with their true pt
    TH1F**                 fHistoTrueDoubleCountingPi0Jet;                       //! Histogram of when a true pi0 is defined to be in multiple jets
    TH2F**                 fHistoTrueEtaJetMotherInvMassPt;                      //! Histogram of true etas in an event with a jet
    TH2F**                 fHistoTrueEtaInJetMotherInvMassPt;                    //! Histogram of true etas in a jet
    TH2F**                 fHistoTruePrimaryEtaJetInvMassPt;                     //! Histogram of true primary etas in an event with a jet
    TH2F**                 fHistoTruePrimaryEtainJetInvMassPt;                   //! Histogram of true primary etas in a jet
    TH1F**                 fHistoTrueDoubleCountingEtaJet;                       //! Histogram of when a true eta is defined to be in multiple jets
    TH2F**                 fHistoTruePi0JetFragmFunc;                            //! Histogram to determine true pi0 fragmentation function
    TH2F**                 fHistoTruePi0JetFragmFuncZInvMass;                    //! Histogram to determine true pi0 Inv Mass distribution with z
    TH2F**                 fHistoTrueEtaJetFragmFunc;                            //! Histogram to determine true eta fragmentation function
    TH2F**                 fHistoTrueEtaJetFragmFuncZInvMass;                    //! Histogram to determine true eta Inv Mass distribution with z
    TH1F**                 fHistoMCPi0JetInAccPt;                                //! Histogram with weighted pi0 in a jet event in acceptance, pT
    TH1F**                 fHistoMCPi0inJetInAccPt;                              //! Histogram with weighted pi0 in a jet in acceptance, pT
    TH1F**                 fHistoMCEtaJetInAccPt;                                //! Histogram with weighted eta in a jet event in acceptance, pT
    TH1F**                 fHistoMCEtainJetInAccPt;                              //! Histogram with weighted eta in a jet in acceptance, pT
    TH1F**                 fHistoMCPi0JetEventGenerated;                         //! Histogram with mesons in a jet event generated, pT
    TH1F**                 fHistoMCPi0inJetGenerated;                            //! Histogram with mesons in a jet generated, pT
    TH1F**                 fHistoMCEtaJetEventGenerated;                         //! Histogram with mesons in a jet event generated, pT
    TH1F**                 fHistoMCEtainJetGenerated;                            //! Histogram with mesons in a jet generated, pT
    TH2F**                 fHistoTrueSecondaryPi0FromK0sJetInvMassPt;            //! Histogram with validated secondary mothers from K0s in an event with a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0FromK0sinJetInvMassPt;          //! Histogram with validated secondary mothers from K0s in a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0FromLambdaJetInvMassPt;         //! Histogram with validated secondary mothers from lambda in an event with a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0FromLambdainJetInvMassPt;       //! Histogram with validated secondary mothers from lambda in a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0FromK0lJetInvMassPt;            //! Histogram with validated secondary mothers from K0l in an event with a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0FromK0linJetInvMassPt;          //! Histogram with validated secondary mothers from K0l in a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0InvJetMassPt;                   //! Histogram with validated secondary mothers in an event with a jet, invMass, pt
    TH2F**                 fHistoTrueSecondaryPi0InvinJetMassPt;                 //! Histogram with validated secondary mothers in a jet, invMass, pt
    TH2F**                 fHistoMotherPi0inJetPtY;                              //! Histogram with the rapidity of the validated pi0s in jets
    TH2F**                 fHistoMotherEtainJetPtY;                              //! Histogram with the rapidity of the validated etas in jets
    TH2F**                 fHistoMotherPi0inJetPtPhi;                            //! Histogram with the phi of the validated pi0s in jets
    TH2F**                 fHistoMotherEtainJetPtPhi;                            //! Histogram with the phi of the validated etas in jets
    TH2F**                 fHistoUnfoldingAsData;                                //! Histogram to use for jet pi0 unfolding
    TH2F**                 fHistoUnfoldingMissed;                                //! Histogram to use for jet pi0 unfolding
    TH2F**                 fHistoUnfoldingReject;                                //! Histogram to use for jet pi0 unfolding
    TH2F**                 fHistoUnfoldingAsDataInvMassZ;                        //! Histogram to use for jet pi0 unfolding for z inmass
    TH2F**                 fHistoUnfoldingMissedInvMassZ;                        //! Histogram to use for jet pi0 unfolding for z inmass
    TH2F**                 fHistoUnfoldingRejectInvMassZ;                        //! Histogram to use for jet pi0 unfolding for z inmass

    vector<Double_t>      fVectorJetPt;                                         //! Vector of JetPt
    vector<Double_t>      fVectorJetPx;                                         //! Vector of JetPx
    vector<Double_t>      fVectorJetPy;                                         //! Vector of JetPy
    vector<Double_t>      fVectorJetPz;                                         //! Vector of JetPz
    vector<Double_t>      fVectorJetEta;                                        //! Vector of JetEta
    vector<Double_t>      fVectorJetPhi;                                        //! Vector of JetPhi
    vector<Double_t>      fVectorJetArea;                                       //! Vector of JetArea
    vector<Double_t>      fTrueVectorJetPt;                                     //! Vector of True JetPt
    vector<Double_t>      fTrueVectorJetPx;                                     //! Vector of True JetPx
    vector<Double_t>      fTrueVectorJetPy;                                     //! Vector of True JetPy
    vector<Double_t>      fTrueVectorJetPz;                                     //! Vector of True JetPz
    vector<Double_t>      fTrueVectorJetEta;                                    //! Vector of True JetEta
    vector<Double_t>      fTrueVectorJetPhi;                                    //! Vector of True JetPhi

    // variable to keep track of multiple & missing reco
    vector<Int_t>           fVectorRecTruePi0s;                                 //! array of strings containing the stack position of the reconstructed validated pi0
    vector<Int_t>           fVectorRecTrueEtas;                                 //! array of strings containing the stack position of the reconstructed validated eta
    TH2F**                  fHistoDoubleCountTruePi0InvMassPt;                  //! array of histos with double counted pi0s, invMass, pT
    TH2F**                  fHistoDoubleCountTrueEtaInvMassPt;                  //! array of histos with double counted etas, invMass, pT
    TH2F**                  fHistoDoubleCountTrueConvGammaRPt;                  //! array of histos with double counted photons, R, pT
    TH2F**                  fHistoDoubleCountTrueClusterGammaPt;                //! array of histos with double counted cluster photons
    vector<Int_t>           fVectorDoubleCountTruePi0s;                         //! vector containing labels of validated pi0
    vector<Int_t>           fVectorDoubleCountTrueEtas;                         //! vector containing labels of validated eta
    vector<Int_t>           fVectorDoubleCountTrueConvGammas;                   //! vector containing labels of validated photons
    vector<Int_t>           fVectorDoubleCountTrueClusterGammas;                //! vector containing labels of validated cluster photons
    TH1F**                  fHistoMultipleCountTruePi0;                         //! array of histos how often TruePi0s are counted
    TH1F**                  fHistoMultipleCountTrueEta;                         //! array of histos how often TrueEtas are counted
    TH1F**                  fHistoMultipleCountTrueConvGamma;                   //! array of histos how often TrueConvGammas are counted
    TH1F**                  fHistoMultipleCountTrueClusterGamma;                //! array of histos how often TrueClusterGammas are counted
    map<Int_t,Int_t>        fMapMultipleCountTruePi0s;                          //! map containing pi0 labels that are counted at least twice
    map<Int_t,Int_t>        fMapMultipleCountTrueEtas;                          //! map containing eta labels that are counted at least twice
    map<Int_t,Int_t>        fMapMultipleCountTrueConvGammas;                    //! map containing photon labels that are counted at least twice
    map<Int_t,Int_t>        fMapMultipleCountTrueClusterGammas;                 //! map containing cluster photon labels that are counted at least twice

    // further QA and calibration hist for Calo clusters
    TH2F**                  fHistoTrueClusGammaEM02;                            //! array of histos with TruePhotons: cluster E vs M02
    TH2F**                  fHistoTrueClusPi0EM02;                              //! array of histos with TruePi0s: cluster E vs M02
    TH2F**                  fHistoTruePi0InvMassECalib;                         //! array of histogram with pure pi0 signal inv Mass, energy of cluster
    TH2F**                  fHistoTruePi0PureGammaInvMassECalib;                //! array of histogram with pure pi0 signal (only pure gammas) inv Mass, energy of cluster
    TH2F**                  fHistoTruePi0InvMassECalibPCM;                         //! array of histogram with pure pi0 signal inv Mass, energy of PCM

    // event histograms
    TH1F**                  fHistoNEvents;                                      //! array of histos with event information
    TH1F**                  fHistoNEventsWOWeight;                              //! array of histos with event information without event weights
    TH1F**                  fHistoNGoodESDTracks;                               //! array of histos with number of good tracks (2010 Standard track cuts)
    TH1F**                  fHistoVertexZ;                                      //! array of histos with vertex z distribution for selected events
    TH1F**                  fHistoVertexX;                                      //! array of histos with vertex x distribution for selected events
    TH1F**                  fHistoVertexY;                                      //! array of histos with vertex y distribution for selected events
    TH1F**                  fHistoNGammaCandidates;                             //! array of histos with number of gamma candidates per event
    TH2F**                  fHistoNGoodESDTracksVsNGammaCandidates;             //! array of histos with number of good tracks vs gamma candidates
    TH2F**                  fHistoSPDClusterTrackletBackground;                 //! array of histos with SPD tracklets vs SPD clusters for background rejection
    TH1F**                  fHistoNV0Tracks;                                    //! array of histos with V0 counts
    TProfile**              fProfileEtaShift;                                   //! array of profiles with eta shift
    TProfile**              fProfileJetJetXSection;                             //! array of profiles with xsection for jetjet
    TH1F**                  fHistoJetJetNTrials;                                //! array of histos with ntrials for jetjet
    TH1F**                  fHistoEventSphericity;                              //! array of histos with event Sphericity

    // additional variables
    Double_t                fEventPlaneAngle;                                   //! EventPlaneAngle
    TTree**                 tESDMesonsInvMassPtDcazMinDcazMaxFlag;              //!
    Float_t                 fInvMass;                                           //!
    Float_t                 fPt;                                                //!
    Float_t                 fDCAzGammaMin;                                      //!
    Float_t                 fDCAzGammaMax;                                      //!
    UChar_t                 iFlag;                                              //!
    UChar_t                 iMesonMCInfo;                                       //!
    TRandom3                fRandom;                                            //! random
    Int_t                   fNGammaCandidates;                                  //!
    Double_t*               fUnsmearedPx;                                       //!
    Double_t*               fUnsmearedPy;                                       //!
    Double_t*               fUnsmearedPz;                                       //!
    Double_t*               fUnsmearedE;                                        //!
    Int_t*                  fMCEventPos;                                        //!
    Int_t*                  fMCEventNeg;                                        //!
    Int_t*                  fESDArrayPos;                                       //!
    Int_t*                  fESDArrayNeg;                                       //!
    Int_t                   fnCuts;                                             // number of cuts to be analysed in parallel
    Int_t                   fiCut;                                              //! current cut
    Bool_t                  fMoveParticleAccordingToVertex;                     // boolean for BG calculation
    Int_t                   fIsHeavyIon;                                        // switch for pp = 0, PbPb = 1, pPb = 2
    Int_t                   fDoLightOutput;                                     // switch for running light output, 0 -> normal mode, 1 -> light mode, 2 -> minimum
    Int_t                   fPi0EtaSwitch;                                     // switch for running only pi0, oly eta or both, 0 -> both, 1 -> pi0, 2 -> eta
    Bool_t                  fDoECalibOutput;                                    // switch for running with E-Calib Histograms in Light Output, kFALSE -> no E-Calib Histograms, kTRUE -> with E-Calib Histograms
    Bool_t                  fDoMesonAnalysis;                                   // flag for meson analysis
    Int_t                   fDoMesonQA;                                         // flag for meson QA
    Int_t                   fDoPhotonQA;                                        // flag for photon QA
    Int_t                   fDoClusterQA;                                       // flag for cluster QA
    Bool_t                  fIsFromDesiredHeader;                               //! flag for MC headers
    Bool_t                  fIsOverlappingWithOtherHeader;                      //! flag for particles in MC overlapping between headers
    Int_t                   fIsMC;                                              // flag for MC information
    Bool_t                  fDoTHnSparse;                                       // flag for using THnSparses for background estimation
    Bool_t                  fSetPlotHistsExtQA;                                 // flag for extended QA hists
    Double_t                fWeightJetJetMC;                                    //! weight for Jet-Jet MC
    Bool_t                  fDoConvGammaShowerShapeTree;                        // flag for tree with conv gamma R vs energy vs shower shape
    Bool_t                  fEnableSortForClusMC;                               // switch on sorting for MC labels in cluster
    Bool_t                  fDoPrimaryTrackMatching;                            // switch for basic track matching for primaries
    Bool_t                  fDoInvMassShowerShapeTree;                          // flag for producing tree tESDInvMassShowerShape
    Bool_t                  fDoMaterialBudgetWeightingOfGammasForTrueMesons;    // flag to set the usage of material budget weights
    TTree*                  tBrokenFiles;                                       //! tree for keeping track of broken files
    TObjString*             fFileNameBroken;                                    //! string object for broken file name
    Bool_t                  fAllowOverlapHeaders;                               // enable overlapping headers for cluster selection
    Int_t                   fTrackMatcherRunningMode;                           // CaloTrackMatcher running mode
    Bool_t                  fDoHBTHistoOutput;                                  // switch for additional HBT output

  private:
    AliAnalysisTaskGammaConvCalo(const AliAnalysisTaskGammaConvCalo&); // Prevent copy-construction
    AliAnalysisTaskGammaConvCalo &operator=(const AliAnalysisTaskGammaConvCalo&); // Prevent assignment

    ClassDef(AliAnalysisTaskGammaConvCalo, 68);
};

#endif
