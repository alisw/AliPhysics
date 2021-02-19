#ifndef ALIANLYSISTASKGAMMACALOMERGED_cxx
#define ALIANLYSISTASKGAMMACALOMERGED_cxx

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
#include <vector>
#include <map>

class AliAnalysisTaskGammaCaloMerged : public AliAnalysisTaskSE {
  public:

    AliAnalysisTaskGammaCaloMerged();
    AliAnalysisTaskGammaCaloMerged(const char *name);
    virtual ~AliAnalysisTaskGammaCaloMerged();

    virtual void   UserCreateOutputObjects();
    virtual Bool_t Notify();
    virtual void   UserExec(Option_t *);
    virtual void   Terminate(const Option_t*);

    void SetV0ReaderName(TString name)                    {
                                                            fV0ReaderName               = name                                                            ;
                                                            return                                                                                        ;
                                                          }
    void SetIsHeavyIon(Int_t flag)                        { fIsHeavyIon                 = flag                                                            ; }

    // Function to set correction task setting
    void SetCorrectionTaskSetting(TString setting)        { fCorrTaskSetting = setting                                                                    ; }
    // base functions for selecting photon and meson candidates in reconstructed data
    void ProcessClusters();

    // MC functions
    void SetIsMC(Int_t isMC)                              { fIsMC                       = isMC                                                            ; }
    void SetSelectedMesonID(Int_t anaMeson)               { fSelectedMesonID            = anaMeson                                                        ; }

    void ProcessMCParticles();
    void ProcessAODMCParticles();

    // Jet functions
    void ProcessJets();

    // determine source according to pdg code of mother
    Int_t GetSourceClassification(Int_t daughter, Int_t pdgCode);

    void ProcessTrueClusterCandidates( AliAODConversionPhoton* TruePhotonCandidate, AliVCluster* cluster, AliAODConversionPhoton *TrueSubClusterCandidate1,
                                    AliAODConversionPhoton *TrueSubClusterCandidate2);
    void ProcessTrueClusterCandidatesAOD( AliAODConversionPhoton* TruePhotonCandidate, AliVCluster* cluster, AliAODConversionPhoton *TrueSubClusterCandidate1,
                                    AliAODConversionPhoton *TrueSubClusterCandidate2);
// //     void ProcessTrueMesonCandidates( AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate0, AliAODConversionPhoton *TrueGammaCandidate1);

    // switches for additional analysis streams or outputs
    void SetLightOutput( Bool_t flag )                    { fDoLightOutput              = flag;}
    void SetDoMesonQA(Int_t flag)                         { fDoMesonQA                  = flag                                                            ; }
    void SetDoClusterQA(Int_t flag)                       { fDoClusterQA                = flag                                                            ; }
    void SetPlotHistsExtQA(Bool_t flag)                   { fSetPlotHistsExtQA          = flag                                                            ; }

    // switch for overlap criterium
    void SetOverlapFromCluster(Bool_t flag)               { fMesonOLFromCluster = flag                                                                    ; }
      // Setting the cut lists for the conversion photons
    void SetEventCutList(Int_t nCuts, TList *CutArray)    {
                                                            fnCuts                      = nCuts                                                           ;
                                                            fEventCutArray              = CutArray                                                        ;
                                                          }

      // Setting the cut lists for the calo photons
    void SetCaloCutList(Int_t nCuts, TList *CutArray)     {
                                                            fnCuts                      = nCuts                                                           ;
                                                            fClusterCutArray            = CutArray                                                        ;
                                                          }

      // Setting the cut lists for the calo photons
    void SetCaloMergedCutList(Int_t nCuts, TList *CutArray){
                                                            fnCuts                      = nCuts                                                           ;
                                                            fClusterMergedCutArray      = CutArray                                                        ;
                                                          }

    // Setting the cut lists for the meson
    void SetMesonCutList(Int_t nCuts, TList *CutArray)      {
                                                              fnCuts                    = nCuts                                                           ;
                                                              fMesonCutArray            = CutArray                                                        ;
                                                            }

    Int_t GetSelectedMesonID()                              { return fSelectedMesonID                                                                     ; }

    // Additional functions for convenience
    void SetLogBinningXTH2(TH2* histoRebin);

    Bool_t CheckVectorForDoubleCount(vector<Int_t> &vec, Int_t tobechecked);
    Bool_t CheckVector1ForEntryAndFillVector2(vector<Int_t> vecCheck, vector<Int_t> &vecFill, Int_t tobechecked);

    void FillMultipleCountMap(map<Int_t,Int_t> &ma, Int_t tobechecked);
    void FillMultipleCountHistoAndClear(map<Int_t,Int_t> &ma, TH1F* hist);

    // Function to enable detailed printouts
    void SetEnableDetailedPrintout(Bool_t enablePO)         { fEnableDetailedPrintOut   = enablePO                                                        ; }

    // Function to enable MC label sorting
    void SetEnableSortingOfMCClusLabels(Bool_t enableSort)  { fEnableSortForClusMC   = enableSort                                                         ; }

    // FUnction to enable very detailed M02 distribution
    void SetEnableDetailedM02Distribtuon(Bool_t enableDetM02) { fDoDetailedM02          = enableDetM02                                                    ; }

    void SetTrackMatcherRunningMode(Int_t mode){fTrackMatcherRunningMode = mode;}

    void SetMinNeutralPionOverlapsMC(Int_t nOverlaps){fMinAllowedPi0OverlapsMC = nOverlaps;}
    void SetMaxNeutralPionOverlapsMC(Int_t nOverlaps){fMaxAllowedPi0OverlapsMC = nOverlaps;}

    Bool_t NumberOfMCEventNeutralPionOverlapInEMCal(AliMCEvent *mcEvent);
    void ProcessNeutralOverlapsMC(AliMCEvent *mcEvent);
    void PrintCaloMCLabelsAndInfoAOD(AliVEvent* event, AliAODConversionPhoton *TrueClusterCandidate, AliVCluster* cluster);

  protected:
    AliV0ReaderV1*          fV0Reader;                                          // basic photon Selection Task
    TString                 fV0ReaderName;
    TString                 fCorrTaskSetting;
    Bool_t                  fDoLightOutput;                                     // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
    AliVEvent*              fInputEvent;                                        // current event
    AliMCEvent*             fMCEvent;                                           // corresponding MC event
    TList**                 fCutFolder;                                         // Array of lists for containers belonging to cut
    TList**                 fESDList;                                           // Array of lists with histograms with reconstructed properties
    TList**                 fTrueList;                                          // Array of lists with histograms with MC validated reconstructed properties
    TList**                 fMCList;                                            // Array of lists with histograms with pure MC information
    TList**                 fHeaderNameList;                                    // Array of lists with header names for MC header selection
    TList*                  fOutputContainer;                                   // Output container
    Int_t                   fNClusterCandidates;                                //! current number of cluster candidates
    Int_t                   fNClusterMergedCandidates;                          //! current number of merged cluster candidates
    TList*                  fEventCutArray;                                     // List with Event Cuts
    AliConvEventCuts*       fEventCuts;                                         // EventCutObject
    TList*                  fClusterCutArray;                                   // List with Cluster Cuts
    TList*                  fClusterMergedCutArray;                             // List with Cluster Cuts for merged clusters
    TList*                  fMesonCutArray;                                     // List with Meson Cuts
    AliConversionMesonCuts* fMesonCuts;                                         // MesonCutObject
    AliAnalysisTaskJetOutlierRemoval*   fOutlierJetReader;                      // JetReader
    AliAnalysisTaskConvJet* fConvJetReader;                                     // JetReader
    Bool_t                  fDoJetAnalysis;                                     // Switch for Jet Analysis
    Bool_t                  fDoJetQA;                                           // Switch for Jet QA
    Int_t                   fDoOutOfJet;                                        // Switch for Jet Analysis out of jet (see AliConversionMesonCuts for detailed description)
    TList**                 fJetHistograms;                                     // List for Jet histograms
    TClonesArray*           fAODMCTrackArray;                                   // MC track array
    TClonesArray*           farrClustersProcess;                                // Cluster array
    std::map<Int_t, Int_t>* fMapNeutralPionOverlap;                             // map with number of neutral pion overlaps for current pion
    Bool_t                  fMesonOLFromCluster;                                // flag to switch between event and cluster classification of overlaps

    //histograms for mesons reconstructed quantities
    TH2F**                  fHistoMotherInvMassPt;                              //! array of histogram with signal + BG for same event photon pairs, inv Mass, pt
    TH2F**                  fHistoMotherPtY;                                    //! array of histograms with signal +BG pt, Y
    TH2F**                  fHistoMotherPtAlpha;                                //! array of histograms with signal +BG  pt, alpha

    // histograms for rec photon clusters/pi0 candidate
    TH1F**                  fHistoClusGammaPt;                                  //! array of histos with cluster, pt
    TH1F**                  fHistoClusGammaE;                                   //! array of histos with cluster, E
    TH1F**                  fHistoClusOverlapHeadersGammaPt;                    //! array of histos with cluster, pt overlapping with other headers
    TH2F**                  fHistoClusNLMPt;                                    //! array of histos with cluster NLM vs Pt
    TH2F**                  fHistoClusMergedPtvsM02;                            //! array of histos with cluster merged, pt vs M02
    TH2F**                  fHistoClusMergedEvsM02;                             //! array of histos with cluster merged, E vs M02
    TH2F**                  fHistoClusMergedPtvsM02Accepted;                    //! array of histos with cluster merged accepted mesons, pt vs M02
    TH2F**                  fHistoClusMergedEvsM02Accepted;                     //! array of histos with cluster merged accepted mesons, E vs M02
    TH2F**                  fHistoClusMergedEvsM20;                             //! array of histos with cluster merged, E vs M20
    TH2F**                  fHistoClusNCellsPt;                                 //! array of histos with cluster NCells vs Pt
    TH2F**                  fHistoClusMergedNCellsPt;                           //! array of histos with merged cluster NCells vs Pt
    TH2F**                  fHistoClusMergedNParticlePt;                        //! array of histos with merged cluster N MC paricles in cluster vs Pt
    TH2F**                  fHistoClusMergedNCellsAroundPt;                     //! array of histos with number of cells surrounding merged cluster vs merged cluster Pt
    TH2F**                  fHistoClusMergedNCellsAroundAndInPt;                //! array of histos with number of cells surrounding merged cluster + Ncells in clus vs merged cluster Pt
    TH2F**                  fHistoClusMergedEAroundE;                           //! array of histos with E surrounding merged cluster vs merged cluster E

    // Jet ESD histograms
    TH1F**                  fHistoPtJet;                                        //! array of histos with Jet pT
    TH1F**                  fHistoJetEta;                                       //! array of histos with Jet eta
    TH1F**                  fHistoJetPhi;                                       //! array of histos with Jet phi
    TH1F**                  fHistoJetArea;                                      //! array of histos with Jet area
    TH2F**                  fClusterEtaPhiJets;                                 //! array of histos with hitmap of clusters in Jet
    TH1F**                  fHistoNJets;                                        //! array of histos with events with Jets
    TH1F**                  fHistoEventwJets;                                   //! array of histos with events with Jets > 0
    TH2F**                  fHistoTruevsRecJetPt;                               //! array of histos Jet rec. vs true pT
    TH2F**                  fHistoClusMergedPtvsRJetAccepted;                   //! array of histos with clusters in Jets vs distance between cluster and Jet axis
    TH2F**                  fHistoJetFragmFunc;                                 //! array of histos with clusters in Jets vs z (z = pJet*pPi0 / |pJet|)

    //histograms for pure MC quantities
    TH1I**                  fHistoMCHeaders;                                    //! array of histos for header names
    TH1F**                  fHistoMCPi0Pt;                                      //! array of histos with weighted pi0, pT
    TH1F**                  fHistoMCPi0ReducedPt;                               //! array of histos with weighted pi0, pT
    TH1F**                  fHistoMCPi0WOWeightPt;                              //! array of histos with unweighted pi0, pT
    TH1F**                  fHistoMCPi0WOEvtWeightPt;                           //! array of histos without event weights pi0, pT
    TH1F**                  fHistoMCEtaPt;                                      //! array of histos with weighted eta, pT
    TH1F**                  fHistoMCEtaWOWeightPt;                              //! array of histos with unweighted eta, pT
    TH1F**                  fHistoMCEtaWOEvtWeightPt;                           //! array of histos without event weights eta, pT
    TH1F**                  fHistoMCPi0DalitzPt;                                //! array of histos with weighted pi0 Dalitz, pT
    TH1F**                  fHistoMCPi0DalitzWOWeightPt;                        //! array of histos with unweighted pi0 Dalitz, pT
    TH1F**                  fHistoMCPi0DalitzWOEvtWeightPt;                     //! array of histos without event weights pi0 Dalitz, pT
    TH1F**                  fHistoMCEtaDalitzPt;                                //! array of histos with weighted eta Dalitz, pT
    TH1F**                  fHistoMCEtaDalitzWOWeightPt;                        //! array of histos with unweighted eta Dalitz, pT
    TH1F**                  fHistoMCEtaDalitzWOEvtWeightPt;                     //! array of histos without event weights eta Dalitz, pT
    TH1F**                  fHistoMCPi0InAccPt;                                 //! array of histos with weighted pi0 in acceptance, pT
    TH1F**                  fHistoMCPi0ReducedInAccPt;                          //! array of histos with weighted pi0 in acceptance, pT
    TH1F**                  fHistoMCEtaInAccPt;                                 //! array of histos with weighted eta in acceptance, pT
    TH1F**                  fHistoMCPi0WOEvtWeightInAccPt;                      //! array of histos without evt weight pi0 in acceptance, pT
    TH1F**                  fHistoMCEtaWOEvtWeightInAccPt;                      //! array of histos without evt weight eta in acceptance, pT
    TH1F**                  fHistoMCPi0DalitzInAccPt;                           //! array of histos with weighted pi0 dalitz in acceptance, pT
    TH1F**                  fHistoMCEtaDalitzInAccPt;                           //! array of histos with weighted eta dalitz in acceptance, pT
    TH1F**                  fHistoMCPi0DalitzWOEvtWeightInAccPt;                //! array of histos without evt weight pi0 in acceptance, pT
    TH1F**                  fHistoMCEtaDalitzWOEvtWeightInAccPt;                //! array of histos without evt weight eta in acceptance, pT

    TH2F**                  fHistoMCSecPi0PtvsSource;                           //! array of histos with weighted pi0 from sec, pT for different sources
    TH2F**                  fHistoMCSecPi0InAccPtvsSource;                      //! array of histos with weighted pi0 from sec in acceptance, pT for different sources
    TH2F**                  fHistoMCPi0PtJetPt;                                 //! array of histos with weighted pi0, pT, hardest jet pt
    TH2F**                  fHistoMCEtaPtJetPt;                                 //! array of histos with weighted eta, pT, hardest jet pt
    TH2F**                  fHistoMCPrimaryPtvsSource;                          //! array of histos with weighted primary particles, pT vs source
    TH2F**                  fHistoMCPrimaryYvsSource;                           //! array of histos with weighted primary particles, Y vs source
    TH1F**                  fHistoMCDecayGammaPt;                               //! array of histos with weighted decay gamma
    TH1F**                  fHistoMCAllGammaPt;                                 //! array of histos with weighted all gamma

    // MC validated cluster histos
    TH2F**                  fHistoTrueClusEFracFirstLabel;                      //!
    TH2F**                  fHistoTrueClusEFracLeadingPi0;                      //!
    TH2F**                  fHistoTrueClusNeutralContamination;                      //!
    TH2F**                  fHistoTrueClusPhotonContamination;                      //!
    TH2F**                  fHistoTrueClusElectronContamination;                      //!
    TH2F**                  fHistoTrueClusOtherContamination;                      //!
    TH2F**                  fHistoTrueClusMergedPtvsM02;                        //!
    TH2F**                  fHistoTrueClusPi0PtvsM02;                           //!
    TH2F**                  fHistoTrueClusMultiplePi0PtvsM02;                   //!
    TH2F**                  fHistoTrueClusPi0DalitzPtvsM02;                     //!
    TH2F**                  fHistoTrueClusPrimPi0PtvsM02;                       //!
    TH2F**                  fHistoTrueClusPrimPi0PtMCPt;                        //!
    TH2F**                  fHistoTruePureMergedClusPrimPi0PtvsM02;             //!
    TH2F**                  fHistoTruePartConvClusPrimPi0PtvsM02;               //!
    TH2F**                  fHistoTrueClusMergedPartConvPi0EVsM20;              //!
    TH2F**                  fHistoTrueClusPrimPi0PureMergedPtMCPt;              //!
    TH2F**                  fHistoTrueClusPrimPi0PartConvMergedPtMCPt;          //!
    TH2F**                  fHistoTrueClusPrimPi01GammaMergedPtMCPt;            //!
    TH2F**                  fHistoTrueClusPrimPi01ElectronMergedPtMCPt;         //!
    TH2F**                  fHistoTrueClusMultiplePrimPi0PtvsM02;               //!
    TH2F**                  fHistoTrueClusSecPi0PtvsM02;                        //!
    TH2F**                  fHistoTrueClusSecPi0FromK0sPtvsM02;                 //!
    TH2F**                  fHistoTrueClusSecPi0FromK0lPtvsM02;                 //!
    TH2F**                  fHistoTrueClusSecPi0FromLambdaPtvsM02;              //!
    TH2F**                  fHistoTrueClusEtaPtvsM02;                           //!
    TH2F**                  fHistoTrueClusEtaDalitzPtvsM02;                     //!
    TH2F**                  fHistoTrueClusMergedPureFromPi0PtvsM02;             //!
    TH2F**                  fHistoTrueClusMergedPureFromEtaPtvsM02;             //!
    TH2F**                  fHistoTrueClusMergedPartConvFromPi0PtvsM02;         //!
    TH2F**                  fHistoTrueClusMergedPartConvFromEtaPtvsM02;         //!
    TH2F**                  fHistoNTrueMultiplePi0vsPt;                         //!
    TH2F**                  fHistoTrueClusGammaFromPi0PtvsM02;                  //!
    TH2F**                  fHistoTrueClusGammaFromEtaPtvsM02;                  //!
    TH2F**                  fHistoTrueClusElectronFromPi0PtvsM02;               //!
    TH2F**                  fHistoTrueClusElectronFromEtaPtvsM02;               //!
    TH2F**                  fHistoTrueSecPi0PtvsDiffReco;                       //!

    TH2F**                  fHistoTrueClusBGPtvsM02;                            //!
    TH2F**                  fHistoTrueClusGammaPtvsM02;                         //!
    TH2F**                  fHistoTrueClusGammaPartConvPtvsM02;                 //!
    TH2F**                  fHistoTrueClusElectronPtvsM02;                      //!
    TH2F**                  fHistoTrueClusElectronFromGammaPtvsM02;             //!
    TH2F**                  fHistoTrueClusMergedInvMassvsPt;                    //!
    TH2F**                  fHistoTrueClusPi0InvMassvsPt;                       //!
    TH2F**                  fHistoTrueClusPrimPi0InvMassvsPt;                   //!
    TH2F**                  fHistoTrueClusEtaInvMassvsPt;                       //!
    TH2F**                  fHistoTrueClusBGInvMassvsPt;                        //!
    TH2F**                  fHistoTrueClusGammaInvMassvsPt;                     //!
    TH2F**                  fHistoTrueClusElectronInvMassvsPt;                  //!
    TH2F**                  fHistoTrueClusBGPtvsSource;                         //!
    TH2F**                  fHistoTrueClusGammaPtvsSource;                      //!
    TH2F**                  fHistoTrueClusElectronPtvsSource;                   //!
//    TH2F**                  fHistoTrueClusElectronPtvsMotherID;                 //!
//    TH2F**                  fHistoTrueClusElectronPtvsTopMotherID;              //!
//    TH2F**                  fHistoTrueClusElectronPtvsConvPhotonTopMotherID;    //!
    TH1F**                  fHistoTrueMergedMissedPDG;                          //!

    TH2F**                  fHistoTrueClusMergedPtvsRJet;                       //! array of histos with clusters from vs distance to Jet axis
    TH2F**                  fHistoTrueClusPi0PtvsRJet;                          //! array of histos with clusters from true pi0s vs distance to Jet axis
    TH2F**                  fHistoTrueClusEtaPtvsRJet;                          //! array of histos with clusters from true etas vs distance to Jet axis

    // MC validated reconstructed quantities mesons
    TH2F**                  fHistoTruePi0PtY;                                   //! array of histos with validated pi0, pt, Y
    TH2F**                  fHistoTrueEtaPtY;                                   //! array of histos with validated eta, pt, Y
    TH2F**                  fHistoTruePi0PtAlpha;                               //! array of histos with validated pi0, pt, alpha
    TH2F**                  fHistoTrueEtaPtAlpha;                               //! array of histos with validated eta, pt, alpha
    TH2F**                  fHistoTrueClusGammaEM02;                            //! array of histos with validated gamma, E, m02
    TH2F**                  fHistoTrueClusElectronEM02;                         //! array of histos with validated electrons, E, m02
    TH2F**                  fHistoTrueClusPi0EM02;                              //! array of histos with validated pi0, E, m02
    TH2F**                  fHistoTrueClusEtaEM02;                              //! array of histos with validated eta, E, m02
    TH2F**                  fHistoTrueClusBGEvsM02;                             //! array of histos with validated hadrons(BG), E, m02

    TH2F**                  fHistoTrueClusGammaEvsM20;                          //! array of histos with validated gamma, E, m02
    TH2F**                  fHistoTrueClusElectronEM20;                         //! array of histos with validated e+-, E, m02
    TH2F**                  fHistoTrueClusMergedPi0EVsM20;                      //! array of histos with validated merged pi0, E, m20
    TH2F**                  fHistoTrueClusMergedPi0EVsM02;                      //! array of histos with validated merged pi0, E, m02
    TH2F**                  fHistoTrueClusMergedEtaEVsM20;                      //! array of histos with validated merged eta, E, m20
    TH2F**                  fHistoTrueClusMergedEtaEVsM02;                      //! array of histos with validated merged eta, E, m02
    TH2F**                  fHistoTrueClusBGEvsM20;                             //! array of histos with validated hadrons(BG), E, m20

    TH2F**                  fHistoTruePrimaryPi0MCPtResolPt;          //! array of histos with validated weighted primary pi0, MCpt, resol pt
    TH2F**                  fHistoTruePrimaryPi0PureMergedMCPtResolPt;          //! array of histos with validated weighted primary pi0, MCpt, resol pt
    TH2F**                  fHistoTruePrimaryPi0MergedPartConvMCPtResolPt;      //! array of histos with validated weighted primary pi0, MCpt, resol pt
    TH2F**                  fHistoTruePrimaryPi01GammaMCPtResolPt;              //! array of histos with validated weighted primary pi0, MCpt, resol pt
    TH2F**                  fHistoTruePrimaryPi01ElectronMCPtResolPt;           //! array of histos with validated weighted primary pi0, MCpt, resol pt
    TH2F**                  fHistoTruePrimaryEtaMCPtResolPt;                    //! array of histos with validated weighted primary eta, MCpt, resol pt
    TH2F**                  fHistoTrueSecondaryPi0MCPtResolPt;                  //! array of histos with validated weighted secondary pi0, MCpt, resol pt

    // MC validated reconstructed quantities photons
    TH2F**                  fHistoDoubleCountTruePi0PtvsM02;                    //! array of histos with double counted pi0s, pT, M02
    TH1F**                  fHistoDoubleCountTrueSecPi0Pt;                      //! array of histos with double counted secondary pi0s, pT, M02
    TH2F**                  fHistoDoubleCountTrueMultiplePi0PtvsM02;            //! array of histos with double counted pi0s, pT, M02
    TH1F**                  fHistoDoubleCountTrueMultipleSecPi0Pt;              //! array of histos with double counted secondary pi0s, pT, M02
    TH2F**                  fHistoDoubleCountTrueEtaPtvsM02;                    //! array of histos with double counted etas, pT, M02
    vector<Int_t>           fVectorLabelsLeadingPi0;                            //! vector containing labels of validated pi0
    vector<Int_t>           fVectorLabelsMultiplePi0;                           //! vector containing labels of validated pi0
    vector<Int_t>           fVectorLabelsMultiplePi0Reduced;                    //! vector containing labels of validated pi0
    vector<Int_t>           fVectorDoubleCountTruePi0s;                         //! vector containing labels of validated pi0
    vector<Int_t>           fVectorDoubleCountTrueMultilePi0s;                  //! vector containing labels of validated pi0
    vector<Int_t>           fVectorDoubleCountTrueEtas;                         //! vector containing labels of validated eta

    // event histograms
    TH1F**                  fHistoNEvents;                                      //! array of histos with event information
    TH1F**                  fHistoNEventsWOWeight;                              //! array of histos with event information without event weights
    TH1F**                  fHistoNGoodESDTracks;                               //! array of histos with number of good tracks (2010 Standard track cuts)
    TH1F**                  fHistoVertexZ;                                      //! array of histos with vertex z distribution for selected events
    TH2F**                  fHistoOverlapsPi0All;                               //! array of histos with number of overlapping particles with generated pi0 within R<0.05
    TH2F**                  fHistoOverlapsPi0Accepted;                          //! array of histos with number of overlapping particles with generated pi0 within R<0.05
    TH1F**                  fHistoNClusterCandidates;                           //! array of histos with number of cluster candidates per event
    TH1F**                  fHistoNClusterMergedCandidates;                     //! array of histos with number of merged cluster candidates per event
    TH2F**                  fHistoNGoodESDTracksVsNClusterCandidates;           //! array of histos with number of good tracks vs gamma candidates
    TH2F**                  fHistoSPDClusterTrackletBackground;                 //! array of histos with SPD tracklets vs SPD clusters for background rejection
    TH1F**                  fHistoNV0Tracks;                                    //! array of histos with V0 counts
    TProfile**              fProfileEtaShift;                                   //! array of profiles with eta shift
    TProfile**              fProfileJetJetXSection;                             //! array of profiles with xsection for jetjet
    TH1F**                  fHistoJetJetNTrials;                                //! array of histos with ntrials for jetjet

    // Jet
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

    // additional variables
    TRandom3                fRandom;                                            // random
    Int_t                   fnCuts;                                             // number of cuts to be analysed in parallel
    Int_t                   fiCut;                                              // current cut
    Int_t                   fIsHeavyIon;                                        // switch for pp = 0, PbPb  = 1, pPb = 2
    Int_t                   fDoMesonQA;                                         // flag for meson QA
    Int_t                   fDoClusterQA;                                       // flag for cluster QA
    Bool_t                  fIsFromMBHeader;                                    // flag for MC headers
    Bool_t                  fIsOverlappingWithOtherHeader;                      // flag for particles in MC overlapping between headers
    Int_t                   fIsMC;                                              // flag for MC information
    Bool_t                  fSetPlotHistsExtQA;                                 // flag for extended QA hists
    Double_t                fWeightJetJetMC;                                    // weight for Jet-Jet MC
    Int_t                   fSelectedMesonID;                                   // switch for meson analysis
    Bool_t                  fEnableDetailedPrintOut;                            // switch on detailed print outs
    Bool_t                  fEnableSortForClusMC;                               // switch on sorting for MC labels in cluster
    TTree*                  tBrokenFiles;                                       // tree for keeping track of broken files
    TObjString*             fFileNameBroken;                                    // string object for broken file name
    Bool_t                  fDoDetailedM02;                                     // detailed M02 distribution
    Int_t                   fTrackMatcherRunningMode;                           // CaloTrackMatcher running mode
    Int_t                   fMinAllowedPi0OverlapsMC;                           // Event rejection based on pi0 overlaps in MC
    Int_t                   fMaxAllowedPi0OverlapsMC;                           // Event rejection based on pi0 overlaps in MC
    TH2F**                  fHistoPi0EvsGammaOverlapE;                          //! array of histos with SPD tracklets vs SPD clusters for background rejection

  private:
    AliAnalysisTaskGammaCaloMerged(const AliAnalysisTaskGammaCaloMerged&); // Prevent copy-construction
    AliAnalysisTaskGammaCaloMerged &operator=(const AliAnalysisTaskGammaCaloMerged&); // Prevent assignment

    ClassDef(AliAnalysisTaskGammaCaloMerged, 42);
};

#endif
