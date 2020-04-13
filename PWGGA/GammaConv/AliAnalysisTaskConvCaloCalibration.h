
#ifndef AliAnalysisTaskConvCaloCalibration_H
#define AliAnalysisTaskConvCaloCalibration_H

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
#include "AliDalitzElectronSelector.h"
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"
#include "THnSparse.h"
#include <vector>
#include <map>

class AliAnalysisTaskConvCaloCalibration : public AliAnalysisTaskSE {
public:

  AliAnalysisTaskConvCaloCalibration();
  AliAnalysisTaskConvCaloCalibration(const char *name);
  virtual ~AliAnalysisTaskConvCaloCalibration();

  virtual void   UserCreateOutputObjects();
  virtual Bool_t Notify();
  virtual void   UserExec(Option_t *);
  virtual void   Terminate(const Option_t*);
  void InitBack();

  void SetV0ReaderName(TString name){fV0ReaderName=name; return;}
  void SetIsHeavyIon(Int_t flag){
      fIsHeavyIon = flag;
  }
  // Function to set correction task setting
  void SetCorrectionTaskSetting(TString setting) {fCorrTaskSetting = setting;}
  void SetDoClusterSelectionForTriggerNorm(Bool_t flag)         { fEnableClusterCutsForTrigger= flag    ;}

  // base functions for selecting photon and meson candidates in reconstructed data
  void ProcessClusters();
  void ProcessPhotonCandidates();
  void CalculateMesonCandidates();
  void SetPhotonVeto();
  void GetV0Electrons();

  // MC functions
  void SetIsMC                        ( Int_t isMC)                                       { fIsMC = isMC                              ;}
  // void ProcessMCParticles             ();
  // void ProcessAODMCParticles          ();
  void RelabelAODPhotonCandidates     ( Bool_t mode);

  // switches for additional analysis streams or outputs
  void SetDoPrimaryTrackMatching      ( Bool_t flag )                                     { fDoPrimaryTrackMatching = flag              ;}
  void SetLightOutput                 ( Bool_t flag )                                     { fDoLightOutput = flag                       ;}
  void SetMesonRecoMode               ( Int_t flag )                                      { fMesonRecoMode = flag                       ;}
  void SetMesonType                   ( Int_t flag )                                      { fMesonType = flag                           ;}
  void SetDoMesonQA                   ( Int_t flag )                                      { fDoMesonQA = flag                           ;}
  void SetDoPhotonQA                  ( Int_t flag )                                      { fDoPhotonQA = flag                          ;}
  void SetDoClusterQA                 ( Int_t flag )                                      { fDoClusterQA = flag                         ;}
  void SetUseTHnSparse                ( Bool_t flag )                                     { fDoTHnSparse = flag                         ;}
  void SetPlotHistsExtQA              ( Bool_t flag )                                     { fSetPlotHistsExtQA = flag                   ;}
  void SetDoTreeConvGammaShowerShape  ( Bool_t flag )                                     { fDoConvGammaShowerShapeTree = flag          ;}
  void SetDoTreeInvMassShowerShape    ( Bool_t flag )                                     { fDoInvMassShowerShapeTree = flag            ;}
  void SetAllowOverlapHeaders         ( Bool_t allowOverlapHeader )                       { fAllowOverlapHeaders = allowOverlapHeader   ;}
  void SetElectronMatchingCalibration ( Int_t flag )                                      { fUseEletronMatchingCalibration = flag       ;}

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

void SetNumOfCaloModules              ( Int_t nModules)                                   {
                                        fnModules = nModules                        ;
                                        if(nModules < 1 || nModules > 20){
                                          fnModules = 20                            ;
                                        }
                                      }

  // BG HandlerSettings
  void CalculateBackground            ();
  void CalculateBackgroundRP          ();
  void RotateParticle                 ( AliAODConversionPhoton *gamma );
  void RotateParticleAccordingToEP    ( AliAODConversionPhoton *gamma,
                                        Double_t previousEventEP,
                                        Double_t thisEventEP );
  void SetMoveParticleAccordingToVertex       ( Bool_t flag )                             { fMoveParticleAccordingToVertex = flag       ;}
  void MoveParticleAccordingToVertex          ( AliAODConversionPhoton* particle,
                                                const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex);
  void UpdateEventByEventData         ();

  // Additional functions for convenience
  void SetLogBinningXTH2              ( TH2* histoRebin );
  Bool_t CheckVectorOnly              ( vector<Int_t> &vec,
                                        Int_t tobechecked );
  Bool_t CheckVectorForDoubleCount    ( vector<Int_t> &vec,
                                        Int_t tobechecked );

  void FillMultipleCountMap           ( map<Int_t,Int_t> &ma,
                                        Int_t tobechecked );
  void FillMultipleCountHistoAndClear ( map<Int_t,
                                        Int_t> &ma,
                                        TH1F* hist );

  // Function to enable MC label sorting
  void SetEnableSortingOfMCClusLabels (Bool_t enableSort) { fEnableSortForClusMC   = enableSort;}

  void SetTrackMatcherRunningMode(Int_t mode){fTrackMatcherRunningMode = mode;}

protected:
  TRandom3                            fRandom;                                // random number
  AliV0ReaderV1*                      fV0Reader;                              // basic photon Selection Task
  AliGammaConversionAODBGHandler**    fBGHandler;                             // BG handler for Conversion
  AliConversionAODBGHandlerRP**       fBGHandlerRP;                           // BG handler for Conversion (possibility to mix with respect to RP)
  AliGammaConversionAODBGHandler**    fBGClusHandler;                         // BG handler for Cluster
  AliConversionAODBGHandlerRP**       fBGClusHandlerRP;                       // BG handler for Cluster (possibility to mix with respect to RP)
  AliVEvent*                          fInputEvent;                            // current event
  AliMCEvent*                         fMCEvent;                               // corresponding MC event
  AliConvEventCuts*                   fEventCuts;                             // EventCutObject
  AliConversionPhotonCuts*            fConversionCuts;                        // ConversionCutObject
  AliCaloPhotonCuts*                  fCaloPhotonCuts;                        // CaloPhotonCutObject
  AliConversionMesonCuts*             fMesonCuts;                             // MesonCutObject
  AliEMCALGeometry*                   fGeomEMCAL;                             // pointer to EMCAL geometry
  AliDalitzElectronSelector*          fElecSelector;					      // basic electron Selection

  TList**                             fCutFolder;                             // Array of lists for containers belonging to cut
  TList**                             fESDList;                               // Array of lists with histograms with reconstructed properties
  TList**                             fBackList;                              // Array of lists with BG THnSparseF
  TList**                             fMotherList;                            // Array of lists with Signal THnSparseF
  TList*                              fOutputContainer;                       // Output container
  TList*                              fGammaCandidates;                       //! current list of photon candidates
  TList*                              fClusterCandidates;                     //! current list of cluster candidates
  TList*                              fEventCutArray;                         // List with Event Cuts
  TList*                              fCutArray;                              // List with Conversion Cuts
  TList*                              fClusterCutArray;                       // List with Cluster Cuts
  TList*                              fMesonCutArray;                         // List with Meson Cuts
  TClonesArray*                       fReaderGammas;                          // Array with conversion photons selected by V0Reader Cut
  std::vector<Int_t>                  fSelectorElectronIndex;                 // Vector with electrons selected with std. cut
  std::vector<Int_t>                  fSelectorPositronIndex;                 // Vector with positrons selected with std. cut
  std::vector<Int_t>                  fV0Electrons;                           // Vector with electrons from V0s

  TString                             fV0ReaderName;                          // V0Reader name to be found in input
  TString                             fCorrTaskSetting;                       // Correction Task Special Name

  TObjString*             fFileNameBroken;                                    // string object for broken file name

  THnSparseF**            fSparseMotherInvMassPtZM;                           //! array of THnSparseF with signal + BG for same event photon pairs, inv Mass, pt
  THnSparseF**            fSparseMotherBackInvMassPtZM;                       //! array of THnSparseF with BG for same event photon pairs, inv Mass, pt

  TTree*                  fTreeBrokenFiles;                                   // tree for keeping track of broken files

  // TProfile2D**            fProfileTruePrimaryMesonWeightsInvMassPt;           //! array of profiles with weights for validated primary mothers, invMass, pt

  TH2F**                  fHistoMotherInvMassPt;                              //! array of histogram with signal + BG for same event photon pairs, inv Mass, pt
  TH2F**                  fHistoMotherMatchedInvMassPt;                       //! array of histogram with signal + BG for same event photon pairs, inv Mass, pt
  TH2F**                  fHistoMotherBackInvMassPt;                          //! array of histogram with BG for mixed event photon pairs, inv Mass, pt
  TH2F**                  fHistoMotherMesonPtY;                               //! array of histograms with invariant mass cut around nominal mass, pt, Y
  TH2F**                  fHistoMotherMesonPtAlpha;                           //! array of histograms with invariant mass cut around nominal mass, pt, alpha
  TH2F**                  fHistoMotherMesonPtOpenAngle;                       //! array of histograms with invariant mass cut around nominal mass, pt, openAngle
  TH2F**                  fHistoMotherMesonConvPhotonEtaPhi;                  //! array of histograms with invariant mass cut around nominal mass ,eta/phi of conversion photon
  TH2F**                  fHistoSPDClusterTrackletBackground;                 //! array of histos for SPD Cluster vs Tracklet plot for pileup monitoring
  TH2F***                 fHistoMotherInvMassECalibSM;                        //! array of histos with signal + background with alpha < 0.1 for NonLin for every Supermodule
  TH2F***                 fHistoMotherBackInvMassECalibSM;                    //! array of histos with mixed event background with alpha < 0.1 for NonLin for every Supermodule
  TH2F**                  fHistoMotherInvMassECalib;                          //! array of histos with signal + background with alpha < 0.1 for NonLin
  TH2F**                  fHistoMotherBackInvMassECalib;                      //! array of histos with mixed event background with alpha < 0.1 for NonLin

  TProfile**              fProfileEtaShift;                                   //! array of profiles with eta shift
  TProfile**              fProfileJetJetXSection;                             //! array of profiles with xsection for jetjet

  // TH1I**                  fHistoMCHeaders;                                    //! array of histos for header names

  TH1F**                  fHistoConvGammaPt;                                  //! array of histogram conversion photon pT
  TH1F**                  fHistoClusGammaPt;                                  //! array of histos with cluster, pt
  TH1F**                  fHistoClusGammaE;                                   //! array of histos with cluster, E
  TH1F***                 fHistoClusGammaPtSM;                                //! array of histos with cluster, pt
  TH1F***                 fHistoClusGammaESM;                                 //! array of histos with cluster, E
  TH1F**                  fHistoMotherInvMassRejected;                        //! array of histos with invariant mass pairs which were rejected
  TH1F**                  fHistoNEvents;                                      //! array of histos with event information
  TH1F**                  fHistoNEventsWOWeight;                              //! array of histos with event information without event weights
  TH1F**                  fHistoNGoodESDTracks;                               //! array of histos with number of good tracks (2010 Standard track cuts)
  TH1F**                  fHistoVertexZ;                                      //! array of histos with vertex z distribution for selected events
  TH1F**                  fHistoVertexX;                                      //! array of histos with vertex x distribution for selected events
  TH1F**                  fHistoVertexY;                                      //! array of histos with vertex y distribution for selected events
  TH1F**                  fHistoNGammaConvCandidates;                         //! array of histos with number of conversion gamma candidates per event
  TH1F**                  fHistoNGammaCaloCandidates;                         //! array of histos with number of calo gamma candidates per event
  TH1F**                  fHistoNV0Tracks;                                    //! array of histos with V0 counts
  TH1F**                  fHistoJetJetNTrials;                                //! array of histos with ntrials for jetjet

  // additional variables
  Double_t*               fUnsmearedPx;                                       //[fNGammaCandidates]
  Double_t*               fUnsmearedPy;                                       //[fNGammaCandidates]
  Double_t*               fUnsmearedPz;                                       //[fNGammaCandidates]
  Double_t*               fUnsmearedE;                                        //[fNGammaCandidates]
  Double_t*               fMesonInvMassWindow;                                // minimum inv mass for histos

  Int_t*                  fMCEventPos;                                        //[fNGammaCandidates]
  Int_t*                  fMCEventNeg;                                        //[fNGammaCandidates]
  Int_t*                  fESDArrayPos;                                       //[fNGammaCandidates]
  Int_t*                  fESDArrayNeg;                                       //[fNGammaCandidates]

  Double_t                fEventPlaneAngle;                                   // EventPlaneAngle
  Double_t                fMesonInvMassMin;                                   // minimum inv mass for histos
  Double_t                fMesonInvMassMax;                                   // maximum inv mass for histos
  Double_t                fMesonInvMassNBins;                                 // Number of bins for inv mass histos
  Double_t                fWeightJetJetMC;                                    // weight for Jet-Jet MC

  Int_t                   fNGammaCandidates;                                  // number of gamma candidates in event
  Int_t                   fnCuts;                                             // number of cuts to be analysed in parallel
  Int_t                   fiCut;                                              // current cut
  Int_t                   fIsHeavyIon;                                        // switch for pp = 0, PbPb = 1, pPb = 2
  Int_t                   fMesonRecoMode;                                     // switch for running with different reconstruction modes: 0 - PCM-PCM, 1 - PCM-Calo, 2 - Calo-Calo
  Int_t                   fMesonType;                                         // selector for meson analysis
  Int_t                   fMesonPDG;                                          // PDG code for selected meson
  Int_t                   fDoMesonQA;                                         // flag for meson QA
  Int_t                   fDoPhotonQA;                                        // flag for photon QA
  Int_t                   fDoClusterQA;                                       // flag for cluster QA
  Int_t                   fIsMC;                                              // flag for MC information
  Int_t                   fnModules;                                          // number of SM of EMCal+DCal

  Bool_t                  fMoveParticleAccordingToVertex;                     // boolean for BG calculation
  Bool_t                  fDoLightOutput;                                     // switch for running light output, kFALSE -> normal mode, kTRUE -> light mode
  Bool_t                  fIsFromDesiredHeader;                               // flag for MC headers
  Bool_t                  fIsOverlappingWithOtherHeader;                      // flag for particles in MC overlapping between headers
  Bool_t                  fDoTHnSparse;                                       // flag for using THnSparses for background estimation
  Bool_t                  fSetPlotHistsExtQA;                                 // flag for extended QA hists
  Bool_t                  fDoConvGammaShowerShapeTree;                        // flag for tree with conv gamma R vs energy vs shower shape
  Bool_t                  fEnableSortForClusMC;                               // switch on sorting for MC labels in cluster
  Bool_t                  fDoPrimaryTrackMatching;                            // switch for basic track matching for primaries
  Bool_t                  fDoInvMassShowerShapeTree;                          // flag for producing tree tESDInvMassShowerShape
  Bool_t                  fAllowOverlapHeaders;                               // enable overlapping headers for cluster selection
  Bool_t                  fEnableClusterCutsForTrigger;                       // enable CLusterCuts output for trigger only
  Int_t                   fTrackMatcherRunningMode;                           // CaloTrackMatcher running mode
  Int_t                   fUseEletronMatchingCalibration;                     // switch for calibration using electron tracks (1) or electrons from V0s (2) to cluster matching

private:
  AliAnalysisTaskConvCaloCalibration(const AliAnalysisTaskConvCaloCalibration&); // Prevent copy-construction
  AliAnalysisTaskConvCaloCalibration &operator=(const AliAnalysisTaskConvCaloCalibration&); // Prevent assignment

  ClassDef(AliAnalysisTaskConvCaloCalibration, 4);
};

#endif // AliAnalysisTaskConvCaloCalibration_H
