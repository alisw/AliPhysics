
#ifndef ALIANALYSISTASKHEAVYNEUTRALMESONTOGG_H
#define ALIANALYSISTASKHEAVYNEUTRALMESONTOGG_H

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
#include "TProfile2D.h"
#include "TH3.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TGenPhaseSpace.h"
#include <vector>
#include <map>

class AliAnalysisTaskHeavyNeutralMesonToGG : public AliAnalysisTaskSE {
public:

  AliAnalysisTaskHeavyNeutralMesonToGG();
  AliAnalysisTaskHeavyNeutralMesonToGG(const char *name);
  virtual ~AliAnalysisTaskHeavyNeutralMesonToGG();

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
  void ProcessTrueMesonCandidatesConvCalo     ( AliAODConversionMother *Pi0Candidate,
                                                AliAODConversionPhoton *TrueGammaCandidate0,
                                                AliAODConversionPhoton *TrueGammaCandidate1,
                                                Bool_t matched);
  void ProcessTrueMesonCandidatesConvCaloAOD  ( AliAODConversionMother *Pi0Candidate,
                                                AliAODConversionPhoton *TrueGammaCandidate0,
                                                AliAODConversionPhoton *TrueGammaCandidate1,
                                                Bool_t matched);
  void ProcessTrueMesonCandidatesCalo     ( AliAODConversionMother *Pi0Candidate,
                                            AliAODConversionPhoton *TrueGammaCandidate0,
                                            AliAODConversionPhoton *TrueGammaCandidate1 );
  void ProcessTrueMesonCandidatesCaloAOD  ( AliAODConversionMother *Pi0Candidate,
                                            AliAODConversionPhoton *TrueGammaCandidate0,
                                            AliAODConversionPhoton *TrueGammaCandidate1 );
  void ProcessTrueMesonCandidatesConv     ( AliAODConversionMother *Pi0Candidate,
                                            AliAODConversionPhoton *TrueGammaCandidate0,
                                            AliAODConversionPhoton *TrueGammaCandidate1 );
  void ProcessTrueMesonCandidatesConvAOD  ( AliAODConversionMother *Pi0Candidate,
                                            AliAODConversionPhoton *TrueGammaCandidate0,
                                            AliAODConversionPhoton *TrueGammaCandidate1 );

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

  TList**                             fCutFolder;                             // Array of lists for containers belonging to cut
  TList**                             fESDList;                               // Array of lists with histograms with reconstructed properties
  TList**                             fBackList;                              // Array of lists with BG THnSparseF
  TList**                             fMotherList;                            // Array of lists with Signal THnSparseF
  TList**                             fTrueList;                              // Array of lists with histograms with MC validated reconstructed properties
  TList**                             fMCList;                                // Array of lists with histograms with pure MC information
  TList*                              fOutputContainer;                       // Output container
  TList*                              fGammaCandidates;                       //! current list of photon candidates
  TList*                              fClusterCandidates;                     //! current list of cluster candidates
  TList*                              fEventCutArray;                         // List with Event Cuts
  TList*                              fCutArray;                              // List with Conversion Cuts
  TList*                              fClusterCutArray;                       // List with Cluster Cuts
  TList*                              fMesonCutArray;                         // List with Meson Cuts
  TClonesArray*                       fReaderGammas;                          // Array with conversion photons selected by V0Reader Cut

  TString                             fV0ReaderName;                          // V0Reader name to be found in input
  TString                             fCorrTaskSetting;                       // Correction Task Special Name

  TObjString*             fFileNameBroken;                                    // string object for broken file name

  THnSparseF**            fSparseMotherInvMassPtZM;                           //! array of THnSparseF with signal + BG for same event photon pairs, inv Mass, pt
  THnSparseF**            fSparseMotherBackInvMassPtZM;                       //! array of THnSparseF with BG for same event photon pairs, inv Mass, pt

  TTree*                  fTreeBrokenFiles;                                   // tree for keeping track of broken files

  TProfile2D**            fProfileTruePrimaryMesonWeightsInvMassPt;           //! array of profiles with weights for validated primary mothers, invMass, pt

  TH2F**                  fHistoMotherInvMassPt;                              //! array of histogram with signal + BG for same event photon pairs, inv Mass, pt
  TH2F**                  fHistoMotherMatchedInvMassPt;                       //! array of histogram with signal + BG for same event photon pairs, inv Mass, pt
  TH2F**                  fHistoMotherBackInvMassPt;                          //! array of histogram with BG for mixed event photon pairs, inv Mass, pt
  TH2F**                  fHistoMotherMesonPtY;                               //! array of histograms with invariant mass cut around nominal mass, pt, Y
  TH2F**                  fHistoMotherMesonPtAlpha;                           //! array of histograms with invariant mass cut around nominal mass, pt, alpha
  TH2F**                  fHistoMotherMesonPtOpenAngle;                       //! array of histograms with invariant mass cut around nominal mass, pt, openAngle
  TH2F**                  fHistoMotherMesonConvPhotonEtaPhi;                  //! array of histograms with invariant mass cut around nominal mass ,eta/phi of conversion photon
  TH2F**                  fHistoTrueMesonInvMassPt;                           //! array of histos with validated meson, invMass, pt
  TH2F**                  fHistoTrueMesonMatchedInvMassPt;                    //! array of histos with rejected meson, invMass, pt
  TH2F**                  fHistoTrueMesonCaloPhotonInvMassPt;                 //! array of histos with validated meson, photon leading, invMass, pt
  TH2F**                  fHistoTrueMesonCaloConvertedPhotonInvMassPt;        //! array of histos with validated meson, converted photon leading, invMass, pt
  TH2F**                  fHistoTrueMesonCaloMixedPhotonConvPhotonInvMassPt;  //! array of histos with validated meson converted photon leading and photon, invMass, pt
  TH2F**                  fHistoTrueMesonCaloConvertedPhotonMatchedInvMassPt; //! array of histos with validated meson matched with conv photon, converted photon leading, invMass, pt
  TH2F**                  fHistoTrueMesonCaloElectronInvMassPt;               //! array of histos with validated mothers, electron leading, invMass, pt
  TH2F**                  fHistoTrueMesonCaloMergedClusterInvMassPt;          //! array of histos with validated mothers, merged cluster invMass, pt
  TH2F**                  fHistoTrueMesonCaloMergedClusterPartConvInvMassPt;  //! array of histos with validated mothers, merged cluster part conv, invMass, pt
  TH2F**                  fHistoTruePrimaryMesonInvMassPt;                    //! array of histos with validated weighted primary mothers, invMass, pt
  TH2F**                  fHistoTruePrimaryMesonW0WeightingInvMassPt;         //! array of histos with validated unweighted primary mothers, invMass, pt
  TH2F**                  fHistoTruePrimaryMesonMCPtResolPt;                  //! array of histos with validated weighted primary meson, MCpt, resol pt
  TH2F**                  fHistoTrueMotherMesonConvPhotonEtaPhi;              //! array of histograms with invariant mass cut around nominal mass ,eta/phi of conversion photon
  TH2F**                  fHistoTrueBckGGInvMassPt;                           //! array of histos with pure gamma gamma combinatorial BG, invMass, pt
  TH2F**                  fHistoTrueBckFullMesonContainedInOneClusterInvMassPt; //! array of histos with meson fully contained in one cluster, invMass, pt
  TH2F**                  fHistoTrueBckAsymEClustersInvMassPt;                //! array of histos with asymmetry energy distributions of clusters, invMass, pt
  TH2F**                  fHistoTrueBckContInvMassPt;                         //! array of histos with contamination BG, invMass, pt
  TH2F**                  fHistoTrueMesonPtY;                                 //! array of histos with validated meson, pt, Y
  TH2F**                  fHistoTrueMesonPtAlpha;                             //! array of histos with validated meson, pt, alpha
  TH2F**                  fHistoTrueMesonPtOpenAngle;                         //! array of histos with validated meson, pt, openAngle
  TH2F**                  fHistoMCMesonPtY;                                   //! array of histos with weighted meson, pT, Y
  TH2F**                  fHistoMCMesonPtAlpha;                               //! array of histos with weighted meson, pT, alpha
  TH2F**                  fHistoMCMesonPtJetPt;                               //! array of histos with weighted meson, pT, hardest jet pt
  TH2F**                  fHistoTrueNLabelsInClusPt;                          //! array of histos with number of labels in cluster
  TH2F**                  fHistoDoubleCountTrueMesonInvMassPt;                //! array of histos with double counted mesons, invMass, pT
  TH2F**                  fHistoDoubleCountTrueConvGammaRPt;                  //! array of histos with double counted photons, R, pT
  TH2F**                  fHistoDoubleCountTrueClusterGammaPt;                //! array of histos with double counted cluster photons
  TH2F**                  fHistoSPDClusterTrackletBackground;                 //! array of histos for SPD Cluster vs Tracklet plot for pileup monitoring

  TProfile**              fProfileEtaShift;                                   //! array of profiles with eta shift
  TProfile**              fProfileJetJetXSection;                             //! array of profiles with xsection for jetjet

  TH1I**                  fHistoMCHeaders;                                    //! array of histos for header names

  TH1F**                  fHistoConvGammaPt;                                  //! array of histogram conversion photon pT
  TH1F**                  fHistoClusGammaPt;                                  //! array of histos with cluster, pt
  TH1F**                  fHistoClusGammaE;                                   //! array of histos with cluster, E
  TH1F**                  fHistoClusOverlapHeadersGammaPt;                    //! array of histos with cluster, pt overlapping with other headers
  TH1F**                  fHistoClusAllHeadersGammaPt;                        //! array of histos with cluster, pt all headers
  TH1F**                  fHistoClusRejectedHeadersGammaPt;                   //! array of histos with cluster, pt rejected headers
  TH1F**                  fHistoMotherInvMassRejected;                        //! array of histos with invariant mass pairs which were rejected
  TH1F**                  fHistoMCMesonPt;                                    //! array of histos with weighted meson, pT
  TH1F**                  fHistoMCMesonWOWeightPt;                            //! array of histos with unweighted meson, pT
  TH1F**                  fHistoMCMesonWOEvtWeightPt;                         //! array of histos without event weights meson, pT
  TH1F**                  fHistoMCMesonInAccPt;                               //! array of histos with weighted meson in acceptance, pT
  TH1F**                  fHistoMCMesonWOWeightInAccPt;                       //! array of histos without weight meson in acceptance, pT
  TH1F**                  fHistoMCMesonWOEvtWeightInAccPt;                    //! array of histos without evt weight meson in acceptance, pT
  TH1F**                  fHistoTrueConvGammaPt;                              //! array of histos with validated conversion photon, pt
  TH1F**                  fHistoTruePrimaryConvGammaPt;                       //! array of histos with validated primary
  TH1F**                  fHistoTrueClusGammaPt;                              //! array of histos with validated cluster (electron or photon), pt
  TH1F**                  fHistoTrueClusConvGammaPt;                          //! array of histos with validated converted photon, pt
  TH1F**                  fHistoTrueClusConvGammaFullyPt;                     //! array of histos with validated converted photon, fully contained, pt
  TH1F**                  fHistoTruePrimaryClusGammaPt;                       //! array of histos with validated primary cluster, photons, pt
  TH1F**                  fHistoTruePrimaryClusConvGammaPt;                   //! array of histos with validated primary cluster, converted photons, pt
  TH1F**                  fHistoMultipleCountTrueMeson;                       //! array of histos how often TrueMesons are counted
  TH1F**                  fHistoMultipleCountTrueConvGamma;                   //! array of histos how often TrueConvGammas are counted
  TH1F**                  fHistoMultipleCountTrueClusterGamma;                //! array of histos how often TrueClusterGammas are counted
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

  // variable to keep track of multiple & missing reco
  map<Int_t,Int_t>        fMapMultipleCountTrueMesons;                        //! map containing meson labels that are counted at least twice
  map<Int_t,Int_t>        fMapMultipleCountTrueConvGammas;                    //! map containing photon labels that are counted at least twice
  map<Int_t,Int_t>        fMapMultipleCountTrueClusterGammas;                 //! map containing cluster photon labels that are counted at least twice
  vector<Int_t>           fVectorRecTrueMesons;                               //! array of strings containing the stack position of the reconstructed validated meson
  vector<Int_t>           fVectorDoubleCountTrueMesons;                       //! vector containing labels of validated meson
  vector<Int_t>           fVectorDoubleCountTrueConvGammas;                   //! vector containing labels of validated photons
  vector<Int_t>           fVectorDoubleCountTrueClusterGammas;                //! vector containing labels of validated cluster photons

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

  TGenPhaseSpace        fGenPhaseSpace;                                       // For generation of decays into two gammas (rotation background)

private:
  AliAnalysisTaskHeavyNeutralMesonToGG(const AliAnalysisTaskHeavyNeutralMesonToGG&); // Prevent copy-construction
  AliAnalysisTaskHeavyNeutralMesonToGG &operator=(const AliAnalysisTaskHeavyNeutralMesonToGG&); // Prevent assignment

  ClassDef(AliAnalysisTaskHeavyNeutralMesonToGG, 4);
};

#endif // ALIANALYSISTASKHEAVYNEUTRALMESONTOGG_H
