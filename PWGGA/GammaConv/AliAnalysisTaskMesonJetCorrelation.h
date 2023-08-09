/**************************************************************************
 * Copyright(c) 1998-2020, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Joshua Koenig <joshua.konig@cern.ch>                                        *
 * Version 1.0                                                            *
 *                                                                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do analysis for light neutral mesons inside Jets
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////

#ifndef ALIANLYSISTASKMesonJetCorrelation_cxx
#define ALIANLYSISTASKMesonJetCorrelation_cxx

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskConvJet.h"
#include "AliAnalysisTaskJetOutlierRemoval.h"
#include "AliGenHepMCEventHeader.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliCaloPhotonCuts.h"
#include "AliConvEventCuts.h"
#include "AliConversionAODBGHandlerRP.h"
#include "AliConversionMesonCuts.h"
#include "AliConversionPhotonCuts.h"
#include "AliESDtrack.h"
#include "AliGammaConversionAODBGHandler.h"
#include "AliKFConversionPhoton.h"
#include "AliV0ReaderV1.h"
#include "TH3.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TGenPhaseSpace.h"
#include "TProfile2D.h"
#include <map>
#include <vector>
#include "AliResponseMatrixHelper.h"
#include "AliGammaConvEventMixing.h"

class AliAnalysisTaskMesonJetCorrelation : public AliAnalysisTaskSE
{
  
 public:
  AliAnalysisTaskMesonJetCorrelation();
  AliAnalysisTaskMesonJetCorrelation(const char* name);
  virtual ~AliAnalysisTaskMesonJetCorrelation();

  virtual void UserCreateOutputObjects();
  virtual Bool_t Notify();
  virtual void UserExec(Option_t*);
  virtual void Terminate(const Option_t*);

  // base functions for selecting photon and meson candidates
  void ProcessClusters();
  void ProcessPhotonCandidates();
  void CalculateMesonCandidates();
  void CalculateBackground();
  void CalculateBackgroundSwapp();
  void CalculateBackgroundMix();
  void UpdateEventMixData();
  void FillInvMassBackHistograms(AliAODConversionMother* backgroundCandidate, const bool isRotBack = true);
  std::array<std::unique_ptr<AliAODConversionPhoton>, 2> GetGammasSwapped(AliAODConversionPhoton* currentEventGoodV0Temp1, AliAODConversionPhoton* currentEventGoodV0Temp2);
  void FillMesonHistograms(AliAODConversionPhoton* gamma0, AliAODConversionPhoton* gamma1, int firstGammaIndex, int secondGammaIndex);
  void ProcessTrueBackgroundCandidatesAOD(AliAODConversionMother* Pi0Candidate, AliAODConversionPhoton* TrueGammaCandidate0, AliAODConversionPhoton* TrueGammaCandidate1, const int matchedJet, const float RJetPi0Cand);
  float GetFrag(AliAODConversionMother* Pi0Candidate, const int matchedJet, int isTrueJet);
  float GetFrag(AliAODMCParticle* Pi0Candidate, const int matchedJet, int isTrueJet);
  float GetRadiusJetPart(AliAODConversionMother* Pi0Candidate, const int matchedJet, int isTrueJet);
  float GetRadiusJetPart(AliAODMCParticle* Pi0Candidate, const int matchedJet, int isTrueJet);
  void FillMesonDCATree(AliAODConversionMother* Pi0Candidate, AliAODConversionPhoton* gamma0, AliAODConversionPhoton* gamma1, const int matchedJet, const bool isTrueMeson);

  // MC functions
  void ProcessAODMCParticles(int isCurrentEventSelected = 0);
  void ProcessTrueClusterCandidatesAOD(AliAODConversionPhoton* TruePhotonCandidate, const int matchedJet = -1);
  void ProcessTruePhotonCandidatesAOD(AliAODConversionPhoton* TruePhotonCandidate);
  bool MCParticleIsSelected(AliAODMCParticle* particle1, AliAODMCParticle* particle2, bool checkConversion);
  bool MCParticleIsSelected(AliAODMCParticle* particle, bool isConv, bool checkConversion);
  int GetPhotonMotherLabel(AliAODConversionPhoton* gammaCand, int& convertedPhotonLabel, bool isCaloPhoton);
  void RelabelAODPhotonCandidates(Bool_t mode);
  bool ProcessTrueMesonCandidatesAOD(AliAODConversionMother* Pi0Candidate, AliAODConversionPhoton* TrueGammaCandidate0, AliAODConversionPhoton* TrueGammaCandidate1, const int matchedJet, const float RJetPi0Cand = 0, const double weightMatBudget = 1.);
  void ProcessTrueMesonCandidatesInTrueJetsAOD(AliAODConversionMother* Pi0Candidate, AliAODConversionPhoton* TrueGammaCandidate0, AliAODConversionPhoton* TrueGammaCandidate1, const int matchedJet, const float RJetPi0Cand = 0, const double weightMatBudget = 1.);
  // void IsTrueParticle(AliAODConversionMother* Pi0Candidate, AliAODConversionPhoton* TrueGammaCandidate0, AliAODConversionPhoton* TrueGammaCandidate1, Bool_t matched);
  bool CheckAcceptance(AliAODMCParticle* gamma0, AliAODMCParticle* gamma1);
  bool IsParticleFromPartonFrag(AliAODMCParticle* particle, int idParton);
  void UnselectStablePi0(AliAODMCParticle* part); // Rleabel Pi0, eta etc. tounstable in case they were labeled stable by the jet framework

  // Jet functions
  void ProcessJets(int isCurrentEventSelected = 0);
  void InitJets();
  bool IsJetInAcc(unsigned int accType, double r, double eta, double phi);
  bool IsInEMCalAcc(double r, double eta, double phi);
  bool IsInDCalAcc(double r, double eta, double phi);

  // Helper functions
  void MakeBinning();
  void CallSumw2ForLists(TList* l);

  // Setters
  void SetIsMC(int isMC) { fIsMC = isMC; }
  void SetLightOutput(int flag) { fDoLightOutput = flag; }
  void SetDoMesonQA(int flag) { fDoMesonQA = flag; }
  void SetDoPhotonQA(int flag) { fDoPhotonQA = flag; }
  void SetDoClusterQA(int flag) { fDoClusterQA = flag; }
  // Function to set correction task setting
  void SetCorrectionTaskSetting(TString setting) { fCorrTaskSetting = setting; }
  void SetDoMaterialBudgetWeightingOfGammasForTrueMesons(Bool_t flag) { fDoMaterialBudgetWeightingOfGammasForTrueMesons = flag; }
  void SetIsCalo(bool isCalo) { fIsCalo = isCalo; }
  void SetIsConv(bool isConv) { fIsConv = isConv; }
  void SetIsConvCalo(bool isConvCalo) { fIsConvCalo = isConvCalo; }
  void SetIsHeavyIon(int flag) { fIsHeavyIon = flag; }
  void SetV0ReaderName(TString name) { fV0ReaderName = name; }
  void SetTrackMatcherRunningMode(int mode) { fTrackMatcherRunningMode = mode; }
  void SetUseTHnSparseForResponse(bool tmp) { fUseThNForResponse = tmp; }
  void SetMesonKind(int meson) { (meson == 0) ? fMesonPDGCode = 111 : fMesonPDGCode = 221; }
  void SetOtherMesons(std::vector<int> vec) { fOtherMesonsPDGCodes = vec; }
  void SetJetContainerAddName(TString name) { fAddNameConvJet = name; }
  void SetFillMesonDCATree(bool tmp) { fFillDCATree = tmp; }
  void SetDoUseCentralEvtSelection(bool tmp) { fUseCentralEventSelection = tmp; }
  void SetUsePtForCalcZ(bool tmp) { fUsePtForZCalc = tmp; }
  void SetForcePi0Unstable(bool tmp) { fUnsetStablePi0 = tmp; }
  void SetUseMixedBackAdd(bool tmp) { fUseMixedBackAdd = tmp; }
  void SetDoRadiusDependence(bool tmp) { fDoRadiusDep = tmp; }
  void SetMesonZPt(int tmp)
  {
    if(tmp == 1){
      fDoAnalysisZ = false;
    } else if(tmp == 2){
      fDoAnalysisPt = false;
    }
  }

  void SetEventCutList(int nCuts,
                       TList* CutArray)
  {
    fnCuts = nCuts;
    fEventCutArray = CutArray;
  }

  // Setting the cut lists for the conversion photons
  void SetConversionCutList(int nCuts,
                            TList* CutArray)
  {
    fnCuts = nCuts;
    fConvCutArray = CutArray;
  }

  // Setting the cut lists for the calo photons
  void SetCaloCutList(int nCuts,
                      TList* CutArray)
  {
    fnCuts = nCuts;
    fClusterCutArray = CutArray;
  }

  // Setting the cut lists for the meson
  void SetMesonCutList(int nCuts,
                       TList* CutArray)
  {
    fnCuts = nCuts;
    fMesonCutArray = CutArray;
  }

 protected:
  //-------------------------------
  // GLobal settings
  //-------------------------------
  AliV0ReaderV1* fV0Reader;       // basic photon Selection Task
  TString fV0ReaderName;          // name of V0Reader
  TClonesArray* fReaderGammas;    // Array with conversion photons selected by V0Reader Cut
  TString fCaloTriggerHelperName; // name of trigger helper for PHOS
  TString fCorrTaskSetting;       // Correction Task Special Name
  AliVEvent* fInputEvent;         // current event
  AliMCEvent* fMCEvent;           // corresponding MC event
  TClonesArray* fAODMCTrackArray; // pointer to track array
  AliEventCuts fAliEventCuts;     ///<  Event cuts (run2 defaults)

  //-------------------------------
  // Lists for cut folders and output containers
  //-------------------------------
  TList** fCutFolder;      // Top level cut folder in which all folders belonging to this cut are stored
  TList** fESDList;        // List for standard histograms for data+MC (Inv mass vs pt etc.)
  TList** fMCList;         // List for MC generated histograms
  TList** fTrueList;       // List for true meson quantities
  TList** fJetList;        // List for jet related observables
  TList** fTrueJetList;    // List of true jet quantities (response matrix)
  TList* fOutputContainer; // Top level output container
  TList* fEventCutArray;   // Event cut output container
  TList* fConvCutArray;    // Conversion cut output container
  TList* fClusterCutArray; // Cluster cut output container
  TList* fMesonCutArray;   // Meson cut output container

  std::vector<AliAODConversionPhoton*> fGammaCandidates;   //! current list of photon candidates
  std::vector<AliAODConversionPhoton*> fClusterCandidates; //! current list of cluster candidates

  EventMixPoolMesonJets* fEventMix; // Eventmixing class used for events with a jet axis (also works with no jet axis but has no z-vertex bins etc.)

  TRandom3 fRandom;              //! random generator needed for the rotation background
  TGenPhaseSpace fGenPhaseSpace; //! For generation of decays into two gammas

  //-------------------------------
  // Jet related variables
  //-------------------------------
  AliAnalysisTaskConvJet* fConvJetReader; // JetReader
  TString fAddNameConvJet;                // Jet container additional name
  TVector3 fHighestJetVector;             // vector px,py,pz of jet with highest momentum in the event. Needed for jet event mixing
  float fMaxPtJet;                        // pt of jet with highest pt in the event. Needed for jet ecvent mixing

  AliAnalysisTaskJetOutlierRemoval* fOutlierJetReader; // Jet outlier Reader

  //-------------------------------
  // global settings and variables
  //-------------------------------
  int fMesonPDGCode;                                    // PDG code of current meson (111 for pi0 etc.)
  std::vector<int> fOtherMesonsPDGCodes;                // PDG code of other mesons (eta code if we are looking for a pi0)
  int fiCut;                                            // index of the current cut
  int fIsMC;                                            // flag for data or MC (JJ MC > 1)
  int fnCuts;                                           // number of cuts
  double fWeightJetJetMC;                               //! weights if jet-jet MC is used
  int fDoLightOutput;                                   // flag if light output should be used
  int fDoMesonQA;                                       // flag if meson QA should be switched on
  int fDoPhotonQA;                                      // flag if photon QA should be switched on
  int fDoClusterQA;                                     // flag if cluster QA should be switched on
  int fDoJetQA;                                         // flag if Jet QA should be switched on
  int fIsHeavyIon;                                      // lag for heavy ion
  bool fIsCalo;                                         // flag if current analysis is calo only
  bool fIsConv;                                         // flag if current analysis is using conversions only
  bool fIsConvCalo;                                     // flag if current analysis is
  bool fIsFromDesiredHeader;                            // flag if particle is from desired header
  bool fDoMaterialBudgetWeightingOfGammasForTrueMesons; // flag if material budget weights should be applied
  double fEventPlaneAngle;                              // event plane angle
  int fTrackMatcherRunningMode;                         // track matcher mode
  bool fUseThNForResponse;                              // flag if THnSparse or TH2 should be used for the 4d response matrices
  bool fEnableSortForClusMC;                            // flag if cluster mc labels should be sorted
  bool fFillDCATree;                                    // flag if DCA tree should be filled or not
  bool fUseCentralEventSelection;                       // flag if central event selection (AliEventSelection.cxx) should be used
  bool fUsePtForZCalc;                                  // flag if z = pt_meson/pt_jet or if its z = (P_{meson}*P_{jet})/|P_{Jet}^{2}|
  bool fUnsetStablePi0;                                 // flag to decide if pi0 need to be reset to Unstable particles
  bool fUseMixedBackAdd;                                // flag to enable a histogram for the mixed jet background in addition to the rotation background. This is to save memory and CPU (As otherwise a completely new wagon would be needed)
  bool fDoRadiusDep;                                    // flag to enable radius dependent histograms
  bool fDoAnalysisPt;                                   // flag to enable filling of pt dependent histograms
  bool fDoAnalysisZ;                                    // flag to enable filling of z dependent histograms
  //-------------------------------
  // conversions
  //-------------------------------
  int* fMCEventPos;  //!
  int* fMCEventNeg;  //!
  int* fESDArrayPos; //!
  int* fESDArrayNeg; //!

  //-------------------------------
  // Materil Budged Weights
  //-------------------------------
  std::vector<double> vecWeightsMatWeightsGammas;    //! vector with weights for photons. Indexing is the same as for fGammaCandidates

  //-------------------------------
  // binning settings
  //-------------------------------
  std::vector<double> fVecBinsMesonInvMass;   //! meson inv. mass binning
  std::vector<double> fVecBinsPhotonPt;       //! photon/cluster pt binning
  std::vector<double> fVecBinsClusterPt;      //! cluster binning until high pT
  std::vector<double> fVecBinsMesonPt;        //! meson pt binning
  std::vector<double> fVecBinsMesonPtCoarse;  //! coarse meson pt binning
  std::vector<double> fVecBinsJetPt;          //! jet pt binning
  std::vector<double> fVecBinsFragment;       //! z (fragmentation function) binning
  std::vector<double> fVecBinsMesonJetRadius; //! radius binning
  std::vector<double> vecEquidistFromMinus05; //! aequdistant binning starting from -0.5

  //-------------------------------
  // Jet related vectors
  //-------------------------------
  vector<double> fVectorJetPt;           //! Vector of JetPt
  vector<double> fVectorJetPx;           //! Vector of JetPx
  vector<double> fVectorJetPy;           //! Vector of JetPy
  vector<double> fVectorJetPz;           //! Vector of JetPz
  vector<double> fVectorJetEta;          //! Vector of JetEta
  vector<double> fVectorJetPhi;          //! Vector of JetPhi
  vector<double> fVectorJetArea;         //! Vector of JetArea
  vector<double> fVectorJetNEF;          //! Vector of jet neutral energy fraction
  vector<double> fVectorJetNch;          //! Vector of jet number of tracks
  vector<double> fVectorJetNclus;        //! Vector of jet number of clusters
  vector<double> fTrueVectorJetPt;       //! Vector of True JetPt
  vector<double> fTrueVectorJetPx;       //! Vector of True JetPx
  vector<double> fTrueVectorJetPy;       //! Vector of True JetPy
  vector<double> fTrueVectorJetPz;       //! Vector of True JetPz
  vector<double> fTrueVectorJetEta;      //! Vector of True JetEta
  vector<double> fTrueVectorJetPhi;      //! Vector of True JetPhi
  vector<double> fTrueVectorJetNPart;    //! Vector of True number of particles in jet
  vector<int> fTrueVectorJetPartonID;    //! Vector of parton id matched to true jet
  vector<double> fTrueVectorJetPartonPt; //! Vector of parton pt matched to true jet
  vector<double> fTrueVectorJetPartonPx; //! Vector of parton pt matched to true jet
  vector<double> fTrueVectorJetPartonPy; //! Vector of parton pt matched to true jet
  vector<double> fTrueVectorJetPartonPz; //! Vector of parton pt matched to true jet

  vector<double> fVectorJetEtaPerp; //! vector of jet -eta (opposite eta to original jet)
  vector<double> fVectorJetPhiPerp; //! vector of jet phi + 90 degree (perpendicular to original jet)

  std::map<int, int> MapRecJetsTrueJets; //! Map containing the reconstructed jet index in vector and mapping it to true Jet index

  //-------------------------------
  // Response Matrix handlers
  //-------------------------------
  std::vector<MatrixHandler4D*> fRespMatrixHandlerMesonPt;                 //! Response matrix for true vs. rec pt for each jet pt true vs. rec. bin
  std::vector<MatrixHandler4D*> fRespMatrixHandlerFrag;                    //! Response matrix for meson z_rec vs z_true for each jet pt true vs. rec. bin
  std::vector<MatrixHandler4D*> fRespMatrixHandlerFragTrueJets;            //! Response matrix for meson z_rec_trueJets (true jet pt taken) vs z_true for each jet pt true vs. rec. bin
  std::vector<MatrixHandler4D*> fRespMatrixHandlerMesonInvMass;            //! Response matrix for meson inv. mass and meson pt for each jet pt true vs. rec. bin
  std::vector<MatrixHandler4D*> fRespMatrixHandlerMesonInvMassVsZ;         //! Response matrix for meson inv. mass and meson z for each jet pt true vs. rec. bin
  std::vector<MatrixHandler4D*> fRespMatrixHandlerMesonBackInvMassVsZ;     //! Response matrix for meson inv. mass and meson z for background candidates (mixed evt/rotation) for each jet pt true vs. rec. bin
  std::vector<MatrixHandler4D*> fRespMatrixHandlerMesonBackInvMassVsPt;    //! Response matrix for meson inv. mass and meson pT for background candidates (mixed evt/rotation) for each jet pt true vs. rec. bin
  std::vector<MatrixHandler4D*> fRespMatrixHandlerMesonBackAddInvMassVsZ;  //! Response matrix for meson inv. mass and meson z for background candidates (mixed evt) for each jet pt true vs. rec. bin
  std::vector<MatrixHandler4D*> fRespMatrixHandlerMesonBackAddInvMassVsPt; //! Response matrix for meson inv. mass and meson pT for background candidates (mixed evt) for each jet pt true vs. rec. bin
  std::vector<MatrixHandler4D*> fRespMatrixHandlerMesonInvMassPerpCone;    //! Same as fRespMatrixHandlerMesonInvMass but in perpendicular cone
  std::vector<MatrixHandler4D*> fRespMatrixHandlerMesonInvMassVsZPerpCone; //! Same as fRespMatrixHandlerMesonInvMassVsZ but in perpendicular cone

  // Ndimensional matrix handler used
  std::vector<MatrixHandlerNDim*> fRespMatrixHandlerMesonPtRadius;            //! Response matrix for true vs. rec pt for each jet pt true vs. rec. bin for different radii
  std::vector<MatrixHandlerNDim*> fRespMatrixHandlerMesonPtTrueRadius;        //! Response matrix for true vs. rec pt for each jet pt true vs. rec. bin for different radii
  std::vector<MatrixHandlerNDim*> fRespMatrixHandlerMesonPtInvMassRadius;     //! Reconstructed jet pt, inv. mass, meson pt and radius of jet and meson
  std::vector<MatrixHandlerNDim*> fRespMatrixHandlerMesonBackPtInvMassRadius; //! Reconstructed jet pt, inv. mass, meson pt and radius of jet and meson
  
  //-------------------------------
  // basic histograms
  //-------------------------------
  std::vector<TH1F*> fHistoNEvents;         //! vector of histos with event information
  std::vector<TH1F*> fHistoNEventsWOWeight; //! vector of histos with event information without event weights in case of JJ MC

  std::vector<TH1F*> fHistoNGoodESDTracks;            //! vector of histos for number of tracks
  std::vector<TH1F*> fHistoNGoodESDTracksEvtWithJets; //! vector of histos for number of tracks in events with jets
  std::vector<TH1F*> fHistoVertexZ;                   //! vector of histos for number of events without weights

  //-------------------------------
  // Jet-Jet MC related
  //-------------------------------
  std::vector<TProfile*> fProfileJetJetXSection; //! vector of profiles for Jet-Jet x-section
  std::vector<TH1F*> fHistoJetJetNTrials;        //! vector of histos for Jet-Jet n-trials
  std::vector<TH2F*> fHistoPtHardJJWeight;       //! vector of histos for Jet-Jet weight

  //-------------------------------
  // cluster related histograms
  //-------------------------------
  std::vector<TH1F*> fHistoClusterPt; //! vector of histos with number of clusters as function of pt
  std::vector<TH1F*> fHistoClusterE;  //! vector of histos with number of clusters as function of E

  std::vector<TH1F*> fHistoClusterPtInJet; //! vector of histos with number of clusters as function of pt inside of jets
  std::vector<TH1F*> fHistoClusterEInJet;  //! vector of histos with number of clusters as function of E inside of jets
  std::vector<TH3F*> fHistoClusterPtResolutionInJet;  //! vector of histos with number of clusters as function of E inside of jets

  std::vector<TH2F*> fHistoClusterPtVsJetPtInJet; //! vector of histos with number of clusters as function of pt vs. jet pt inside of jets

  // perpendicular cone
  std::vector<TH1F*> fHistoClusterPtPerpCone; //! vector of histos with number of clusters as function of pt in perpendicular cone

  //-------------------------------
  // conversion photon histograms
  //-------------------------------
  std::vector<TH1F*> fHistoConvGammaPt;             //! vector of histos conversion photons vs. pt
  std::vector<TH1F*> fHistoConvGammaPtInJet;        //! vector of histos conversion photons vs. pt inside of jet
  std::vector<TH2F*> fHistoConvGammaPtVsJetPtInJet; //! vector of histos conversion photons vs. pt  vs Jet Pt inside of jet
  

  // perpendicular cone
  std::vector<TH1F*> fHistoConvGammaPtPerpCone; //! vector of histos conversion photons vs. pt in perpendicular cone

  //-------------------------------
  // Inv. Mass histograms
  //-------------------------------
  std::vector<TH2F*> fHistoInvMassVsPt;          //! vector of histos with inv. mass vs pt
  std::vector<TH2F*> fHistoInvMassVsPt_Incl;     //! vector of histos with inv. mass vs pt for all mesons (no in-jet criterium)
  std::vector<TH2F*> fHistoInvMassVsZ;           //! vector of histos with inv. mass vs z
  std::vector<TH2F*> fHistoJetPtVsFrag;          //! vector of histos for jet pt vs meson-z
  std::vector<TH2F*> fHistoInvMassVsPtMassCut;   //! vector of histos with inv. mass vs. pT after cut
  std::vector<TH2F*> fHistoInvMassVsPtMassCutSB; //! vector of histos with inv. mass vs. pT after cut in Sideband

  std::vector<TH2F*> fHistoInvMassVsPtBack; //! vector of histos  with inv. mass vs pt for background distribution

  std::vector<TH2F*> fHistoInvMassVsPtPerpCone;    //! same as fHistoInvMassVsPt but in perpendicular cone
  std::vector<TH2F*> fHistoJetPtVsFragPerpCone;    //! same as fHistoJetPtVsFrag but in perp cone

  //-------------------------------
  // Jet related histograms
  //-------------------------------
  std::vector<TH1F*> fHistoEventwJets;                //! vector of histos for events with jets
  std::vector<TH1F*> fHistoNJets;                     //! vector of histos with number of jets per event
  std::vector<TH1F*> fHistoPtJet;                     //! vector of histos with pt of jets
  std::vector<TH1F*> fHistoJetEta;                    //! vector of histos with eta of jets
  std::vector<TH1F*> fHistoJetPhi;                    //! vector of histos with phi of jets
  std::vector<TH1F*> fHistoJetArea;                   //! vector of histos with jet area
  std::vector<TH2F*> fHistoTruevsRecJetPt;            //! vector of histos response matrix for jets
  std::vector<TH2F*> fHistoTruevsRecJetPtForTrueJets; //! vector of histos response matrix for true jets
  std::vector<TH2F*> fHistoTrueJetPtVsPartonPt;       //! vector of histos true jet pt vs. parton pt
  std::vector<TH1F*> fHistoMatchedPtJet;              //! vector of histos with pt of jets for jets that got matched with a true jet
  std::vector<TH1F*> fHistoUnMatchedPtJet;            //! vector of histos with pt of jets for jets that did not get matched with a true jet
  std::vector<TH1F*> fHistoTruePtJet;                 //! vector of histos with pt of true jets
  std::vector<TH1F*> fHistoTruePtJetInAcc;            //! vector of histos with pt of true jets for generated jets in rec. acceptance
  std::vector<TH1F*> fHistoTruePtJetNotTriggered;     //! vector of histos with pt of true jets for events that did not trigger
  std::vector<TH1F*> fHistoTrueJetEta;                //! vector of histos with eta of jets
  std::vector<TH1F*> fHistoTrueJetPhi;                //! vector of histos with phi of jets
  std::vector<TH1F*> fHistoTrueMatchedPtJet;          //! vector of histos with pt of true jets that are matched to a rec jet
  std::vector<TH1F*> fHistoTrueUnMatchedPtJet;        //! vector of histos with pt of jets that are not matched to a rec jet
  std::vector<TH2F*> fHistoNEFVsPtJet;                //! vector of histos with pt of jets vs neutral energy fraction
  std::vector<TH2F*> fHistoNchVsPtJet;                //! vector of histos with pt of jets vs number of charged tracks in jet
  std::vector<TH2F*> fHistoNclusVsPtJet;              //! vector of histos with pt of jets vs neutral clusters in jet
  std::vector<TH2F*> fHistoNPartInTrueJetVsJetPt;     //! vector of histos with pt jet vs number of true jet particles
  
  //-------------------------------
  // True meson histograms
  //-------------------------------
  std::vector<MatrixHandler4D*> fRespMatrixHandlerTrueMesonInvMassVsPt;          //! vector of histos inv. mass vs. pT for true mesons
  std::vector<MatrixHandler4D*> fRespMatrixHandlerTrueMesonInvMassVsZ;           //! vector of histos inv. mass vs. Z for true mesons
  std::vector<MatrixHandler4D*> fRespMatrixHandlerTrueOtherMesonInvMassVsPt;     //! vector of histos inv. mass vs. pT for true other mesons (if selected meson is pi0, other mesons are etas and eta prime etc.)
  std::vector<MatrixHandler4D*> fRespMatrixHandlerTrueOtherMesonInvMassVsZ;      //! vector of histos inv. mass vs. Z for true other mesons (if selected meson is pi0, other mesons are etas and eta prime etc.)
  std::vector<MatrixHandler4D*> fRespMatrixHandlerTrueSecondaryMesonInvMassVsPt; //! vector of histos inv. mass vs. pT for true secondary mesons
  std::vector<MatrixHandler4D*> fRespMatrixHandlerTrueSecondaryMesonInvMassVsZ;  //! vector of histos inv. mass vs. pT for true secondary mesons
  std::vector<MatrixHandlerNDim*> fRespMatrixHandlerTrueMesonPtRadius;           //! vector of histos inv. mass vs. pT vs. radius fo true mesons
  std::vector<MatrixHandlerNDim*> fRespMatrixHandlerTrueMesonTruePtRadius;       //! vector of histos inv. mass vs. pT vs. radius fo true mesons
  std::vector<MatrixHandlerNDim*> fRespMatrixHandlerTrueSecondaryMesonPtRadius;  //! vector of histos inv. mass vs. pT vs. radius fo true secondary mesons
  std::vector<TH2F*> fHistoTrueMesonInvMassVsTruePt;                             //! vector of histos inv. mass vs. true pT for true mesons
  std::vector<TH2F*> fHistoTruePrimaryMesonInvMassPt;                            //! vector of histos inv. mass vs. pT for true primary mesons
  std::vector<TH2F*> fHistoTrueSecondaryMesonInvMassPt;                          //! vector of histos inv. mass vs. pT for true secondary mesons
  std::vector<TH2F*> fHistoTrueMesonJetPtVsTruePt;                               //! vector of histos true meson pt vs true jet pt
  std::vector<TH2F*> fHistoTrueMesonJetPtVsTrueZ;                                //! vector of histos true meson z vs true jet pt
  std::vector<TH2F*> fHistoTrueMesonJetPtVsRecPt;                                //! vector of histos rec meson pt vs rec jet pt
  std::vector<TH2F*> fHistoTrueMesonJetPtVsRecZ;                                 //! vector of histos rec meson z vs rec jet pt
  std::vector<TH2F*> fHistoTrueSecMesonJetPtVsRecPt;                             //! vector of histos rec meson pt vs rec jet pt
  std::vector<TH2F*> fHistoTrueSecMesonJetPtVsRecZ;                              //! vector of histos rec meson z vs rec jet pt
  std::vector<TH2F*> fHistoTrueMesonInTrueJet_JetPtVsTruePt;                     //! vector of histos true meson pt vs true jet pt inside true jets
  std::vector<TH2F*> fHistoTrueMesonInTrueJet_JetPtVsTrueZ;                      //! vector of histos true meson z vs true jet pt inside true jets
  std::vector<TH2F*> fHistoMesonResponse;                                        //! vector of histos with meson response matrix
  std::vector<TH3F*> fHistoMesonResolutionJetPt;                                 //! vector of histos with meson resolution as function of jet momentum

  std::vector<TH2F*> fHistoTrueMesonBothDaughtersInJet;                          //! vector of histos with meson pt vs. jet pt for number of mesons with both decay daughters inside jet cone
  std::vector<TH2F*> fHistoTrueMesonOneDaughtersInJet;                           //! vector of histos with meson pt vs. jet pt for number of mesons with only one decay daughters inside jet cone
  std::vector<TH2F*> fHistoTrueMesonNoDaughtersInJet;                            //! vector of histos with meson pt vs. jet pt for number of mesons with no decay daughters inside jet cone
  std::vector<TH2F*> fHistoTrueMesonDaughtersInOtherJet;                         //! vector of histos with meson pt vs. jet pt for number of mesons with decay photons in jet but not in same as pi0
  
  //-------------------------------
  // Meson double counting
  //-------------------------------
  std::vector<int> fMesonDoubleCount;                                              //! keeps track of MC IDs of true mesons and check if one double counts
  std::vector<MatrixHandler4D*> fRespMatrixHandlerTrueMesonInvMassVsPtDoubleCount; //! vector of histos inv. mass vs. pT for true mesons which are counted more than once
  std::vector<MatrixHandler4D*> fRespMatrixHandlerTrueMesonInvMassVsZDoubleCount;  //! vector of histos inv. mass vs. Z for true mesons which are counted more than once

  //-------------------------------
  // True conversion photon histograms
  //-------------------------------
  std::vector<TH1F*> fHistoTrueConvGammaPt;               //! vector of histos true conversion pt
  std::vector<TH1F*> fHistoTruePrimaryConvGammaPt;        //! vector of histos true primary conversion pt
  std::vector<TH2F*> fHistoTruePrimaryConvGammaESDPtMCPt; //! vector of histos true vs. rec. for conversions

  //-------------------------------
  // True cluster histograms
  //-------------------------------
  std::vector<TH1F*> fHistoTrueClusGammaPt; //! vector of histos true EM clusters

  //-------------------------------
  // mc generated histograms
  //-------------------------------
  std::vector<TH1F*> fHistoMCGammaPtNotTriggered;      //! vector of histos with photons in events which are not triggered vs pT
  std::vector<TH1F*> fHistoMCGammaPtNoVertex;          //! vector of histos with photons in events which have no vertex vs pT
  std::vector<TH1F*> fHistoMCAllGammaPt;               //! vector of histos with all photons vs pT
  std::vector<TH1F*> fHistoMCDecayGammaPi0Pt;          //! vector of histos gammas pt from pi0s
  std::vector<TH1F*> fHistoMCDecayGammaRhoPt;          //! vector of histos gammas pt from rho
  std::vector<TH1F*> fHistoMCDecayGammaEtaPt;          //! vector of histos gammas pt from eta
  std::vector<TH1F*> fHistoMCDecayGammaOmegaPt;        //! vector of histos  gammas pt from omega
  std::vector<TH1F*> fHistoMCDecayGammaEtapPt;         //! vector of histos gammas pt from eta prime
  std::vector<TH1F*> fHistoMCDecayGammaPhiPt;          //! vector of histos gammas pt from phi
  std::vector<TH1F*> fHistoMCDecayGammaSigmaPt;        //! vector of histos  gammas pt from sigma
  std::vector<TH2F*> fHistoMCPrimaryPtvsSource;        //! vector of histos primary gammas for different particles
  std::vector<TH1F*> fHistoMCMesonPtNotTriggered;      //! vector of histos mesons which are in events that are not triggered
  std::vector<TH1F*> fHistoMCMesonPtNoVertex;          //! vector of histos mesons which are in events that have no vertex
  std::vector<TH1F*> fHistoMCMesonPt;                  //! vector of histos meson pt
  std::vector<TH1F*> fHistoInclusiveMCMesonPt;         //! vector of histos meson pt for in+outside
  std::vector<TH2F*> fHistoMCMesonPtVsEta;             //! vector of histos meson pt vs. rapidity
  std::vector<TH1F*> fHistoMCMesonWOEvtWeightPt;       //! vector of histos meson pt without event weights
  std::vector<TH1F*> fHistoMCMesonInAccPt;             //! vector of histos mesons in acceptance
  std::vector<TH1F*> fHistoMCMesonInAccPtNotTriggered; //! vector of histos mesons in acceptance which are in events that are not triggered
  std::vector<TH1F*> fHistoMCMesonWOWeightInAccPt;     //! vector of histos mesons in acceptance without event weight
  std::vector<TH1F*> fHistoMCMesonWOEvtWeightInAccPt;  //! vector of histos mesons in acceptance without event weight
  std::vector<TH2F*> fHistoMCSecMesonPtvsSource;       //! vector of histos secondary mesons from different sources vs. pt
  std::vector<TH1F*> fHistoMCSecMesonSource;           //! vector of histos secondary mesons from different sources
  std::vector<TH2F*> fHistoMCSecMesonInAccPtvsSource;  //! vector of histos accepted secondary mesons from different sources vs. pt
  
  //-------------------------------
  // mc generated histograms jet-meson corr
  //-------------------------------
  std::vector<TH2F*> fHistoMCJetPtVsMesonPt;              //! vector of histos True Jet pT vs. true Meson Pt (mc particle based distribution)
  std::vector<TH2F*> fHistoMCJetPtVsMesonPtInAcc;         //! vector of histos True Jet pT vs. true Meson Pt (for mesons in detector acceptance)
  std::vector<TH2F*> fHistoMCJetPtVsFrag;                 //! vector of histos True Jet pT vs. true Frag (mc particle based distribution)
  std::vector<TH2F*> fHistoMCJetPtVsFragInAcc;            //! vector of histos True Jet pT vs. true Frag (for mesons in detector acceptance)
  std::vector<TH2F*> fHistoMCJetPtVsFrag_Sec;             //! vector of histos True Jet pT vs. true Frag (mc particle based distribution for secondaries)
  std::vector<TH2F*> fHistoMCRecJetPtVsFrag;              //! vector of histos with True Jet pT vs true Frag inside rec. jets (to find the meson reconstruction effi)
  std::vector<TH2F*> fHistoMCRecJetPtVsMesonPt;           //! vector of histos with True Jet pT vs true meson pt inside rec. jets (to find the meson reconstruction effi)
  std::vector<TH2F*> fHistoMCJetPtVsMesonPt_Sec;          //! vector of histos True Jet pT vs. true meson pt (mc particle based distribution for secondaries)
  std::vector<TH2F*> fHistoMCPartonPtVsFrag;              //! vector of histos True parton pT vs. true Frag (mc particle based distribution)
  std::vector<TH2F*> fHistoMCJetPtVsFragTrueParton;       //! vector of histos True Jet pT vs. true Frag (mc particle based distribution) for particles originating from hard parton from Jet
  std::vector<TH2F*> fHistoMCPartonPtVsFragTrueParton;    //! vector of histos True parton pT vs. true Frag (mc particle based distribution) for particles originating from hard parton from Jet
  std::vector<TH3F*> fHistoMCJetPtVsMesonPtVsRadius;      //! vector of histos True jet pT vs. true meson pt vs jet-meson radius for particles originating from hard parton from Jet
  std::vector<TH3F*> fHistoMCJetPtVsMesonPtVsRadiusInAcc; //! vector of histos True jet pT vs. true meson pt vs jet-meson radius for particles originating from hard parton from Jet

  //-------------------------------
  // DCA tree for PCM pile-up estimation
  //-------------------------------
  std::vector<TTree*> fDCATree;                       //! vector of trees with the following info: DCA for both photons, pt of meson, mass of meson, corresponding jet pT
  unsigned short fDCATree_InvMass;                    //! inv mass of meson for DCA tree
  unsigned short fDCATree_Pt;                         //! pt of meson for DCA tree
  short fDCATree_DCAzGammaMin;                        //! DCA of photon with smaller DCA for DCA tree
  short fDCATree_DCAzGammaMax;                        //! DCA of photon with larger DCA for DCA tree
  short fDCATree_QualityFlag;                         //! Quality flag (categorie 1 - 6) of meson for DCA tree
  unsigned short fDCATree_JetPt;                      //! jet pt meson belongs to for DCA tree
  bool fDCATree_isTrueMeson;                          //! flag if meson is true meson or not
  float fDCATree_EvtWeight;                           //! event weight for the tree in case of MC

 private:
  static constexpr bool fLocalDebugFlag = false;
  AliAnalysisTaskMesonJetCorrelation(const AliAnalysisTaskMesonJetCorrelation&);            // Prevent copy-construction
  AliAnalysisTaskMesonJetCorrelation& operator=(const AliAnalysisTaskMesonJetCorrelation&); // Prevent assignment

  ClassDef(AliAnalysisTaskMesonJetCorrelation, 21);
};

#endif