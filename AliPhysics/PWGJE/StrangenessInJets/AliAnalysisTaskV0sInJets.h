#ifndef AliAnalysisTaskV0sInJets_cxx
#define AliAnalysisTaskV0sInJets_cxx

// task for analysis of V0s (K0S, (anti-)Lambda) in charged jets
// Author: Vit Kucera (vit.kucera@cern.ch)

class TH1D;
class TH2D;
class THnSparse;
class TRandom;
class TClonesArray;

class AliAODv0;
class AliAODVertex;
class AliAODJet;

#include "AliAnalysisTaskSE.h"
#include "THnSparse.h"
//#include "AuxFunctions.h"

class AliAnalysisTaskV0sInJets : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskV0sInJets(); // Default constructor
  AliAnalysisTaskV0sInJets(const char* name); // Constructor
  virtual ~AliAnalysisTaskV0sInJets(); // Destructor
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t*) {}

  void SetTypeAOD(Int_t type = 1) {fiAODAnalysis = type;}
  void SetIsPbPb(Bool_t val = 1) {fbIsPbPb = val;}
  void SetJetBranchName(char* line) {fsJetBranchName = line;}
  void SetJetBgBranchName(char* line) {fsJetBgBranchName = line;}
  void SetCuts(Double_t z = 10, Double_t r = 1, Double_t cL = 0, Double_t cH = 80) {fdCutVertexZ = z; fdCutVertexR2 = r * r; fdCutCentLow = cL; fdCutCentHigh = cH;}
  void SetPtJetMin(Double_t ptMin = 0) {fdCutPtJetMin = ptMin;}
  void SetPtTrackMin(Double_t ptMin = 0) {fdCutPtTrackMin = ptMin;}
  void SetJetRadius(Double_t r = 0.4) {fdRadiusJet = r;}
  void SetJetRadiusBg(Double_t r = 0.4) {fdRadiusJetBg = r;}
  void SetJetSelection(Bool_t select = kTRUE) {fbJetSelection = select;}
  void SetMCAnalysis(Bool_t select = kTRUE) {fbMCAnalysis = select;}
//  void SetTreeOutput(Bool_t select = kTRUE){fbTreeOutput = select;}
  void FillQAHistogramV0(AliAODVertex* vtx, const AliAODv0* vZero, Int_t iIndexHisto, Bool_t IsCandK0s, Bool_t IsCandLambda, Bool_t IsInPeakK0s, Bool_t IsInPeakLambda);
//  virtual Double_t MassPeakSigma(Double_t pt, Int_t particle);
//  virtual Double_t MassPeakSigma(Int_t iCent, Double_t pt, Int_t particle);
  void FillCandidates(Double_t mK, Double_t mL, Double_t mAL, Bool_t isK, Bool_t isL, Bool_t isAL, Int_t iCut, Int_t iCent);
  Bool_t IsParticleInCone(const AliVParticle* part1, const AliVParticle* part2, Double_t dRMax) const; // decides whether a particle is inside a jet cone
  Bool_t OverlapWithJets(const TClonesArray* array, const AliVParticle* cone, Double_t dDistance) const; // decides whether a cone overlaps with other jets
  AliAODJet* GetRandomCone(const TClonesArray* array, Double_t dEtaConeMax, Double_t dDistance) const; // generate a random cone which does not overlap with selected jets
  AliAODJet* GetMedianCluster(const TClonesArray* array, Double_t dEtaConeMax) const; // get median kt cluster
  Double_t AreaCircSegment(Double_t dRadius, Double_t dDistance) const; // area of circular segment

  void SetCutDCAToPrimVtxMin(Double_t val = 0.1) {fdCutDCAToPrimVtxMin = val;}
  void SetCutDCADaughtersMax(Double_t val = 1.) {fdCutDCADaughtersMax = val;}
  void SetCutNSigmadEdxMax(Double_t val = 3.) {fdCutNSigmadEdxMax = val;}
  void SetCutCPAMin(Double_t val = 0.998) {fdCutCPAMin = val;}
  void SetCutNTauMax(Double_t val = 5.) {fdCutNTauMax = val;}

  Bool_t IsSelectedForJets(AliAODEvent* fAOD, Double_t dVtxZCut, Double_t dVtxR2Cut, Double_t dCentCutLo, Double_t dCentCutUp, Bool_t bCutDeltaZ = kFALSE, Double_t dDeltaZMax = 100.);
  Int_t GetCentralityBinIndex(Double_t centrality);
  Int_t GetCentralityBinEdge(Int_t index);
  TString GetCentBinLabel(Int_t index);
  Double_t MassPeakSigmaOld(Double_t pt, Int_t particle);
  static bool CompareClusters(const std::vector<Double_t> cluster1, const std::vector<Double_t> cluster2); // compare clusters by their pt/area

  // upper edges of centrality bins
  static const Int_t fgkiNBinsCent = 1; // number of centrality bins
  static const Int_t fgkiCentBinRanges[fgkiNBinsCent]; // upper edges of centrality bins
  // axis: pT of V0
  static const Double_t fgkdBinsPtV0[2]; // [GeV/c] minimum and maximum or desired binning of the axis (intended for the rebinned axis)
  static const Int_t fgkiNBinsPtV0; // number of bins (intended for the rebinned axis)
  static const Int_t fgkiNBinsPtV0Init; // initial number of bins (uniform binning)
  // axis: pT of jets
  static const Double_t fgkdBinsPtJet[2]; // [GeV/c] minimum and maximum or desired binning of the axis (intended for the rebinned axis)
  static const Int_t fgkiNBinsPtJet; // number of bins (intended for the rebinned axis)
  static const Int_t fgkiNBinsPtJetInit; // initial number of bins (uniform binning)
  // axis: K0S invariant mass
  static const Int_t fgkiNBinsMassK0s; // number of bins (uniform binning)
  static const Double_t fgkdMassK0sMin; // minimum
  static const Double_t fgkdMassK0sMax; // maximum
  // axis: Lambda invariant mass
  static const Int_t fgkiNBinsMassLambda; // number of bins (uniform binning)
  static const Double_t fgkdMassLambdaMin; // minimum
  static const Double_t fgkdMassLambdaMax; // maximum

private:
  AliAODEvent* fAODIn; //! Input AOD event
  AliAODEvent* fAODOut; //! Output AOD event
  TList* fOutputListStd; //! Output list for standard analysis results
  TList* fOutputListQA; //! Output list for quality assurance
  TList* fOutputListCuts; //! Output list for checking cuts
  TList* fOutputListMC; //! Output list for MC related results
//  TTree* ftreeOut; //! output tree

  Int_t fiAODAnalysis; // switch for input AOD/ESD
  Bool_t fbIsPbPb; // switch Pb-Pb / p-p collisions

  // V0 selection
  Double_t fdCutDCAToPrimVtxMin; // [cm] min DCA of daughters to the prim vtx
  Double_t fdCutDCADaughtersMax; // [sigma of TPC tracking] max DCA between daughters
  Double_t fdCutNSigmadEdxMax; // [sigma dE/dx] max difference between measured and expected signal of dE/dx in the TPC
  Double_t fdCutCPAMin; // min cosine of the pointing angle
  Double_t fdCutNTauMax; // [tau] max proper lifetime in multiples of the mean lifetime
  // jet selection
  TString fsJetBranchName; // name of the branch with jets
  TString fsJetBgBranchName; // name of the branch with kt clusters used for the rho calculation
  Double_t fdCutPtJetMin; // [GeV/c] minimum jet pt
  Double_t fdCutPtTrackMin; // [GeV/c] minimum pt of leading jet-track
  Double_t fdRadiusJet; // R of jet finder used for finding V0s in the jet cone
  Double_t fdRadiusJetBg; // R of kt jet finder used for reconstruction of bg clusters
  Bool_t fbJetSelection; // switch for the analysis of V0s in jets

  Bool_t fbMCAnalysis; // switch for the analysis of simulated data
//  Bool_t fbTreeOutput; // switch for the output tree
  TRandom* fRandom; //! random-number generator

  // event cuts
  Double_t fdCutVertexZ; // [cm] maximum |z| of primary vertex
  Double_t fdCutVertexR2; // [cm^2] maximum r^2 of primary vertex
  Double_t fdCutCentLow; // [%] minimum centrality
  Double_t fdCutCentHigh; // [%] maximum centrality
  /*
  // output branches
  TClonesArray* fBranchV0Rec; //! output branch for reconstructed V0s
  TClonesArray* fBranchV0Gen; //! output branch for generated V0s
  TClonesArray* fBranchJet; //! output branch for selected jets
  AliEventInfoObject* fEventInfo; //! class to store info about events
  */
  Double_t fdCentrality; //!

  // event histograms
  TH1D* fh1EventCounterCut; //! number of events for different selection steps
  TH1D* fh1EventCounterCutCent[fgkiNBinsCent]; //! number of events for different selection steps and different centralities
  TH1D* fh1EventCent; //! number of events for different centralities
  TH1D* fh1EventCent2; //! number of events for different centralities
  TH1D* fh1EventCent2Jets; //! number of events for different centralities
  TH1D* fh1EventCent2NoJets; //! number of events for different centralities
  TH2D* fh2EventCentTracks; //! number of tracks vs centrality
  TH1D* fh1VtxZ[fgkiNBinsCent]; //! z coordinate of the primary vertex
  TH2D* fh2VtxXY[fgkiNBinsCent]; //! xy coordinates of the primary vertex
  TH1D* fh1V0CandPerEvent; //! number of V0 cand per event

  // jet histograms
  TH1D* fh1PtJet[fgkiNBinsCent]; //! pt spectra of jets for normalisation of in-jet V0 spectra
  TH1D* fh1EtaJet[fgkiNBinsCent]; //! jet eta
  TH2D* fh2EtaPtJet[fgkiNBinsCent]; //! jet eta-pT
  TH1D* fh1PhiJet[fgkiNBinsCent]; //! jet phi
  TH1D* fh1NJetPerEvent[fgkiNBinsCent]; //! number of jets per event
  TH1D* fh1NRndConeCent; //! number of generated random cones in centrality bins
  TH2D* fh2EtaPhiRndCone[fgkiNBinsCent]; //! random cone eta-pT
  TH1D* fh1NMedConeCent; //! number of found median-cluster cones in centrality bins
  TH2D* fh2EtaPhiMedCone[fgkiNBinsCent]; //! median-cluster cone eta-phi
  TH1D* fh1AreaExcluded; //! area of excluded cones for outside-cones V0s

  static const Int_t fgkiNCategV0 = 17; // number of V0 selection steps

  // QA histograms
  static const Int_t fgkiNQAIndeces = 2; // 0 - before cuts, 1 - after cuts
  TH1D* fh1QAV0Status[fgkiNQAIndeces]; //! online vs offline reconstructed V0 candidates
  TH1D* fh1QAV0TPCRefit[fgkiNQAIndeces]; //! TPC refit on vs off
  TH1D* fh1QAV0TPCRows[fgkiNQAIndeces]; //! crossed TPC pad rows
  TH1D* fh1QAV0TPCFindable[fgkiNQAIndeces]; //! findable clusters
  TH1D* fh1QAV0TPCRowsFind[fgkiNQAIndeces]; //! ratio rows/clusters
  TH1D* fh1QAV0Eta[fgkiNQAIndeces]; //! pseudorapidity
  TH2D* fh2QAV0EtaRows[fgkiNQAIndeces]; //! pseudorapidity vs TPC rows
  TH2D* fh2QAV0PtRows[fgkiNQAIndeces]; //! pt vs TPC rows
  TH2D* fh2QAV0PhiRows[fgkiNQAIndeces]; //! azimuth vs TPC rows
  TH2D* fh2QAV0NClRows[fgkiNQAIndeces]; //! clusters vs TPC rows
  TH2D* fh2QAV0EtaNCl[fgkiNQAIndeces]; //! pseudorapidity vs clusters

  // K0s
  TH1D* fh1V0CounterCentK0s[fgkiNBinsCent]; //! number of K0s candidates after various cuts
  TH1D* fh1V0InvMassK0sAll[fgkiNCategV0]; //! V0 invariant mass, selection steps
  TH2D* fh2QAV0EtaPtK0sPeak[fgkiNQAIndeces]; //! daughters pseudorapidity vs V0 pt, in mass peak
  TH2D* fh2QAV0EtaEtaK0s[fgkiNQAIndeces]; //! daughters pseudorapidity vs pseudorapidity
  TH2D* fh2QAV0PhiPhiK0s[fgkiNQAIndeces]; //! daughters azimuth vs azimuth
  TH1D* fh1QAV0RapK0s[fgkiNQAIndeces]; //! V0 rapidity
  TH2D* fh2QAV0PtPtK0sPeak[fgkiNQAIndeces]; //! daughters pt vs pt, in mass peak
  TH2D* fh2ArmPodK0s[fgkiNQAIndeces]; //! Armenteros-Podolanski
  TH1D* fh1V0CandPerEventCentK0s[fgkiNBinsCent]; //! number of K0s candidates per event, in centrality bins
  TH1D* fh1V0InvMassK0sCent[fgkiNBinsCent]; //! V0 invariant mass, in centrality bins
  // K0s Inclusive
  THnSparse* fhnV0InclusiveK0s[fgkiNBinsCent]; //! V0 inv mass vs pt before and after cuts, in centrality bins
  // K0s Cones
  THnSparse* fhnV0InJetK0s[fgkiNBinsCent]; //! V0 invariant mass vs V0 pt vs jet pt, in centrality bins
  THnSparse* fhnV0InPerpK0s[fgkiNBinsCent]; //! V0 invariant mass vs V0 pt vs jet pt, in centrality bins
  THnSparse* fhnV0InRndK0s[fgkiNBinsCent]; //! V0 invariant mass vs V0 pt vs jet pt, in centrality bins
  THnSparse* fhnV0InMedK0s[fgkiNBinsCent]; //! V0 invariant mass vs V0 pt vs jet pt, in centrality bins
  THnSparse* fhnV0OutJetK0s[fgkiNBinsCent]; //! V0 invariant mass vs V0 pt, in centrality bins
  THnSparse* fhnV0NoJetK0s[fgkiNBinsCent]; //! V0 invariant mass vs V0 pt, in centrality bins

  TH2D* fh2V0PtJetAngleK0s[fgkiNBinsCent]; //! pt jet vs angle V0-jet, in centrality bins
  TH1D* fh1DCAInK0s[fgkiNBinsCent]; //! DCA between daughters of V0 inside jets, in centrality bins
  TH1D* fh1DCAOutK0s[fgkiNBinsCent]; //! DCA between daughters of V0 outside jets, in centrality bins
//  TH1D* fh1DeltaZK0s[fgkiNBinsCent]; //! z-distance between V0 vertex and primary vertex, in centrality bins
  // MC histograms
  // inclusive
  TH1D* fh1V0K0sPtMCGen[fgkiNBinsCent]; //! pt spectrum of all generated K0s in event
  TH2D* fh2V0K0sPtMassMCRec[fgkiNBinsCent]; //! pt-mass spectrum of successfully reconstructed K0s in event
  TH1D* fh1V0K0sPtMCRecFalse[fgkiNBinsCent]; //! pt spectrum of false reconstructed K0s in event
  // inclusive eta-pT efficiency
  TH2D* fh2V0K0sEtaPtMCGen[fgkiNBinsCent]; //! eta-pt spectrum of all generated K0s in event
  THnSparse* fh3V0K0sEtaPtMassMCRec[fgkiNBinsCent]; //! eta-pt-mass spectrum of successfully reconstructed K0s in event
  // MC daughter eta inclusive
//  THnSparse* fhnV0K0sInclDaughterEtaPtPtMCGen[fgkiNBinsCent]; //! eta_daughter-pt_daughter-pt_V0 generated
  THnSparse* fhnV0K0sInclDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! eta_daughter-pt_daughter-pt_V0 reconstructed
  // in jets
  TH2D* fh2V0K0sInJetPtMCGen[fgkiNBinsCent]; //! pt spectrum of generated K0s in jet
  THnSparse* fh3V0K0sInJetPtMassMCRec[fgkiNBinsCent]; //! mass-pt spectrum of successfully reconstructed K0s in jet
  // in jets eta-pT efficiency
  THnSparse* fh3V0K0sInJetEtaPtMCGen[fgkiNBinsCent]; //! eta-pt spectrum of generated K0s in jet
  THnSparse* fh4V0K0sInJetEtaPtMassMCRec[fgkiNBinsCent]; //! mass-eta-pt spectrum of successfully reconstructed K0s in jet
  // MC daughter eta in JC
//  THnSparse* fhnV0K0sInJetsDaughterEtaPtPtMCGen[fgkiNBinsCent]; //! eta_daughter-pt_daughter-pt_V0 generated
  THnSparse* fhnV0K0sInJetsDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! eta_daughter-pt_daughter-pt_V0 reconstructed

  // resolution
  TH2D* fh2V0K0sMCResolMPt[fgkiNBinsCent]; //! K0s mass resolution vs pt
  TH2D* fh2V0K0sMCPtGenPtRec[fgkiNBinsCent]; //! K0s generated pt vs reconstructed pt

  // Lambda
  TH1D* fh1V0CounterCentLambda[fgkiNBinsCent]; //! number of Lambda candidates after various cuts
  TH1D* fh1V0InvMassLambdaAll[fgkiNCategV0]; //!
  TH2D* fh2QAV0EtaPtLambdaPeak[fgkiNQAIndeces]; //!
  TH2D* fh2QAV0EtaEtaLambda[fgkiNQAIndeces]; //!
  TH2D* fh2QAV0PhiPhiLambda[fgkiNQAIndeces]; //!
  TH1D* fh1QAV0RapLambda[fgkiNQAIndeces]; //!
  TH2D* fh2QAV0PtPtLambdaPeak[fgkiNQAIndeces]; //!
  TH2D* fh2ArmPodLambda[fgkiNQAIndeces]; //!
  TH1D* fh1V0CandPerEventCentLambda[fgkiNBinsCent]; //!
  TH1D* fh1V0InvMassLambdaCent[fgkiNBinsCent]; //!
  // Lambda Inclusive
  THnSparse* fhnV0InclusiveLambda[fgkiNBinsCent]; //!
  // Lambda Cones
  THnSparse* fhnV0InJetLambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0InPerpLambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0InRndLambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0InMedLambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0OutJetLambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0NoJetLambda[fgkiNBinsCent]; //!

  TH2D* fh2V0PtJetAngleLambda[fgkiNBinsCent]; //!
  TH1D* fh1DCAInLambda[fgkiNBinsCent]; //!
  TH1D* fh1DCAOutLambda[fgkiNBinsCent]; //!
//  TH1D* fh1DeltaZLambda[fgkiNBinsCent]; //!
  // MC histograms
  // inclusive
  TH1D* fh1V0LambdaPtMCGen[fgkiNBinsCent]; //!
  TH2D* fh2V0LambdaPtMassMCRec[fgkiNBinsCent]; //!
  TH1D* fh1V0LambdaPtMCRecFalse[fgkiNBinsCent]; //!
  // inclusive eta-pT efficiency
  TH2D* fh2V0LambdaEtaPtMCGen[fgkiNBinsCent]; //!
  THnSparse* fh3V0LambdaEtaPtMassMCRec[fgkiNBinsCent]; //!
  // MC daughter eta inclusive
//  THnSparse* fhnV0LambdaInclDaughterEtaPtPtMCGen[fgkiNBinsCent]; //! eta_daughter-pt_daughter-pt_V0 generated
  THnSparse* fhnV0LambdaInclDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! eta_daughter-pt_daughter-pt_V0 reconstructed
  // in jets
  TH2D* fh2V0LambdaInJetPtMCGen[fgkiNBinsCent]; //!
  THnSparse* fh3V0LambdaInJetPtMassMCRec[fgkiNBinsCent]; //!
  // in jets eta-pT efficiency
  THnSparse* fh3V0LambdaInJetEtaPtMCGen[fgkiNBinsCent]; //!
  THnSparse* fh4V0LambdaInJetEtaPtMassMCRec[fgkiNBinsCent]; //!
  // MC daughter eta in JC
//  THnSparse* fhnV0LambdaInJetsDaughterEtaPtPtMCGen[fgkiNBinsCent]; //! eta_daughter-pt_daughter-pt_V0 generated
  THnSparse* fhnV0LambdaInJetsDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! eta_daughter-pt_daughter-pt_V0 reconstructed

  // resolution
  TH2D* fh2V0LambdaMCResolMPt[fgkiNBinsCent]; //!
  TH2D* fh2V0LambdaMCPtGenPtRec[fgkiNBinsCent]; //!
  // feed-down
  THnSparseD* fhnV0LambdaInclMCFD[fgkiNBinsCent]; //!
  THnSparseD* fhnV0LambdaInJetsMCFD[fgkiNBinsCent]; //!
  THnSparseD* fhnV0LambdaBulkMCFD[fgkiNBinsCent]; //!
  TH1D* fh1V0XiPtMCGen[fgkiNBinsCent]; //!

  // ALambda
  TH1D* fh1V0CounterCentALambda[fgkiNBinsCent]; //! number of ALambda candidates after various cuts
  TH1D* fh1V0InvMassALambdaAll[fgkiNCategV0]; //!
  TH2D* fh2QAV0EtaPtALambdaPeak[fgkiNQAIndeces]; //!
  TH2D* fh2QAV0EtaEtaALambda[fgkiNQAIndeces]; //!
  TH2D* fh2QAV0PhiPhiALambda[fgkiNQAIndeces]; //!
  TH1D* fh1QAV0RapALambda[fgkiNQAIndeces]; //!
  TH2D* fh2QAV0PtPtALambdaPeak[fgkiNQAIndeces]; //!
  TH2D* fh2ArmPodALambda[fgkiNQAIndeces]; //!
  TH1D* fh1V0CandPerEventCentALambda[fgkiNBinsCent]; //!
  TH1D* fh1V0InvMassALambdaCent[fgkiNBinsCent]; //!
  TH1D* fh1V0ALambdaPt[fgkiNBinsCent]; //!
  // ALambda Inclusive
  THnSparse* fhnV0InclusiveALambda[fgkiNBinsCent]; //!
  // ALambda Cones
  THnSparse* fhnV0InJetALambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0InPerpALambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0InRndALambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0InMedALambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0OutJetALambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0NoJetALambda[fgkiNBinsCent]; //!

  TH2D* fh2V0PtJetAngleALambda[fgkiNBinsCent]; //!
  TH1D* fh1DCAInALambda[fgkiNBinsCent]; //!
  TH1D* fh1DCAOutALambda[fgkiNBinsCent]; //!
//  TH1D* fh1DeltaZALambda[fgkiNBinsCent]; //!
  // MC histograms
  // inclusive
  TH1D* fh1V0ALambdaPtMCGen[fgkiNBinsCent]; //!
  TH1D* fh1V0ALambdaPtMCRec[fgkiNBinsCent]; //!
  TH2D* fh2V0ALambdaPtMassMCRec[fgkiNBinsCent]; //!
  TH1D* fh1V0ALambdaPtMCRecFalse[fgkiNBinsCent]; //!
  // inclusive eta-pT efficiency
  TH2D* fh2V0ALambdaEtaPtMCGen[fgkiNBinsCent]; //!
  THnSparse* fh3V0ALambdaEtaPtMassMCRec[fgkiNBinsCent]; //!
  // MC daughter eta inclusive
//  THnSparse* fhnV0ALambdaInclDaughterEtaPtPtMCGen[fgkiNBinsCent]; //! eta_daughter-pt_daughter-pt_V0 generated
  THnSparse* fhnV0ALambdaInclDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! eta_daughter-pt_daughter-pt_V0 reconstructed
  // in jets
  TH2D* fh2V0ALambdaInJetPtMCGen[fgkiNBinsCent]; //!
  TH2D* fh2V0ALambdaInJetPtMCRec[fgkiNBinsCent]; //!
  THnSparse* fh3V0ALambdaInJetPtMassMCRec[fgkiNBinsCent]; //!
  // in jets eta-pT efficiency
  THnSparse* fh3V0ALambdaInJetEtaPtMCGen[fgkiNBinsCent]; //!
  THnSparse* fh4V0ALambdaInJetEtaPtMassMCRec[fgkiNBinsCent]; //!
  // MC daughter eta in JC
//  THnSparse* fhnV0ALambdaInJetsDaughterEtaPtPtMCGen[fgkiNBinsCent]; //! eta_daughter-pt_daughter-pt_V0 generated
  THnSparse* fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! eta_daughter-pt_daughter-pt_V0 reconstructed

  // resolution
  TH2D* fh2V0ALambdaMCResolMPt[fgkiNBinsCent]; //!
  TH2D* fh2V0ALambdaMCPtGenPtRec[fgkiNBinsCent]; //!
  // feed-down
  THnSparseD* fhnV0ALambdaInclMCFD[fgkiNBinsCent]; //!
  THnSparseD* fhnV0ALambdaInJetsMCFD[fgkiNBinsCent]; //!
  THnSparseD* fhnV0ALambdaBulkMCFD[fgkiNBinsCent]; //!
  TH1D* fh1V0AXiPtMCGen[fgkiNBinsCent]; //!

  TH1D* fh1QAV0Pt[fgkiNQAIndeces]; //! pt
  TH1D* fh1QAV0Charge[fgkiNQAIndeces]; //! charge
  TH1D* fh1QAV0DCAVtx[fgkiNQAIndeces]; //! DCA of daughters to prim vtx
  TH1D* fh1QAV0DCAV0[fgkiNQAIndeces]; //! DCA between daughters
  TH1D* fh1QAV0Cos[fgkiNQAIndeces]; //! cosine of pointing angle (CPA)
  TH1D* fh1QAV0R[fgkiNQAIndeces]; //! radial distance between prim vtx and decay vertex
  TH1D* fh1QACTau2D[fgkiNQAIndeces]; //! lifetime calculated in xy
  TH1D* fh1QACTau3D[fgkiNQAIndeces]; //! lifetime calculated in xyz
  TH2D* fh2ArmPod[fgkiNQAIndeces]; //! Armenteros-Podolanski
  TH2D* fh2CCK0s; //! K0s candidates in Lambda peak
  TH2D* fh2CCLambda; //! Lambda candidates in K0s peak
  THnSparse* fh3CCMassCorrelBoth; //! mass correlation of candidates
  THnSparse* fh3CCMassCorrelKNotL; //! mass correlation of candidates
  THnSparse* fh3CCMassCorrelLNotK; //! mass correlation of candidates

  // Cut tuning
  // crossed/findable, daughter pt, dca, cpa, r, pseudorapidity, y, decay length, PID sigma
  /*
  TH2D* fh2CutTPCRowsK0s[fgkiNQAIndeces]; //! inv mass vs TPC rows
  TH2D* fh2CutTPCRowsLambda[fgkiNQAIndeces]; //!
  TH2D* fh2CutPtPosK0s[fgkiNQAIndeces]; //! inv mass vs pt of positive daughter
  TH2D* fh2CutPtNegK0s[fgkiNQAIndeces]; //! inv mass vs pt of negative daughter
  TH2D* fh2CutPtPosLambda[fgkiNQAIndeces]; //!
  TH2D* fh2CutPtNegLambda[fgkiNQAIndeces]; //!
  TH2D* fh2CutDCAVtx[fgkiNQAIndeces]; //! inv mass vs DCA of daughters to prim vtx
  TH2D* fh2CutDCAV0[fgkiNQAIndeces]; //! inv mass vs DCA between daughters
  TH2D* fh2CutCos[fgkiNQAIndeces]; //! inv mass vs CPA
  TH2D* fh2CutR[fgkiNQAIndeces]; //! inv mass vs R
  TH2D* fh2CutEtaK0s[fgkiNQAIndeces]; //! inv mass vs pseudorapidity
  TH2D* fh2CutEtaLambda[fgkiNQAIndeces]; //!
  TH2D* fh2CutRapK0s[fgkiNQAIndeces]; //! inv mass vs rapidity
  TH2D* fh2CutRapLambda[fgkiNQAIndeces]; //!
  TH2D* fh2CutCTauK0s[fgkiNQAIndeces]; //! inv mass vs lifetime
  TH2D* fh2CutCTauLambda[fgkiNQAIndeces]; //!
  TH2D* fh2CutPIDPosK0s[fgkiNQAIndeces]; //! inv mass vs number of dE/dx sigmas for positive daughter
  TH2D* fh2CutPIDNegK0s[fgkiNQAIndeces]; //! inv mass vs number of dE/dx sigmas for negative daughter
  TH2D* fh2CutPIDPosLambda[fgkiNQAIndeces]; //!
  TH2D* fh2CutPIDNegLambda[fgkiNQAIndeces]; //!

  TH2D* fh2Tau3DVs2D[fgkiNQAIndeces]; //! pt vs ratio 3D lifetime / 2D lifetime
  */

  AliAnalysisTaskV0sInJets(const AliAnalysisTaskV0sInJets&); // not implemented
  AliAnalysisTaskV0sInJets& operator=(const AliAnalysisTaskV0sInJets&); // not implemented

  ClassDef(AliAnalysisTaskV0sInJets, 3) // example of analysis
};

#endif
