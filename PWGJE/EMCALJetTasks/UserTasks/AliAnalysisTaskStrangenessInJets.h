#ifndef AliAnalysisTaskStrangenessInJets_cxx
#define AliAnalysisTaskStrangenessInJets_cxx

//-------------------------------------------------------------------------
// Task for V0 analysis in charged jets with 
// the strange particles (instead of daughters) added to the  jet finder
// Author: Ekaterina Grecka (ermeeka@fjfi.cvut.cz)
// Modification of the AlianalysisTaskV0sInJetsEmcal task (author Vit Kucera) 
//-------------------------------------------------------------------------


class TH1D;
class TH2D;
class THnSparse;
class TClonesArray;
class TRandom;

class AliAODv0;
class AliAODVertex;
class AliAODcascade;

class AliParticleContainer;
class AliClusterContainer;

class AliEventPoolManager;
class AliEventCuts;

#include "AliFJWrapper.h"
#include "FJ_includes.h"

#include "AliAnalysisTaskEmcal.h"

namespace fastjet {
  class PseudoJet;
}

class AliAnalysisTaskStrangenessInJets : public AliAnalysisTaskEmcal
{
public:
  
#if !defined(__CINT__) && !defined(__MAKECINT__) 
  typedef fastjet::JetAlgorithm FJJetAlgo;
  typedef fastjet::RecombinationScheme FJRecoScheme;
#endif

  AliAnalysisTaskStrangenessInJets(); // Default constructor
  AliAnalysisTaskStrangenessInJets(const char* name); // Constructor
  virtual ~AliAnalysisTaskStrangenessInJets(); // Destructor
  void UserCreateOutputObjects();
  void Terminate(Option_t*) {}

  void SetIsPbPb(Bool_t val = 1) {fbIsPbPb = val;}
  void SetMCAnalysis(Bool_t select = kTRUE) {fbMCAnalysis = select;}
  void SetGeneratorName(TString name) {fsGeneratorName = name;}
  Bool_t IsFromGoodGenerator(Int_t index); // True if the MC particle with the given index comes from the selected generator


  // Event selection setters 
  void SetEventCuts(Double_t z = 10, Double_t r = 1, Double_t cL = 0, Double_t cH = 80, Double_t dZ = 0.1, Int_t iNC = 1) {fdCutVertexZ = z; fdCutVertexR2 = r * r; fdCutCentLow = cL; fdCutCentHigh = cH; fdCutDeltaZMax = dZ; fiNContribMin = iNC;} 
  void SetUseMultiplicity(Bool_t val = kTRUE) {fbUseMultiplicity = val;} 
  void SetUseIonutCut(Bool_t val = kTRUE) {fbUseIonutCut = val;}
  // V0 selection setters 
  void SetCutTPCRefit(Bool_t val = kTRUE) {fbTPCRefit = val;}
  void SetCutRejectKinks(Bool_t val = kTRUE) {fbRejectKinks = val;}
  void SetCutFindableClusters(Bool_t val = kTRUE) {fbFindableClusters = val;}
  void SetCutV0PtMin(Double_t val = 1.) {fdCutV0PtMin = val;}
  void SetCutNCrossedRowsTPCMin(Double_t val = 70.) {fdCutNCrossedRowsTPCMin = val;}
  void SetCutCrossedRowsOverFindMin(Double_t val = 0.8) {fdCutCrossedRowsOverFindMin = val;}
  void SetCutCrossedRowsOverFindMax(Double_t val = 1e3) {fdCutCrossedRowsOverFindMax = val;}
  void SetCutPtDaughterMin(Double_t val = 0.150) {fdCutPtDaughterMin = val;} 
  void SetCutDCAToPrimVtxMin(Double_t val = 0.1) {fdCutDCAToPrimVtxMin = val;}  
  void SetCutDCADaughtersMax(Double_t val = 1.) {fdCutDCADaughtersMax = val;}
  void SetCutEtaDaughterMax(Double_t val = 0.8) {fdCutEtaDaughterMax = val;}
  void SetCutNSigmadEdxMax(Double_t val = 3.) {fdCutNSigmadEdxMax = val;}
  void SetPtProtonPIDMax(Double_t val = 1.) {fdPtProtonPIDMax = val;}
  void SetOnFly(Bool_t val = 0) {fbOnFly = val;}
  void SetCutCPAKMin(Double_t val = 0.998) {fdCutCPAKMin = val;}
  void SetCutCPALMin(Double_t val = 0.998) {fdCutCPALMin = val;}
  void SetCutRadiusDecayMin(Double_t val = 5.) {fdCutRadiusDecayMin = val;}
  void SetCutRadiusDecayMax(Double_t val = 100.) {fdCutRadiusDecayMax = val;}
  void SetCutEtaV0Max(Double_t val = 0.7) {fdCutEtaV0Max = val;}
  void SetCutRapV0Max(Double_t val = 0.75) {fdCutRapV0Max = val;}
  void SetCutNTauKMax(Double_t val = 5.0) {fdCutNTauKMax = val;}
  void SetCutNTauLMax(Double_t val = 5.0) {fdCutNTauLMax = val;}
  void SetCutArmPod(Bool_t val = kTRUE) {fbCutArmPod = val;}
  void SetCutCross(Bool_t val = kTRUE) {fbCutCross = val;}

  //jet analysis selection: 
  void SetGhostArea(Double_t gharea) { fdGhostArea = gharea; }
  void SetRadius(Double_t r) { fdRadius = r; } 
  void SetMinJetArea(Double_t a)                  { fdMinJetArea       = a     ; }
  void SetMinJetPt(Double_t j)                    { fdMinJetPt         = j     ; }
  void SetJetEtaRange(Double_t emi, Double_t ema) { fdJetEtaMin        = emi   ; fdJetEtaMax = ema; }
  void SetJetPhiRange(Double_t pmi, Double_t pma) { fdJetPhiMin        = pmi   ; fdJetPhiMax = pma; }
  void SetMinJetTrackPt(Double_t j)               { fdJetTrackPtMin        = j     ; }
  void SetMaxJetTrackEta(Double_t j)               { fdJetTrackEtaMax      = j     ; }
  void SetMCDistPrimaryMax(Double_t val = 0.01) {fdDistPrimaryMax = val;}
  void SetDistanceV0JetMax(Double_t val = 0.4) {fdDistanceV0JetMax = val;}
  void SetPtJetMin(Double_t ptMin = 0) {fdCutPtJetMin = ptMin;}
  void SetPtTrackJetMin(Double_t ptMin = 0) {fdCutPtTrackJetMin = ptMin;}
  void SetAreaPercJetMin(Double_t area = 0) {fdCutAreaPercJetMin = area;}
  
  //getters
  Bool_t GetIsPbPb() const         { return fbIsPbPb; }
  Bool_t GetMCAnalysis() const     { return fbMCAnalysis; }
  Bool_t GetUseMultiplicity() const   { return fbUseMultiplicity; }
  Bool_t GetUseIonutCut() const   { return fbUseIonutCut; }  
  Bool_t GetCutTPCRefit() const   { return fbTPCRefit; }  
  Bool_t GetCutRejectKinks() const   { return fbRejectKinks; }
  Bool_t GetCutFindableClusters() const   { return fbFindableClusters; }
  Double_t GetCutNCrossedRowsTPCMin() const   { return fdCutNCrossedRowsTPCMin; }
  Double_t GetCutCrossedRowsOverFindMin() const   { return fdCutCrossedRowsOverFindMin; }  
  Double_t GetCutCrossedRowsOverFindMax() const   { return fdCutCrossedRowsOverFindMax; }
  Double_t GetCutPtDaughterMin() const   { return fdCutPtDaughterMin; }
  Double_t GetCutDCAToPrimVtxMin() const   { return fdCutDCAToPrimVtxMin; }
  Double_t GetCutDCADaughtersMax() const   { return fdCutDCADaughtersMax; }
  Double_t GetCutEtaDaughterMax() const   { return fdCutEtaDaughterMax; }
  Double_t GetCutNSigmadEdxMax() const   { return fdCutNSigmadEdxMax; }
  Double_t GetPtProtonPIDMax() const   { return fdPtProtonPIDMax; } 
  Bool_t GetOnFly() const   { return fbOnFly; } 
  Double_t GetCutCPAKMin() const   { return fdCutCPAKMin; }
  Double_t GetCutCPALMin() const   { return fdCutCPALMin; }
  Double_t GetCutRadiusDecayMin() const   { return fdCutRadiusDecayMin; }
  Double_t GetCutRadiusDecayMax() const   { return fdCutRadiusDecayMax; }
  Double_t GetCutEtaV0Max() const   { return fdCutEtaV0Max; }
  Double_t GetCutRapV0Max() const   { return fdCutRapV0Max; }
  Double_t GetCutNTauKMax() const   { return fdCutNTauKMax; }
  Double_t GetCutNTauLMax() const   { return fdCutNTauLMax; }
  Bool_t GetCutArmPod() const   { return fbCutArmPod; } 
  Bool_t GetCutCross() const   { return fbCutCross; } 
  Double_t GetGhostArea() { return fdGhostArea; }
  Double_t GetRadius() { return fdRadius; }

  // upper edges of centrality bins
  static const Int_t fgkiNBinsCent = 1; // number of centrality bins
  static const Int_t fgkiCentBinRanges[fgkiNBinsCent]; // upper edges of centrality bins
  // axis: pT of V0
  static const Double_t fgkdBinsPtV0[2]; // [GeV/c] minimum and maximum or desired binning of the axis (intended for the rebinned axis)
  static const Int_t fgkiNBinsPtV0; // number of bins (intended for the rebinned axis)
  static const Int_t fgkiNBinsPtV0Init; // initial number of bins (uniform binning)
  static const Int_t fgkiNBinsPtV0InitInJet; // initial number of bins for V0s in jets (uniform binning)
  // axis: K0S invariant mass
  static const Int_t fgkiNBinsMassK0s; // number of bins (uniform binning)
  static const Double_t fgkdMassK0sMin; // minimum K0S mass
  static const Double_t fgkdMassK0sMax; // maximum K0S mass
  // axis: Lambda invariant mass
  static const Int_t fgkiNBinsMassLambda; // number of bins (uniform binning)
  static const Double_t fgkdMassLambdaMin; // minimum Lambda mass
  static const Double_t fgkdMassLambdaMax; // maximum Lambda mass
  // axis: pT of jets
  static const Double_t fgkdBinsPtJet[2]; // [GeV/c] minimum and maximum or desired binning of the axis (intended for the rebinned axis)
  static const Int_t fgkiNBinsPtJet; // number of bins (intended for the rebinned axis)
  static const Int_t fgkiNBinsPtJetInit; // initial number of bins (uniform binning)
  // PDG codes of used particles
  static const Int_t iPdgCodePion;
  static const Int_t iPdgCodeProton;
  static const Int_t iPdgCodeK0s;
  static const Int_t iPdgCodeLambda;

  //Index codes to distinguish particles in the jet 
  static const Int_t iK0Id;
  static const Int_t iLambdaId;
  static const Int_t iALambdaId;
  static const Int_t iK0LId;
  static const Int_t iK0ALId;

  void FillCandidates(Double_t mK, Double_t mL, Double_t mAL, Bool_t isK, Bool_t isL, Bool_t isAL, Int_t iCut, Int_t iCent); //Fills histograms according to the V0 type 
  Bool_t IsParticleInCone(const AliVParticle* part1, const AliVParticle* part2, Double_t dRMax) const; // decides whether a particle is inside a jet cone
  Bool_t OverlapWithJets(const TClonesArray* array, const AliVParticle* cone, Double_t dDistance) const; // decides whether a cone overlaps with other jets
  AliAODJet* GetRandomCone(const TClonesArray* array, Double_t dEtaConeMax, Double_t dDistance) const; // generate a random cone which does not overlap with selected jets
  Bool_t IsSelectedForAnalysis(); // Function for the event selection
  Int_t GetCentralityBinIndex(Double_t centrality);
  Int_t GetCentralityBinEdge(Int_t index);
  TString GetCentBinLabel(Int_t index);
  static Double_t AddDaughters(AliAODRecoDecay* cand, TObjArray& daughters);
  Bool_t AssociateRecV0withMC( AliAODv0* v, AliEmcalJet *xjet, Bool_t bIsK, Bool_t bIsL, Bool_t bIsAL, Int_t iCent);
  Bool_t GeneratedMCParticles( TClonesArray* track, Int_t iCent );
  Double_t MassPeakSigma(Double_t pt, Int_t particle);
  Double_t MassPeakMean(Double_t pt, Int_t particle);
  
protected:
  void ExecOnce();
  Bool_t FillHistograms();
  Bool_t Run();
  void AddEventTracks(TClonesArray* coll, TClonesArray* tracks, std::vector<fastjet::PseudoJet>& VectorBgPart);  
  void AddEventTracksMC(TClonesArray* coll, TClonesArray* tracks, std::vector<fastjet::PseudoJet>& VectorBgPartMC);  
  Bool_t GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const;

  TList* fOutputListStd; //! Output list for standard analysis results
  TList* fOutputListStdJets; //! Output list for jet analysis results
  TList* fOutputListMC;  //! Output list for MC analysis results
  TClonesArray   *fV0CandidateArray;          //! contains selected V0 candidates
  TClonesArray   *fGenMCV0;                    //! contains MC generated V0s
 
  Int_t           fNCand;                   //! number of selected V0 candidates already added to fCandidateArray

  TClonesArray          *fJets;                   //!<!jet collection
  
  AliFJWrapper           fFastJetWrapper;         //!<!fastjet wrapper 
  //AliFJWrapper           fFastJetWrapperBG;       //!<!fastjet wrapper for the bg jets
  AliFJWrapper           fFastJetWrapperMCGen;    //!<!fastjet wrapper for the bg jets

  //std::vector<fastjet::PseudoJet> InputBgParticles;
  fastjet::Selector selectorBG;

private:
  AliAODEvent* fAODIn; //! Input AOD event
  AliAODEvent* fAODOut; //! Output AOD event
  AliMCEvent* fEventMC; //! MC event
  TRandom* fRandom; //! random-number generator
  AliEventCuts fEventCutsStrictAntipileup; //! Event cuts class
  
  static const Int_t fgkiNCategV0 = 18; // number of V0 selection steps
 
 // Data selection
  Bool_t fbIsPbPb; // switch: Pb+Pb / p+p collisions
  Bool_t fbMCAnalysis; // switch: simulated / real data
  TString fsGeneratorName; // pattern for selecting only V0s from a specific MC generator
  // Event selection
  Double_t fdCutVertexZ; // [cm] maximum |z| of primary vertex
  Double_t fdCutVertexR2; // [cm^2] maximum r^2 of primary vertex
  Double_t fdCutCentLow; // [%] minimum centrality
  Double_t fdCutCentHigh; // [%] maximum centrality
  Double_t fdCutDeltaZMax; // [cm] maximum |Delta z| between nominal prim vtx and SPD vtx
  Int_t fiNContribMin; // minimum number of prim vtx contributors
  Double_t fdCentrality; //! [%] centrality
  Bool_t fbUseMultiplicity; // switch for getting centrality from AliMultSelection instead of from AliCentrality
  Bool_t fbUseIonutCut; // Ionut's cut on the correlation between event variables
  // V0 selection
  // Daughter tracks
  Bool_t fbTPCRefit; // (yes) TPC refit for daughter tracks
  Bool_t fbRejectKinks; // (no) reject kink-like production vertices of daughter tracks
  Bool_t fbFindableClusters; // (no) require positive number of findable clusters
  Double_t fdCutNCrossedRowsTPCMin; // (70.) min number of crossed TPC rows
  Double_t fdCutCrossedRowsOverFindMin; // (0.8) min ratio crossed rows / findable clusters
  Double_t fdCutCrossedRowsOverFindMax; // (1e3) max ratio crossed rows / findable clusters
  Double_t fdCutPtDaughterMin; // (0.150) [GeV/c] min transverse momentum of daughter tracks, to reject primaries which do not make it to the TPC
  Double_t fdCutDCAToPrimVtxMin; // (0.1) [cm] min DCA of daughters to the prim vtx
  Double_t fdCutDCADaughtersMax; // (1.) [sigma of TPC tracking] max DCA between daughters
  Double_t fdCutEtaDaughterMax; // (0.8) max |pseudorapidity| of daughter tracks, historical reasons: tracking in MC for 2010 was restricted to 0.7
  Double_t fdCutNSigmadEdxMax; // (3.) [sigma dE/dx] max difference between measured and expected signal of dE/dx in the TPC
  Double_t fdPtProtonPIDMax; // (1.) [GeV/c] maxium pT of proton for applying PID cut in Pb-Pb

  Bool_t fbCascadeOn; // switch: cascade on/off

  // V0 candidate
  Bool_t fbOnFly; // (0) on-the-fly (yes) or offline (no) reconstructed
  Double_t fdCutV0PtMin; // (1 GeV) min pT of the selected v0 particles
  Double_t fdCutCPAKMin; // (0.998) min cosine of the pointing angle, K0S
  Double_t fdCutCPALMin; // (0.998) min cosine of the pointing angle, Lambda
  Double_t fdCutRadiusDecayMin; // (5.) [cm] min radial distance of the decay vertex
  Double_t fdCutRadiusDecayMax; // (100.) [cm] max radial distance of the decay vertex
  Double_t fdCutEtaV0Max; // (0.7) max |pseudorapidity| of V0
  Double_t fdCutRapV0Max; // (0.75) max |rapidity| of V0 (turned off)
  Double_t fdCutNTauKMax; // (5.0) [tau] max proper lifetime in multiples of the mean lifetime, K0S
  Double_t fdCutNTauLMax; // (5.0) [tau] max proper lifetime in multiples of the mean lifetime, Lambda
  Bool_t fbCutArmPod; // (yes) Armenteros-Podolanski for K0S
  Bool_t fbCutCross; // (no) cross-contamination

//jet analysis variables 
  Double_t fdGhostArea;              ///< ghost area  
  Double_t fdRadius;                 ///< jet radius
  Double_t fdMinJetArea;             ///< min area to keep jet in output
  Double_t fdMinJetPt;               ///< min jet pt to keep jet in output
  Double_t fdJetPhiMin;              ///< minimum phi to keep jet in output
  Double_t fdJetPhiMax;              ///< maximum phi to keep jet in output
  Double_t fdJetEtaMin;              ///< minimum eta to keep jet in output
  Double_t fdJetEtaMax;              ///< maximum eta to keep jet in output
  Double_t fdJetTrackPtMin;          ///< minimum jet constituent track pt
  Double_t fdJetTrackEtaMax;          ///< max jet constituent track eta
  Double_t fdMaxEtaJetBG;            ///< Maximum jet eta cut for the bacground estimation
  Double_t fdBgRadius;               ///< Background jet radius

  Double_t fdCutPtJetMin; // [GeV/c] minimum jet pt
  Double_t fdCutPtTrackJetMin; // [GeV/c] minimum pt of leading jet-track
  Double_t fdCutAreaPercJetMin; // [pi*R^2] minimum jet area with respect to the expected value
  Double_t fdDistanceV0JetMax; // (R) D - maximum distance between V0 and jet axis used for finding V0s in the perp, rnd and median jet cone
  
  //MC var
  Double_t fdDistPrimaryMax;          ///< [cm] max distance of production point to the primary vertex (criterion for choice of MC particles considered as primary) 

  // Event histograms
  TH1D* fh1EventCounterCut; //! number of events for different selection steps
  TH1D* fh1EventCounterCutCent[fgkiNBinsCent]; //! number of events for different selection steps and different centralities
  TH1D* fh1EventCent; //! number of events for different centralities
  TH1D* fh1EventCent2; //! number of events for different centralities
  TH1D* fh1EventCent2Jets; //! number of events for different centralities
  TH1D* fh1EventCent2NoJets; //! number of events for different centralities
  TH2D* fh2EventCentTracks; //! number of tracks vs centrality
  TH2D* fh2EventCentMult; //! reference multiplicity vs centrality
  TH1D* fh1VtxZ[fgkiNBinsCent]; //! z coordinate of the primary vertex
  TH2D* fh2VtxXY[fgkiNBinsCent]; //! xy coordinates of the primary vertex
  TH1D* fh1V0CandPerEvent; //! number of V0 cand per event

  //jet histograms 
  TH1D* fh1PtJet[fgkiNBinsCent]; //! pt spectra of jets for normalisation of in-jet V0 spectra
  TH1D* fh1EtaJet[fgkiNBinsCent]; //! jet eta
  TH2D* fh2EtaPtJet[fgkiNBinsCent]; //! jet eta-pT
  TH1D* fh1PhiJet[fgkiNBinsCent]; //! jet phi
  TH2D* fh2PtJetPtTrackLeading[fgkiNBinsCent]; //! pt_jet; pt of leading jet track
  TH1D* fh1NJetPerEvent[fgkiNBinsCent]; //! number of jets per event
  TH1D* fh1NV0JetPerEvent[fgkiNBinsCent]; //! number of V0 jets per event
  TH1D* fh1NV0sInJetStats; //! V0s in jets statistics

  TH1D* fh1NRndConeCent; //! number of generated random cones in centrality bins
  TH2D* fh2EtaPhiRndCone[fgkiNBinsCent]; //! random cone eta-pT
  //TH1D* fh1NMedConeCent; //! number of found median-cluster cones in centrality bins
  //TH2D* fh2EtaPhiMedCone[fgkiNBinsCent]; //! median-cluster cone eta-phi
  
  // K0s
  TH1D* fh1V0CounterCentK0s[fgkiNBinsCent]; //! number of K0s candidates after various cuts
  TH1D* fh1V0InvMassK0sAll[fgkiNCategV0]; //! V0 invariant mass for each selection steps
  TH1D* fh1V0CandPerEventCentK0s[fgkiNBinsCent]; //! number of K0s candidates per event, in centrality bins
  TH1D* fh1V0InvMassK0sCent[fgkiNBinsCent]; //! V0 invariant mass, in centrality bins
  THnSparse* fhnV0InclusiveK0s[fgkiNBinsCent]; //! V0 inclusive, in a centrality bin, m_V0; pt_V0; eta_V0
  THnSparse* fhnV0InvMassCutK0s[fgkiNBinsCent]; //! V0 after invariant mass signal window cut, in a centrality bin, m_V0; pt_V0; eta_V0
  THnSparse* fhnV0InJetK0s[fgkiNBinsCent]; //! V0 in jet cones, in a centrality bin, m_V0; pt_V0; eta_V0; pt_jet
  //TH2D* fh2ArmPodK0s[fgkiNQAIndeces]; //! Armenteros-Podolanski
  
  THnSparse* fhnV0InPerpK0s[fgkiNBinsCent]; //! V0 in perpendicular cones, in a centrality bin, m_V0; pt_V0; eta_V0; pt_jet
  THnSparse* fhnV0InRndK0s[fgkiNBinsCent]; //! V0 in random cones, in a centrality bin, m_V0; pt_V0; eta_V0
  THnSparse* fhnV0InMedK0s[fgkiNBinsCent]; //! V0 in medium cones, in a centrality bin, m_V0; pt_V0; eta_V0
  THnSparse* fhnV0OutJetK0s[fgkiNBinsCent]; //! V0 outside jet cones, in a centrality bin, m_V0; pt_V0; eta_V0
  THnSparse* fhnV0NoJetK0s[fgkiNBinsCent]; //! V0 in no-jet events, in a centrality bin, m_V0; pt_V0; eta_V0

  // Lambda
  TH1D* fh1V0CounterCentLambda[fgkiNBinsCent]; //! number of Lambda candidates after various cuts
  TH1D* fh1V0InvMassLambdaAll[fgkiNCategV0]; //!
  TH1D* fh1V0CandPerEventCentLambda[fgkiNBinsCent]; //!
  TH1D* fh1V0InvMassLambdaCent[fgkiNBinsCent]; //!
  THnSparse* fhnV0InclusiveLambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0InvMassCutLambda[fgkiNBinsCent]; //! V0 after invariant mass signal window cut, in a centrality bin, m_V0; pt_V0; eta_V0
  THnSparse* fhnV0InJetLambda[fgkiNBinsCent]; //!
  //TH2D* fh2ArmPodLambda[fgkiNQAIndeces]; //!

  THnSparse* fhnV0InPerpLambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0InRndLambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0InMedLambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0OutJetLambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0NoJetLambda[fgkiNBinsCent]; //!

  // ALambda
  TH1D* fh1V0CounterCentALambda[fgkiNBinsCent]; //! number of ALambda candidates after various cuts
  TH1D* fh1V0InvMassALambdaAll[fgkiNCategV0]; //!
  TH1D* fh1V0CandPerEventCentALambda[fgkiNBinsCent]; //!
  TH1D* fh1V0InvMassALambdaCent[fgkiNBinsCent]; //!
  THnSparse* fhnV0InclusiveALambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0InvMassCutALambda[fgkiNBinsCent]; //! V0 after invariant mass signal window cut, in a centrality bin, m_V0; pt_V0; eta_V0
  THnSparse* fhnV0InJetALambda[fgkiNBinsCent]; //!
  //TH2D* fh2ArmPodALambda[fgkiNQAIndeces]; //!

  THnSparse* fhnV0InPerpALambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0InRndALambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0InMedALambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0OutJetALambda[fgkiNBinsCent]; //!
  THnSparse* fhnV0NoJetALambda[fgkiNBinsCent]; //!

  // MC histograms
  TH1D* fh1MCStats; //! V0s in jets statistics
  // K0 inclusive
  TH1D* fh1V0K0sPtMCGen[fgkiNBinsCent]; //! pt spectrum of all generated K0s in event
  TH2D* fh2V0K0sPtMassMCRec[fgkiNBinsCent]; //! pt-mass spectrum of successfully reconstructed K0s in event
  TH1D* fh1V0K0sPtMCRecFalse[fgkiNBinsCent]; //! pt spectrum of false reconstructed K0s in event
  TH2D* fh2V0K0sEtaPtMCGen[fgkiNBinsCent]; //! eta-pt spectrum of all generated K0s in event
  THnSparse* fh3V0K0sEtaPtMassMCRec[fgkiNBinsCent]; //! eta-pt-mass spectrum of successfully reconstructed K0s in event
  THnSparse* fhnV0K0sInclDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! V0 inclusive, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_V0; pt_V0; pt_jet
  //K0 in jets
  TH2D* fh2V0K0sInJetPtMCGen[fgkiNBinsCent]; //! pt spectrum of generated K0s in jet
  THnSparse* fh3V0K0sInJetPtMassMCRec[fgkiNBinsCent]; //! mass-pt spectrum of successfully reconstructed K0s in jet
  THnSparse* fh3V0K0sInJetEtaPtMCGen[fgkiNBinsCent]; //! eta-pt spectrum of generated K0s in jet
  THnSparse* fh4V0K0sInJetEtaPtMassMCRec[fgkiNBinsCent]; //! mass-eta-pt spectrum of successfully reconstructed K0s in jet
  THnSparse* fhnV0K0sInJetsDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! V0 in jets, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_V0; pt_V0; pt_jet
  // resolution
  TH2D* fh2V0K0sMCResolMPt[fgkiNBinsCent]; //! K0s mass resolution vs pt
  TH2D* fh2V0K0sMCPtGenPtRec[fgkiNBinsCent]; //! K0s generated pt vs reconstructed pt

  // Lambda inclusive
  TH1D* fh1V0LambdaPtMCGen[fgkiNBinsCent]; //!
  TH2D* fh2V0LambdaPtMassMCRec[fgkiNBinsCent]; //!
  TH1D* fh1V0LambdaPtMCRecFalse[fgkiNBinsCent]; //!
  TH2D* fh2V0LambdaEtaPtMCGen[fgkiNBinsCent]; //!
  THnSparse* fh3V0LambdaEtaPtMassMCRec[fgkiNBinsCent]; //!
  THnSparse* fhnV0LambdaInclDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! V0 inclusive, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_V0; pt_V0; pt_jet
  // Lambda in jets
  TH2D* fh2V0LambdaInJetPtMCGen[fgkiNBinsCent]; //!
  THnSparse* fh3V0LambdaInJetPtMassMCRec[fgkiNBinsCent]; //!
  THnSparse* fh3V0LambdaInJetEtaPtMCGen[fgkiNBinsCent]; //!
  THnSparse* fh4V0LambdaInJetEtaPtMassMCRec[fgkiNBinsCent]; //!
  THnSparse* fhnV0LambdaInJetsDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! V0 in jets, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_V0; pt_V0; pt_jet
  // resolution
  TH2D* fh2V0LambdaMCResolMPt[fgkiNBinsCent]; //!
  TH2D* fh2V0LambdaMCPtGenPtRec[fgkiNBinsCent]; //!

  TH1D* fh1V0XiPtMCGen[fgkiNBinsCent]; //!
  TH1D* fh1V0AXiPtMCGen[fgkiNBinsCent]; //!

  // ALambda inclusive
  TH1D* fh1V0ALambdaPtMCGen[fgkiNBinsCent]; //!
  TH2D* fh2V0ALambdaPtMassMCRec[fgkiNBinsCent]; //!
  TH1D* fh1V0ALambdaPtMCRecFalse[fgkiNBinsCent]; //!
  TH2D* fh2V0ALambdaEtaPtMCGen[fgkiNBinsCent]; //!
  THnSparse* fh3V0ALambdaEtaPtMassMCRec[fgkiNBinsCent]; //!
  THnSparse* fhnV0ALambdaInclDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! V0 inclusive, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_V0; pt_V0; pt_jet
  // ALambda in jets
  TH2D* fh2V0ALambdaInJetPtMCGen[fgkiNBinsCent]; //!
  THnSparse* fh3V0ALambdaInJetPtMassMCRec[fgkiNBinsCent]; //!
  THnSparse* fh3V0ALambdaInJetEtaPtMCGen[fgkiNBinsCent]; //!
  THnSparse* fh4V0ALambdaInJetEtaPtMassMCRec[fgkiNBinsCent]; //!
  THnSparse* fhnV0ALambdaInJetsDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! V0 in jets, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_V0; pt_V0; pt_jet

  // resolution
  TH2D* fh2V0ALambdaMCResolMPt[fgkiNBinsCent]; //!
  TH2D* fh2V0ALambdaMCPtGenPtRec[fgkiNBinsCent]; //!

  //TH1D* fh1V0AXiPtMCGen[fgkiNBinsCent]; //!

  AliAnalysisTaskStrangenessInJets(const AliAnalysisTaskStrangenessInJets&); // not implemented
  AliAnalysisTaskStrangenessInJets& operator=(const AliAnalysisTaskStrangenessInJets&); // not implemented

  ClassDef(AliAnalysisTaskStrangenessInJets, 2) // task for analysis of V0s (K0S, (anti-)Lambda) in charged jets
};

#endif
