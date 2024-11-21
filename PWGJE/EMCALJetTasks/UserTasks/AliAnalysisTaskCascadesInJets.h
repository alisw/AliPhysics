#ifndef AliAnalysisTaskCascadesInJets_cxx
#define AliAnalysisTaskCascadesInJets_cxx

//-------------------------------------------------------------------------
// Task for Cascade analysis in charged jets with 
// the strange particles (instead of daughters) added to the  jet finder
// Author: Ekaterina Grecka (ermeeka@fjfi.cvut.cz)
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
#include "AliAODRecoDecay.h"
#include "AliAODJet.h"
#include "AliAnalysisTaskEmcal.h"

namespace fastjet {
  class PseudoJet;
}

class AliAnalysisTaskCascadesInJets : public AliAnalysisTaskEmcal
{
public:
  
#if !defined(__CINT__) && !defined(__MAKECINT__) 
  typedef fastjet::JetAlgorithm FJJetAlgo;
  typedef fastjet::RecombinationScheme FJRecoScheme;
#endif

  AliAnalysisTaskCascadesInJets(); // Default constructor
  AliAnalysisTaskCascadesInJets(const char* name); // Constructor
  virtual ~AliAnalysisTaskCascadesInJets(); // Destructor
  void UserCreateOutputObjects();
  void Terminate(Option_t*) {}

  void SetIsPbPb(Bool_t val = 1) {fbIsPbPb = val;}
  void SetMCAnalysis(Bool_t select = kTRUE) {fbMCAnalysis = select;}
  void SetGeneratorName(TString name) {fsGeneratorName = name;}
  Bool_t IsFromGoodGenerator(Int_t index); // True if the MC particle with the given index comes from the selected generator

  void SetSignalInBG(Bool_t val = 0) {fbSignalInBG = val;}
  void SetNSigmas(Double_t val = 9) {fdNSigmas = val;}

  // Event selection setters 
  void SetEventCuts(Double_t z = 10, Double_t r = 1, Double_t cL = 0, Double_t cH = 80, Double_t dZ = 0.1, Int_t iNC = 1) {fdCutVertexZ = z; fdCutVertexR2 = r * r; fdCutCentLow = cL; fdCutCentHigh = cH; fdCutDeltaZMax = dZ; fiNContribMin = iNC;} 
  void SetUseMultiplicity(Bool_t val = kTRUE) {fbUseMultiplicity = val;} 
  void SetUseIonutCut(Bool_t val = kTRUE) {fbUseIonutCut = val;}
  // Cascades selection setters 

  void SetCutTPCRefit(Bool_t val = kTRUE) {fbTPCRefit = val;}
  void SetCutRejectKinks(Bool_t val = kTRUE) {fbRejectKinks = val;}
  void SetCutFindableClusters(Bool_t val = kTRUE) {fbFindableClusters = val;}
  void SetCutCascadePtMin(Double_t val = 1.) {fdCutCascadePtMin = val;}
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
  void SetCutRadiusDecayMin(Double_t val = 5.) {fdCutRadiusDecayMin = val;}
  void SetCutRadiusDecayMax(Double_t val = 100.) {fdCutRadiusDecayMax = val;}
  void SetCutEtaCascadeMax(Double_t val = 0.7) {fdCutEtaCascadeMax = val;}
  void SetCutRapCascadeMax(Double_t val = 0.75) {fdCutRapCascadeMax = val;}

  void SetCutDCACascadeBachToPrimVtxMin(Double_t val = 0.04) {fdCutDCACascadeBachToPrimVtxMin = val;}
  void SetCutDCACascadeV0ToPrimVtxMin(Double_t val = 0.1) {fdCutDCACascadeV0ToPrimVtxMin = val;}
  void SetCutDCACascadeDaughtersToPrimVtxMin(Double_t val = 0.03) {fdCutDCACascadeDaughtersToPrimVtxMin = val;}
  void SetCutDCACascadeV0DaughtersMax(Double_t val = 1.) {fdCutDCACascadeV0DaughtersMax = val;}
  void SetCutDCACascadeBachToV0Max(Double_t val = 1.3) {fdCutDCACascadeBachToV0Max = val;}
  void SetCutCPACascadeMin(Double_t val = 0.97) {fdCutCPACascadeMin = val;} 
  void SetCutCPACascadeV0Min(Double_t val = 0.998) {fdCutCPACascadeV0Min = val;} 

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
  void SetDistanceCascadeJetMax(Double_t val = 0.4) {fdDistanceCascadeJetMax = val;}
  void SetPtJetMin(Double_t ptMin = 0) {fdCutPtJetMin = ptMin;}
  void SetPtTrackJetMin(Double_t ptMin = 0) {fdCutPtTrackJetMin = ptMin;}
  void SetAreaPercJetMin(Double_t area = 0) {fdCutAreaPercJetMin = area;}
  void SetLeadingV0(Bool_t b = 0) {bdLeadingV0 = b;}
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
  Double_t GetCutRadiusDecayMin() const   { return fdCutRadiusDecayMin; }
  Double_t GetCutRadiusDecayMax() const   { return fdCutRadiusDecayMax; }
  Double_t GetCutEtaCascadeMax() const   { return fdCutEtaCascadeMax; }
  Double_t GetCutRapCascadeMax() const   { return fdCutRapCascadeMax; }
  Double_t GetGhostArea() { return fdGhostArea; }
  Double_t GetRadius() { return fdRadius; }

  // upper edges of centrality bins
  static const Int_t fgkiNBinsCent = 1; // number of centrality bins
  static const Int_t fgkiCentBinRanges[fgkiNBinsCent]; // upper edges of centrality bins
  // axis: pT of Cascade
  static const Double_t fgkdBinsPtCascade[2]; // [GeV/c] minimum and maximum or desired binning of the axis (intended for the rebinned axis)
  static const Int_t fgkiNBinsPtCascade; // number of bins (intended for the rebinned axis)
  static const Int_t fgkiNBinsPtCascadeInit; // initial number of bins (uniform binning)
  static const Int_t fgkiNBinsPtCascadeInitInJet; // initial number of bins for Cascades in jets (uniform binning)
  // axis: Xi invariant mass
  static const Int_t fgkiNBinsMassXi; // number of bins (uniform binning)
  static const Double_t fgkdMassXiMin; // minimum
  static const Double_t fgkdMassXiMax; // maximum  -
  // axis: Omega invariant mass
  static const Int_t fgkiNBinsMassOmega; // number of bins (uniform binning)
  static const Double_t fgkdMassOmegaMin; // minimum
  static const Double_t fgkdMassOmegaMax; // maximum  
  // axis: pT of jets
  static const Double_t fgkdBinsPtJet[2]; // [GeV/c] minimum and maximum or desired binning of the axis (intended for the rebinned axis)
  static const Int_t fgkiNBinsPtJet; // number of bins (intended for the rebinned axis)
  static const Int_t fgkiNBinsPtJetInit; // initial number of bins (uniform binning)
  // PDG codes of used particles
  static const Int_t iPdgCodePion;
  static const Int_t iPdgCodeProton;
  static const Int_t iPdgCodeK0s;
  static const Int_t iPdgCodeLambda;
  static const Int_t iPdgCodeKaon;
  static const Int_t iPdgCodeXi;
  static const Int_t iPdgCodeOmega;

  //Index codes to distinguish particles in the jet 
  static const Int_t iXiMinusId;
  static const Int_t iXiPlusId;
  static const Int_t iOmegaMinusId;
  static const Int_t iOmegaPlusId;

  void FillCandidates(Double_t mXi, Double_t mOmega, Bool_t isXiMinus, Bool_t isXiPlus, Bool_t isOmegaMinus, Bool_t isOmegaPlus, Int_t iCut, Int_t iCent); //Fills histograms according to the Cascade type 
  Bool_t IsParticleInCone(const AliVParticle* part1, const AliVParticle* part2, Double_t dRMax) const; // decides whether a particle is inside a jet cone
  Bool_t OverlapWithJets(const TClonesArray* array, const AliVParticle* cone, Double_t dDistance) const; // decides whether a cone overlaps with other jets
  AliAODJet* GetRandomCone(const TClonesArray* array, Double_t dEtaConeMax, Double_t dDistance) const; // generate a random cone which does not overlap with selected jets
  Double_t AreaCircSegment(Double_t dRadius, Double_t dDistance) const; // area of circular segment
  Bool_t IsSelectedForAnalysis(); // Function for the event selection
  Int_t GetCentralityBinIndex(Double_t centrality);
  Int_t GetCentralityBinEdge(Int_t index);
  TString GetCentBinLabel(Int_t index);
  static Double_t AddDaughters(AliAODRecoDecay* cand, TObjArray& daughters);
  Bool_t AssociateRecCascadeWithMC( AliAODcascade* cascpart, AliEmcalJet *xjet, Bool_t bIsXiMinus, Bool_t bIsXiPlus, Bool_t bIsOmegaMinus, Bool_t bIsOmegaPlus, Int_t iCent);
  Bool_t GeneratedMCParticles( TClonesArray* track, Int_t iCent );
  Double_t MassPeakSigma(Double_t pt, Int_t particle);
  Double_t MassPeakMean(Double_t pt, Int_t particle);
  
protected:
  void ExecOnce();
  Bool_t FillHistograms();
  Bool_t Run();
  void AddEventTracks(TClonesArray* coll, std::vector<fastjet::PseudoJet>& VectorBgPart);  
  void AddEventTracksMC(TClonesArray* coll, std::vector<fastjet::PseudoJet>& VectorBgPartMC);  
  Bool_t GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const;

  TList* fOutputListStd; //! Output list for standard analysis results
  TList* fOutputListStdJets; //! Output list for jet analysis results
  TList* fOutputListMC;  //! Output list for MC analysis results
  TClonesArray   *fCascadeCandidateArray;          //! contains selected Cascade candidates
  TClonesArray   *fGenMCCascade;                    //! contains MC generated Cascades
 
  Int_t           fNCand;                   //! number of selected Cascade candidates already added to fCandidateArray

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
  
  static const Int_t fgkiNCategCascade = 21; // number of Cascade selection steps
 
 // Data selection
  Bool_t fbIsPbPb; // switch: Pb+Pb / p+p collisions
  Bool_t fbMCAnalysis; // switch: simulated / real data
  TString fsGeneratorName; // pattern for selecting only Cascades from a specific MC generator
  
  Bool_t fbSignalInBG; //switch: takes Cascades from BG region insted of signal for the jet analysis
  Double_t fdNSigmas; //multiple of sigmas for the bg estimation

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
  // Cascade selection
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

  Double_t fdCutDCACascadeBachToPrimVtxMin; // [cm] min DCA of bachelor track to the prim vtx
  Double_t fdCutDCACascadeV0ToPrimVtxMin;   // [cm] min DCA of V0 (from cascade decay) to the prim vtx
  Double_t fdCutDCACascadeDaughtersToPrimVtxMin; // [cm] min DCA of daughters (from secondary V0 decay) to the prim vtx
  Double_t fdCutDCACascadeV0DaughtersMax; // 1.5; // [sigma of TPC tracking] max DCA between Cascade V0 daughters
  Double_t fdCutDCACascadeBachToV0Max; // 1.3; // [cm] max DCA between bachelor track to V0   
  Double_t fdCutCPACascadeMin;   // min cosine of the pointing angle of the cascade
  Double_t fdCutCPACascadeV0Min; // min cosine of the pointing angle of the V0 in the cascade

  Double_t fdCutNTauXMax; // (5.0) [tau] max proper lifetime in multiples of the mean lifetime, Xi

  // Cascade candidate
  Bool_t fbOnFly; // (0) on-the-fly (yes) or offline (no) reconstructed
  Double_t fdCutCascadePtMin; // (1 GeV) min pT of the selected v0 particles
  Double_t fdCutRadiusDecayMin; // (5.) [cm] min radial distance of the decay vertex
  Double_t fdCutRadiusDecayMax; // (100.) [cm] max radial distance of the decay vertex
  Double_t fdCutEtaCascadeMax; // (0.7) max |pseudorapidity| of Cascade
  Double_t fdCutRapCascadeMax; // (0.75) max |rapidity| of Cascade (turned off)


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
  Double_t fdDistanceCascadeJetMax; // (R) D - maximum distance between Cascade and jet axis used for finding Cascades in the perp, rnd and median jet cone
  Bool_t bdLeadingV0; ///< 0 - leading track pt cut on all jets, 1 - if leading is V0 do not apply cut  
  //MC var
  Double_t fdDistPrimaryMax;          ///< [cm] max distance of production point to the primary vertex (criterion for choice of MC particles considered as primary) 
  
  AliTrackContainer* fTracksCont; //! Tracks

  // Event histograms
  TH1D* fh1EventCounterCut; //! number of events for different selection steps
  TH1D* fh1EventCounterCutCent[fgkiNBinsCent]; //! number of events for different selection steps and different centralities
  TH1D* fh1EventCent; //! number of events for different centralities
  TH1D* fh1EventCent2; //! number of events for different centralities
  TH2D* fh2EventCentTracks; //! number of tracks vs centrality
  TH2D* fh2EventCentMult; //! reference multiplicity vs centrality
  TH1D* fh1VtxZ[fgkiNBinsCent]; //! z coordinate of the primary vertex
  TH2D* fh2VtxXY[fgkiNBinsCent]; //! xy coordinates of the primary vertex
  TH1D* fh1CascadeCandPerEvent; //! number of Cascade cand per event

  //jet histograms 
  TH1D* fh1PtJet[fgkiNBinsCent]; //! pt spectra of jets for normalisation of in-jet Cascade spectra
  TH1D* fh1EtaJet[fgkiNBinsCent]; //! jet eta
  TH2D* fh2EtaPtJet[fgkiNBinsCent]; //! jet eta-pT
  TH1D* fh1PhiJet[fgkiNBinsCent]; //! jet phi
  TH2D* fh2PtJetPtTrackLeading[fgkiNBinsCent]; //! pt_jet; pt of leading jet track
  TH1D* fh1NJetPerEvent[fgkiNBinsCent]; //! number of jets per event
  TH1D* fh1NCascadeJetPerEvent[fgkiNBinsCent]; //! number of Cascade jets per event
  TH1D* fh1NCascadesInJetStats; //! Cascades in jets statistics
  TH1D* fh1AreaExcluded; //! area of excluded cones for outside-cones V0s


  TH1D* fh1NRndConeCent; //! number of generated random cones in centrality bins
  TH2D* fh2EtaPhiRndCone[fgkiNBinsCent]; //! random cone eta-pT
  //TH1D* fh1NMedConeCent; //! number of found median-cluster cones in centrality bins
  //TH2D* fh2EtaPhiMedCone[fgkiNBinsCent]; //! median-cluster cone eta-phi
  
  // XiMinus
  TH1D* fh1CascadeCounterCentXiMinus[fgkiNBinsCent]; //! number of XiMinus candidates after various cuts
  TH1D* fh1CascadeInvMassXiMinusAll[fgkiNCategCascade]; //! Cascade invariant mass for each selection steps
  TH1D* fh1CascadeCandPerEventCentXiMinus[fgkiNBinsCent]; //! number ofXiMinus candidates per event, in centrality bins
  TH1D* fh1CascadeInvMassXiMinusCent[fgkiNBinsCent]; //! Cascade invariant mass, in centrality bins
  THnSparse* fhnCascadeInclusiveXiMinus[fgkiNBinsCent]; //! Cascade inclusive, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeInvMassCutXiMinus[fgkiNBinsCent]; //! Cascade after invariant mass signal window cut, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeInJetXiMinus[fgkiNBinsCent]; //! Cascade in jet cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade; pt_jet

  THnSparse* fhnCascadeInPerpXiMinus[fgkiNBinsCent]; //! Cascade in perpendicular cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade; pt_jet
  THnSparse* fhnCascadeInRndXiMinus[fgkiNBinsCent]; //! Cascade in random cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeInMedXiMinus[fgkiNBinsCent]; //! Cascade in medium cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeOutJetXiMinus[fgkiNBinsCent]; //! Cascade outside jet cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeNoJetXiMinus[fgkiNBinsCent]; //! Cascade in no-jet events, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  // XiPlus
  TH1D* fh1CascadeCounterCentXiPlus[fgkiNBinsCent]; //! number of XiPlus candidates after various cuts
  TH1D* fh1CascadeInvMassXiPlusAll[fgkiNCategCascade]; //! Cascade invariant mass for each selection steps
  TH1D* fh1CascadeCandPerEventCentXiPlus[fgkiNBinsCent]; //! number ofXiMinus candidates per event, in centrality bins
  TH1D* fh1CascadeInvMassXiPlusCent[fgkiNBinsCent]; //! Cascade invariant mass, in centrality bins
  THnSparse* fhnCascadeInclusiveXiPlus[fgkiNBinsCent]; //! Cascade inclusive, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeInvMassCutXiPlus[fgkiNBinsCent]; //! Cascade after invariant mass signal window cut, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeInJetXiPlus[fgkiNBinsCent]; //! Cascade in jet cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade; pt_jet

  THnSparse* fhnCascadeInPerpXiPlus[fgkiNBinsCent]; //! Cascade in perpendicular cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade; pt_jet
  THnSparse* fhnCascadeInRndXiPlus[fgkiNBinsCent]; //! Cascade in random cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeInMedXiPlus[fgkiNBinsCent]; //! Cascade in medium cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeOutJetXiPlus[fgkiNBinsCent]; //! Cascade outside jet cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeNoJetXiPlus[fgkiNBinsCent]; //! Cascade in no-jet events, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  // OmegaMinus
  TH1D* fh1CascadeCounterCentOmegaMinus[fgkiNBinsCent]; //! number of OmegaMinus candidates after various cuts
  TH1D* fh1CascadeInvMassOmegaMinusAll[fgkiNCategCascade]; //! Cascade invariant mass for each selection steps
  TH1D* fh1CascadeCandPerEventCentOmegaMinus[fgkiNBinsCent]; //! number ofOmegaMinus candidates per event, in centrality bins
  TH1D* fh1CascadeInvMassOmegaMinusCent[fgkiNBinsCent]; //! Cascade invariant mass, in centrality bins
  THnSparse* fhnCascadeInclusiveOmegaMinus[fgkiNBinsCent]; //! Cascade inclusive, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeInvMassCutOmegaMinus[fgkiNBinsCent]; //! Cascade after invariant mass signal window cut, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeInJetOmegaMinus[fgkiNBinsCent]; //! Cascade in jet cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade; pt_jet

  THnSparse* fhnCascadeInPerpOmegaMinus[fgkiNBinsCent]; //! Cascade in perpendicular cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade; pt_jet
  THnSparse* fhnCascadeInRndOmegaMinus[fgkiNBinsCent]; //! Cascade in random cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeInMedOmegaMinus[fgkiNBinsCent]; //! Cascade in medium cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeOutJetOmegaMinus[fgkiNBinsCent]; //! Cascade outside jet cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeNoJetOmegaMinus[fgkiNBinsCent]; //! Cascade in no-jet events, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  // OmegaPlus
  TH1D* fh1CascadeCounterCentOmegaPlus[fgkiNBinsCent]; //! number of OmegaPlus candidates after various cuts
  TH1D* fh1CascadeInvMassOmegaPlusAll[fgkiNCategCascade]; //! Cascade invariant mass for each selection steps
  TH1D* fh1CascadeCandPerEventCentOmegaPlus[fgkiNBinsCent]; //! number ofOmegaMinus candidates per event, in centrality bins
  TH1D* fh1CascadeInvMassOmegaPlusCent[fgkiNBinsCent]; //! Cascade invariant mass, in centrality bins
  THnSparse* fhnCascadeInclusiveOmegaPlus[fgkiNBinsCent]; //! Cascade inclusive, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeInvMassCutOmegaPlus[fgkiNBinsCent]; //! Cascade after invariant mass signal window cut, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeInJetOmegaPlus[fgkiNBinsCent]; //! Cascade in jet cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade; pt_jet

  THnSparse* fhnCascadeInPerpOmegaPlus[fgkiNBinsCent]; //! Cascade in perpendicular cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade; pt_jet
  THnSparse* fhnCascadeInRndOmegaPlus[fgkiNBinsCent]; //! Cascade in random cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeInMedOmegaPlus[fgkiNBinsCent]; //! Cascade in medium cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeOutJetOmegaPlus[fgkiNBinsCent]; //! Cascade outside jet cones, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade
  THnSparse* fhnCascadeNoJetOmegaPlus[fgkiNBinsCent]; //! Cascade in no-jet events, in a centrality bin, m_Cascade; pt_Cascade; eta_Cascade

  // MC histograms
  TH1D* fh1MCStats; //! Cascades in jets statistics
  // XiMinus inclusive
  TH1D* fh1CascadeXiMinusPtMCGen[fgkiNBinsCent]; //! pt spectrum of all generated XiMinus in event
  TH2D* fh2CascadeXiMinusPtMassMCRec[fgkiNBinsCent]; //! pt-mass spectrum of successfully reconstructed XiMinus in event
  TH1D* fh1CascadeXiMinusPtMCRecFalse[fgkiNBinsCent]; //! pt spectrum of false reconstructed XiMinus in event
  TH2D* fh2CascadeXiMinusEtaPtMCGen[fgkiNBinsCent]; //! eta-pt spectrum of all generated XiMinus in event
  THnSparse* fh3CascadeXiMinusEtaPtMassMCRec[fgkiNBinsCent]; //! eta-pt-mass spectrum of successfully reconstructed XiMinus in event
  THnSparse* fhnCascadeXiMinusInclDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! Cascade inclusive, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_Cascade; pt_Cascade; pt_jet
  //XiMinus in jets
  TH2D* fh2CascadeXiMinusInJetPtMCGen[fgkiNBinsCent]; //! pt spectrum of generated XiMinus in jet
  THnSparse* fh3CascadeXiMinusInJetPtMassMCRec[fgkiNBinsCent]; //! mass-pt spectrum of successfully reconstructed XiMinus in jet
  THnSparse* fh3CascadeXiMinusInJetEtaPtMCGen[fgkiNBinsCent]; //! eta-pt spectrum of generated XiMinus in jet
  THnSparse* fh4CascadeXiMinusInJetEtaPtMassMCRec[fgkiNBinsCent]; //! mass-eta-pt spectrum of successfully reconstructed XiMinus in jet
  THnSparse* fhnCascadeXiMinusInJetsDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! Cascade in jets, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_Cascade; pt_Cascade; pt_jet
  // resolution
  TH2D* fh2CascadeXiMinusMCResolMPt[fgkiNBinsCent]; //! XiMinus mass resolution vs pt
  TH2D* fh2CascadeXiMinusMCPtGenPtRec[fgkiNBinsCent]; //! K0s generated pt vs reconstructed pt

  // XiPlus inclusive
  TH1D* fh1CascadeXiPlusPtMCGen[fgkiNBinsCent]; //! pt spectrum of all generated XiPlus in event
  TH2D* fh2CascadeXiPlusPtMassMCRec[fgkiNBinsCent]; //! pt-mass spectrum of successfully reconstructed XiPlus in event
  TH1D* fh1CascadeXiPlusPtMCRecFalse[fgkiNBinsCent]; //! pt spectrum of false reconstructed XiPlus in event
  TH2D* fh2CascadeXiPlusEtaPtMCGen[fgkiNBinsCent]; //! eta-pt spectrum of all generated XiPlus in event
  THnSparse* fh3CascadeXiPlusEtaPtMassMCRec[fgkiNBinsCent]; //! eta-pt-mass spectrum of successfully reconstructed XiPlus in event
  THnSparse* fhnCascadeXiPlusInclDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! Cascade inclusive, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_Cascade; pt_Cascade; pt_jet
  //XiPlus in jets
  TH2D* fh2CascadeXiPlusInJetPtMCGen[fgkiNBinsCent]; //! pt spectrum of generated XiPlus in jet
  THnSparse* fh3CascadeXiPlusInJetPtMassMCRec[fgkiNBinsCent]; //! mass-pt spectrum of successfully reconstructed XiPlus in jet
  THnSparse* fh3CascadeXiPlusInJetEtaPtMCGen[fgkiNBinsCent]; //! eta-pt spectrum of generated XiPlus in jet
  THnSparse* fh4CascadeXiPlusInJetEtaPtMassMCRec[fgkiNBinsCent]; //! mass-eta-pt spectrum of successfully reconstructed XiPlus in jet
  THnSparse* fhnCascadeXiPlusInJetsDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! Cascade in jets, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_Cascade; pt_Cascade; pt_jet
  // resolution
  TH2D* fh2CascadeXiPlusMCResolMPt[fgkiNBinsCent]; //! XiPlus mass resolution vs pt
  TH2D* fh2CascadeXiPlusMCPtGenPtRec[fgkiNBinsCent]; //! K0s generated pt vs reconstructed pt
  
  // OmegaMinus inclusive
  TH1D* fh1CascadeOmegaMinusPtMCGen[fgkiNBinsCent]; //! pt spectrum of all generated OmegaMinus in event
  TH2D* fh2CascadeOmegaMinusPtMassMCRec[fgkiNBinsCent]; //! pt-mass spectrum of successfully reconstructed OmegaMinus in event
  TH1D* fh1CascadeOmegaMinusPtMCRecFalse[fgkiNBinsCent]; //! pt spectrum of false reconstructed OmegaMinus in event
  TH2D* fh2CascadeOmegaMinusEtaPtMCGen[fgkiNBinsCent]; //! eta-pt spectrum of all generated OmegaMinus in event
  THnSparse* fh3CascadeOmegaMinusEtaPtMassMCRec[fgkiNBinsCent]; //! eta-pt-mass spectrum of successfully reconstructed OmegaMinus in event
  THnSparse* fhnCascadeOmegaMinusInclDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! Cascade inclusive, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_Cascade; pt_Cascade; pt_jet
  //OmegaMinus in jets
  TH2D* fh2CascadeOmegaMinusInJetPtMCGen[fgkiNBinsCent]; //! pt spectrum of generated OmegaMinus in jet
  THnSparse* fh3CascadeOmegaMinusInJetPtMassMCRec[fgkiNBinsCent]; //! mass-pt spectrum of successfully reconstructed OmegaMinus in jet
  THnSparse* fh3CascadeOmegaMinusInJetEtaPtMCGen[fgkiNBinsCent]; //! eta-pt spectrum of generated OmegaMinus in jet
  THnSparse* fh4CascadeOmegaMinusInJetEtaPtMassMCRec[fgkiNBinsCent]; //! mass-eta-pt spectrum of successfully reconstructed OmegaMinus in jet
  THnSparse* fhnCascadeOmegaMinusInJetsDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! Cascade in jets, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_Cascade; pt_Cascade; pt_jet
  // resolution
  TH2D* fh2CascadeOmegaMinusMCResolMPt[fgkiNBinsCent]; //! OmegaMinus mass resolution vs pt
  TH2D* fh2CascadeOmegaMinusMCPtGenPtRec[fgkiNBinsCent]; //! K0s generated pt vs reconstructed pt
  
  // OmegaPlus inclusive
  TH1D* fh1CascadeOmegaPlusPtMCGen[fgkiNBinsCent]; //! pt spectrum of all generated OmegaPlus in event
  TH2D* fh2CascadeOmegaPlusPtMassMCRec[fgkiNBinsCent]; //! pt-mass spectrum of successfully reconstructed OmegaPlus in event
  TH1D* fh1CascadeOmegaPlusPtMCRecFalse[fgkiNBinsCent]; //! pt spectrum of false reconstructed OmegaPlus in event
  TH2D* fh2CascadeOmegaPlusEtaPtMCGen[fgkiNBinsCent]; //! eta-pt spectrum of all generated OmegaPlus in event
  THnSparse* fh3CascadeOmegaPlusEtaPtMassMCRec[fgkiNBinsCent]; //! eta-pt-mass spectrum of successfully reconstructed OmegaPlus in event
  THnSparse* fhnCascadeOmegaPlusInclDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! Cascade inclusive, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_Cascade; pt_Cascade; pt_jet
  //OmegaPlus in jets
  TH2D* fh2CascadeOmegaPlusInJetPtMCGen[fgkiNBinsCent]; //! pt spectrum of generated OmegaPlus in jet
  THnSparse* fh3CascadeOmegaPlusInJetPtMassMCRec[fgkiNBinsCent]; //! mass-pt spectrum of successfully reconstructed OmegaPlus in jet
  THnSparse* fh3CascadeOmegaPlusInJetEtaPtMCGen[fgkiNBinsCent]; //! eta-pt spectrum of generated OmegaPlus in jet
  THnSparse* fh4CascadeOmegaPlusInJetEtaPtMassMCRec[fgkiNBinsCent]; //! mass-eta-pt spectrum of successfully reconstructed OmegaPlus in jet
  THnSparse* fhnCascadeOmegaPlusInJetsDaughterEtaPtPtMCRec[fgkiNBinsCent]; //! Cascade in jets, reconstructed: charge_daughter; eta_daughter; pt_daughter; eta_Cascade; pt_Cascade; pt_jet
  // resolution
  TH2D* fh2CascadeOmegaPlusMCResolMPt[fgkiNBinsCent]; //! OmegaPlus mass resolution vs pt
  TH2D* fh2CascadeOmegaPlusMCPtGenPtRec[fgkiNBinsCent]; //! K0s generated pt vs reconstructed pt

  AliAnalysisTaskCascadesInJets(const AliAnalysisTaskCascadesInJets&); // not implemented
  AliAnalysisTaskCascadesInJets& operator=(const AliAnalysisTaskCascadesInJets&); // not implemented

  ClassDef(AliAnalysisTaskCascadesInJets, 7) // task for analysis of Cascades (Xi+-, Omega+-) in charged jets
};

#endif
