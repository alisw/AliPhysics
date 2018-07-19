#ifndef AliAnalysisTaskEmcalJetHCorrelations_H
#define AliAnalysisTaskEmcalJetHCorrelations_H

/**
 * @class AliAnalysisTaskEmcalJetHCorrelations
 * @brief Jet-hadron correlations analysis task for central Pb-Pb and pp
 *
 * %Analysis task for jet-hadron correlations in central Pb-Pb and pp.
 * Includes the ability to weight entries by the efficiency, as well
 * as utilize a jet energy scale correction in the form of a histogram.
 *
 * This code has been cross checked against the code used for the event
 * plane dependent analysis (AliAnalysisTaskEmcalJetHadEPpid).
 *
 * @author Raymond Ehlers <raymond.ehlers@cern.ch>, Yale University
 * @author Megan Connors, Georgia State University
 * @date 1 Jan 2017
 */

#include <vector>

class TH1;
class TH2;
class TH3;
class THnSparse;
class AliEmcalJet;
class AliEventPoolManager;
class AliTLorentzVector;

class AliEmcalJet;
class AliJetContainer;
class AliParticleContainer;

#include "THistManager.h"
#include "AliAnalysisTaskEmcalJet.h"

namespace PWGJE {
namespace EMCALJetTasks {

class AliAnalysisTaskEmcalJetHCorrelations : public AliAnalysisTaskEmcalJet {
 public:
  /**
   * @enum jetBias_t
   * @brief Default value used to disable constituent bias
   */
  enum jetBias_t {
    kDisableBias = 10000            //!<! Arbitrarily large value which can be used to disable the constituent bias. Can be used for either tracks or clusters.
  };

  /**
   * @enum ESingleTrackEfficiency_t
   * @brief Define the single track efficiency to apply
   */
  enum ESingleTrackEfficiency_t {
    kEffDisable = 0,                //!<! Disable single track efficiency
    kEffAutomaticConfiguration = 1, //!<! Auto configure the single track efficiency based on the beam type, centrality, and run quality (number)
    kEffPP = 10,                    //!<! Explicitly select pp single track efficiency
  };

  AliAnalysisTaskEmcalJetHCorrelations();
  AliAnalysisTaskEmcalJetHCorrelations(const char *name);
  virtual ~AliAnalysisTaskEmcalJetHCorrelations() {}

  Double_t                GetTrackBias() const                       { return fTrackBias; }
  Double_t                GetClusterBias() const                     { return fClusterBias; }

  // Jet bias - setters
  /// Require a track with pt > b in jet
  virtual void            SetTrackBias(Double_t b)                   { fTrackBias    = b; }
  /// Require a cluster with pt > b in jet
  virtual void            SetClusterBias(Double_t b)                 { fClusterBias = b; }

  // Event trigger/mixed selection - setters
  /// Set the trigger event trigger selection
  virtual void            SetTriggerType(UInt_t te)                  { fTriggerType = te; }
  /// Set the mixed event trigger selection
  virtual void            SetMixedEventTriggerType(UInt_t me)        { fMixingEventType = me; }
  /// True if the task should be disabled for the fast partition
  void                    SetDisableFastPartition(Bool_t b = kTRUE)  { fDisableFastPartition = b; }
  Bool_t                  GetDisableFastPartition() const            { return fDisableFastPartition; }
  /// Require jet to be matched when embedding
  Bool_t                  GetRequireMatchedJetWhenEmbedding() const   { return fRequireMatchedJetWhenEmbedding; }
  void                    SetRequireMatchedJetWhenEmbedding(Bool_t b) { fRequireMatchedJetWhenEmbedding = b; }
  /// Mimimum shared momentum fraction for matched jet
  double                  GetMinimumSharedMomentumFraction() const    { return fMinSharedMomentumFraction; }
  void                    SetMinimumSharedMomentumFraction(double d)  { fMinSharedMomentumFraction = d; }
  double                  GetMaximumMatchedJetDistance() const        { return fMaxMatchedJetDistance; }
  void                    SetMaximumMatchedJetDistance(double d)      { fMaxMatchedJetDistance = d; }
  bool                    GetRequireMatchedPartLevelJet() const       { return fRequireMatchedPartLevelJet; }
  void                    SetRequireMatchedPartLevelJet(bool b)       { fRequireMatchedPartLevelJet = b; }

  // Mixed events
  virtual void            SetEventMixing(Bool_t enable)              { fDoEventMixing = enable;}
  virtual void            SetNumberOfMixingTracks(Int_t tracks)      { fNMixingTracks = tracks; }
  virtual void            SetMinNTracksForMixedEvents(Int_t nmt)     { fMinNTracksMixedEvents = nmt; }
  virtual void            SetMinNEventsForMixedEvents(Int_t nme)     { fMinNEventsMixedEvents = nme; }
  virtual void            SetNCentBinsMixedEvent(Bool_t centbins)    { fNCentBinsMixedEvent = centbins; }
  // Switch to cut out some unneeded sparse axis
  void                    SetDoLessSparseAxes(Bool_t dlsa)           { fDoLessSparseAxes = dlsa; }
  void                    SetDoWiderTrackBin(Bool_t wtrbin)          { fDoWiderTrackBin = wtrbin; }
  // Set efficiency correction
  void                    SetSingleTrackEfficiencyType(ESingleTrackEfficiency_t trackEffType) { fSingleTrackEfficiencyCorrectionType = trackEffType; }
  /// Artificial tracking inefficiency from 0 to 1. 1.0 (default) will disable it.
  void                    SetArtificialTrackingInefficiency(double eff) { fArtificialTrackInefficiency = eff; }
  // Setup JES correction
  void                    SetJESCorrectionHist(TH2D * hist)          { fJESCorrectionHist = hist; }
  void                    SetNoMixedEventJESCorrection(Bool_t b) { fNoMixedEventJESCorrection = b; }

  Bool_t                  RetrieveAndInitializeJESCorrectionHist(TString filename, TString histName, Double_t trackBias = AliAnalysisTaskEmcalJetHCorrelations::kDisableBias, Double_t clusterBias = AliAnalysisTaskEmcalJetHCorrelations::kDisableBias);

  // Corrections
  // Public so that it can be tested externally
  Double_t       EffCorrection(Double_t trkETA, Double_t trkPT, AliAnalysisTaskEmcal::BeamType beamType) const;

  virtual void            UserCreateOutputObjects();

  // AddTask
  static AliAnalysisTaskEmcalJetHCorrelations * AddTaskEmcalJetHCorrelations(
     const char *nTracks              = "usedefault",
     const char *nCaloClusters        = "usedefault",
     // Jet options
     const Double_t trackBias         = 5,
     const Double_t clusterBias       = 5,
     // Mixed event options
     const Int_t nTracksMixedEvent    = 0,  // Additionally acts as a switch for enabling mixed events
     const Int_t minNTracksMixedEvent = 5000,
     const Int_t minNEventsMixedEvent = 5,
     const UInt_t nCentBinsMixedEvent = 10,
     // Triggers
     UInt_t trigEvent                 = AliVEvent::kAny,
     UInt_t mixEvent                  = AliVEvent::kAny,
     // Options
     const Bool_t lessSparseAxes      = kFALSE,
     const Bool_t widerTrackBin       = kFALSE,
     // Corrections
     const AliAnalysisTaskEmcalJetHCorrelations::ESingleTrackEfficiency_t singleTrackEfficiency = AliAnalysisTaskEmcalJetHCorrelations::kEffDisable,
     const Bool_t JESCorrection = kFALSE,
     const char * JESCorrectionFilename = "alien:///alice/cern.ch/user/r/rehlersi/JESCorrection.root",
     const char * JESCorrectionHistName = "JESCorrection",
     const char *suffix               = "biased"
   );

  bool ConfigureForStandardAnalysis(std::string trackName = "usedefault",
    std::string clusName = "usedefault",
    const double jetConstituentPtCut = 3,
    const double trackEta = 0.8,
    const double jetRadius = 0.2);

  bool ConfigureForEmbeddingAnalysis(std::string trackName = "usedefault",
    std::string clusName = "caloClustersCombined",
    const double jetConstituentPtCut = 3,
    const double trackEta = 0.8,
    const double jetRadius = 0.2,
    const std::string & jetTag = "hybridLevelJets",
    const std::string & correlationsTracksCutsPeriod = "lhc11a");

 protected:

  // NOTE: This is not an ideal way to resolve the size of histogram initialization.
  //       Will be resolved when we move fully to the THnSparse
  /**
   * @enum binArrayLimits_t
   * @brief Define the number of elements in various arrays
   */
  enum binArrayLimits_t {
    kMaxTrackPtBins = 7,             //!<! Number of elements in track pt binned arrays
    kMaxCentralityBins = 5,          //!<! Number of elements in centrality binned arrays
    kMixedEventMultiplicityBins = 8, //!<! Number of elements in mixed event multiplicity binned arrays
  };

  // EMCal framework functions
  void                   ExecOnce();
  Bool_t                 Run();

  // Utility functions
  AliParticleContainer * CreateParticleOrTrackContainer(const std::string & collectionName) const;

  // Determine if a jet has been matched
  bool CheckForMatchedJet(AliJetContainer * jets, AliEmcalJet * jet, const std::string & histName);
  // Apply artificial tracking inefficiency
  bool CheckArtificialTrackEfficiency(unsigned int trackIndex, std::vector<unsigned int> & rejectedTrackIndices, bool useRejectedList);
  // Reduce event mixing memory usage
  TObjArray*             CloneAndReduceTrackList(std::vector<unsigned int> & rejectedTrackIndices, const bool useRejectedList);
  // Histogram helper functions
  virtual THnSparse*     NewTHnSparseF(const char* name, UInt_t entries);
  virtual void           GetDimParams(Int_t iEntry,TString &label, Int_t &nbins, Double_t &xmin, Double_t &xmax);
  // Binning helper functions
  Int_t                  GetTrackPtBin(Double_t pt) const;
  UInt_t                 RetrieveTriggerMask() const;
  // Helper functions
  void                   InitializeArraysToZero();
  void                   GetDeltaEtaDeltaPhiDeltaR(AliTLorentzVector & particleOne, AliVParticle * particleTwo, Double_t & deltaEta, Double_t & deltaPhi, Double_t & deltaR);
  Double_t               GetRelativeEPAngle(Double_t jetAngle, Double_t epAngle) const;
  // Test for biased jet
  Bool_t                 BiasedJet(AliEmcalJet * jet);
  // Corrections
  Double_t               EffCorrection(Double_t trkETA, Double_t trkPT) const;

  // Fill methods which allow for the JES correction
  void                   FillHist(TH1 * hist, Double_t fillValue, Double_t weight = 1.0, Bool_t noCorrection = kFALSE);
  void                   FillHist(THnSparse * hist, Double_t *fillValue, Double_t weight = 1.0, Bool_t noCorrection = kFALSE);
  void                   AccessSetOfYBinValues(TH2D * hist, Int_t xBin, std::vector <Double_t> & yBinsContent, Double_t scaleFactor = -1.0);
  
  // Jet bias
  Double_t               fTrackBias;               ///< Jet track bias
  Double_t               fClusterBias;             ///< Jet cluster bias
  // Event Mixing
  Bool_t                 fDoEventMixing;           ///< flag to do event mixing
  Int_t                  fNMixingTracks;           ///< size of track buffer for event mixing
  Int_t                  fMinNTracksMixedEvents;   ///< threshold to use event pool # tracks
  Int_t                  fMinNEventsMixedEvents;   ///< threshold to use event pool # events
  UInt_t                 fNCentBinsMixedEvent;     ///< N cent bins for the event mixing pool
  AliEventPoolManager   *fPoolMgr;                 //!<! Event pool manager
  // Event selection types
  UInt_t                 fTriggerType;             ///< Event selection for jets (ie triggered events).
  UInt_t                 fMixingEventType;         ///< Event selection for mixed events
  Bool_t                 fDisableFastPartition;    ///< True if task should be disabled for the fast partition, where the EMCal is not included.
  // Efficiency correction
  ESingleTrackEfficiency_t fSingleTrackEfficiencyCorrectionType; ///< Control the efficiency correction. See EffCorrection() for meaning of values.
  Double_t               fArtificialTrackInefficiency; ///< Artificial track inefficiency. Enabled if < 1.0
  // JES correction
  Bool_t                 fNoMixedEventJESCorrection; ///< True if the jet energy scale correction should be applied to mixed event histograms
  TH2D                  *fJESCorrectionHist;       ///< Histogram containing the jet energy scale correction
  // Histogram binning variables
  Bool_t                 fDoLessSparseAxes;        ///< True if there should be fewer THnSparse axes
  Bool_t                 fDoWiderTrackBin;         ///< True if the track pt bins in the THnSparse should be wider
  Bool_t                 fRequireMatchedJetWhenEmbedding; ///< True if jets are required to be matched (ie. jet->MatchedJet() != nullptr)
  Double_t               fMinSharedMomentumFraction; ///< Minimum shared momentum with matched jet
  bool                   fRequireMatchedPartLevelJet; ///< True if matched jets are required to be matched to a particle level jet
  Double_t               fMaxMatchedJetDistance;  ///< Maximum distance between two matched jets

  // Histograms
  THistManager           fHistManager;             ///<  Histogram manager
  TH1                   *fHistTrackPt;             //!<! Track pt spectrum
  TH2                   *fHistJetEtaPhi;           //!<! Jet eta-phi distribution
  TH2                   *fHistTrackEtaPhi[7];      //!<! Track eta-phi distribution (the array corresponds to track pt)
  TH2                   *fHistJetHEtaPhi;          //!<! Eta-phi distribution of jets which are in jet-hadron correlations

  TH1                   *fHistJetPt[6];            //!<! Jet pt spectrum (the array corresponds to centrality bins)
  TH1                   *fHistJetPtBias[6];        //!<! Jet pt spectrum of jets which meet the constituent bias criteria (the array corresponds to centrality bins)
  THnSparse             *fhnMixedEvents;           //!<! Mixed events THnSparse
  THnSparse             *fhnJH;                    //!<! JetH THnSparse
  THnSparse             *fhnTrigger;               //!<! JetH trigger sparse

  // Pb-Pb Efficiency correction coefficients
  static Double_t p0_10SG[17];                    ///< 0-10% centrality semi-good runs
  static Double_t p10_30SG[17];                   ///< 10-30% centrality semi-good runs
  static Double_t p30_50SG[17];                   ///< 30-50% centrality semi-good runs
  static Double_t p50_90SG[17];                   ///< 50-90% centrality semi-good runs
  // Good Runs
  static Double_t p0_10G[17];                     ///< 0-10% centrality good runs
  static Double_t p10_30G[17];                    ///< 10-30% centrality good runs
  static Double_t p30_50G[17];                    ///< 30-50% centrality good runs
  static Double_t p50_90G[17];                    ///< 50-90% centrality good runs

 private:

  AliAnalysisTaskEmcalJetHCorrelations(const AliAnalysisTaskEmcalJetHCorrelations&); // not implemented
  AliAnalysisTaskEmcalJetHCorrelations& operator=(const AliAnalysisTaskEmcalJetHCorrelations&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetHCorrelations, 16);
};

} /* namespace EMCALJetTasks */
} /* namespace PWGJE */

#endif
