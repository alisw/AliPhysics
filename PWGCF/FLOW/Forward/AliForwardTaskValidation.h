#ifndef AliForwardTaskValidation_cxx
#define AliForwardTaskValidation_cxx

#include <string>
#include <vector>

#include "TClonesArray.h"

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliForwardSettings.h"
#include "AliForwardFlowUtil.h"

class TH2;

class AliForwardTaskValidation : public AliAnalysisTaskSE {
 public:
  AliForwardTaskValidation();
  /// `is_reconstructed` is used to toggle some event selections
  AliForwardTaskValidation(const char *name);
  AliForwardTaskValidation(const AliForwardTaskValidation&);
  /// Set up this task. This function acts as the AddTask macro
  /// `is_reconstructed` is passed on to the constructor of this task
  static AliForwardTaskValidation* ConnectTask(TString name, TString suffix);
  /// The Exchange container which is to be accessed by other classes
  AliAnalysisDataContainer* GetExchangeContainter();
  virtual ~AliForwardTaskValidation() {};

  // Enums describing each event validator. These can be pushed into
  // fEventValidators by the user when configuring their task
  enum EventValidation {
      kNoEventCut,
      kPassesAliEventCuts,
      kTrigger,
      kHasFMD,
      kHasEntriesFMD,
      kHasValidFMD,
      kPassesFMD_V0CorrelatioCut
  };

  enum EventValidationMC {
    kNoEventCutMC,
    kHasEntriesFMDMC,
    kHasValidFMDMC,
    kHasPrimariesMC
  };

  // Enums describing each event validator. These can be pushed into
  // fEventValidators by the user when configuring their task
  enum TrackValidation {
    kNoTrackCut,
    kTPCOnly,
    kEtaCut,
    kPtCut
  };

  // A simple struct that combines track/tracklet/hit information such that it can be
  // used as a type for every detector
  struct Track {
    Float_t eta;
    Float_t phi;
    Float_t pt;
    Float_t weight;

    //Constructor
    Track(Float_t _eta, Float_t _phi, Float_t _pt, Float_t _weight)
      :eta(_eta), phi(_phi), pt(_pt), weight(_weight) {};
  };
  typedef std::vector<AliForwardTaskValidation::Track> Tracks;

  /// Check if the given event is Valid. Return true if so
  Bool_t IsValidEvent() {return fIsValidEvent;};

  // Get FMD hits. The eta values are calculated with respect to the
  // primary vertex of the collision
  AliForwardTaskValidation::Tracks GetFMDhits() const;

  // Get V0 hits. The eta values are calculated with respect to the
  // primary vertex of the collision(?)
  AliForwardTaskValidation::Tracks GetV0hits() const;

  // Get SPD tracklets. Note that tracklets have no pt resolution; pt is set to 0!
  AliForwardTaskValidation::Tracks GetTracklets() const;

  // Get central barrel tracks
  AliForwardTaskValidation::Tracks GetTracks();


  AliForwardSettings fSettings; 

 protected:
  /// The Holy Grail: Is this a valid event? To be read be following tasks
  Bool_t fIsValidEvent;
  /// Vector with all the event validators as enums. Can be set by the
  /// user when setting up the task
  std::vector<AliForwardTaskValidation::EventValidation> fEventValidators;
  std::vector<AliForwardTaskValidation::EventValidationMC> fEventValidatorsMC;

  /// Vector with all the _track_ validators as enums. Can be set by the
  /// user when setting up the task
  std::vector<AliForwardTaskValidation::TrackValidation> fTrackValidators;

  void UserCreateOutputObjects();
  TList *fOutputList;  //!

  void UserExec(Option_t *);
  //Bool_t UserNotify();
  /// Create QA histograms based on the set validators
  void CreateQAHistograms(TList* outlist);
  /// Histogram showing why an even got discarded to be read from left to right
  TH1F *fQA_event_discard_flow;   //!
  TH1F *fQA_event_discard_flow_MC;//!

  /// Histogram showing why a _Track_ was discarded to be read from left to right
  TH1F *fQA_track_discard_flow;//!

  // A class applying the recommended event cuts
  AliEventCuts fEventCuts;

 private:
  /// Get __ALL TRACKS__ in the central barrel
  /// No checks are done on these tracks - thus this is a private helper function.
  TClonesArray* GetAllCentralBarrelTracks();

  /// Extra cut on the FMD, rejects events with hot spots in FMD
  Bool_t HasValidFMD();

  /// Utils class used by some of the cuts
  AliAnalysisUtils fUtils;

  /// Returns alwasy true. Used to have the "all event bin" in the qa histograms
  Bool_t NoCut() {return true;};
  Bool_t HasPrimaries();

  // Check if the given event is an aod event
  Bool_t IsAODEvent();
  Bool_t AcceptTrigger(AliVEvent::EOfflineTriggerTypes TriggerType);
  /// Check if the current event has a FMD object
  Bool_t HasFMD();
  /// Check if the current event has tracklets
  Bool_t HasTracklets();
  /// Check if the current event has any counts in the FMD
  Bool_t HasEntriesFMD();
  /// Check if there is a least one channel with a signal in the V0
  Bool_t HasEntriesV0();
  /// Passes the default cuts in AliEventCuts class
  Bool_t PassesAliEventCuts();
  /// Event is not an outlier in the FMD-V0 correlation
  Bool_t PassesFMDV0CorrelatioCut(Bool_t fill_qa=false);
  /// Event is not an outlier in the FMD-V0C correlation
  /// Is not out of bunch pileup
  Bool_t NotOutOfBunchPU() {return !fUtils.IsOutOfBunchPileUp(this->InputEvent());};
  /// Is not multi-vertex pileup
  Bool_t NotMultiVertexPU() {return !fUtils.IsPileUpMV(this->InputEvent());};
  /// Is not SPD pile-up
  Bool_t NotSPDPU() {return !fUtils.IsPileUpSPD(this->InputEvent());};
  /// Is not SPD Cluster vs Tracklet background
  Bool_t NotSPDClusterVsTrackletBG() {return !fUtils.IsSPDClusterVsTrackletBG(this->InputEvent());};

  // Correlation between FMD and V0s
  TH2 *fFMDV0;       //!
  TH2 *fFMDV0_post;  //!
  TH2 *fFMDV0A;      //!
  TH2 *fFMDV0A_post; //!
  TH2 *fFMDV0C;      //!
  TH2 *fFMDV0C_post; //!
  TH2 *fOutliers;    //!

  TH1D* fCentrality;//!
  TH1D* fVertex;//!

  TH2D* centralDist;//!
  TH2D* refDist;    //!
  TH2D* forwardDist;//!

  AliForwardFlowUtil fUtil;

  ClassDef(AliForwardTaskValidation, 1);
};
#endif
