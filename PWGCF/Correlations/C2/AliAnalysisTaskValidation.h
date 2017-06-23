#ifndef AliAnalysisTaskValidation_cxx
#define AliAnalysisTaskValidation_cxx

#include <string>
#include <vector>

#include "TClonesArray.h"

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

class TH2;

class AliAnalysisTaskValidation : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskValidation();
  AliAnalysisTaskValidation(const char *name);
  virtual ~AliAnalysisTaskValidation() {};

  // Enums describing each event validator. These can be pushed into
  // fValidators by the user when configuring their task
  enum EventValidation {
      kNoCut,
      kIsAODEvent,
      kHasFMD,
      kHasEntriesFMD,
      kHasEntriesV0,
      kPassesAliEventCuts,
      kPassesFMD_V0CorrelatioCut,
      kHasValidVertex,
      kHasMultSelection,
      kNotOutOfBunchPU,
      kNotMultiVertexPU,
      kNotSPDPU,
      kNotSPDClusterVsTrackletBG,
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
  typedef std::vector<AliAnalysisTaskValidation::Track> Tracks;
  
  /// Check if the given event is Valid. Return true if so
  Bool_t IsValidEvent() {return fIsValidEvent;};

  // Get FMD hits. The eta values are calculated with respect to the
  // primary vertex of the collision
  AliAnalysisTaskValidation::Tracks GetFMDhits() const;

  // Get V0 hits. The eta values are calculated with respect to the
  // primary vertex of the collision(?)
  AliAnalysisTaskValidation::Tracks GetV0hits() const;

  // Get SPD tracklets. Note that tracklets have no pt resolution; pt is set to 0!
  AliAnalysisTaskValidation::Tracks GetSPDtracklets() const;

  // Get SPD clusters. Note that closters have no pt resolution; pt is set to 0!
  AliAnalysisTaskValidation::Tracks GetSPDclusters() const;

 protected:
  /// The Holy Grail: Is this a valid event? To be read be following tasks
  Bool_t fIsValidEvent;
  /// Vector with all the event validators as enums. Can be set by the
  /// user when setting up the task
  std::vector<AliAnalysisTaskValidation::EventValidation> fValidators;

  void UserCreateOutputObjects();
  TList *fOutputList;  //!

  void UserExec(Option_t *);

  /// Create QA histograms based on the set validators
  void CreateQAHistograms(TList* outlist);
  TH1F *fQADiscard_flow;

  // A class applying the recommended event cuts
  AliEventCuts fEventCuts;

 private:
  /// Get __ALL TRACKS__ in the central barrel
  /// No checks are done on these tracks - thus this is a private helper function.
  TClonesArray* GetAllCentralBarrelTracks();

  /// Get __ALL MC TRUTH TRACKS__.
  /// No checks are done on these tracks and they could be anywhere in the detector!
  TClonesArray* GetAllMCTruthTracks();

  /// Utils class used by some of the cuts
  AliAnalysisUtils fUtils;
  /// Returns alwasy true. Used to have the "all event bin" in the qa histograms
  Bool_t NoCut() {return true;};
  // Check if the given event is an aod event
  Bool_t IsAODEvent();
  /// Check if the current event has a FMD object
  Bool_t HasFMD();
  /// Check if the current event has any counts in the FMD
  Bool_t HasEntriesFMD();
  /// Check if there is a least one channel with a signal in the V0
  Bool_t HasEntriesV0();
  /// A primary vertex exists
  Bool_t HasValidVertex();
  /// Passes the default cuts in AliEventCuts class
  Bool_t PassesAliEventCuts();
  /// Event is not an outlier in the FMD-V0 correlation
  Bool_t PassesFMDV0CorrelatioCut(Bool_t fill_qa=false);
  /// Event is not an outlier in the FMD-V0C correlation
  Bool_t HasMultSelection();
  /// Is not out of bunch pileup
  Bool_t NotOutOfBunchPU() {return !fUtils.IsOutOfBunchPileUp(this->InputEvent());};
  /// Is not multi-vertex pileup 
  Bool_t NotMultiVertexPU() {return !fUtils.IsPileUpMV(this->InputEvent());};
  /// Is not SPD pile-up
  Bool_t NotSPDPU() {return !fUtils.IsPileUpSPD(this->InputEvent());};
  /// Is not SPD Cluster vs Tracklet background
  Bool_t NotSPDClusterVsTrackletBG() {return !fUtils.IsSPDClusterVsTrackletBG(this->InputEvent());};

  // Correlation between FMD and V0s
  TH2 *fFMDV0;  //!
  TH2 *fFMDV0_post;  //!
  TH2 *fFMDV0A; //!
  TH2 *fFMDV0A_post; //!
  TH2 *fFMDV0C; //!
  TH2 *fFMDV0C_post; //!

  ClassDef(AliAnalysisTaskValidation, 1);
};
#endif


