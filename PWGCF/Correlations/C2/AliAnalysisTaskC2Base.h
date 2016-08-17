#ifndef AliAnalysisTaskC2Base_cxx
#define AliAnalysisTaskC2Base_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisC2Settings.h"
#include "AliAnalysisC2NanoTrack.h"

class AliAODITSsaTrackCuts;
class TH1F;
class TH2F;

class AliAnalysisTaskC2Base : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskC2Base();
  AliAnalysisTaskC2Base(const char *name);
  virtual ~AliAnalysisTaskC2Base() {};

  AliAnalysisC2Settings fSettings;
  virtual void UserCreateOutputObjects();
  // Check if the current event is valid
  Bool_t IsValidEvent();
  // Check if the given particle fulfilles the track and quality cuts
  Bool_t IsValidParticle(AliVParticle *particle);

  // Return a value classifying the event. This might be multiplicity
  // percentile, absolute values or something like sphericity.
  // Make sure to only call this function
  // on valid events!!!
  Float_t GetEventClassifierValue();
  // Get an array of all (including invalid) tracks in this event
  TClonesArray* GetAllTracks();

  // Get all tracks that pass the track cut including FMD hits
  std::vector< AliAnalysisC2NanoTrack > &GetValidTracks();

 protected:
  // Get FMD hits
  std::vector< AliAnalysisC2NanoTrack > GetFMDhits();

  // Load and set up things that need to be updated for each event
  void SetupEventForBase();
  // Evaluate asymmetry in V0; Taken from Leonardo / Katarina; TODO:Cross-check
  Bool_t IsAsymmetricV0();
  // Taken from Leonardo / Katarina; TODO:Cross-check
  Bool_t IsOutOfBunchPileup();

  // Vector holding all tracks that pass the track selection
  std::vector< AliAnalysisC2NanoTrack > fValidTracks;

  TList *fOutputList;  //!
  AliAODITSsaTrackCuts* fitssatrackcuts;  //!
  TH1 *fDiscardedEvents;  //!
  TH1 *fDiscardedTracks;  //!
  TH1F *fmultDistribution; //!
  TH2F *fetaVsZvtx;        //!

  /// Possible reasons to discard an event as use in QA histogram
  struct cDiscardEventReasons {
    enum type {
      _eventIsValid,
      invalidxVertex,
      isIncomplete,
      isOutOfBunchPileup,
      multEstimatorNotAvailable,
      noMultSelectionObject,
      MeanMult0,
      noTracks,
      noTracksInPtRegion,
      noTrigger,
      noForwardMultObj,
      noEntriesInFMD,
      notV0AND,
      spdFastOr,
      spdPipeup,
      spdVertexContributors,
      tklClusterCut,
      v0asymmetryCut,
      zvtxPosition,
      // not a reason, just the number of reasons
      nDiscardEventReasons
    };
  };

  /// Possible reasons to discard a track as used in QA histogram
  struct cDiscardTrackReasons{
    enum type {
      _trackIsValid,
      dca,
      etaAcceptance,
      failedFilterBits,
      failedITSCut,
      neutralCharge,
      notAODPrimary,
      notHybridGCG,  // filter BIT(20)
      notMCPrimary,
      // not a reason, just the number of reasons
      nDiscardTrackReasons
    };
  };  
  struct cCachedValues{
    enum type {
      validEvent,
      validTracks,
      multiplicity,
      nCachedValues
    };
  };
  Bool_t fCachedValues[cCachedValues::nCachedValues];

  // Declaring these shuts up warnings from Weffc++
  AliAnalysisTaskC2Base(const AliAnalysisTaskC2Base&); // not implemented
  AliAnalysisTaskC2Base& operator=(const AliAnalysisTaskC2Base&); // not implemented

  ClassDef(AliAnalysisTaskC2Base, 1);
};

#endif
