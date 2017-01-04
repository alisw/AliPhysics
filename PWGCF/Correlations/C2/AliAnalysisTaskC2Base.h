#ifndef AliAnalysisTaskC2Base_cxx
#define AliAnalysisTaskC2Base_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisC2Settings.h"
#include "AliAnalysisC2NanoTrack.h"

class TH1F;
class THn;

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

  // Get all tracks that pass the track cut including FMD hits
  std::vector< AliAnalysisC2NanoTrack > &GetValidTracks();

  // Return true if we are currently running over an AOD dataset, else false.
  Bool_t IsAODdataset();

 protected:
  // Get an array of all (including invalid) __tracks__ in this event. Note that
  // this does not include FMD hits in reconstructed data
  TClonesArray* GetAllTracks();

  // Get FMD hits; pT is set to the minimum value of analysis!
  std::vector< AliAnalysisC2NanoTrack > GetFMDhits();

  // Get SPD hits based on clusters; pT is set to the minimum value of analysis!
  std::vector< AliAnalysisC2NanoTrack > GetSPDhits();

  // Get NanoTracks based on *tracks* from the central region, ie, they have a proper pT.
  // Only tracks passing the track selection criteria are returned.
  std::vector< AliAnalysisC2NanoTrack > GetValidCentralTracks();

  // Load and set up things that need to be updated for each event
  void SetupEventForBase();
  // Evaluate asymmetry in V0; Taken from Leonardo / Katarina; TODO:Cross-check
  Bool_t IsAsymmetricV0();

  // Vector holding all tracks that pass the track selection
  std::vector< AliAnalysisC2NanoTrack > fValidTracks;

  TList *fOutputList;  //!
  TH1 *fDiscardedEvents;  //!
  TH1 *fDiscardedTracks;  //!
  TH1F *fmultDistribution; //!
  THn *fEtaPhiZvtx_max_res;        //!

  /// Possible reasons to discard an event as use in QA histogram
  struct cDiscardEventReasons {
    enum type {
      _eventIsValid,
      isMB,
      invalidxVertex,
      isIncomplete,
      isOutOfBunchPileup,
      SPDClusterVsTrackletBG,
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
      neutralCharge,
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
