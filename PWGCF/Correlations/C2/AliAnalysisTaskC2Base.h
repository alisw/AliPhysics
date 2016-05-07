#ifndef AliAnalysisTaskC2Base_cxx
#define AliAnalysisTaskC2Base_cxx

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisC2Settings.h"

class AliAODITSsaTrackCuts;
class TH1;

class AliAnalysisTaskC2Base : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskC2Base();
  AliAnalysisTaskC2Base(const char *name);
  virtual ~AliAnalysisTaskC2Base() {};

  AliAnalysisC2Settings fSettings;
  virtual void UserCreateOutputObjects();
  Bool_t IsValidEvent();
  Bool_t IsValidParticle(AliVParticle *particle);

 protected:
  // Load and set up things that need to be updated for each event
  void SetupEventForBase();
  // Evaluate asymmetry in V0; Taken from Leonardo / Katarina; TODO:Cross-check
  Bool_t IsAsymmetricV0();
  // Taken from Leonardo / Katarina; TODO:Cross-check
  Bool_t IsOutOfBunchPileup();

  TList *fOutputList;  //!
  AliAODITSsaTrackCuts* fitssatrackcuts;  //!
  TH1 *fDiscardedEvents;  //!
  TH1 *fDiscardedTracks;  //!
  /// Possible reasons to discard an event as use in QA histogram
  struct cDiscardEventReasons {
    enum type {
      _eventIsValid,
      badRun,
      invalidxVertex,
      isIncomplete,
      isOutOfBunchPileup,
      multEstimatorNotAvailable,
      noMultSelectionObject,
      noTracks,
      noTracksInPtRegion,
      noTrigger,
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

  // Declaring these shuts up warnings from Weffc++
  AliAnalysisTaskC2Base(const AliAnalysisTaskC2Base&); // not implemented
  AliAnalysisTaskC2Base& operator=(const AliAnalysisTaskC2Base&); // not implemented

  ClassDef(AliAnalysisTaskC2Base, 1);
};

#endif
