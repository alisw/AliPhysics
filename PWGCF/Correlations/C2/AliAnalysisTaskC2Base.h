#ifndef AliAnalysisTaskC2Base_cxx
#define AliAnalysisTaskC2Base_cxx

#include "TH1F.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskC2Base : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskC2Base();
  AliAnalysisTaskC2Base(const char *name, Int_t mode);
  virtual ~AliAnalysisTaskC2Base() {};

  virtual void UserCreateOutputObjects();
  Bool_t IsValidEvent();
  Bool_t IsValidParticle(AliVParticle *particle);

  enum {kMCTRUTH, kRECON};  /// Types of data this analysis can access
  const Int_t fMode;        /// The type of data this task is accessing

 protected:
  TList      *fOutputList;
  /// Possible reasons to discard an event as use in QA histogram
  struct cDiscardEventReasons {
    enum type {invalidxVertex, zvtxPosition, noTracks,
	       noTracksInPtRegion, nDiscardEventReasons}; };

  /// Possible reasons to discard a track as used in QA histogram
  struct cDiscardTrackReasons{
    enum type {neutralCharge, etaAcceptance, notHybridGCG, notAODPrimary,
	       notMCPrimary, dca, nDiscardTrackReasons};
  };

  /// Struct to bundle basic settings for this analysis
  /// The values are set in the constructor of this class 
  struct {
    Double_t etaAcceptanceLowEdge;
    Double_t etaAcceptanceUpEdge;
  } fSettings;

  TH1F *fDiscardedEvents;
  TH1F *fDiscardedTracks;

  // Declaring these shuts up warnings from Weffc++
  AliAnalysisTaskC2Base(const AliAnalysisTaskC2Base&); // not implemented
  AliAnalysisTaskC2Base& operator=(const AliAnalysisTaskC2Base&); // not implemented

  ClassDef(AliAnalysisTaskC2Base, 1); // example of analysis
};

#endif
