#ifndef AliAnalysisTaskC2_cxx
#define AliAnalysisTaskC2_cxx

#include "AliAnalysisC2Settings.h"
#include "AliAnalysisTaskValidation.h"

class TH1;
class THn;

class AliAODTrack;
class AliVEvent;

class AliAnalysisTaskC2 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskC2();
  AliAnalysisTaskC2(const char *name);
  virtual ~AliAnalysisTaskC2() {};

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  // A class combining all the settings for this analysis
  AliAnalysisC2Settings fSettings;

 private:
  // Helper function to load all tracks/hits used in the analysis
  AliAnalysisTaskValidation::Tracks GetValidTracks();
  // Helper function to find the event classifier value (e.g. the
  // multiplicity in the estimator set in fSettings
  Float_t GetEventClassifierValue();
  TList *fOutputList;  //!

  // Histograms to construct C2
  THn *fEventCounter;  //!

  // Eta, phi, zvtx track count at max resolution; used to identifie dead regions in post processing
  THn *fEtaPhiZvtx_max_res; //!

  std::vector< THn* > fsingleHists; //!
  std::vector< THn* > fpairHists; //!
  // Get the index in the `fpairsHists` for this pair of tracks eta1 <
  // eta2. The point here is that the call to filling a histogram is
  // very slow and is done in a nested loop. At the same time, Fill()
  // does not show if it was an over or underflow. This is an attempt
  // to make things smarter.
  UInt_t GetPairHistIndex(Float_t eta1, Float_t eta2);
  struct cEventCounterDims {
    enum type {kMult, kZvtx, kNdimensions};
  };
  struct cPairsDims {
    enum type {kEta1, kEta2, kPhi1, kPhi2, kPtPair, kMult, kZvtx, kNdimensions};
  };
  struct cSinglesDims {
    enum type {kEta, kPhi, kPt, kMult, kZvtx, kNdimensions};
  };

  // Declaring these shuts up warnings from Weffc++
  AliAnalysisTaskC2(const AliAnalysisTaskC2&); // not implemented
  AliAnalysisTaskC2& operator=(const AliAnalysisTaskC2&); // not implemented

  ClassDef(AliAnalysisTaskC2, 1); // example of analysis
};

#endif
