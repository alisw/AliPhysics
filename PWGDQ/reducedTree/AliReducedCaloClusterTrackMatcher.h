//
// Author: Lucas Altenkamper, 18/04/2019
// email: lucas.altenkamper@cern.ch
//

#ifndef ALIREDUCEDCALOCLUSTERTRACKMATCHER
#define ALIREDUCEDCALOCLUSTERTRACKMATCHER

#include <vector>
#include "TList.h"
#include "TH1.h"
#include "TVector3.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedCaloClusterInfo.h"

//_____________________________________________________________________
class AliReducedCaloClusterTrackMatcher : public TObject {

public:
  AliReducedCaloClusterTrackMatcher();
  AliReducedCaloClusterTrackMatcher(const Char_t* name);
  virtual ~AliReducedCaloClusterTrackMatcher();
  
  // setters
  void SetMaximumMatchingDistance(Float_t max);
  void SetMaximumMatchingDeltaPhi(Float_t min, Float_t max);
  void SetMaximumMatchingDeltaPhi(Float_t max);
  void SetMaximumMatchingDeltaEta(Float_t min, Float_t max);
  void SetMaximumMatchingDeltaEta(Float_t max);

  // getters
  Float_t             GetMaximumMatchingDistance() const { return fMaxMatchingDistance; }
  void                GetMaximumMatchingDeltaPhi(Float_t &min, Float_t &max) const { min=fMinMatchingDeltaPhi; max=fMaxMatchingDeltaPhi; }
  void                GetMaximumMatchingDeltaEta(Float_t &min, Float_t &max) const { min=fMinMatchingDeltaEta; max=fMaxMatchingDeltaEta; }
  std::vector<Int_t>  GetMatchedClusterIDsBefore() const { return fMatchedClusterIDsBefore; }
  std::vector<Int_t>  GetMatchedClusterIDsAfter() const { return fMatchedClusterIDsAfter; }

  void ClearMatchedClusterIDsBefore() { fMatchedClusterIDsBefore.clear(); } // needs to be called in analysis task after each event!
  void ClearMatchedClusterIDsAfter() { fMatchedClusterIDsAfter.clear(); }   // needs to be called in analysis task after each event!
  void FillMultipleMatchesHistogram(TH1I* hist, std::vector<Int_t> clusterIDs);
  Bool_t IsClusterMatchedToTrack(AliReducedTrackInfo* track, AliReducedCaloClusterInfo* cluster, Float_t &deltaPhi, Float_t &deltaEta, Float_t &dist);
  
private:
  
  Bool_t              fUseDistance;               // true: distance used for cluster-track matching
  Float_t             fMaxMatchingDistance;       // maximum distance for cluster-track matching
  Bool_t              fUseDeltaEtaDeltaPhi;       // true: delta eta, delta phi used for cluster-track matching
  Float_t             fMinMatchingDeltaPhi;       // minimum distance in phi for cluster-track matching
  Float_t             fMaxMatchingDeltaPhi;       // maximum distance in phi for cluster-track matching
  Float_t             fMinMatchingDeltaEta;       // minimum distance in eta for cluster-track matching
  Float_t             fMaxMatchingDeltaEta;       // maximum distance in eta for cluster-track matching
  std::vector<Int_t>  fMatchedClusterIDsBefore;   // cluster IDs before matching
  std::vector<Int_t>  fMatchedClusterIDsAfter;    // cluster IDs after matching

  ClassDef(AliReducedCaloClusterTrackMatcher, 2)
};

#endif
