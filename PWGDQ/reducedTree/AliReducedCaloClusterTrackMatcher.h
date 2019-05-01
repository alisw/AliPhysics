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
  void SetMaximumMatchingDistance(Float_t max) { fMaxMatchingDistance=max; }
  
  // getters
  Float_t             GetMaximumMatchingDistance() const { return fMaxMatchingDistance; }
  std::vector<Int_t>  GetMatchedClusterIDsBefore() const { return fMatchedClusterIDsBefore; }
  std::vector<Int_t>  GetMatchedClusterIDsAfter() const { return fMatchedClusterIDsAfter; }

  void ClearMatchedClusterIDsBefore() { fMatchedClusterIDsBefore.clear(); } // needs to be called in analysis task after each event!
  void ClearMatchedClusterIDsAfter() { fMatchedClusterIDsAfter.clear(); }   // needs to be called in analysis task after each event!
  void FillMultipleMatchesHistogram(TH1I* hist, std::vector<Int_t> clusterIDs);
  Bool_t IsClusterMatchedToTrack(AliReducedTrackInfo* track, AliReducedCaloClusterInfo* cluster, Float_t &deltaPhi, Float_t &deltaEta, Float_t &dist);
  
private:
  
  Float_t             fMaxMatchingDistance;       // maximum distance for cluster-track matching
  std::vector<Int_t>  fMatchedClusterIDsBefore;   // cluster IDs before matching
  std::vector<Int_t>  fMatchedClusterIDsAfter;    // cluster IDs after matching

  ClassDef(AliReducedCaloClusterTrackMatcher, 1)
};

#endif
