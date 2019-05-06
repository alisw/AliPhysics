/*
 ***********************************************************
 Implementation of the AliReducedCaloClusterTrackMatcher class
 Contact: lucas.altenkamper@cern.ch
 18/04/2019
 *********************************************************
 */

#include "AliReducedCaloClusterTrackMatcher.h"

#include <iostream>
#include <fstream>
#include <map>
using namespace std;

ClassImp(AliReducedCaloClusterTrackMatcher)

//_______________________________________________________________________________
AliReducedCaloClusterTrackMatcher::AliReducedCaloClusterTrackMatcher() :
  fMaxMatchingDistance(-999.),
  fMatchedClusterIDsBefore(),
  fMatchedClusterIDsAfter()
{
  //
  // default constructor
  //
}

//_______________________________________________________________________________
AliReducedCaloClusterTrackMatcher::AliReducedCaloClusterTrackMatcher(const Char_t* name) :
  fMaxMatchingDistance(-999.),
  fMatchedClusterIDsBefore(),
  fMatchedClusterIDsAfter()
{
  //
  // named constructor
  //
}

//_______________________________________________________________________________
AliReducedCaloClusterTrackMatcher::~AliReducedCaloClusterTrackMatcher() {
  //
  // destructor
  //
}

//_______________________________________________________________________________
void AliReducedCaloClusterTrackMatcher::FillMultipleMatchesHistogram(TH1I* hist, std::vector<Int_t> clusterIDs) {
  //
  // fill an return multiple track matches histogram
  //
  if (!hist) return;
  if (!clusterIDs.size()) return;
  
  std::map<Int_t, Int_t> countMap;
  for (auto & elem : clusterIDs) {
    auto result = countMap.insert(std::pair<Int_t, Int_t>(elem, 1));
    if (result.second==kFALSE) result.first->second++;
  }
  
  for (auto & elem : countMap) {
    if (elem.second) {
      hist->Fill(elem.second);
    }
  }
}

//_______________________________________________________________________________
Bool_t AliReducedCaloClusterTrackMatcher::IsClusterMatchedToTrack(AliReducedTrackInfo* track,
                                                                  AliReducedCaloClusterInfo* cluster,
                                                                  Float_t &deltaPhi,
                                                                  Float_t &deltaEta,
                                                                  Float_t &dist) {
  //
  // check for cluster-track match below defined maximum distance
  //
  if (!track) return kFALSE;
  if (!cluster) return kFALSE;
  
  Float_t trackPhi = track->PhiOnCalo();
  Float_t trackEta = track->EtaOnCalo();
  if (trackPhi==-999. || trackEta==-999.) return kFALSE;
  
  TVector3 clusterVector(cluster->X(), cluster->Y(), cluster->Z());
  Float_t clusterPhi = clusterVector.Phi();
  if (clusterPhi<0) clusterPhi += 2*TMath::Pi();
  Float_t clusterEta = clusterVector.Eta();

  deltaPhi  = trackPhi-clusterPhi;
  deltaEta  = trackEta-clusterEta;
  dist      = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

  fMatchedClusterIDsBefore.push_back(track->CaloClusterId());
  if (dist>fMaxMatchingDistance) return kFALSE;
  fMatchedClusterIDsAfter.push_back(track->CaloClusterId());
  return kTRUE;
}
