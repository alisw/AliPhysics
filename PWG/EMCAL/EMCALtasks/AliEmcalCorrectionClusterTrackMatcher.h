#ifndef ALIEMCALCORRECTIONCLUSTERTRACKMATCHER_H
#define ALIEMCALCORRECTIONCLUSTERTRACKMATCHER_H

#include "AliEmcalCorrectionComponent.h"

class TH1;
class TClonesArray;

class AliVParticle;

/**
 * @class AliEmcalCorrectionClusterTrackMatcher
 * @ingroup EMCALCOREFW
 * @brief Cluster-track matcher component in the EMCal correction framework.
 *
 * Tracks and clusters are matched using a simple geometrical algorithm. Multiple tracks can be matched to a single cluster; however only one cluster can be matched to a track. The default configuration of the task is such that it will attempt track propagation to the EMCal surface (440 cm) if the track is not already propagated. This means that the OCDB has to be loaded beforehand (e.g. using the CDBConnect task), as well as the geometry (handled automatically by AliEmcalCorrectionTask). This should usually work in both AOD and ESD events.
 
 The number of tracks matched to a cluster can be retrieved using `cluster->GetNTracksMatched()`. Unfortunately the method to access the tracks matched to a cluster depend on the data format. For ESD clusters (AliESDCaloClusters):
 ~~~{.cxx}
 Int_t iTrack = cluster->GetTrackMatchedIndex(i);
 ~~~
 will return the position of the track in the array. The integer i is a number from 0 to `cluster->GetNTracksMatched() -1`. A pointer to the track object can be retrieved using:
 ~~~{.cxx}
 AliVTrack* track = static_cast<AliVTrack*>(GetParticleContainer(0)->GetParticle(iTrack));
 ~~~
 (assuming that the task is derived from AliAnalysisTaskEmcal or AliAnalysisTaskEmcalJet).
 
 For AOD clusters (AliAODCaloClusters) the method:
 ~~~{.cxx}
 AliVTrack* track = static_cast<AliVTrack*>(cluster->GetTrackMatched(i));
 ~~~
 will directly return a pointer to the matched track.
 
 To get the cluster matched to a track one can use (both ESD and AOD):
 ~~~{.cxx}
 Int_t iCluster = track->GetEMCALcluster();
 ~~~
 This will return the index of the cluster. To get a pointer to the cluster object:
 ~~~{.cxx}
 AliVCluster *cluster = GetClusterContainer(0)->GetCluster(iCluster);
 ~~~
 (again assuming that the task is derived from AliAnalysisTaskEmcal or AliAnalysisTaskEmcalJet).
 *
 * Based on code in AliEmcalClusTrackMatcherTask. 
 *
 * @author Constantin Loizides, LBNL, AliEmcalClusTrackMatcherTask
 * @author Salvatore Aiola, LBNL, AliEmcalClusTrackMatcherTask
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University, centralize EMCal corrections using components
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University, centralize EMCal corrections using components
 * @date Jul 8, 2016
 */

class AliEmcalCorrectionClusterTrackMatcher : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionClusterTrackMatcher();
  virtual ~AliEmcalCorrectionClusterTrackMatcher();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  Bool_t Run();
  
 protected:
  Int_t         GetMomBin(Double_t p) const;
  void          GenerateEmcalParticles();
  void          DoMatching();
  void          UpdateTracks();
  void          UpdateClusters();
  Bool_t        IsTrackInEmcalAcceptance(AliVParticle* part, Double_t edges=0.9) const;
  
  Double_t      fPropDist;              ///< distance to surface (440cm default)
  Bool_t        fDoPropagation;         ///< if true then propagate all hybrid tracks to EMCal surface
  Bool_t        fAttemptProp;           ///< if true then attempt to propagate if not done yet
  Bool_t        fAttemptPropMatch;      ///< if true then attempt to propagate if not done yet but IsEMCAL is true
  Double_t      fMaxDistance;           ///< maximum distance to match clusters and tracks
  Bool_t        fUpdateTracks;          ///< update tracks with matching info
  Bool_t        fUpdateClusters;        ///< update clusters with matching info
  
  TClonesArray *fEmcalTracks;           //!<!emcal tracks
  TClonesArray *fEmcalClusters;         //!<!emcal clusters
  Int_t         fNEmcalTracks;          //!<!number of emcal tracks
  Int_t         fNEmcalClusters;        //!<!number of emcal clusters
  TH1          *fHistMatchEtaAll;       //!<!deta distribution
  TH1          *fHistMatchPhiAll;       //!<!dphi distribution
  TH1          *fHistMatchEta[8][9][2]; //!<!deta distribution
  TH1          *fHistMatchPhi[8][9][2]; //!<!dphi distribution

 private:
  AliEmcalCorrectionClusterTrackMatcher(const AliEmcalCorrectionClusterTrackMatcher &);               // Not implemented
  AliEmcalCorrectionClusterTrackMatcher &operator=(const AliEmcalCorrectionClusterTrackMatcher &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionClusterTrackMatcher> reg;

  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionClusterTrackMatcher, 1); // EMCal cluster track matcher correction component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCLUSTERTRACKMATCHER_H */
