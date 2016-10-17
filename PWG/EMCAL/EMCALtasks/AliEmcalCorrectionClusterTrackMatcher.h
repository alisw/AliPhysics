#ifndef ALIEMCALCORRECTIONCLUSTERTRACKMATCHER_H
#define ALIEMCALCORRECTIONCLUSTERTRACKMATCHER_H

#include "AliEmcalCorrectionComponent.h"

class TH1;
class TClonesArray;

class AliVParticle;

/**
 * @class AliEmcalCorrectionClusterTrackMatcher
 * @brief Cluster-track matcher component in the EMCal correction framework
 * @ingroup EMCALCOREFW
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
