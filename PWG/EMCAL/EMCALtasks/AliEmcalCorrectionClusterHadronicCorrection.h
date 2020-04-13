#ifndef ALIEMCALCORRECTIONCLUSTERHADRONICCORRECTION_H
#define ALIEMCALCORRECTIONCLUSTERHADRONICCORRECTION_H

#include "AliEmcalCorrectionComponent.h"

#if !(defined(__CINT__) || defined(__MAKECINT__))
#include "AliEmcalContainerIndexMap.h"
#endif

class TH1;
class TH2;

/**
 * @class AliEmcalCorrectionClusterHadronicCorrection
 * @ingroup EMCALCORRECTIONFW
 * @brief Hadronic correction component in the EMCal correction framework.
 *
 * Charged particles deposit some energy in the calorimeter. Most of the charged particle are hadrons, such as pions, kaons and protons. The hadronic response of the calorimeter has been studied in some details. Most of the high energetic particles ( > 1 GeV) only release a small amount of energy. These are usually called "minimum ionizing particles" (MIP). Occasionally hadrons may interact strongly with the nuclei of the material in the calorimeter and start a hadronic shower. In this case the energy deposition is much higher. High momentum muons are also MIP, but they never shower in the calorimeter. Finally electrons do shower in the calorimeter, in a way that is quite similar to a shower initiated by a photon.
 
 All charged particles are measured in the tracking detectors (ITS+TPC+TOF). To avoid double counting their contribution to the jet energy flow, they need to be subtracted cluster-by-cluster. This is done by matching EMCal/DCal clusters with charged tracks and then subtracting a certain fraction of the matched tracks from the cluster energy. The most common choice is to subtract 100% of the momentum of the sum of the matched tracks.
 
 The energy of the cluster **after** the hadronic correction can be obtained using the method `cluster->GetHadCorrEnergy()`.
 *
 * Based on code in AliHadCorrTask.
 *
 * @author Rosi Reed, AliHadCorrTask
 * @author Constantin Loizides, LBNL, AliHadCorrTask
 * @author Salvatore Aiola, LBNL, AliHadCorrTask
 * @author James Mulligan <james.mulligan@yale.edu>, Yale University, centralize EMCal corrections using components
 * @author Raymond Ehlers <raymond.ehlers@yale.edu>, Yale University, centralize EMCal corrections using components
 * @date Jul 8, 2016
 */

class AliEmcalCorrectionClusterHadronicCorrection : public AliEmcalCorrectionComponent {
 public:
  AliEmcalCorrectionClusterHadronicCorrection();
  virtual ~AliEmcalCorrectionClusterHadronicCorrection();

  // Sets up and runs the task
  Bool_t Initialize();
  void UserCreateOutputObjects();
  void ExecOnce();
  Bool_t Run();
  
protected:
  void                   GenerateTrackToContainerMap();
  Double_t               ApplyHadCorrOneTrack(Int_t icluster, Double_t hadCorr);
  Double_t               ApplyHadCorrAllTracks(Int_t icluster, Double_t hadCorr);
  void                   DoMatchedTracksLoop(Int_t icluster, Double_t &totalTrkP, Int_t &Nmatches, Double_t &trkPMCfrac, Int_t &NMCmatches);
  void                   DoTrackLoop();
  Double_t               GetEtaSigma(Int_t pbin)                   const;
  UInt_t                 GetMomBin(Double_t pt)                    const;
  Double_t               GetPhiMean(Int_t pbin, Int_t centbin)     const;
  Double_t               GetPhiSigma(Int_t pbin, Int_t centbin)    const;
  Double_t               ComputeM02Subtraction(const AliVCluster* cluster, Double_t energyclus, Int_t Nmatches, Double_t totalTrkP, Double_t hadCorr);
  
  // Task configuration
  Double_t               fPhiMatch;                  ///< phi match value (pp=0.050)
  Double_t               fEtaMatch;                  ///< eta match value (pp=0.025)
  Bool_t                 fDoTrackClus;               ///< loop over tracks first
  Bool_t                 fDoMomDepMatching;          ///< set the values fPhiMatch and fEtaMatch depending on the track momentum
  Double_t               fHadCorr;                   ///< hadronic correction (fraction)
  Double_t               fEexclCell;                 ///< energy/cell that we cannot subtract from the clusters
  Bool_t                 fPlotOversubtractionHistograms; ///< compute and plot oversubtracted energy from embedded/signal matches (embedding only)
  Bool_t                 fDoNotOversubtract;         ///< do not oversubtract energy from cluster (embedding only)
  Bool_t                 fUseM02SubtractionScheme;   ///< Flag to enable hadronic correction scheme using cluster M02 value
  Bool_t                 fUseConstantSubtraction;    ///< Flag to perform constant rather than fractional subtract (only applicable if using M02 scheme)
  Double_t               fConstantSubtractionValue;  ///< Value to be used for constant subtraction (only applicable if using constant subtraction in M02 scheme)

#if !(defined(__CINT__) || defined(__MAKECINT__))
  // Handle mapping between index and containers
  AliEmcalContainerIndexMap <AliClusterContainer, AliVCluster> fClusterContainerIndexMap;    //!<! Mapping between index and cluster containers
  AliEmcalContainerIndexMap <AliParticleContainer, AliVParticle> fParticleContainerIndexMap; //!<! Mapping between index and particle containers
#endif
  std::map <AliVTrack *, AliParticleContainer *> fTrackToContainerMap;                       //!<! Mapping between AliVTracks and their respective particle containers. Needed for AODs only. See GenerateTrackToContainerMap() for more information.
  
  // QA plots
  TH2                   *fHistMatchEtaPhi[10][9][2];  //!<!deta vs. dphi of matched cluster-track pairs
  TH2                   *fHistMatchEtaPhiAll;        //!<!deta vs. dphi of matched cluster-track pairs
  TH2                   *fHistMatchEtaPhiAllCl;      //!<!deta vs. dphi of all cluster-track pairs (cl loop)
  TH2                   *fHistMatchEvsP[5];          //!<!cluster energy vs. track momentum of matched pairs
  TH2                   *fHistNMatchEnergy[5];       //!<!n matches vs. cluster energy
  TH2                   *fHistNCellsEnergy[5][4];    //!<!n cells vs. cluster energy
  TH2                   *fHistMatchdRvsEP[5];        //!<!matching distance vs. E/P
  TH1                   *fHistNclusvsCent;           //!<!n clusters vs. centrality
  TH1                   *fHistNclusMatchvsCent;      //!<!n clusters matched to some track vs. centrality
  TH1                   *fHistEbefore;               //!<!average energy of clusters before correction vs. centrality
  TH1                   *fHistEafter;                //!<!average energy of clusters after correction vs. centrality
  TH2                   *fHistEoPCent;               //!<!E/P vs. centrality
  TH2                   *fHistNMatchCent;            //!<!n matches vs. centraity
  TH2                   *fHistNClusMatchCent;        //!<!n clusters macthed to some track (tracks allowed to match more than one cluster)
  TH1                   *fHistEsubPch[10];           //!<!Esub vs. total momentum of matched tracks (only 1 match)
  TH2                   *fHistEsubPchRat[10];        //!<!Esub/momentum of matched tracks vs. total momentum of matched tracks (only 1 match)
  TH2                   *fHistEsubPchRatAll[5];      //!<!Esub/momentum of matched tracks vs. total momentum of matched tracks (all number of matches)
  TH2                   *fHistEmbTrackMatchesOversub[5];    //!<!Over-subtracted energy / cluster energy with embedded track matches (non-embedded matches < 5%)
  TH2                   *fHistNonEmbTrackMatchesOversub[5]; //!<!Over-subtracted energy / cluster energy with non-embedded track matches (embedded matches < 5%)
  TH2                   *fHistOversubMCClusters[5];         //!<!Over-subtracted energy / cluster energy (cluster MC energy fraction > 95%)
  TH2                   *fHistOversubNonMCClusters[5];      //!<!Over-subtracted energy / cluster energy (cluster MC energy fraction < 5%)
  TH2                   *fHistOversub[5];                   //!<!Over-subtracted energy / cluster energy

 private:
  AliEmcalCorrectionClusterHadronicCorrection(const AliEmcalCorrectionClusterHadronicCorrection &);               // Not implemented
  AliEmcalCorrectionClusterHadronicCorrection &operator=(const AliEmcalCorrectionClusterHadronicCorrection &);    // Not implemented

  // Allows the registration of the class so that it is available to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionClusterHadronicCorrection> reg;
  
  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionClusterHadronicCorrection, 5); // EMCal cluster hadronic correction component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCLUSTERHADRONICCORRECTION_H */
