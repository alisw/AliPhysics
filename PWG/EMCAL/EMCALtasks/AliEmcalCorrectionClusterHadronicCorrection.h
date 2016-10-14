#ifndef ALIEMCALCORRECTIONCLUSTERHADRONICCORRECTION_H
#define ALIEMCALCORRECTIONCLUSTERHADRONICCORRECTION_H

#include "AliEmcalCorrectionComponent.h"

class TH1;
class TH2;

/**
 * @class AliEmcalCorrectionClusterHadronicCorrection
 * @brief Hadronic correction component in the EMCal correction framework
 * @ingroup EMCALCOREFW
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
  Bool_t Run();
  
protected:
  Double_t               ApplyHadCorrOneTrack(Int_t icluster, Double_t hadCorr);
  Double_t               ApplyHadCorrAllTracks(Int_t icluster, Double_t hadCorr);
  void                   DoMatchedTracksLoop(Int_t icluster, Double_t &totalTrkP, Int_t &Nmatches, Double_t &trkPMCfrac, Int_t &NMCmatches);
  void                   DoTrackLoop();
  Double_t               GetEtaSigma(Int_t pbin)                   const;
  UInt_t                 GetMomBin(Double_t pt)                    const;
  Double_t               GetPhiMean(Int_t pbin, Int_t centbin)     const;
  Double_t               GetPhiSigma(Int_t pbin, Int_t centbin)    const;
  
  // Task configuration
  Double_t               fPhiMatch;                  ///< phi match value (pp=0.050)
  Double_t               fEtaMatch;                  ///< eta match value (pp=0.025)
  Bool_t                 fDoTrackClus;               ///< loop over tracks first
  Double_t               fHadCorr;                   ///< hadronic correction (fraction)
  Double_t               fEexclCell;                 ///< energy/cell that we cannot subtract from the clusters
  Bool_t                 fDoExact;                   ///< do exact correction (embedding only)
  
  // QA plots
  TH2                   *fHistMatchEtaPhi[8][9][2];  //!<!deta vs. dphi of matched cluster-track pairs
  TH2                   *fHistMatchEtaPhiAll;        //!<!deta vs. dphi of matched cluster-track pairs
  TH2                   *fHistMatchEtaPhiAllTr;      //!<!deta vs. dphi of all cluster-track pairs (tr loop)
  TH2                   *fHistMatchEtaPhiAllCl;      //!<!deta vs. dphi of all cluster-track pairs (cl loop)
  TH2                   *fHistMatchEvsP[4];          //!<!cluster energy vs. track momentum of matched pairs
  TH2                   *fHistNMatchEnergy[4];       //!<!n matches vs. cluster energy
  TH2                   *fHistNCellsEnergy[4][4];    //!<!n cells vs. cluster energy
  TH2                   *fHistMatchdRvsEP[4];        //!<!matching distance vs. E/P
  TH1                   *fHistNclusvsCent;           //!<!n clusters vs. centrality
  TH1                   *fHistNclusMatchvsCent;      //!<!n clusters matched to some track vs. centrality
  TH1                   *fHistEbefore;               //!<!average energy of clusters before correction vs. centrality
  TH1                   *fHistEafter;                //!<!average energy of clusters after correction vs. centrality
  TH2                   *fHistEoPCent;               //!<!E/P vs. centrality
  TH2                   *fHistNMatchCent;            //!<!n matches vs. centraity
  TH2                   *fHistNClusMatchCent;        //!<!n clusters macthed to some track (tracks allowed to match more than one cluster)
  TH1                   *fHistEsubPch[8];            //!<!Esub vs. total momentum of matched tracks (only 1 match)
  TH2                   *fHistEsubPchRat[8];         //!<!Esub/momentum of matched tracks vs. total momentum of matched tracks (only 1 match)
  TH2                   *fHistEsubPchRatAll[8];      //!<!Esub/momentum of matched tracks vs. total momentum of matched tracks (all number of matches)
  TH2                   *fHistEmbTrackMatchesOversub[4];    //!<!Over-subtracted energy / cluster energy with embedded track matches (non-embedded matches < 5%)
  TH2                   *fHistNonEmbTrackMatchesOversub[4]; //!<!Over-subtracted energy / cluster energy with non-embedded track matches (embedded matches < 5%)
  TH2                   *fHistOversubMCClusters[4];         //!<!Over-subtracted energy / cluster energy (cluster MC energy fraction > 95%)
  TH2                   *fHistOversubNonMCClusters[4];      //!<!Over-subtracted energy / cluster energy (cluster MC energy fraction < 5%)
  TH2                   *fHistOversub[4];                   //!<!Over-subtracted energy / cluster energy

 private:
  AliEmcalCorrectionClusterHadronicCorrection(const AliEmcalCorrectionClusterHadronicCorrection &);               // Not implemented
  AliEmcalCorrectionClusterHadronicCorrection &operator=(const AliEmcalCorrectionClusterHadronicCorrection &);    // Not implemented

  // Allows the registration of the class so that it is availble to be used by the correction task.
  static RegisterCorrectionComponent<AliEmcalCorrectionClusterHadronicCorrection> reg;
  
  /// \cond CLASSIMP
  ClassDef(AliEmcalCorrectionClusterHadronicCorrection, 1); // EMCal cluster hadronic correction component
  /// \endcond
};

#endif /* ALIEMCALCORRECTIONCLUSTERHADRONICCORRECTION_H */
