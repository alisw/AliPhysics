#ifndef ALIANALYSISETRECONSTRUCTED_H
#define ALIANALYSISETRECONSTRUCTED_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD analysis
//  - reconstruction output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEt.h"
class TH2F;
class TH3F;
class TH2D;
class TH2I;
class AliVParticle;
class AliESDEvent;
class AliAnalysisHadEtCorrections;

class AliAnalysisEtReconstructed : public AliAnalysisEt
{

public:

    AliAnalysisEtReconstructed();
    virtual ~AliAnalysisEtReconstructed();

    virtual Int_t AnalyseEvent(AliVEvent* event);

    virtual void Init();

    /** Fill the objects you want to output, classes which add new histograms should overload this. */
    virtual void FillOutputList(TList *list);
    void SetCorrections(AliAnalysisHadEtCorrections *corr) {
        fCorrections = corr;
    }

    /** Create the histograms, must be overloaded if you want to add your own */
    virtual void CreateHistograms();

    Float_t GetNumberOfChargedHadronsMatched(){return nChargedHadronsMeasured;}
    Float_t GetTotalNumberOfChargedHadrons(){return nChargedHadronsTotal;}
    
    void SetEMinCorrection(const Double_t factor) { fEMinCorrection = factor; }
    void MakeQATree(){fMakeQATree = kTRUE;}

protected:

    bool CheckGoodVertex(AliVParticle *track);
    virtual bool TrackHitsCalorimeter(AliVParticle *track, Double_t magField);

    virtual Double_t GetCorrectionModification(const AliESDCaloCluster& cluster,Int_t nonLinCorr, Int_t effCorr, Int_t cent);//nonLinCorr 0 = nominal 1 = high -1 = low, effCorr  0 = nominal 1 = high -1 = low

    TTree *fQATree; //! Tree for holding information about funny events in PHOS
    Bool_t fMakeQATree;//boolean for whether or not to make it
    Int_t fClusterMultiplicity;
    Int_t fTrackMultiplicity;
    Int_t fEventID;

    AliAnalysisHadEtCorrections *fCorrections;//!//corrections needed for hadronic et

    Double_t fPidCut; // cut on the pid probability
    Float_t nChargedHadronsMeasured;
    Float_t nChargedHadronsTotal;

    TH2F *fHistChargedPionEnergyDeposit;//! /** Energy deposited in calorimeter by charged pions */
    TH2F *fHistProtonEnergyDeposit;//! /** Energy deposited in calorimeter by protons */
    TH2F *fHistAntiProtonEnergyDeposit;//! /** Energy deposited in calorimeter by anti-protons */
    TH2F *fHistChargedKaonEnergyDeposit;//! /** Energy deposited in calorimeter by charged kaons */
    TH2F *fHistMuonEnergyDeposit;//! /** Energy deposited in calorimeter by muons */

    TH1F *fHistRemovedEnergy;//! // removed energy
    
    Double_t fGeomCorrection; // geometry correction
    Double_t fEMinCorrection; // Emin correction
    
    Double_t fRecEffCorrection; // Eff correction
    
    TH2D *fClusterPositionAccepted;//! // Position of clusters
    TH2D *fClusterPositionAll;//! // Position of clusters
    TH2D *fClusterPositionAcceptedEnergy;//! // Position of clusters
    TH2D *fClusterPositionAllEnergy;//! // Position of clusters
    TH1F *fClusterEnergy;//! // Distribution of cluster energies
    TH2F *fClusterEnergyCent;//! // Distribution of cluster energies vs centrality bin
    TH2F *fClusterEnergyModifiedTrackMatchesCent;//! // Distribution of cluster energies vs centrality bin
    TH2F *fClusterEnergyCentMatched;//! // Distribution of cluster energies vs centrality bin
    TH2F *fClusterEnergyCentNotMatched;//! // Distribution of cluster energies vs centrality bin
    TH2F *fClusterEt;//! // Distribution of cluster energies
    
    TH2D *fHistChargedEnergyRemoved;//! // Charged energy removed
    TH2D *fHistNeutralEnergyRemoved;//! // Neutral energy removed
    TH2D *fHistGammaEnergyAdded;//! // gamma energy added

    TH3F *fHistMatchedTracksEvspTvsCent;//!   //For measuring hadron deposits
    TH3F *fHistMatchedTracksEvspTvsCentEffCorr;//!   //For measuring hadron deposits
    TH3F *fHistMatchedTracksEvspTvsCentEffTMCorr;//!   //For measuring hadron deposits
    TH3F *fHistPeripheralMatchedTracksEvspTvsCentEffTMCorr;//!   //For measuring hadron deposits - uses peripheral bins and different centralities' efficiences
    TH3F *fHistMatchedTracksEvspTvsCentEffTMCorr500MeV;//!   //For measuring hadron deposits
    TH2F *fHistFoundHadronsvsCent;//!   //For measuring hadron deposits
    TH2F *fHistNotFoundHadronsvsCent;//!   //For measuring hadron deposits
    TH2F *fHistFoundHadronsEtvsCent;//!   //For measuring hadron deposits
    TH2F *fHistNotFoundHadronsEtvsCent;//!   //For measuring hadron deposits
    TH2F *fHistFoundHadronsvsCent500MeV;//!   //For measuring hadron deposits
    TH2F *fHistNotFoundHadronsvsCent500MeV;//!   //For measuring hadron deposits
    TH2F *fHistFoundHadronsEtvsCent500MeV;//!   //For measuring hadron deposits
    TH2F *fHistNotFoundHadronsEtvsCent500MeV;//!   //For measuring hadron deposits
    TH2D *fHistNominalRawEt;//!//Total ET from clusters with nominal reconstruction efficiency and nonlinearity correction vs centrality
    TH2D *fHistNominalNonLinHighEt;//!//Total ET from clusters with nominal reconstruction efficiency and high bound of nonlinearity correction vs centrality
    TH2D *fHistNominalNonLinLowEt;//!//Total ET from clusters with nominal reconstruction efficiency and low bound of nonlinearity correction vs centrality
    TH2D *fHistNominalEffHighEt;//!//Total ET from clusters with high bound on reconstruction efficiency and nominal nonlinearity correction vs centrality
    TH2D *fHistNominalEffLowEt;//!//Total ET from clusters with low bound on reconstruction efficiency and nominal nonlinearity correction vs centrality

    TH2F *fHistTotRawEtEffCorr;//! // gamma efficiency applied
    TH2F *fHistTotRawEt;//! //no  gamma efficiency applied
    TH2F *fHistTotRawEtEffCorr500MeV;//!//Total ET from clusters with nominal reconstruction efficiency and nonlinearity correction vs centrality
    TH2F *fHistTotAllRawEt;//! // all clusters no reco eff
    TH2F *fHistTotAllRawEtEffCorr;//! // all clusters reco eff applied
    Double_t ApplyModifiedCorrections(const AliESDCaloCluster& cluster,Int_t nonLinCorr, Int_t effCorr, Int_t cent);//nonLinCorr 0 = nominal 1 = high -1 = low, effCorr  0 = nominal 1 = high -1 = low

    TH2F *fHistNUsedClusters;//! // Distribution of cluster energies vs centrality bin
    TH3F *fHistNClustersPhosVsEmcal;//! // all clusters no reco eff
    TH2F *fHistClusterSizeVsCent;//! // all clusters no reco eff
    TH2F *fHistMatchedClusterSizeVsCent;//! // all clusters no reco eff
    TH2F *fHistTotAllRawEtVsTotalPt;//! // all clusters no reco eff
    //fHistTotAllRawEtVsTotalPtCent
    TH3F *fHistTotAllRawEtVsTotalPtVsCent;//! // all clusters no reco eff
    TH3F *fHistTotMatchedRawEtVsTotalPtVsCent;//! // all clusters no reco eff
    TH2F *fHistPIDProtonsTrackMatchedDepositedVsNch;//!
    TH2F *fHistPIDAntiProtonsTrackMatchedDepositedVsNch;//!
    TH2F *fHistPIDProtonsTrackMatchedDepositedVsNcl;//!
    TH2F *fHistPIDAntiProtonsTrackMatchedDepositedVsNcl;//!
    TH2F *fHistPiKPTrackMatchedDepositedVsNch;//!
    TH2F *fHistPIDProtonsTrackMatchedDepositedVsNchNoEff;//!
    TH2F *fHistPIDAntiProtonsTrackMatchedDepositedVsNchNoEff;//!
    TH2F *fHistPIDProtonsTrackMatchedDepositedVsNclNoEff;//!
    TH2F *fHistPIDAntiProtonsTrackMatchedDepositedVsNclNoEff;//!
    TH2F *fHistPiKPTrackMatchedDepositedVsNchNoEff;//!
    TH3F *fHistCentVsNchVsNclReco;//!

    TH1F *fHistRawSignalReco;//!
    TH1F *fHistEffCorrSignalReco;//!
    TH3F *fHistRecoRCorrVsPtVsCent;//! // enter comment here

private:

    AliAnalysisEtReconstructed(const AliAnalysisEtReconstructed& g);
    AliAnalysisEtReconstructed & operator=(const AliAnalysisEtReconstructed&);



    ClassDef(AliAnalysisEtReconstructed, 1);
};

#endif // ALIANALYSISETRECONSTRUCTED_H
