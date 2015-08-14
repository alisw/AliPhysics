#ifndef ALIANALYSISETMONTECARLO_H
#define ALIANALYSISETMONTECARLO_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis
//  - MC output
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________
#include "AliAnalysisEt.h"
class TParticle;
class TH3F;
class TH2I;
class AliPHOSGeoUtils;
class AliPHOSGeometry;
class AliStack;
class AliMCEvent;
class AliGenEventHeader;
class AliStack;
//class AliMCEvent;
//class AliESDEvent;

class AliAnalysisEtMonteCarlo : public AliAnalysisEt
{

public:

    AliAnalysisEtMonteCarlo();
    virtual ~AliAnalysisEtMonteCarlo();

    virtual Int_t AnalyseEvent(AliVEvent* event);
    virtual Int_t AnalyseEvent(AliVEvent* event, AliVEvent* event2);
    //virtual Int_t AnalyseEvent(AliMCEvent* event, AliESDEvent* event2);

    virtual void Init();
    virtual void ResetEventValues();
    virtual void CreateHistograms();
    virtual void FillOutputList(TList* list);

    virtual void FillHistograms();

    void CalcTrackMatchVsMult(){fCalcTrackMatchVsMult = kTRUE;}
    void CalcForKaonCorrection(){fCalcForKaonCorrection = kTRUE;}
    void IsData(){fIsMC = kFALSE;}

    void SetNumberOfChargedHadronsMatched(Float_t val){nChargedHadronsMeasured = val;}
    void SetTotalNumberOfChargedHadrons(Float_t val){nChargedHadronsTotal = val;}
    Float_t GetNumberOfChargedHadronsMatched(){return nChargedHadronsMeasured;}
    Float_t GetTotalNumberOfChargedHadrons(){return nChargedHadronsTotal;}


    void SetGeneratorMinMaxParticles(AliMCEvent *eventMC);
    AliGenEventHeader* GetGenEventHeader(AliMCEvent *eventMC) const;
    Bool_t IsHIJINGLabel(Int_t label,AliMCEvent *eventMC,AliStack *stack);
    void CheckForHIJINGLabel(){checkLabelForHIJING = kTRUE;}

protected:

    virtual bool TrackHitsCalorimeter(TParticle *part, Double_t magField=0.5);


    Int_t GetPrimMother(Int_t partIdx, AliStack *stack);

    Int_t GetK0InFamily(Int_t partIdx, AliStack *stack);

    Int_t PrintFamilyTree(Int_t partIdx, AliStack *stack);
    Int_t PrintMothers(Int_t partIdx, AliStack *stack, Int_t gen);
    Int_t PrintFamilyTreeShort(Int_t partIdx, AliStack *stack);
    Int_t PrintMothersShort(Int_t partIdx, AliStack *stack, Int_t gen);



protected:

    Float_t nChargedHadronsMeasured;
    Float_t nChargedHadronsTotal;

    Bool_t fIsMC;//if we are running over data, we still need this object to exist but we don't want to do anything.
    Bool_t checkLabelForHIJING;//Boolean - do we need to check this data set for HIJING labels?  We don't need to do this unless there were signals embedded

    Double_t fImpactParameter; // b(fm), for Hijing; 0 otherwise
    Int_t fNcoll; // Ncoll, for Hijing; 1 otherwise
    Int_t fNpart; // Ncoll, for Hijing; 2 otherwise

    TTree *fPrimaryTree; // Tree holding info on primaries

    Double_t fTotEtWithSecondaryRemoved; // enter comment here
    Double_t fTotEtSecondaryFromEmEtPrimary; // enter comment here
    Double_t fTotEtSecondary; // enter comment here
    
    Int_t fPrimaryCode; // enter comment here
    Int_t fPrimaryCharge; // enter comment here

    Double_t fPrimaryE; // enter comment here
    Double_t fPrimaryEt; // enter comment here

    Double_t fPrimaryPx; // enter comment here
    Double_t fPrimaryPy; // enter comment here
    Double_t fPrimaryPz; // enter comment here
    
    Double_t fPrimaryVx; // enter comment here
    Double_t fPrimaryVy; // enter comment here
    Double_t fPrimaryVz; // enter comment here
    
    Bool_t fPrimaryAccepted; // enter comment here
    Bool_t fPrimaryMatched;
    Int_t fDepositedCode; // enter comment here Double_t fDepositedEt; // enter comment here
    Double_t fDepositedE;
    Double_t fDepositedEt;
    Int_t fDepositedCharge; // enter comment here

    Double_t fDepositedVx; // enter comment here
    Double_t fDepositedVy; // enter comment here
    Double_t fDepositedVz; // enter comment here

    Bool_t fSecondary;

    Double_t fReconstructedE;
    Double_t fReconstructedEt;
    
    Double_t fTotPx;
    Double_t fTotPy;
    Double_t fTotPz;
    

    Int_t fClusterMult;

    TH3F *fHistDecayVertexNonRemovedCharged; //! Decay vertex for non-removed charged particles
    TH3F *fHistDecayVertexRemovedCharged; //! Decay vertex for non-removed charged particles
    TH3F *fHistDecayVertexNonRemovedNeutral; //! Decay vertex for non-removed charged particles
    TH3F *fHistDecayVertexRemovedNeutral; //! Decay vertex for non-removed charged particles

    TH2F *fHistRemovedOrNot; //! If charged/neutral particles were removed or not

    TH2F *fHistEtNonRemovedProtons; //! enter comment here
    TH2F *fHistEtNonRemovedAntiProtons; //! enter comment here
    TH2F *fHistEtNonRemovedPiPlus; //! enter comment here
    TH2F *fHistEtNonRemovedPiMinus; //! enter comment here
    TH2F *fHistEtNonRemovedKaonPlus; //! enter comment here
    TH2F *fHistEtNonRemovedKaonMinus; //! enter comment here
    TH2F *fHistEtNonRemovedK0s; //! enter comment here
    TH2F *fHistEtNonRemovedK0L; //! enter comment here
    TH2F *fHistEtNonRemovedLambdas; //! enter comment here
    TH2F *fHistEtNonRemovedElectrons; //! enter comment here
    TH2F *fHistEtNonRemovedPositrons; //! enter comment here
    TH2F *fHistEtNonRemovedMuPlus; //! enter comment here
    TH2F *fHistEtNonRemovedMuMinus; //! enter comment here
    TH2F *fHistEtNonRemovedNeutrons; //! enter comment here
    TH2F *fHistEtNonRemovedAntiNeutrons; //! enter comment here
    TH2F *fHistEtNonRemovedGammas; //! enter comment here
    TH2F *fHistEtNonRemovedGammasFromPi0; //! enter comment here

    TH2F *fHistEtRemovedGammas; //! enter comment here
    TH2F *fHistEtRemovedNeutrons; //! enter comment here
    TH2F *fHistEtRemovedAntiNeutrons; //! enter comment here

    TH2F *fHistEtRemovedCharged; //! enter comment here
    TH2F *fHistEtRemovedNeutrals; //! enter comment here

    TH2F *fHistEtNonRemovedCharged; //! enter comment here
    TH2F *fHistEtNonRemovedNeutrals; //! enter comment here

    TH2F *fHistMultNonRemovedProtons; //! enter comment here
    TH2F *fHistMultNonRemovedAntiProtons; //! enter comment here
    TH2F *fHistMultNonRemovedPiPlus; //! enter comment here
    TH2F *fHistMultNonRemovedPiMinus; //! enter comment here
    TH2F *fHistMultNonRemovedKaonPlus; //! enter comment here
    TH2F *fHistMultNonRemovedKaonMinus; //! enter comment here
    TH2F *fHistMultNonRemovedK0s; //! enter comment here
    TH2F *fHistMultNonRemovedK0L; //! enter comment here
    TH2F *fHistMultNonRemovedLambdas; //! enter comment here
    TH2F *fHistMultNonRemovedElectrons; //! enter comment here
    TH2F *fHistMultNonRemovedPositrons; //! enter comment here
    TH2F *fHistMultNonRemovedMuPlus; //! enter comment here
    TH2F *fHistMultNonRemovedMuMinus; //! enter comment here
    TH2F *fHistMultNonRemovedNeutrons; //! enter comment here
    TH2F *fHistMultNonRemovedAntiNeutrons; //! enter comment here
    TH2F *fHistMultNonRemovedGammas; //! enter comment here

    TH2F *fHistMultRemovedGammas; //! enter comment here
    TH2F *fHistMultRemovedNeutrons; //! enter comment here
    TH2F *fHistMultRemovedAntiNeutrons; //! enter comment here

    TH2F *fHistMultRemovedCharged; //! enter comment here
    TH2F *fHistMultRemovedNeutrals; //! enter comment here

    TH2F *fHistMultNonRemovedCharged; //! enter comment here
    TH2F *fHistMultNonRemovedNeutrals; //! enter comment here

    TH2F *fHistTrackMultvsNonRemovedCharged; //! enter comment here
    TH2F *fHistTrackMultvsNonRemovedNeutral; //! enter comment here
    TH2F *fHistTrackMultvsRemovedGamma; //! enter comment here

    TH2F *fHistClusterMultvsNonRemovedCharged; //! enter comment here
    TH2F *fHistClusterMultvsNonRemovedNeutral; //! enter comment here
    TH2F *fHistClusterMultvsRemovedGamma; //! enter comment here

    TH2F *fHistMultvsNonRemovedChargedE; //! enter comment here
    TH2F *fHistMultvsNonRemovedNeutralE; //! enter comment here
    TH2F *fHistMultvsRemovedGammaE; //! enter comment here

    Bool_t fCalcForKaonCorrection;//turns on and off creation of kaon correction histograms
    TH3F *fHistK0EDepositsVsPtInAcceptance; //! enter comment here
    TH3F *fHistK0EGammaVsPtInAcceptance; //! enter comment here
    TH3F *fHistK0EDepositsVsPtOutOfAcceptance; //! enter comment here
    TH3F *fHistK0EGammaVsPtOutOfAcceptance; //! enter comment here
    TH1F *fHistSimKaonsInAcceptance;//! enter comment here
    TH1F *fHistSimK0SInAcceptance;//! enter comment here
    TH1F *fHistSimKPlusInAcceptance;//! enter comment here
    TH1F *fHistSimKMinusInAcceptance;//! enter comment here
    TH1F *fHistSimK0LInAcceptance;//! enter comment here
    TH1F *fHistSimKaonsOutOfAcceptance;//! enter comment here
    TH1F *fHistSimKaonsInAcceptanceWithDepositsPrimaries;//! enter comment here
    TH1F *fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries;//! enter comment here
    TH1F *fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries;//! enter comment here

    Float_t fEtNonRemovedProtons; // enter comment here
    Float_t fEtNonRemovedAntiProtons; // enter comment here
    Float_t fEtNonRemovedPiPlus; // enter comment here
    Float_t fEtNonRemovedPiMinus; // enter comment here
    Float_t fEtNonRemovedKaonPlus; // enter comment here
    Float_t fEtNonRemovedKaonMinus; // enter comment here
    Float_t fEtNonRemovedK0S; // enter comment here
    Float_t fEtNonRemovedK0L; // enter comment here
    Float_t fEtNonRemovedLambdas; // enter comment here
    Float_t fEtNonRemovedElectrons; // enter comment here
    Float_t fEtNonRemovedPositrons; // enter comment here
    Float_t fEtNonRemovedMuMinus; // enter comment here
    Float_t fEtNonRemovedMuPlus; // enter comment here
    Float_t fEtNonRemovedGammas; // enter comment here
    Float_t fEtNonRemovedGammasFromPi0; // enter comment here
    Float_t fEtNonRemovedNeutrons; // enter comment here
    Float_t fEtNonRemovedAntiNeutrons; // enter comment here

    Float_t fEtRemovedProtons; // enter comment here
    Float_t fEtRemovedAntiProtons; // enter comment here
    Float_t fEtRemovedPiPlus; // enter comment here
    Float_t fEtRemovedPiMinus; // enter comment here
    Float_t fEtRemovedKaonPlus; // enter comment here
    Float_t fEtRemovedKaonMinus; // enter comment here
    Float_t fEtRemovedK0s; // enter comment here
    Float_t fEtRemovedK0L; // enter comment here
    Float_t fEtRemovedLambdas; // enter comment here
    Float_t fEtRemovedElectrons; // enter comment here
    Float_t fEtRemovedPositrons; // enter comment here
    Float_t fEtRemovedMuMinus; // enter comment here
    Float_t fEtRemovedMuPlus; // enter comment here

    Float_t fEtRemovedGammasFromPi0; // enter comment here
    Float_t fEtRemovedGammas; // enter comment here
    Float_t fEtRemovedNeutrons; // enter comment here
    Float_t fEtRemovedAntiNeutrons; // enter comment here

    Int_t fMultNonRemovedProtons; // enter comment here
    Int_t fMultNonRemovedAntiProtons; // enter comment here
    Int_t fMultNonRemovedPiPlus; // enter comment here
    Int_t fMultNonRemovedPiMinus; // enter comment here
    Int_t fMultNonRemovedKaonPlus; // enter comment here
    Int_t fMultNonRemovedKaonMinus; // enter comment here
    Int_t fMultNonRemovedK0s; // enter comment here
    Int_t fMultNonRemovedK0L; // enter comment here
    Int_t fMultNonRemovedLambdas; // enter comment here
    Int_t fMultNonRemovedElectrons; // enter comment here
    Int_t fMultNonRemovedPositrons; // enter comment here
    Int_t fMultNonRemovedMuMinus; // enter comment here
    Int_t fMultNonRemovedMuPlus; // enter comment here
    Int_t fMultNonRemovedGammas; // enter comment here
    Int_t fMultNonRemovedNeutrons; // enter comment here
    Int_t fMultNonRemovedAntiNeutrons; // enter comment here

    Int_t fMultRemovedProtons; // enter comment here
    Int_t fMultRemovedAntiProtons; // enter comment here
    Int_t fMultRemovedPiPlus; // enter comment here
    Int_t fMultRemovedPiMinus; // enter comment here
    Int_t fMultRemovedKaonPlus; // enter comment here
    Int_t fMultRemovedKaonMinus; // enter comment here
    Int_t fMultRemovedK0s; // enter comment here
    Int_t fMultRemovedK0L; // enter comment here

    Int_t fMultRemovedLambdas; // enter comment here
    Int_t fMultRemovedElectrons; // enter comment here
    Int_t fMultRemovedPositrons; // enter comment here
    Int_t fMultRemovedMuMinus; // enter comment here
    Int_t fMultRemovedMuPlus; // enter comment here

    Int_t fMultRemovedGammas; // enter comment here
    Int_t fMultRemovedNeutrons; // enter comment here
    Int_t fMultRemovedAntiNeutrons; // enter comment here

    Int_t fTrackMultInAcc; // enter comment here


    TH2F *fHistDxDzNonRemovedCharged; //! enter comment here
    TH2F *fHistDxDzRemovedCharged; //! enter comment here
    TH2F *fHistDxDzNonRemovedNeutral; //! enter comment here
    TH2F *fHistDxDzRemovedNeutral; //! enter comment here

    TH1F *fHistPiPlusMult; //! enter comment here
    TH1F *fHistPiMinusMult; //! enter comment here
    TH1F *fHistPiZeroMult; //! enter comment here

    TH1F *fHistPiPlusMultAcc; //! enter comment here
    TH1F *fHistPiMinusMultAcc; //! enter comment here
    TH1F *fHistPiZeroMultAcc; //! enter comment here

   // Int_t fPiPlusMult; // enter comment here
   // Int_t fPiMinusMult; // enter comment here
    
    Int_t fPiZeroMult; // enter comment here

    Int_t fPiPlusMultAcc; // enter comment here
    Int_t fPiMinusMultAcc; // enter comment here
    Int_t fPiZeroMultAcc; // enter comment here


    Int_t fNeutralRemoved; // number of neutral particles that where removed by track matching
    Int_t fChargedRemoved; // number of charged particles that where removed by track matching
    Int_t fChargedNotRemoved; // number of charged particles that were not removed
    Int_t fNeutralNotRemoved; // number of neutral particles that were not removed
    Int_t fGammaRemoved; // number of gammas removed

    Int_t fSecondaryNotRemoved;

    Double_t fEnergyNeutralRemoved; // energy of neutral particles that where removed by track matching
    Double_t fEnergyChargedRemoved; // energy of charged particles that where removed by track matching
    Double_t fEnergyChargedNotRemoved; // energy of charged particles that were not removed
    Double_t fEnergyNeutralNotRemoved; // energy of neutral particles that were not removed
    Double_t fEnergyGammaRemoved; // energy of gammas that were removed

    Int_t fNClusters; // Number of clusters in event

    Double_t fTotNeutralEtAfterMinEnergyCut; // enter comment here
    
    TH2F *fHistSimEmEtCent;//!

    Bool_t fCalcTrackMatchVsMult;
    TH1F *fHistGammasFound;//!
    TH1F *fHistGammasGenerated;//!
    TH2F *fHistGammasFoundCent;//!
    TH2F *fHistGammasFoundOutOfAccCent;//!
    TH2F *fHistGammasFoundAltCent;//!
    TH2F *fHistGammasFoundOutOfAccAltCent;//!
    TH2F *fHistGammasGeneratedCent;//!
    TH2F *fHistGammasFoundRecoEnergyCent;//!
    TH3F *fHistGammasFoundRecoEnergyTrueEnergyCent;//!
    TH2F *fHistAllGammasFoundRecoEnergyCent;//!
    TH2F *fHistGammasFoundOutOfAccRecoEnergyCent;//!
    TH2F *fHistAllGammasFoundOutOfAccRecoEnergyCent;//!
    TH1F *fHistChargedTracksCut;//!
    TH1F *fHistChargedTracksAccepted;//!
    TH1F *fHistGammasCut;//!
    TH1F *fHistGammasAccepted;//!
    TH2F *fHistChargedTrackDepositsAcceptedVsPt;//!
    TH2F *fHistChargedTrackDepositsAllVsPt;//!
    TH2F *fHistChargedTrackDepositsAcceptedVsPtEffCorr;//!
    TH2F *fHistChargedTrackDepositsAllVsPtEffCorr;//!
    TH2F *fHistChargedTracksCutMult;//!
    TH2F *fHistChargedTracksAcceptedMult;//!
    TH2F *fHistChargedTracksAcceptedLowPtCentEffCorr;//!
    TH2F *fHistChargedTracksAcceptedLowPtCent;//!
    TH2F *fHistChargedTracksAcceptedLowPtCent500MeV;//!
    TH2F *fHistChargedTracksAcceptedLowPtCentNoAntiProtons;//!
    TH2F *fHistChargedTracksAcceptedLowPtCentAntiProtons;//!
    TH2F *fHistGammasCutMult;//!
    TH2F *fHistGammasAcceptedMult;//!
    TH1F *fHistBadTrackMatches;//!
    TH2F *fHistMatchedTracksEvspTBkgd;//!
    TH2F *fHistMatchedTracksEvspTSignal;//!
    TH2F *fHistMatchedTracksEvspTBkgdPeripheral;//!
    TH2F *fHistMatchedTracksEvspTSignalPeripheral;//!
    TH3F *fHistMatchedTracksEvspTBkgdvsCent;//!
    TH3F *fHistMatchedTracksEvspTSignalvsCent;//!
    TH3F *fHistMatchedTracksEvspTBkgdvsCentEffCorr;//!
    TH3F *fHistMatchedTracksEvspTSignalvsCentEffCorr;//!
    TH1F *fHistChargedTracksCutPeripheral;//!
    TH1F *fHistChargedTracksAcceptedPeripheral;//!
    TH1F *fHistGammasCutPeripheral;//!
    TH1F *fHistGammasAcceptedPeripheral;//!
    TH2F *fHistBadTrackMatchesdPhidEta;//!
    TH2F *fHistGoodTrackMatchesdPhidEta;//!
    TH1F *fHistHadronDepositsAll;//!
    TH1F *fHistHadronDepositsReco;//!
    TH2F *fHistHadronDepositsAllCent;//!
    TH2F *fHistHadronDepositsAllCent500MeV;//!
    TH2F *fHistHadronDepositsRecoCent;//!
    TH2F *fHistHadronDepositsAllvsECent;//!
    TH2F *fHistHadronDepositsRecovsECent;//!
    TH2F *fHistHadronsAllCent;//!
    TH3F *fHistMultChVsSignalVsMult;//!
    TH2F *fHistNeutralRemovedSecondaryEtVsCent;//!
    TH2F *fHistChargedRemovedSecondaryEtVsCent;//!
    TH2F *fHistNeutralNotRemovedSecondaryEtVsCent;//!
    TH2F *fHistChargedNotRemovedSecondaryEtVsCent;//!
    TH2F *fHistNeutralRemovedSecondaryNumVsNCluster;//!
    TH2F *fHistChargedRemovedSecondaryNumVsNCluster;//!
    TH2F *fHistNeutralNotRemovedSecondaryNumVsNCluster;//!
    TH2F *fHistChargedNotRemovedSecondaryNumVsNCluster;//!
    TH2F *fHistNeutralRemovedSecondaryNumVsCent;//!
    TH2F *fHistChargedRemovedSecondaryNumVsCent;//!
    TH2F *fHistNeutralNotRemovedSecondaryNumVsCent;//!
    TH2F *fHistChargedNotRemovedSecondaryNumVsCent;//!
    TH2F *fHistNeutronsEtVsCent;//!
    TH2F *fHistNeutronsNumVsCent;//!
    TH2F *fHistNotNeutronsNumVsCent;//!
    TH2F *fHistPiKPDepositedVsNch;//!
    TH2F *fHistPiKPNotTrackMatchedDepositedVsNch;//!

    TH2F *fHistNeutronsDepositedVsNch;//!
    TH2F *fHistAntiNeutronsDepositedVsNch;//!
    TH2F *fHistProtonsDepositedVsNch;//!
    TH2F *fHistAntiProtonsDepositedVsNch;//!
    TH2F *fHistProtonsNotTrackMatchedDepositedVsNch;//!
    TH2F *fHistAntiProtonsNotTrackMatchedDepositedVsNch;//!//
    TH2F *fHistNeutronsDepositedVsNcl;//!
    TH2F *fHistAntiNeutronsDepositedVsNcl;//!
    TH2F *fHistProtonsDepositedVsNcl;//!
    TH2F *fHistAntiProtonsDepositedVsNcl;//!
    TH2F *fHistProtonsNotTrackMatchedDepositedVsNcl;//!
    TH2F *fHistAntiProtonsNotTrackMatchedDepositedVsNcl;//!
    TH2F *fHistSecondariesVsNch;//!
    TH2F *fHistSecondariesVsNcl;//!
    TH2F *fHistSecondariesEffCorrVsNch;//!
    TH2F *fHistSecondariesEffCorrVsNcl;//!
    TH2F *fHistSecondariesOutOfAccEffCorrVsNch;//!
    TH2F *fHistSecondariesDetectorCoverEffCorrVsNch;//!


    TH2F *fHistNeutronsDepositedVsNchNoEffCorr;//!//filled
    TH2F *fHistAntiNeutronsDepositedVsNchNoEffCorr;//!//filled
    TH2F *fHistProtonsDepositedVsNchNoEffCorr;//!//filled
    TH2F *fHistAntiProtonsDepositedVsNchNoEffCorr;//!//filled
    TH2F *fHistProtonsNotTrackMatchedDepositedVsNchNoEffCorr;//!//filled
    TH2F *fHistAntiProtonsNotTrackMatchedDepositedVsNchNoEffCorr;//!//filled
    TH2F *fHistNeutronsDepositedVsNclNoEffCorr;//!//filled
    TH2F *fHistAntiNeutronsDepositedVsNclNoEffCorr;//!//filled
    TH2F *fHistProtonsDepositedVsNclNoEffCorr;//!//filled
    TH2F *fHistAntiProtonsDepositedVsNclNoEffCorr;//!//filled
    TH2F *fHistProtonsNotTrackMatchedDepositedVsNclNoEffCorr;//!//filled
    TH2F *fHistAntiProtonsNotTrackMatchedDepositedVsNclNoEffCorr;//!//filled

    TH3F *fHistCentVsNchVsNcl;//!
    TH3F *fHistSecondaryPositionInDetector;//!
    //TH2F *fHistSecondaryPositionInDetector2D;//!
    TH2F *fClusterPositionWeird;//!
    //TH3F *fHistSecondaryPositionInDetectorMultiple;//!
    TH1F *fSecondaryClusterEnergy;//! // Distribution of cluster energies

    TH2F *fHistGammaCrossCheck;//!
    TH2F *fHistGammaCrossCheckAlt;//!//includes gammas that did not hit as gammas
    TH2F *fHistGammaEnergyCrossCheck;//!
    TH3F *fHistGammaEnergyCrossCheckCent;//!
    TH2F *fHistGammaEnergyCrossCheckAlt;//!//includes gammas that did not hit as gammas
    TH2F *fHistNeutronCrossCheck;//!
    TH2F *fHistSecondaryCrossCheck;//!
    TH2F *fHistHadronCrossCheck;//!
    TH2F *fHistKaonCrossCheck;//!
    TH2F *fHistNeutronCorrection;//!
    TH2F *fHistSecondaryCorrection;//!
    TH2F *fHistHadronCorrection;//!
    TH2F *fHistKaonCorrection;//!

    TH1F *fHistAllEnergy;//!//energy of all clusters passing cuts vs centrality
    TH1F *fHistSignalEnergy;//!//signal of signal clusters passing cuts vs centrality
    TH1F *fHistNeutronEnergy;//!//signal of neutron clusters passing cuts vs centrality
    TH1F *fHistKaonEnergy;//!//signal of kaon clusters passing cuts vs centrality
    TH1F *fHistHadronEnergy;//!//signal of hadron clusters passing cuts vs centrality
    TH1F *fHistSecondaryEnergy;//!//signal of secondary clusters passing cuts vs centrality
    TH1F *fHistSecondaryChargedEnergy;//!//signal of secondary clusters passing cuts vs centrality
    TH1F *fHistSecondaryNeutronEnergy;//!//signal of secondary clusters passing cuts vs centrality
    TH1F *fHistSecondaryGammaEnergy;//!//signal of secondary clusters passing cuts vs centrality
    TH1F *fHistSecondaryElectronEnergy;//!//signal of secondary clusters passing cuts vs centrality
    TH1F *fHistSecondaryOtherEnergy;//!//signal of secondary clusters passing cuts vs centrality
    TH1F *fHistSimulatedGammaEnergy;//!//signal of signal clusters passing cuts vs centrality
    TH1F *fHistReconstructedGammaEnergy;//!//signal of signal clusters passing cuts vs centrality
    TH1F *fHistSimulatedGammaEnergyAboveThreshold;//!//signal of signal clusters passing cuts vs centrality
    TH1F *fHistReconstructedSignalEnergy;//!//signal of signal clusters passing cuts vs centrality

    TH2F *fHistFracSignalVsNClusters;//!//Fraction of signal vs number of clusters
    TH2F *fHistFracHadronsVsNClusters;//!//Fraction of hadrons vs number of clusters
    TH2F *fHistFracNeutronsVsNClusters;//!//Fraction of signal vs number of clusters
    TH2F *fHistFracKaonsVsNClusters;//!//Fraction of signal vs number of clusters
    TH2F *fHistFracSecondariesVsNClusters;//!//Fraction of signal vs number of clusters
    TH2F *fHistFracSignalVsNMultiplicity;//!//Fraction of signal vs charged track multiplicity
    TH2F *fHistFracHadronsVsNMultiplicity;//!//Fraction of hadrons vs charged track multiplicity
    TH2F *fHistFracNeutronsVsNMultiplicity;//!//Fraction of signal vs charged track multiplicity
    TH2F *fHistFracKaonsVsNMultiplicity;//!//Fraction of signal vs charged track multiplicity
    TH2F *fHistFracSecondariesVsNMultiplicity;//!//Fraction of signal vs charged track multiplicity
    TH2F *fHistFracSignalVsNMatchedTracks;//!//Fraction of signal vs number of matched tracks
    TH2F *fHistFracHadronsVsNMatchedTracks;//!//Fraction of hadrons vs number of matched tracks
    TH2F *fHistFracNeutronsVsNMatchedTracks;//!//Fraction of signal vs number of matched tracks
    TH2F *fHistFracKaonsVsNMatchedTracks;//!//Fraction of signal vs number of matched tracks
    TH2F *fHistFracSecondariesVsNMatchedTracks;//!//Fraction of signal vs number of matched tracks
    TH2F *fHistFracSignalVsNTotalTracks;//!//Fraction of signal vs number of total tracks
    TH2F *fHistFracHadronsVsNTotalTracks;//!//Fraction of hadrons vs number of total tracks
    TH2F *fHistFracNeutronsVsNTotalTracks;//!//Fraction of signal vs number of total tracks
    TH2F *fHistFracKaonsVsNTotalTracks;//!//Fraction of signal vs number of total tracks
    TH2F *fHistFracSecondariesVsNTotalTracks;//!//Fraction of signal vs number of total tracks
    TH3F *fHistRCorrVsPtVsCent;//! // enter comment here

    Int_t            fNMCProducedMin;             // In case of cocktail, select particles in the list with label from this value
    Int_t            fNMCProducedMax;             // In case of cocktail, select particles in the list with label up to this value

private:

    //Declare it private to avoid compilation warning
    AliAnalysisEtMonteCarlo & operator = (const AliAnalysisEtMonteCarlo & g) ;//cpy assignment
    AliAnalysisEtMonteCarlo(const AliAnalysisEtMonteCarlo & g) ; // cpy ctor
    ClassDef(AliAnalysisEtMonteCarlo, 2);
};

#endif // ALIANALYSISETMONTECARLO_H
