//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis
//  - MC output
//  implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "AliAnalysisEtMonteCarlo.h"
#include "AliAnalysisEtCuts.h"
#include "AliAnalysisEtSelector.h"
#include "AliAnalysisEtSelector.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliMCEvent.h"
#include "AliESDEvent.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "TList.h"
#include "AliESDCaloCluster.h"
#include "AliLog.h"
#include <iostream>
#include <AliCentrality.h>
#include "AliPHOSGeoUtils.h"
#include "AliPHOSGeometry.h"
#include "TFile.h"
#include "AliESDtrackCuts.h"
#include "AliGenEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
using namespace std;

ClassImp(AliAnalysisEtMonteCarlo);


// ctor
AliAnalysisEtMonteCarlo::AliAnalysisEtMonteCarlo():AliAnalysisEt()
						  ,nChargedHadronsMeasured(0)
						  ,nChargedHadronsTotal(0)
						  ,fIsMC(kTRUE)
						  ,checkLabelForHIJING(kFALSE)
    ,fImpactParameter(0)
    ,fNcoll(0)
    ,fNpart(0)
    ,fPrimaryTree(0)
    ,fTotEtWithSecondaryRemoved(0)
    ,fTotEtSecondaryFromEmEtPrimary(0)
    ,fTotEtSecondary(0)
    ,fPrimaryCode(0)
    ,fPrimaryCharge(0)
    ,fPrimaryE(0)
    ,fPrimaryEt(0)
    ,fPrimaryPx(0)
    ,fPrimaryPy(0)
    ,fPrimaryPz(0)
    ,fPrimaryVx(0)
    ,fPrimaryVy(0)
    ,fPrimaryVz(0)
    ,fPrimaryAccepted(0)
    ,fPrimaryMatched(0)
    ,fDepositedCode(0)
    ,fDepositedE(0)
    ,fDepositedEt(0)
    ,fDepositedCharge(0)
    ,fDepositedVx(0)
    ,fDepositedVy(0)
    ,fDepositedVz(0)
    ,fSecondary(kFALSE)
    ,fReconstructedE(0)
    ,fReconstructedEt(0)
    ,fTotPx(0)
    ,fTotPy(0)
    ,fTotPz(0)
    ,fClusterMult(0)
    ,fHistDecayVertexNonRemovedCharged(0)
    ,fHistDecayVertexRemovedCharged(0)
    ,fHistDecayVertexNonRemovedNeutral(0)
    ,fHistDecayVertexRemovedNeutral(0)

    ,fHistRemovedOrNot(0)
    ,fHistEtNonRemovedProtons(0)
    ,fHistEtNonRemovedAntiProtons(0)
    ,fHistEtNonRemovedPiPlus(0)
    ,fHistEtNonRemovedPiMinus(0)
    ,fHistEtNonRemovedKaonPlus(0)
    ,fHistEtNonRemovedKaonMinus(0)
    ,fHistEtNonRemovedK0s(0)
    ,fHistEtNonRemovedK0L(0)
    ,fHistEtNonRemovedLambdas(0)
    ,fHistEtNonRemovedElectrons(0)
    ,fHistEtNonRemovedPositrons(0)
    ,fHistEtNonRemovedMuPlus(0)
    ,fHistEtNonRemovedMuMinus(0)
    ,fHistEtNonRemovedNeutrons(0)
    ,fHistEtNonRemovedAntiNeutrons(0)
    ,fHistEtNonRemovedGammas(0)
    ,fHistEtNonRemovedGammasFromPi0(0)
    ,fHistEtRemovedGammas(0)
    ,fHistEtRemovedNeutrons(0)
    ,fHistEtRemovedAntiNeutrons(0)
    ,fHistEtRemovedCharged(0)
    ,fHistEtRemovedNeutrals(0)
    ,fHistEtNonRemovedCharged(0)
    ,fHistEtNonRemovedNeutrals(0)
    ,fHistMultNonRemovedProtons(0)
    ,fHistMultNonRemovedAntiProtons(0)
    ,fHistMultNonRemovedPiPlus(0)
    ,fHistMultNonRemovedPiMinus(0)
    ,fHistMultNonRemovedKaonPlus(0)
    ,fHistMultNonRemovedKaonMinus(0)
    ,fHistMultNonRemovedK0s(0)
    ,fHistMultNonRemovedK0L(0)
    ,fHistMultNonRemovedLambdas(0)
    ,fHistMultNonRemovedElectrons(0)
    ,fHistMultNonRemovedPositrons(0)
    ,fHistMultNonRemovedMuPlus(0)
    ,fHistMultNonRemovedMuMinus(0)
    ,fHistMultNonRemovedNeutrons(0)
    ,fHistMultNonRemovedAntiNeutrons(0)
    ,fHistMultNonRemovedGammas(0)
    ,fHistMultRemovedGammas(0)
    ,fHistMultRemovedNeutrons(0)
    ,fHistMultRemovedAntiNeutrons(0)
    ,fHistMultRemovedCharged(0)
    ,fHistMultRemovedNeutrals(0)
    ,fHistMultNonRemovedCharged(0)
    ,fHistMultNonRemovedNeutrals(0)
    ,fHistTrackMultvsNonRemovedCharged(0)
    ,fHistTrackMultvsNonRemovedNeutral(0)
    ,fHistTrackMultvsRemovedGamma(0)
    ,fHistClusterMultvsNonRemovedCharged(0)
    ,fHistClusterMultvsNonRemovedNeutral(0)
    ,fHistClusterMultvsRemovedGamma(0)
    ,fHistMultvsNonRemovedChargedE(0)
    ,fHistMultvsNonRemovedNeutralE(0)
    ,fHistMultvsRemovedGammaE(0)
						  ,fCalcForKaonCorrection(kFALSE)
						  ,fHistK0EDepositsVsPtInAcceptance(0)
						  ,fHistK0EGammaVsPtInAcceptance(0)
						  ,fHistK0EDepositsVsPtOutOfAcceptance(0)
						  ,fHistK0EGammaVsPtOutOfAcceptance(0)
						  ,fHistSimKaonsInAcceptance(0)
						  ,fHistSimK0SInAcceptance(0)
						  ,fHistSimKPlusInAcceptance(0)
						  ,fHistSimKMinusInAcceptance(0)
						  ,fHistSimK0LInAcceptance(0)
						  ,fHistSimKaonsOutOfAcceptance(0)
						  ,fHistSimKaonsInAcceptanceWithDepositsPrimaries(0)
						  ,fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries(0)
						  ,fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries(0)
    ,fEtNonRemovedProtons(0)
    ,fEtNonRemovedAntiProtons(0)
    ,fEtNonRemovedPiPlus(0)
    ,fEtNonRemovedPiMinus(0)
    ,fEtNonRemovedKaonPlus(0)
    ,fEtNonRemovedKaonMinus(0)
    ,fEtNonRemovedK0S(0)
    ,fEtNonRemovedK0L(0)
    ,fEtNonRemovedLambdas(0)
    ,fEtNonRemovedElectrons(0)
    ,fEtNonRemovedPositrons(0)
    ,fEtNonRemovedMuMinus(0)
    ,fEtNonRemovedMuPlus(0)
    ,fEtNonRemovedGammas(0)
    ,fEtNonRemovedGammasFromPi0(0)
    ,fEtNonRemovedNeutrons(0)
    ,fEtNonRemovedAntiNeutrons(0)
    ,fEtRemovedProtons(0)
    ,fEtRemovedAntiProtons(0)
    ,fEtRemovedPiPlus(0)
    ,fEtRemovedPiMinus(0)
    ,fEtRemovedKaonPlus(0)
    ,fEtRemovedKaonMinus(0)
    ,fEtRemovedK0s(0)
    ,fEtRemovedK0L(0)
    ,fEtRemovedLambdas(0)
    ,fEtRemovedElectrons(0)
    ,fEtRemovedPositrons(0)
    ,fEtRemovedMuMinus(0)
    ,fEtRemovedMuPlus(0)
    ,fEtRemovedGammasFromPi0(0)
    ,fEtRemovedGammas(0)
    ,fEtRemovedNeutrons(0)
    ,fEtRemovedAntiNeutrons(0)
    ,fMultNonRemovedProtons(0)
    ,fMultNonRemovedAntiProtons(0)
    ,fMultNonRemovedPiPlus(0)
    ,fMultNonRemovedPiMinus(0)
    ,fMultNonRemovedKaonPlus(0)
    ,fMultNonRemovedKaonMinus(0)
    ,fMultNonRemovedK0s(0)
    ,fMultNonRemovedK0L(0)
    ,fMultNonRemovedLambdas(0)
    ,fMultNonRemovedElectrons(0)
    ,fMultNonRemovedPositrons(0)
    ,fMultNonRemovedMuMinus(0)
    ,fMultNonRemovedMuPlus(0)
    ,fMultNonRemovedGammas(0)
    ,fMultNonRemovedNeutrons(0)
    ,fMultNonRemovedAntiNeutrons(0)
    ,fMultRemovedProtons(0)
    ,fMultRemovedAntiProtons(0)
    ,fMultRemovedPiPlus(0)
    ,fMultRemovedPiMinus(0)
    ,fMultRemovedKaonPlus(0)
    ,fMultRemovedKaonMinus(0)
    ,fMultRemovedK0s(0)
    ,fMultRemovedK0L(0)
    ,fMultRemovedLambdas(0)
    ,fMultRemovedElectrons(0)
    ,fMultRemovedPositrons(0)
    ,fMultRemovedMuMinus(0)
    ,fMultRemovedMuPlus(0)
    ,fMultRemovedGammas(0)
    ,fMultRemovedNeutrons(0)
    ,fMultRemovedAntiNeutrons(0)
    ,fTrackMultInAcc(0)
    ,fHistDxDzNonRemovedCharged(0)
    ,fHistDxDzRemovedCharged(0)
    ,fHistDxDzNonRemovedNeutral(0)
    ,fHistDxDzRemovedNeutral(0)
    ,fHistPiPlusMult(0)
    ,fHistPiMinusMult(0)
    ,fHistPiZeroMult(0)
    ,fHistPiPlusMultAcc(0)
    ,fHistPiMinusMultAcc(0)
    ,fHistPiZeroMultAcc(0)
//     ,fPiPlusMult(0)
//     ,fPiMinusMult(0)
    ,fPiZeroMult(0)
    ,fPiPlusMultAcc(0)
    ,fPiMinusMultAcc(0)
    ,fPiZeroMultAcc(0)
    ,fNeutralRemoved(0)
    ,fChargedRemoved(0)
    ,fChargedNotRemoved(0)
    ,fNeutralNotRemoved(0)
    ,fGammaRemoved(0)
    ,fSecondaryNotRemoved(0)
    ,fEnergyNeutralRemoved(0)
    ,fEnergyChargedRemoved(0)
    ,fEnergyChargedNotRemoved(0)
    ,fEnergyNeutralNotRemoved(0)
    ,fEnergyGammaRemoved(0)
    ,fNClusters(0)
    ,fTotNeutralEtAfterMinEnergyCut(0)
						  ,fHistSimEmEtCent(0)
						  ,fCalcTrackMatchVsMult(kFALSE)
						  ,fHistGammasFound(0)
						  ,fHistGammasGenerated(0)
						  ,fHistGammasFoundCent(0)
						  ,fHistGammasFoundOutOfAccCent(0)
						  ,fHistGammasFoundAltCent(0)
						  ,fHistGammasFoundOutOfAccAltCent(0)
						  ,fHistGammasGeneratedCent(0)
						  ,fHistGammasFoundRecoEnergyCent(0)
						  ,fHistGammasFoundRecoEnergyTrueEnergyCent(0)
						  ,fHistAllGammasFoundRecoEnergyCent(0)
						  ,fHistGammasFoundOutOfAccRecoEnergyCent(0)
						  ,fHistAllGammasFoundOutOfAccRecoEnergyCent(0)
						  ,fHistChargedTracksCut(0)
						  ,fHistChargedTracksAccepted(0)
						  ,fHistGammasCut(0)
						  ,fHistGammasAccepted(0)
						  ,fHistChargedTrackDepositsAcceptedVsPt(0)
						  ,fHistChargedTrackDepositsAllVsPt(0)
						  ,fHistChargedTrackDepositsAcceptedVsPtEffCorr(0)
						  ,fHistChargedTrackDepositsAllVsPtEffCorr(0)
						  ,fHistChargedTracksCutMult(0)
						  ,fHistChargedTracksAcceptedMult(0)
						  ,fHistChargedTracksAcceptedLowPtCentEffCorr(0)
						  ,fHistChargedTracksAcceptedLowPtCent(0)
						  ,fHistChargedTracksAcceptedLowPtCent500MeV(0)
						  ,fHistChargedTracksAcceptedLowPtCentNoAntiProtons(0)
						  ,fHistChargedTracksAcceptedLowPtCentAntiProtons(0)
						  ,fHistGammasCutMult(0)
						  ,fHistGammasAcceptedMult(0)
						  ,fHistBadTrackMatches(0)
						  ,fHistMatchedTracksEvspTBkgd(0)
						  ,fHistMatchedTracksEvspTSignal(0)
						  ,fHistMatchedTracksEvspTBkgdPeripheral(0)
						  ,fHistMatchedTracksEvspTSignalPeripheral(0)
						  ,fHistMatchedTracksEvspTBkgdvsCent(0)
						  ,fHistMatchedTracksEvspTSignalvsCent(0)
						  ,fHistMatchedTracksEvspTBkgdvsCentEffCorr(0)
						  ,fHistMatchedTracksEvspTSignalvsCentEffCorr(0)

						  ,fHistChargedTracksCutPeripheral(0)
						  ,fHistChargedTracksAcceptedPeripheral(0)
						  ,fHistGammasCutPeripheral(0)
						  ,fHistGammasAcceptedPeripheral(0)
						  ,fHistBadTrackMatchesdPhidEta(0)
						  ,fHistGoodTrackMatchesdPhidEta(0)
						  ,fHistHadronDepositsAll(0)
						  ,fHistHadronDepositsReco(0)
						  ,fHistHadronDepositsAllCent(0)
						  ,fHistHadronDepositsAllCent500MeV(0)
						  ,fHistHadronDepositsRecoCent(0)
						  ,fHistHadronDepositsAllvsECent(0)
						  ,fHistHadronDepositsRecovsECent(0)
						  ,fHistHadronsAllCent(0)
						  ,fHistMultChVsSignalVsMult(0)
						  ,fHistNeutralRemovedSecondaryEtVsCent(0)
						  ,fHistChargedRemovedSecondaryEtVsCent(0)
						  ,fHistNeutralNotRemovedSecondaryEtVsCent(0)
						  ,fHistChargedNotRemovedSecondaryEtVsCent(0)
						  ,fHistNeutralRemovedSecondaryNumVsNCluster(0)
						  ,fHistChargedRemovedSecondaryNumVsNCluster(0)
						  ,fHistNeutralNotRemovedSecondaryNumVsNCluster(0)
						  ,fHistChargedNotRemovedSecondaryNumVsNCluster(0)
						  ,fHistNeutralRemovedSecondaryNumVsCent(0)
						  ,fHistChargedRemovedSecondaryNumVsCent(0)
						  ,fHistNeutralNotRemovedSecondaryNumVsCent(0)
						  ,fHistChargedNotRemovedSecondaryNumVsCent(0)
						  ,fHistNeutronsEtVsCent(0)
						  ,fHistNeutronsNumVsCent(0)
						  ,fHistNotNeutronsNumVsCent(0)
						  ,fHistPiKPDepositedVsNch(0)
						  ,fHistPiKPNotTrackMatchedDepositedVsNch(0)
						  ,fHistNeutronsDepositedVsNch(0)//start
						  ,fHistAntiNeutronsDepositedVsNch(0)
						  ,fHistProtonsDepositedVsNch(0)
						  ,fHistAntiProtonsDepositedVsNch(0)
						  ,fHistProtonsNotTrackMatchedDepositedVsNch(0)
						  ,fHistAntiProtonsNotTrackMatchedDepositedVsNch(0)
						  ,fHistNeutronsDepositedVsNcl(0)
						  ,fHistAntiNeutronsDepositedVsNcl(0)
						  ,fHistProtonsDepositedVsNcl(0)
						  ,fHistAntiProtonsDepositedVsNcl(0)
						  ,fHistProtonsNotTrackMatchedDepositedVsNcl(0)
						  ,fHistAntiProtonsNotTrackMatchedDepositedVsNcl(0)
						  ,fHistSecondariesVsNch(0)
						  ,fHistSecondariesVsNcl(0)
						  ,fHistSecondariesEffCorrVsNch(0)
						  ,fHistSecondariesEffCorrVsNcl(0)
						  ,fHistSecondariesOutOfAccEffCorrVsNch(0)
						  ,fHistSecondariesDetectorCoverEffCorrVsNch(0)//end
						  ,fHistNeutronsDepositedVsNchNoEffCorr(0)//start
						  ,fHistAntiNeutronsDepositedVsNchNoEffCorr(0)
						  ,fHistProtonsDepositedVsNchNoEffCorr(0)
						  ,fHistAntiProtonsDepositedVsNchNoEffCorr(0)
						  ,fHistProtonsNotTrackMatchedDepositedVsNchNoEffCorr(0)
						  ,fHistAntiProtonsNotTrackMatchedDepositedVsNchNoEffCorr(0)
						  ,fHistNeutronsDepositedVsNclNoEffCorr(0)
						  ,fHistAntiNeutronsDepositedVsNclNoEffCorr(0)
						  ,fHistProtonsDepositedVsNclNoEffCorr(0)
						  ,fHistAntiProtonsDepositedVsNclNoEffCorr(0)
						  ,fHistProtonsNotTrackMatchedDepositedVsNclNoEffCorr(0)
						  ,fHistAntiProtonsNotTrackMatchedDepositedVsNclNoEffCorr(0)//end
						  ,fHistCentVsNchVsNcl(0)
							,fHistSecondaryPositionInDetector(0)
						  ,fClusterPositionWeird(0)
						//,fHistSecondaryPositionInDetectorMultiple(0)
						  ,fSecondaryClusterEnergy(0)
						  ,fHistGammaCrossCheck(0)
						  ,fHistGammaCrossCheckAlt(0)
						  ,fHistGammaEnergyCrossCheck(0)
						  ,fHistGammaEnergyCrossCheckCent(0)
						  ,fHistGammaEnergyCrossCheckAlt(0)
    ,fHistNeutronCrossCheck(0)
    ,fHistSecondaryCrossCheck(0)
    ,fHistHadronCrossCheck(0)
    ,fHistKaonCrossCheck(0)
    ,fHistNeutronCorrection(0)
    ,fHistSecondaryCorrection(0)
    ,fHistHadronCorrection(0)
    ,fHistKaonCorrection(0)
    ,fHistAllEnergy(0)
    ,fHistSignalEnergy(0)
    ,fHistNeutronEnergy(0)
    ,fHistKaonEnergy(0)
    ,fHistHadronEnergy(0)
    ,fHistSecondaryEnergy(0)
    ,fHistSecondaryChargedEnergy(0)
    ,fHistSecondaryNeutronEnergy(0)
    ,fHistSecondaryGammaEnergy(0)
    ,fHistSecondaryElectronEnergy(0)
    ,fHistSecondaryOtherEnergy(0)
    ,fHistSimulatedGammaEnergy(0)
    ,fHistReconstructedGammaEnergy(0)
						  ,fHistSimulatedGammaEnergyAboveThreshold(0)
						  ,fHistReconstructedSignalEnergy(0)
    ,fHistFracSignalVsNClusters(0)
    ,fHistFracHadronsVsNClusters(0)
    ,fHistFracNeutronsVsNClusters(0)
    ,fHistFracKaonsVsNClusters(0)
    ,fHistFracSecondariesVsNClusters(0)
    ,fHistFracSignalVsNMultiplicity(0)
    ,fHistFracHadronsVsNMultiplicity(0)
    ,fHistFracNeutronsVsNMultiplicity(0)
    ,fHistFracKaonsVsNMultiplicity(0)
    ,fHistFracSecondariesVsNMultiplicity(0)
    ,fHistFracSignalVsNMatchedTracks(0)
    ,fHistFracHadronsVsNMatchedTracks(0)
    ,fHistFracNeutronsVsNMatchedTracks(0)
    ,fHistFracKaonsVsNMatchedTracks(0)
    ,fHistFracSecondariesVsNMatchedTracks(0)
    ,fHistFracSignalVsNTotalTracks(0)
    ,fHistFracHadronsVsNTotalTracks(0)
    ,fHistFracNeutronsVsNTotalTracks(0)
    ,fHistFracKaonsVsNTotalTracks(0)
    ,fHistFracSecondariesVsNTotalTracks(0)
    ,fHistRCorrVsPtVsCent(0)
						  ,fNMCProducedMin(0)
						  ,fNMCProducedMax(0)
{
}

// dtor
AliAnalysisEtMonteCarlo::~AliAnalysisEtMonteCarlo()
{   //Destructor

  if(fPrimaryTree){
    fPrimaryTree->Clear();
    delete fPrimaryTree;
  }
    delete fHistDecayVertexNonRemovedCharged; // Decay vertex for non-removed charged particles
    delete fHistDecayVertexRemovedCharged; // Decay vertex for non-removed charged particles
    delete fHistDecayVertexNonRemovedNeutral; // Decay vertex for non-removed charged particles
    delete fHistDecayVertexRemovedNeutral; // Decay vertex for non-removed charged particles

    delete fHistRemovedOrNot; // If charged/neutral particles were removed or not

    delete fHistEtNonRemovedProtons; // enter comment here
    delete fHistEtNonRemovedAntiProtons; // enter comment here
    delete fHistEtNonRemovedPiPlus; // enter comment here
    delete fHistEtNonRemovedPiMinus; // enter comment here
    delete fHistEtNonRemovedKaonPlus; // enter comment here
    delete fHistEtNonRemovedKaonMinus; // enter comment here
    delete fHistEtNonRemovedK0s; // enter comment here
    delete fHistEtNonRemovedK0L; // enter comment here
    delete fHistEtNonRemovedLambdas; // enter comment here
    delete fHistEtNonRemovedElectrons; // enter comment here
    delete fHistEtNonRemovedPositrons; // enter comment here
    delete fHistEtNonRemovedMuPlus; // enter comment here
    delete fHistEtNonRemovedMuMinus; // enter comment here
    delete fHistEtNonRemovedNeutrons; // enter comment here
    delete fHistEtNonRemovedAntiNeutrons; // enter comment here
    delete fHistEtNonRemovedGammas; // enter comment here
    delete fHistEtNonRemovedGammasFromPi0; // enter comment here

    delete fHistEtRemovedGammas; // enter comment here
    delete fHistEtRemovedNeutrons; // enter comment here
    delete fHistEtRemovedAntiNeutrons; // enter comment here


    delete fHistMultNonRemovedProtons; // enter comment here
    delete fHistMultNonRemovedAntiProtons; // enter comment here
    delete fHistMultNonRemovedPiPlus; // enter comment here
    delete fHistMultNonRemovedPiMinus; // enter comment here
    delete fHistMultNonRemovedKaonPlus; // enter comment here
    delete fHistMultNonRemovedKaonMinus; // enter comment here
    delete fHistMultNonRemovedK0s; // enter comment here
    delete fHistMultNonRemovedK0L; // enter comment here
    delete fHistMultNonRemovedLambdas; // enter comment here
    delete fHistMultNonRemovedElectrons; // enter comment here
    delete fHistMultNonRemovedPositrons; // enter comment here
    delete fHistMultNonRemovedMuPlus; // enter comment here
    delete fHistMultNonRemovedMuMinus; // enter comment here
    delete fHistMultNonRemovedNeutrons; // enter comment here
    delete fHistMultNonRemovedAntiNeutrons; // enter comment here
    delete fHistMultNonRemovedGammas; // enter comment here

    delete fHistMultRemovedGammas; // enter comment here
    delete fHistMultRemovedNeutrons; // enter comment here
    delete fHistMultRemovedAntiNeutrons; // enter comment here

    delete fHistTrackMultvsNonRemovedCharged; // enter comment here
    delete fHistTrackMultvsNonRemovedNeutral; // enter comment here
    delete fHistTrackMultvsRemovedGamma; // enter comment here

    delete fHistClusterMultvsNonRemovedCharged; // enter comment here
    delete fHistClusterMultvsNonRemovedNeutral; // enter comment here
    delete fHistClusterMultvsRemovedGamma; // enter comment here

    delete fHistMultvsNonRemovedChargedE; // enter comment here
    delete fHistMultvsNonRemovedNeutralE; // enter comment here
    delete fHistMultvsRemovedGammaE; // enter comment here
    delete fHistK0EDepositsVsPtInAcceptance;//enter comment here
    delete fHistK0EGammaVsPtInAcceptance;//enter comment here
    delete fHistK0EDepositsVsPtOutOfAcceptance;
    delete fHistK0EGammaVsPtOutOfAcceptance;

    delete fHistSimKaonsInAcceptance;
    delete fHistSimK0SInAcceptance;
    delete fHistSimKPlusInAcceptance;
    delete fHistSimKMinusInAcceptance;
    delete fHistSimK0LInAcceptance;
    delete fHistSimKaonsOutOfAcceptance;
    delete fHistSimKaonsInAcceptanceWithDepositsPrimaries;
    delete fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries;
    delete fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries;

    delete fHistDxDzNonRemovedCharged; // enter comment here
    delete fHistDxDzRemovedCharged; // enter comment here
    delete fHistDxDzNonRemovedNeutral; // enter comment here
    delete fHistDxDzRemovedNeutral; // enter comment here

    delete fHistPiPlusMult; // enter comment here
    delete fHistPiMinusMult; // enter comment here
    delete fHistPiZeroMult; // enter comment here

    delete fHistPiPlusMultAcc; // enter comment here
    delete fHistPiMinusMultAcc; // enter comment here
    delete fHistPiZeroMultAcc; // enter comment here
    delete fHistSimEmEtCent;
    delete fHistGammasFound; // enter comment here
    delete fHistGammasGenerated; // enter comment here
    delete fHistGammasFoundOutOfAccCent; // enter comment here
    delete fHistGammasFoundCent; // enter comment here
    delete fHistGammasFoundOutOfAccAltCent; // enter comment here
    delete fHistGammasFoundAltCent; // enter comment here
    delete fHistGammasGeneratedCent; // enter comment here
    delete fHistGammasFoundRecoEnergyCent;
    delete fHistGammasFoundRecoEnergyTrueEnergyCent;
    delete fHistAllGammasFoundRecoEnergyCent;
    delete fHistGammasFoundOutOfAccRecoEnergyCent;
    delete fHistAllGammasFoundOutOfAccRecoEnergyCent;
    delete fHistChargedTracksCut;
    delete fHistChargedTracksAccepted;
    delete fHistGammasCut;
    delete fHistGammasAccepted;
    delete fHistChargedTracksCutMult;
    delete fHistChargedTrackDepositsAcceptedVsPt;
    delete fHistChargedTrackDepositsAllVsPt;
    delete fHistChargedTrackDepositsAcceptedVsPtEffCorr;
    delete fHistChargedTrackDepositsAllVsPtEffCorr;
    delete fHistChargedTracksAcceptedMult;
    delete fHistChargedTracksAcceptedLowPtCentEffCorr;
    delete fHistChargedTracksAcceptedLowPtCent;
    delete fHistChargedTracksAcceptedLowPtCent500MeV;
    delete fHistChargedTracksAcceptedLowPtCentNoAntiProtons;
    delete fHistChargedTracksAcceptedLowPtCentAntiProtons;
    delete fHistGammasCutMult;
    delete fHistGammasAcceptedMult;
    delete fHistBadTrackMatches;
    delete fHistMatchedTracksEvspTBkgd;
    delete fHistMatchedTracksEvspTSignal;
    delete fHistMatchedTracksEvspTBkgdPeripheral;
    delete fHistMatchedTracksEvspTSignalPeripheral;
    delete fHistMatchedTracksEvspTBkgdvsCent;
    delete fHistMatchedTracksEvspTSignalvsCent;
    delete fHistMatchedTracksEvspTBkgdvsCentEffCorr;
    delete fHistMatchedTracksEvspTSignalvsCentEffCorr;
    delete fHistChargedTracksCutPeripheral;
    delete fHistChargedTracksAcceptedPeripheral;
    delete fHistGammasCutPeripheral;
    delete fHistGammasAcceptedPeripheral;
    delete fHistBadTrackMatchesdPhidEta;
    delete fHistGoodTrackMatchesdPhidEta;
    delete fHistHadronDepositsAll;
    delete fHistHadronDepositsReco;
    delete fHistHadronDepositsAllCent;
    delete fHistHadronDepositsAllCent500MeV;
    delete fHistHadronDepositsRecoCent;
    delete fHistHadronDepositsAllvsECent;
    delete fHistHadronDepositsRecovsECent;
    delete fHistHadronsAllCent;
    delete fHistMultChVsSignalVsMult;
    delete fHistNeutralRemovedSecondaryEtVsCent;
    delete fHistChargedRemovedSecondaryEtVsCent;
    delete fHistNeutralNotRemovedSecondaryEtVsCent;
    delete fHistChargedNotRemovedSecondaryEtVsCent;
    delete fHistNeutralRemovedSecondaryNumVsNCluster;
    delete fHistChargedRemovedSecondaryNumVsNCluster;
    delete fHistNeutralNotRemovedSecondaryNumVsNCluster;
    delete fHistChargedNotRemovedSecondaryNumVsNCluster;
    delete fHistNeutralRemovedSecondaryNumVsCent;
    delete fHistChargedRemovedSecondaryNumVsCent;
    delete fHistNeutralNotRemovedSecondaryNumVsCent;
    delete fHistChargedNotRemovedSecondaryNumVsCent;
    delete fHistNeutronsEtVsCent;
    delete fHistNeutronsNumVsCent;
    delete fHistNotNeutronsNumVsCent;
    delete fHistPiKPDepositedVsNch;
    delete fHistPiKPNotTrackMatchedDepositedVsNch;

    delete fHistNeutronsDepositedVsNch;
    delete fHistAntiNeutronsDepositedVsNch;
    delete fHistProtonsDepositedVsNch;
    delete fHistAntiProtonsDepositedVsNch;
    delete fHistProtonsNotTrackMatchedDepositedVsNch;
    delete fHistAntiProtonsNotTrackMatchedDepositedVsNch;
    delete fHistNeutronsDepositedVsNcl;
    delete fHistAntiNeutronsDepositedVsNcl;
    delete fHistProtonsDepositedVsNcl;
    delete fHistAntiProtonsDepositedVsNcl;
    delete fHistProtonsNotTrackMatchedDepositedVsNcl;
    delete fHistAntiProtonsNotTrackMatchedDepositedVsNcl;
    delete fHistSecondariesVsNch;
    delete fHistSecondariesVsNcl;
    delete fHistSecondariesEffCorrVsNch;
    delete fHistSecondariesEffCorrVsNcl;
    delete fHistSecondariesOutOfAccEffCorrVsNch;
    delete fHistSecondariesDetectorCoverEffCorrVsNch;


    delete fHistNeutronsDepositedVsNchNoEffCorr;
    delete fHistAntiNeutronsDepositedVsNchNoEffCorr;
    delete fHistProtonsDepositedVsNchNoEffCorr;
    delete fHistAntiProtonsDepositedVsNchNoEffCorr;
    delete fHistProtonsNotTrackMatchedDepositedVsNchNoEffCorr;
    delete fHistAntiProtonsNotTrackMatchedDepositedVsNchNoEffCorr;
    delete fHistNeutronsDepositedVsNclNoEffCorr;
    delete fHistAntiNeutronsDepositedVsNclNoEffCorr;
    delete fHistProtonsDepositedVsNclNoEffCorr;
    delete fHistAntiProtonsDepositedVsNclNoEffCorr;
    delete fHistProtonsNotTrackMatchedDepositedVsNclNoEffCorr;
    delete fHistAntiProtonsNotTrackMatchedDepositedVsNclNoEffCorr;

    delete fHistCentVsNchVsNcl;
    delete fHistSecondaryPositionInDetector;
    delete fClusterPositionWeird;
    //delete fHistSecondaryPositionInDetectorMultiple;
    delete fSecondaryClusterEnergy;
    delete fHistGammaCrossCheck;
    delete fHistGammaCrossCheckAlt;
    delete fHistGammaEnergyCrossCheck;
    delete fHistGammaEnergyCrossCheckCent;
    delete fHistGammaEnergyCrossCheckAlt;
    delete fHistNeutronCrossCheck;
    delete fHistSecondaryCrossCheck;
    delete fHistHadronCrossCheck;
    delete fHistKaonCrossCheck;
    delete fHistNeutronCorrection;
    delete fHistSecondaryCorrection;
    delete fHistHadronCorrection;
    delete fHistKaonCorrection;

    delete fHistAllEnergy;
    delete fHistSignalEnergy;
    delete fHistNeutronEnergy;
    delete fHistKaonEnergy;
    delete fHistHadronEnergy;
    delete fHistSecondaryEnergy;
    delete fHistSecondaryChargedEnergy;
    delete fHistSecondaryNeutronEnergy;
    delete fHistSecondaryGammaEnergy;
    delete fHistSecondaryElectronEnergy;
    delete fHistSecondaryOtherEnergy;
    delete fHistSimulatedGammaEnergy;
    delete fHistReconstructedGammaEnergy;
    delete fHistSimulatedGammaEnergyAboveThreshold;
    delete fHistReconstructedSignalEnergy;
    delete fHistFracSignalVsNClusters;
    delete fHistFracHadronsVsNClusters;
    delete fHistFracNeutronsVsNClusters;
    delete fHistFracKaonsVsNClusters;
    delete fHistFracSecondariesVsNClusters;
    delete fHistFracSignalVsNMultiplicity;
    delete fHistFracHadronsVsNMultiplicity;
    delete fHistFracNeutronsVsNMultiplicity;
    delete fHistFracKaonsVsNMultiplicity;
    delete fHistFracSecondariesVsNMultiplicity;
    delete fHistFracSignalVsNMatchedTracks;
    delete fHistFracHadronsVsNMatchedTracks;
    delete fHistFracNeutronsVsNMatchedTracks;
    delete fHistFracKaonsVsNMatchedTracks;
    delete fHistFracSecondariesVsNMatchedTracks;
    delete fHistFracSignalVsNTotalTracks;
    delete fHistFracHadronsVsNTotalTracks;
    delete fHistFracNeutronsVsNTotalTracks;
    delete fHistFracKaonsVsNTotalTracks;
    delete fHistFracSecondariesVsNTotalTracks;
    delete fHistRCorrVsPtVsCent;
}

Int_t AliAnalysisEtMonteCarlo::AnalyseEvent(AliVEvent* ev)
{   // analyse MC event
  if(!fIsMC) return 0;
    ResetEventValues();

    
    fPiPlusMult = 0;
    fPiMinusMult = 0;
    fPiZeroMult = 0;
    if (fCentrality)
    {
        fCentClass = fCentrality->GetCentralityClass5(fCentralityMethod);

    }

    // Get us an mc event
    if (!ev) {
        AliFatal("ERROR: Event does not exist");
        return 0;
    }
    AliMCEvent *event = dynamic_cast<AliMCEvent*>(ev);
    if (!event) {
        AliFatal("ERROR: MC Event does not exist");
        return 0;
    }

    Double_t protonMass =fgProtonMass;

    // Hijing header
    AliGenEventHeader* genHeader = event->GenEventHeader();
    if (!genHeader) {
        Printf("ERROR: Event generation header does not exist");
        return 0;
    }
    AliGenHijingEventHeader* hijingGenHeader = dynamic_cast<AliGenHijingEventHeader*>(genHeader);
    if (hijingGenHeader) {
        fImpactParameter = hijingGenHeader->ImpactParameter();
        fNcoll = hijingGenHeader->HardScatters(); // or should this be some combination of NN() NNw() NwN() NwNw() ?
        fNpart = hijingGenHeader->ProjectileParticipants() + hijingGenHeader->TargetParticipants();
    }

    // Let's play with the stack!
    AliStack *stack = event->Stack();

    Int_t nPrim = stack->GetNtrack();
    //cout<<endl<<endl<<"new event simulated nPrim "<<nPrim<<endl;

    Int_t partCount = 0;
    for (Int_t iPart = 0; iPart < nPrim; iPart++)
    {

	//Some productions have signals added in.  This switch allows the explicit exclusion of these added in signals when running over the data.
	if(checkLabelForHIJING && !IsHIJINGLabel(iPart,event,stack) ) continue;
        TParticle *part = stack->Particle(iPart);
	
	

        if (!part)
        {
	  Printf("ERROR: Could not get particle %d", iPart);
            continue;
        }
        TParticlePDG *pdg = part->GetPDG(0);

        if (!pdg)
        {
	  //Printf("ERROR: Could not get particle PDG %d", iPart);
            continue;
        }

        Double_t particleMassPart = 0; //The mass part in the Et calculation for this particle
        Int_t code = pdg->PdgCode();
	
	if(stack->IsPhysicalPrimary(iPart))
	{
	  fTotPx += part->Px();
	  fTotPy += part->Py();
	  fTotPz += part->Pz();
	}
	
        // Check for reasonable (for now neutral and singly charged) charge on the particle
        if (fSelector->IsNeutralMcParticle(iPart,*stack,*pdg))
        {

            fMultiplicity++;
             //PrintFamilyTreeShort(iPart, stack);
//
            if (TMath::Abs(part->Eta()) < fCuts->GetCommonEtaCut())
            {
                //Printf("Particle with eta: %f, pid: %d", part->Eta(), code);
                // calculate E_T
         
	      //fHistGammasFoundRecoEnergyCent->Fill(fReconstructedE,fCentClass);      
	      //cout<<"Messed up filling would have filled with "<<fReconstructedE<<" cent "<<fCentClass<<endl;
		  if (
                    TMath::Abs(code) == fgProtonCode ||
                    TMath::Abs(code) == fgNeutronCode ||
                    TMath::Abs(code) == fgLambdaCode ||
                    TMath::Abs(code) == fgXiCode ||
                    TMath::Abs(code) == fgXi0Code ||
                    TMath::Abs(code) == fgOmegaCode
                )
                {
                    if (code > 0) {
                        particleMassPart = - protonMass;
                    }
                    if (code < 0) {
                        particleMassPart = protonMass;
                    }
                }
                Double_t et = part->Energy() * TMath::Sin(part->Theta()) + particleMassPart;


                if(code == fgGammaCode || code == fgPi0Code || code == fgEtaCode || code==fgOmega0Code || code==fgEPlusCode || code==fgEMinusCode)
                {

                    if(fSelector->CutGeometricalAcceptance(*part) )
                    {
                        fNeutralMultiplicity++;
                        fTotNeutralEt += et;
                        if(fSelector->PassMinEnergyCut(*part))
                        {
                            fTotNeutralEtAfterMinEnergyCut += et;
                        }
                        if (part->Energy() > 0.05) partCount++;
                    }
                }
                //Charged particles
                else if (TMath::Abs( pdg->Charge() - fCuts->GetMonteCarloNeutralParticle())>1e-3 )
                {

                    // inside EMCal acceptance
                    if (fSelector->CutGeometricalAcceptance(*part))
                    {

                        fChargedMultiplicity++;

                        fTotChargedEt += et;

                        if (code == fgProtonCode || code == fgAntiProtonCode)
                        {
                        }
                        if (code == fgPiPlusCode || code == fgPiMinusCode)
                        {
                        }
                        if (code == fgKPlusCode || code == fgKMinusCode)
                        {
                        }
                        if (code == fgMuPlusCode || code == fgMuMinusCode)
                        {
                        }

                        if (code == fgEPlusCode || code == fgEMinusCode)
                        {
                            fTotNeutralEt += et; // calling electrons neutral
                            fTotChargedEt -= et;
                        }
                    } // inside EMCal acceptance

                    if (TrackHitsCalorimeter(part)) // magnetic field info not filled?
                    {
                        if (pdg->Charge() > 0) fHistPhivsPtPos->Fill(part->Phi(),part->Pt());
                        else fHistPhivsPtNeg->Fill(part->Phi(), part->Pt());
                    }
                }
            }

        }
    }
   // std::cout << "Total: p_x = " << fTotPx << ", p_y = " << fTotPy << ", p_z = " << fTotPz << std::endl;
    fTotEt = fTotChargedEt + fTotNeutralEt;
    //fTotEtAcc = fTotChargedEtAcc + fTotNeutralEtAcc;//
    //std::cout << "Event done! # of particles: " << partCount << std::endl;
    return 0;
}
Int_t AliAnalysisEtMonteCarlo::AnalyseEvent(AliVEvent* ev,AliVEvent* ev2)
{   // analyse MC and real event info
    //if(!mcEvent || !realEvent){
  if(!fIsMC) return 0;
    if (!ev || !ev2) {
        AliFatal("ERROR: Event does not exist");
        return 0;
    }
    AliAnalysisEt::AnalyseEvent(ev);
    AliMCEvent *mcEvent = dynamic_cast<AliMCEvent*>(ev);
    AliESDEvent *realEvent = dynamic_cast<AliESDEvent*>(ev2);
    if (!mcEvent || !realEvent) {
        AliFatal("ERROR: mcEvent or realEvent does not exist");
       
    }
    if(checkLabelForHIJING) SetGeneratorMinMaxParticles(mcEvent);
    std::vector<Int_t> foundGammas;
    
    fSelector->SetEvent(realEvent);

    AnalyseEvent(ev);

    AliStack *stack = mcEvent->Stack();

    TObjArray* list = fEsdtrackCutsTPC->GetAcceptedTracks(realEvent);
    Int_t nGoodTracks = list->GetEntries();

    //cout<<"fcuts max phi "<<fCuts->GetGeometryEmcalPhiAccMaxCut()<<endl;
    //Note that this only returns clusters for the selected detector.  fSelector actually calls the right GetClusters... for the detector
    //It does not apply any cuts on these clusters
    TRefArray *caloClusters = fSelector->GetClusters();

    Int_t nCluster = caloClusters->GetEntries();
    fClusterMult = nCluster;

    //cout<<endl<<"new event reconstructed nclusters "<<nCluster<<endl;
    fNClusters = 0;
    Int_t fClusterMultChargedTracks = 0;
    Int_t fClusterMultGammas = 0;
    Int_t nChargedSecondariesRemoved = 0;
    Int_t nChargedSecondariesNotRemoved = 0;
    Int_t nNeutralSecondariesRemoved = 0;
    Int_t nNeutralSecondariesNotRemoved = 0;
    Int_t nNeutrons = 0;
    Int_t nNotNeutrons = 0;
    //All efficiency corrected except etSecondaries -->
    Float_t etPiKPDeposited = 0.0;//
    Float_t etPiKPDepositedNotTrackMatched = 0.0;//
    Float_t etProtonDepositedNotTrackMatched = 0.0;//
    Float_t etAntiProtonDepositedNotTrackMatched = 0.0;//
    Float_t etProtonDepositedNotTrackMatchedNoEffCorr = 0.0;//
    Float_t etAntiProtonDepositedNotTrackMatchedNoEffCorr = 0.0;//
    //Float_t etPIDProtonDepositedNotTrackMatched = 0.0;//Still has to be filled!
    //1Float_t etPIDAntiProtonDepositedNotTrackMatched = 0.0;//Still has to be filled
    Float_t etProtonDeposited = 0.0;//
    Float_t etAntiProtonDeposited = 0.0;//
    Float_t etNeutronDeposited = 0.0;
    Float_t etAntiNeutronDeposited = 0.0;
    Float_t etSecondaries = 0.0;
    Float_t etSecondariesEffCorr = 0.0;
    Float_t etSecondariesOutOfAccEffCorr = 0.0;
    Float_t etSecondariesDetectorCoverEffCorr = 0.0;
    Float_t etProtonDepositedNoEffCorr = 0.0;//
    Float_t etAntiProtonDepositedNoEffCorr = 0.0;//
    Float_t etNeutronDepositedNoEffCorr = 0.0;
    Float_t etAntiNeutronDepositedNoEffCorr = 0.0;
    Float_t etGammaCrossCheck = 0.0;
    Float_t etGammaCrossCheckAlt = 0.0;
//     Float_t etNeutronCrossCheck = 0.0;
//     Float_t etSecondaryCrossCheck = 0.0;
//     Float_t etHadronCrossCheck = 0.0;
//     Float_t etKaonCrossCheck = 0.0;
    Float_t etGammaCrossCheckTrue = 0.0;//--
    Float_t etNeutronCrossCheckTrue = 0.0;//--
    Float_t etSecondaryCrossCheckTrue = 0.0;//--
    Float_t etHadronCrossCheckTrue = 0.0;//--
    Float_t etKaonCrossCheckTrue = 0.0;//
    Float_t multiplicity = fEsdtrackCutsTPC->GetReferenceMultiplicity(realEvent,kTRUE);

    Float_t subtotalAllEnergy = 0.0;//energy of all clusters passing cuts vs centrality
    Float_t subtotalSignalEnergy = 0.0;//signal of signal clusters passing cuts vs centrality
    Float_t subtotalNeutronEnergy = 0.0;//signal of neutron clusters passing cuts vs centrality
    Float_t subtotalKaonEnergy = 0.0;//signal of kaon clusters passing cuts vs centrality
    Float_t subtotalHadronEnergy = 0.0;//signal of hadron clusters passing cuts vs centrality
    Float_t subtotalSecondaryEnergy = 0.0;//signal of secondary clusters passing cuts vs centrality
    Float_t subtotalSecondaryChargedEnergy = 0.0;//signal of secondary clusters passing cuts vs centrality
    Float_t subtotalSecondaryNeutronEnergy = 0.0;//signal of secondary clusters passing cuts vs centrality
    Float_t subtotalSecondaryGammaEnergy = 0.0;//signal of secondary clusters passing cuts vs centrality
    Float_t subtotalSecondaryElectronEnergy = 0.0;//signal of secondary clusters passing cuts vs centrality
    Float_t subtotalSecondaryOtherEnergy = 0.0;//signal of secondary clusters passing cuts vs centrality
    Float_t subtotalSimulatedGammaEnergy = 0.0;//signal of signal clusters passing cuts vs centrality
    Float_t subtotalSimulatedGammaEnergyAboveThreshold = 0.0;//signal of signal clusters passing cuts vs centrality
    Float_t subtotalReconstructedGammaEnergy = 0.0;//signal of signal clusters passing cuts vs centrality
    Float_t subtotalReconstructedSignalEnergy = 0.0;//signal of signal clusters passing cuts vs centrality
    // cout<<"fsub "<<fsub<<endl;

    // loop the clusters
    for (int iCluster = 0; iCluster < nCluster; iCluster++ )
    {
        Int_t cf = 0;
        AliESDCaloCluster* caloCluster = ( AliESDCaloCluster* )caloClusters->At( iCluster );
        //Float_t caloE = caloCluster->E()
        if (!fSelector->CutGeometricalAcceptance(*caloCluster)) continue;
        fNClusters++;
        //const UInt_t iPart = (UInt_t)TMath::Abs(fSelector->GetLabel(caloCluster));//->GetLabel());
	const UInt_t iPart = fSelector->GetLabel(caloCluster,stack);
	//if(checkLabelForHIJING) cerr<<"I am checking the label"<<endl;
	//Some productions have signals added in.  This switch allows the explicit exclusion of these added in signals when running over the data.
	if(checkLabelForHIJING && !IsHIJINGLabel(iPart,mcEvent,stack) ) continue;
	//const UInt_t iPart = (UInt_t)TMath::Abs(caloCluster->GetLabel());
        TParticle *part  =  stack->Particle(iPart);

        if (!part)
        {
            Printf("No MC particle %d", iCluster);
            continue;
        }

        int primIdx = iPart;
        if (!stack->IsPhysicalPrimary(iPart)) // check whether particle is primary. we keep secondary electron and gamma for testing.
        {
            primIdx = GetPrimMother(iPart, stack);
	    if(primIdx != stack->GetPrimary(iPart))
	    {
	      //std::cout << primIdx << " = " << stack->GetPrimary(iPart) << std::endl;
	      //PrintFamilyTree(iPart, stack);
	    }
	    //if it is from a K0S
            if(primIdx < 0)
            {
	      //std::cout << "What!? No primary?" << std::endl;
	      //PrintFamilyTree(iPart, stack);
	      //continue;
	      //This is a work around to fix a bug.  For the EMCal when you use the tender supply, the parent particle ID gets messed up.
	      primIdx = iPart;
            }

        } // end of primary particle check
        //const int primCode = stack->Particle(primIdx)->GetPdgCode();
        TParticlePDG *pdg = part->GetPDG();
        if (!pdg)
        {
            Printf("ERROR: Could not get particle PDG %d", iPart);
            continue;
        }

        //Int_t code = pdg->PdgCode();
// 	if(primCode == fgGammaCode) 
// 	{
	  
   
	
	Bool_t nottrackmatched = kTRUE;//default to no track matched
	Bool_t countasmatched = kFALSE;
	Float_t matchedTrackp = 0.0;
	Float_t matchedTrackpt = 0.0;
        fDepositedCode = part->GetPdgCode();
	fReconstructedE = caloCluster->E();
        Float_t pos[3];
	//PrintFamilyTree(
        caloCluster->GetPosition(pos);
        TVector3 cp(pos);
	fReconstructedEt = caloCluster->E()*TMath::Sin(cp.Theta());
	nottrackmatched = fSelector->PassTrackMatchingCut(*caloCluster);
	//by default ALL matched tracks are accepted, whether or not the match is good.  So we check to see if the track is good.
	if(!nottrackmatched){//if the track is trackmatched
	  Int_t trackMatchedIndex = caloCluster->GetTrackMatchedIndex();
	  if(trackMatchedIndex < 0) nottrackmatched=kTRUE;
	  AliESDtrack *track = realEvent->GetTrack(trackMatchedIndex);
	  // cout<<"track code "<<fTrackDepositedCode<<" cluster code "<<fDepositedCode<<" track label "<<trackLabel<<" cluster label "<<iPart<<endl;

	  //if this is a good track, accept track will return true.  The track matched is good, so not track matched is false
	  nottrackmatched = !(fEsdtrackCutsTPC->AcceptTrack(track));//if the track is bad, this is not track matched
	  if(!nottrackmatched){//if the track that was matched is a good track
	    matchedTrackp = track->P();
	    matchedTrackpt = track->Pt();
	    UInt_t trackLabel = (UInt_t)TMath::Abs(track->GetLabel());
	    if(checkLabelForHIJING && !IsHIJINGLabel(trackLabel,mcEvent,stack) ) nottrackmatched = kTRUE;//if this was matched to something we don't want to count
	    TParticle  *trackSimPart  = stack->Particle(trackLabel);
	    Int_t fTrackDepositedCode = trackSimPart->GetPdgCode();
	    //if(!nottrackmatched) cout<<"Matched track p: "<<matchedTrackp<<" sim "<<part->P()<<endl;
	    //now we want to fill our Rcorr histos
	    Float_t rcorr = fReconstructedE/track->P();
	    fHistRCorrVsPtVsCent->Fill(rcorr,matchedTrackpt, fCentClass);
	    //cout<<"rcorr "<<rcorr<<endl;
	    Int_t n=caloCluster->GetNLabels() ;
	    //if(fReconstructedE - fsub* track->P() > 0.0){//if more energy was deposited than the momentum of the track  and more than one particle led to the cluster
	    //cout<<"was matched"<<endl;
	    if(fSelector->PassMinEnergyCut( (fReconstructedE - fsub* track->P())*TMath::Sin(cp.Theta()) ) ){//if more energy was deposited than the momentum of the track  and more than one particle led to the cluster
	      //then we say the cluster was not track matched but correct the energy
	      nottrackmatched = kTRUE;
	      //but we don't want to double count the matched tracks...
	      countasmatched = kTRUE;
	      //cout<<" fsub "<<fsub;
	      //cout<<" Reassigning energy "<<fReconstructedEt;
	      fReconstructedE = fReconstructedE - fsub* track->P();
	      fReconstructedEt = fReconstructedE*TMath::Sin(cp.Theta());
	      //cout<<" to "<<fReconstructedEt<<endl;
	      if(fDepositedCode==fTrackDepositedCode && n>1){
		//the energy deposited was more than the momentum of the track but the cluster was assigned to the particle that created the track.  We therefore need to reassign the cluster label
		//this is a sub-optimal way to re-assign the particle but it is rare that a particle needs to be reassigned
		//cout<<"Particle was "<<part->GetName()<<" reassigning label from "<<fDepositedCode<<", testing ";
		Int_t iMax=-1;
		Double_t*  Ekin=  new  Double_t[n] ;
		for(Int_t i=0;  i<n;  i++){
		  TParticle*  p=  stack->Particle(caloCluster->GetLabelAt(i)) ;
		  if(p->GetPdgCode()==fgPiPlusCode || p->GetPdgCode()==fgKPlusCode || p->GetPdgCode()==fgProtonCode || p->GetPdgCode()==fgPiMinusCode || p->GetPdgCode()==fgKMinusCode || p->GetPdgCode()==fgAntiProtonCode){ 
		    Ekin[i]=0.3;//p->P()/3.0 ;  // estimate of MIP peak, more likely to be what was deposited if the MIP peak overlapped with another peak
		  }
		  else{
		    Ekin[i]=p->Energy() ;  // what's deposited by electromagnetic particles
		  }
		  if(p->GetPdgCode()==fgAntiProtonCode){
		    Ekin[i]+=1.8  ;  //due to annihilation
		  }
		}
		Double_t eMax=0.;//eSubMax=0. ;
		//cout<<"n="<<n<<", ";
		for(Int_t i=0;  i<n;  i++){
		    Int_t label = caloCluster->GetLabelAt(i);
		    if(label!=trackMatchedIndex){
		      //TParticle*  p=  stack->Particle(caloCluster->GetLabelAt(i)) ;
		      //cout<<i<<" "<<p->GetName()<<" with E="<<Ekin[i]<<",";
		      if(Ekin[i]>eMax){
			//      eSubMax=eMax;
			eMax=Ekin[i];
			iMax=i;
		      }
		    }
		}
		delete [] Ekin;
		UInt_t correctLabel = caloCluster->GetLabelAt(iMax);
		if(iMax>0){
		  TParticle *newPart  =  stack->Particle(correctLabel);
		  if(newPart){
		    part = newPart;
		    fDepositedCode = part->GetPdgCode();
		    //cout<<", to "<<fDepositedCode<<" and particle is now "<<part->GetName();
		  }
		}
		//cout<<endl;
	      }
	      // 	      if(fDepositedCode==fgProtonCode ||  fDepositedCode==fgAntiProtonCode || fDepositedCode==fgPiPlusCode || fDepositedCode==fgPiMinusCode || fDepositedCode==fgKPlusCode || fDepositedCode==fgKMinusCode){
	      // 	for(UInt_t i = 0; i < caloCluster->GetNLabels(); i++)
	      // 	{
	      // 	  Int_t pIdx = caloCluster->GetLabelAt(i);
	      // 	const UInt_t iPart = fSelector->GetLabel(caloCluster,stack);
	      // 	//const UInt_t iPart = (UInt_t)TMath::Abs(caloCluster->GetLabel());
	      //         TParticle *part  =  stack->Particle(iPart);


	    }

	  }
	}
     
	Bool_t containsGamma = kFALSE;
	//Int_t gammaCode = -1;
	Float_t gammaEnergy = 0.0;
	for(UInt_t i = 0; i < caloCluster->GetNLabels(); i++)
	{
	  Int_t pIdx = caloCluster->GetLabelAt(i);
	  //Some productions have signals added in.  This switch allows the explicit exclusion of these added in signals when running over the data.
	  if(checkLabelForHIJING && !IsHIJINGLabel(pIdx,mcEvent,stack) ) continue;
	  //Int_t initialLabel = pIdx;
	  //TParticle *p = stack->Particle(pIdx);
	  
	  if(!stack->IsPhysicalPrimary(pIdx))
	    {//This should count gammas as reconstructed if we found even part of them
	      pIdx = GetPrimMother(pIdx, stack);
	      //TParticle *part2  =  stack->Particle(initialLabel);
// 	      TParticle *part2  =  stack->Particle(pIdx);
// 	      if(part2){
// 		TParticlePDG *pdg2 = part2->GetPDG();
// 		if(pdg2){
// 		  Int_t code2 = pdg2->PdgCode();
// 		  if(code2 == fgGammaCode){
// 		    cout<<"Relabeling primary index initial label "<<initialLabel<<" final label "<<pIdx<<endl;
// 		    PrintFamilyTree(initialLabel, stack);
// 		  }
// 		}
// 	      }
	    }
 	  if(fSelector->PassDistanceToBadChannelCut(*caloCluster))//&&fSelector->CutGeometricalAcceptance(*(stack->Particle(primIdx))))
	  {
	    //if this contained a gamma...
	    if(pIdx>0){
	      TParticle *part2  =  stack->Particle(pIdx);
	      if(part2){
		TParticlePDG *pdg2 = part2->GetPDG();
		if(pdg2){
		  Int_t code2 = pdg2->PdgCode();
		  if(code2 == fgGammaCode){//if this is counted as a reconstructed gamma, we want to know what was reconstructed about it
		    //gammaCode = pIdx;
		    containsGamma = kTRUE;
		    gammaEnergy += part2->Energy();//do += because if we have multiple gammas we want to add that to our reconstructed gamma energy, even if we didn't separat them
		  }
		}
	      }
	      //	    std::cout << "Gamma primary: " << primIdx << std::endl;
	      // 	    foundGammas.push_back(primIdx); 
	      if(nottrackmatched){
		foundGammas.push_back(pIdx); 
	      }
	    }
	  }
	}
	fCutFlow->Fill(cf++);
        if(!fSelector->PassDistanceToBadChannelCut(*caloCluster)) continue;
        Double_t clEt = CorrectForReconstructionEfficiency(*caloCluster,fReconstructedE,fCentClass);
//	if(code == fgK0SCode) std::cout << "K0 energy: " << caloCluster->E() << std::endl;
	//if(!fSelector->PassMinEnergyCut(*caloCluster)) continue;
	 	 if(!fSelector->PassMinEnergyCut(fReconstructedE)) continue;

	
        fCutFlow->Fill(cf++);

        TParticle *primPart = stack->Particle(primIdx);
        fPrimaryCode = primPart->GetPdgCode();
        if(primPart->GetPDG()) fPrimaryCharge = (Int_t) primPart->GetPDG()->Charge();

        fPrimaryE = primPart->Energy();
        fPrimaryEt = primPart->Energy()*TMath::Sin(primPart->Theta());
        fPrimaryPx = primPart->Px();
        fPrimaryPy = primPart->Py();
        fPrimaryPz = primPart->Pz();
	//cout<<"I have a cluster and it's good energy "<<caloCluster->E()<<" simulated "<<fPrimaryE<<endl;
        fPrimaryVx = primPart->Vx();
        fPrimaryVy = primPart->Vy();
        fPrimaryVz = primPart->Vz();

        fPrimaryAccepted = false;
        fPrimaryMatched = false;

	fDepositedE = part->Energy();
        fDepositedEt = part->Energy()*TMath::Sin(part->Theta());
        if(part->GetPDG()) fDepositedCharge = (Int_t) part->GetPDG()->Charge();

        fDepositedVx = part->Vx();
        fDepositedVy = part->Vy();
        fDepositedVz = part->Vz();

	if(nottrackmatched) subtotalAllEnergy += fReconstructedEt;

	//fSecondary = fSelector->FromSecondaryInteraction(*primPart, *stack);
	fSecondary =fSelector->FromSecondaryInteraction(iPart, *stack);
	//if(fSecondary && fReconstructedEt<0.3) fSecondary = kFALSE;//patch to do quick cross check THIS SHOULD NOT GET COMMITTED!!!  IF IT DID YELL AT CHRISTINE ASAP
// 	if(fSecondary) 
// 	{
// 	  //std::cout << "Have secondary!" << std::endl;
// 	  //PrintFamilyTree(iPart, stack);
// 	}
	
	pdg = primPart->GetPDG(0);
        if (!pdg)
        {
	  //Printf("ERROR: Could not get particle PDG %d", iPart);
            continue;
        }
	//Int_t code = primPart->GetPdgCode();

	Bool_t written = kFALSE;

// 	Bool_t nottrackmatched = kTRUE;//default to no track matched
// 	nottrackmatched = fSelector->PassTrackMatchingCut(*caloCluster);
// 	//by default ALL matched tracks are accepted, whether or not the match is good.  So we check to see if the track is good.
// 	if(!nottrackmatched){
// 	  Int_t trackMatchedIndex = caloCluster->GetTrackMatchedIndex();
// 	  if(trackMatchedIndex < 0) nottrackmatched=kTRUE;
// 	  AliESDtrack *track = realEvent->GetTrack(trackMatchedIndex);
// 	  //if this is a good track, accept track will return true.  The track matched is good, so not track matched is false
// 	  nottrackmatched = !(fEsdtrackCutsTPC->AcceptTrack(track));
// 	}

	if(fSecondary){//all particles from secondary interactions 
	  //===================================BEGIN SECONDARIES==================================================
	  written = kTRUE;

	  if (fDepositedCharge != 0 && fDepositedCode!=fgEMinusCode && fDepositedCode!=fgEPlusCode){//if the particle hitting the calorimeter is pi/k/p/mu
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BEGIN CHARGED SECONDARIES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    fHistChargedTrackDepositsAllVsPt->Fill(part->Pt(), fCentClass,fReconstructedEt);
	    fHistChargedTrackDepositsAllVsPtEffCorr->Fill(part->Pt(), fCentClass,clEt);
	    if(!nottrackmatched){//if it is track matched
// 	      cout<<"Charged Secondary ";
// 	      PrintFamilyTreeShort(iPart, stack);
	      fHistChargedTrackDepositsAcceptedVsPt->Fill(part->Pt(), fCentClass,fReconstructedEt);
	      fHistChargedTrackDepositsAcceptedVsPtEffCorr->Fill(part->Pt(), fCentClass,fReconstructedEt);
	      fHistHadronDepositsReco->Fill(part->Pt());
	      fHistHadronDepositsRecoCent->Fill(part->Pt(), fCentClass);
	      fHistHadronDepositsRecovsECent->Fill(fReconstructedE, fCentClass);
	      fChargedRemoved++;
	      fEnergyChargedRemoved += clEt;
	      fHistRemovedOrNot->Fill(0.0, fCentClass);
	      fHistDxDzRemovedCharged->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());
	      fHistChargedRemovedSecondaryEtVsCent->Fill(fReconstructedEt,fCentClass);
	      nChargedSecondariesRemoved++;
	      fHistChargedTracksCut->Fill(fReconstructedEt);
	      if(fCalcTrackMatchVsMult){
		fHistChargedTracksCutMult->Fill(fReconstructedEt,fClusterMult);
		if(fClusterMult<25){fHistChargedTracksCutPeripheral->Fill(fReconstructedEt);}
	      }
	      fHistMatchedTracksEvspTBkgd->Fill(matchedTrackp,fReconstructedE);
	      if(fCalcTrackMatchVsMult){
		  if(fClusterMult<25){fHistMatchedTracksEvspTBkgdPeripheral->Fill(matchedTrackp,fReconstructedEt);}
		  fHistMatchedTracksEvspTBkgdvsCent->Fill(matchedTrackp,fReconstructedEt, fCentClass);
		  fHistMatchedTracksEvspTBkgdvsCentEffCorr->Fill(matchedTrackp,clEt, fCentClass);//Fill with the efficiency corrected energy
	      }
	      //Int_t trackindex = (caloCluster->GetLabelsArray())->At(1);
	      UInt_t trackindex = fSelector->GetLabel(caloCluster,stack);//(caloCluster->GetLabelsArray())->At(1);
	      if(((UInt_t)caloCluster->GetLabel())!=trackindex){
		//if(fSelector->GetLabel(caloCluster,stack) !=trackindex){
		fHistBadTrackMatches->Fill(part->Pt(),fReconstructedE);
		fHistBadTrackMatchesdPhidEta->Fill(caloCluster->GetTrackDx(),caloCluster->GetTrackDz());
		//cout<<"Track matched, label cluster "<<caloCluster->GetLabel()<<" track "<<trackindex<<endl;
	      }
	      else{
		fHistGoodTrackMatchesdPhidEta->Fill(caloCluster->GetTrackDx(),caloCluster->GetTrackDz());
	      }
	    }
	    else{
	      fChargedNotRemoved++;
	      fEnergyChargedNotRemoved += clEt;
	      fHistRemovedOrNot->Fill(2.0, fCentClass);
	      fHistDxDzNonRemovedCharged->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());
	      fHistChargedNotRemovedSecondaryEtVsCent->Fill(fReconstructedEt,fCentClass);
	      nChargedSecondariesNotRemoved++;
	    }
	    fHistHadronDepositsAll->Fill(part->Pt());
	    fHistHadronDepositsAllCent->Fill(part->Pt(), fCentClass);
	    fHistHadronDepositsAllvsECent->Fill(fReconstructedE, fCentClass);
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END CHARGED SECONDARIES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	  }
	  else{
	  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BEGIN NEUTRAL SECONDARIES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    if(nottrackmatched  && !countasmatched){//secondaries not removed
// 	      cout<<"Neutral Secondary ";
// 	      PrintFamilyTreeShort(iPart, stack);

	      subtotalSecondaryEnergy += fReconstructedEt;
	      Bool_t tmpIsWritten = kFALSE;
	      if(fDepositedCode==fgProtonCode ||  fDepositedCode==fgAntiProtonCode || fDepositedCode==fgPiPlusCode || fDepositedCode==fgPiMinusCode || fDepositedCode==fgKPlusCode || fDepositedCode==fgKMinusCode){
		subtotalSecondaryChargedEnergy += fReconstructedEt;
		tmpIsWritten = kTRUE;
	      }
	      if(fDepositedCode==fgNeutronCode ||  fDepositedCode==fgAntiNeutronCode){
		subtotalSecondaryNeutronEnergy += fReconstructedEt;
		tmpIsWritten = kTRUE;
	      }
	      if(fDepositedCode==fgGammaCode){
		subtotalSecondaryNeutronEnergy += fReconstructedEt;
		tmpIsWritten = kTRUE;
	      }
	      if(fDepositedCode==fgEPlusCode || fDepositedCode==fgEMinusCode){
		subtotalSecondaryElectronEnergy += fReconstructedEt;
		tmpIsWritten = kTRUE;
	      }
	      if(!tmpIsWritten){
		subtotalSecondaryOtherEnergy += fReconstructedEt;
	      }
	      
	      if(!fSelector->CutGeometricalAcceptance(*part)){
		if(TMath::Sqrt(part->Vx() * part->Vx() + part->Vy() * part->Vy() + part->Vz()* part->Vz())>100){
		  etSecondariesOutOfAccEffCorr += clEt;
		}
	      }
	      else{
		if(TMath::Sqrt(part->Vx() * part->Vx() + part->Vy() * part->Vy() + part->Vz()* part->Vz())>430){
		  //clusters which are in the cover of the PHOS/EMCal
		  etSecondariesDetectorCoverEffCorr += clEt;
		}
	      }
	      fSecondaryClusterEnergy->Fill(fReconstructedEt);
 	      fHistSecondaryPositionInDetector->Fill(part->Vx(),part->Vy(),part->Vz());
	      fSecondaryNotRemoved++;
	      etSecondaries += fReconstructedEt;
	      etSecondariesEffCorr += clEt;
	      fHistNeutralNotRemovedSecondaryEtVsCent->Fill(fReconstructedEt,fCentClass);
	      nNeutralSecondariesNotRemoved++;
	    }
	    else{//secondaries removed
		fNeutralRemoved++;
		fEnergyNeutralRemoved += clEt;
		fHistRemovedOrNot->Fill(1.0, fCentClass);
		fHistDxDzRemovedNeutral->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());
		fHistNeutralRemovedSecondaryEtVsCent->Fill(fReconstructedEt,fCentClass);
		
		nNeutralSecondariesRemoved++;
	    }
	  }//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END NEUTRAL SECONDARIES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	}//===================================END SECONDARIES==================================================
	else{//not a secondary
	  //===================================BEGIN CHARGED TRACKS==================================================
	  if (fDepositedCharge != 0 && fDepositedCode!=fgEMinusCode && fDepositedCode!=fgEPlusCode){//if the particle hitting the calorimeter is pi/k/p/mu
// 	      cout<<"Hadron ";
// 	      PrintFamilyTreeShort(iPart, stack);
	    fHistChargedTrackDepositsAllVsPt->Fill(part->Pt(), fCentClass,fReconstructedEt);
	    fHistChargedTrackDepositsAllVsPtEffCorr->Fill(part->Pt(), fCentClass,clEt);

	    written = kTRUE;
	    fClusterMultChargedTracks++;
	    etPiKPDeposited += clEt;
	    if(fDepositedCode==fgProtonCode){
	      etProtonDeposited += clEt;
	      etProtonDepositedNoEffCorr += fReconstructedEt;
	    }
	    if(fDepositedCode==fgAntiProtonCode){
	      etAntiProtonDeposited += clEt;
	      etAntiProtonDepositedNoEffCorr += fReconstructedEt;
	    }
	    if(nottrackmatched && !countasmatched){//not removed but should be
	      subtotalHadronEnergy += fReconstructedEt;
	      etHadronCrossCheckTrue += clEt;
	      etPiKPDepositedNotTrackMatched += clEt;
	      if(fDepositedCode==fgProtonCode){
		etProtonDepositedNotTrackMatched += clEt;
		etProtonDepositedNotTrackMatchedNoEffCorr += fReconstructedEt;
	      }
	      if(fDepositedCode==fgAntiProtonCode){
		etAntiProtonDepositedNotTrackMatched += clEt;
		etAntiProtonDepositedNotTrackMatchedNoEffCorr += fReconstructedEt;
	      }
	      fHistHadronDepositsAll->Fill(part->Pt());
	      fHistHadronDepositsAllCent->Fill(part->Pt(), fCentClass);//denominator in track matching efficiency
	      //we want this to be all hadrons where we actually reconstructed and would accept the track

	      for (Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++){
		AliESDtrack *track = dynamic_cast<AliESDtrack*> (list->At(iTrack));
		if (!track){
		  Printf("ERROR: Could not get track %d", iTrack);
		  continue;
		}
		else{
		  if((Int_t)fSelector->GetLabel(caloCluster,stack)==track->GetLabel()){//then we found the track from the particle that created this
		    fHistHadronDepositsAllvsECent->Fill(fReconstructedE, fCentClass);
		  }
		}
	      }
	      if(fReconstructedEt>0.5) fHistHadronDepositsAllCent500MeV->Fill(part->Pt(), fCentClass);
	      fChargedNotRemoved++;
	      fEnergyChargedNotRemoved += clEt;
	      fHistRemovedOrNot->Fill(2.0, fCentClass);
	      fHistDxDzNonRemovedCharged->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());
	      fHistChargedTracksAccepted->Fill(fReconstructedEt);
	      if(fCalcTrackMatchVsMult){
		if(matchedTrackpt<0.5){//if we could never have matched this because of its pt, how much energy did it deposit?
		  fHistChargedTracksAcceptedLowPtCent->Fill(fReconstructedEt, fCentClass);
		  fHistChargedTracksAcceptedLowPtCentEffCorr->Fill(clEt, fCentClass);
		  if(fReconstructedEt>=0.5) fHistChargedTracksAcceptedLowPtCent500MeV->Fill(fReconstructedEt, fCentClass);
		  if(pdg && pdg->PdgCode()!=fgAntiProtonCode){
		    fHistChargedTracksAcceptedLowPtCentNoAntiProtons->Fill(fReconstructedEt, fCentClass);
		  }
		  else{
		    fHistChargedTracksAcceptedLowPtCentAntiProtons->Fill(fReconstructedEt, fCentClass);
		  }
		}
		fHistChargedTracksAcceptedMult->Fill(fReconstructedEt,fClusterMult);
		if(fClusterMult<25){fHistChargedTracksAcceptedPeripheral->Fill(fReconstructedEt);}
	      }
	    }
	    else{//removed and should have been
	      //if(countasmatched) cout<<" I was counted as matched even though some of my energy might have been saved."<<endl;
	      //cout<<" t.m. primary"<<endl;
	      Int_t trackindex =  fSelector->GetLabel(caloCluster,stack);// (caloCluster->GetLabelsArray())->At(0);
	      fHistChargedTrackDepositsAcceptedVsPt->Fill(part->Pt(), fCentClass,fReconstructedEt);
	      fHistChargedTrackDepositsAcceptedVsPtEffCorr->Fill(part->Pt(), fCentClass,clEt);
	      fHistHadronDepositsReco->Fill(part->Pt());
	      fHistHadronDepositsRecoCent->Fill(part->Pt(), fCentClass);
	      fHistHadronDepositsRecovsECent->Fill(fReconstructedE, fCentClass);
	      fHistHadronDepositsAll->Fill(part->Pt());
	      fHistHadronDepositsAllCent->Fill(part->Pt(), fCentClass);
	      for (Int_t iTrack = 0; iTrack < nGoodTracks; iTrack++){
		AliESDtrack *track = dynamic_cast<AliESDtrack*> (list->At(iTrack));
		if (!track){
		  Printf("ERROR: Could not get track %d", iTrack);
		  continue;
		}
		else{
		  if((Int_t) fSelector->GetLabel(caloCluster,stack)==track->GetLabel()){//then we found the track from the particle that created this
		    fHistHadronDepositsAllvsECent->Fill(fReconstructedE, fCentClass);
		  }
		}
	      }
	      if(fReconstructedEt>0.5) fHistHadronDepositsAllCent500MeV->Fill(part->Pt(), fCentClass);
	      if((Int_t)fSelector->GetLabel(caloCluster,stack)!= trackindex){
		fHistBadTrackMatches->Fill(part->Pt(),fReconstructedE);
		fHistBadTrackMatchesdPhidEta->Fill(caloCluster->GetTrackDx(),caloCluster->GetTrackDz());
		//cout<<"Track matched, label cluster "<<caloCluster->GetLabel()<<" track "<<trackindex<<endl;
	      }
	      else{
		fHistGoodTrackMatchesdPhidEta->Fill(caloCluster->GetTrackDx(),caloCluster->GetTrackDz());
	      }
	      fChargedRemoved++;
	      fEnergyChargedRemoved += clEt;
	      fHistRemovedOrNot->Fill(0.0, fCentClass);
	      fHistDxDzRemovedCharged->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());
	      fHistChargedTracksCut->Fill(fReconstructedEt);
	      if(fCalcTrackMatchVsMult){
		fHistChargedTracksCutMult->Fill(fReconstructedEt,fClusterMult);
		if(fClusterMult<25){fHistChargedTracksCutPeripheral->Fill(fReconstructedEt);}
	      }
	      fHistMatchedTracksEvspTBkgd->Fill(matchedTrackp,fReconstructedE);
	      if(fCalcTrackMatchVsMult){
		if(fClusterMult<25){fHistMatchedTracksEvspTBkgdPeripheral->Fill(matchedTrackp,fReconstructedEt);}
		fHistMatchedTracksEvspTBkgdvsCent->Fill(matchedTrackp,fReconstructedEt, fCentClass);
		fHistMatchedTracksEvspTBkgdvsCentEffCorr->Fill(matchedTrackp,clEt, fCentClass);//fill with the efficiency corrected energy
	      }
	    }
	  }
	  //===================================END CHARGED TRACKS==================================================
	  //===================================BEGIN KAONS==================================================
	  //K0L and any neutral particles from the decay of K+/- or K0S
	  if(!written && (fPrimaryCode==fgKPlusCode || fPrimaryCode==fgKMinusCode || fPrimaryCode==fgK0SCode ||fPrimaryCode==fgK0LCode) && nottrackmatched){
// 	    cout<<"Kaon ";
// 	    PrintFamilyTreeShort(iPart, stack);

	    etKaonCrossCheckTrue += clEt;
	    subtotalKaonEnergy += fReconstructedEt;
	    //etKaonCrossCheckTrue += fReconstructedEt;
	    written = kTRUE;//At this point we are not tracking them but we don't count them as neutrals accidentally removed
	  }
	  //===================================END KAONS==================================================
	  //===================================BEGIN SIGNAL==================================================
	  if(!written && (fDepositedCode==fgGammaCode || fDepositedCode==fgEMinusCode || fDepositedCode ==fgEPlusCode)){//if the particle hitting the calorimeter is gamma, electron and not from a kaon
// 	    cout<<"Signal ";
// 	    PrintFamilyTreeShort(iPart, stack);
//efficiencies without cutting track matched gammas
	    if(containsGamma){
	      if(!fSelector->CutGeometricalAcceptance(*primPart)){
		fHistAllGammasFoundOutOfAccRecoEnergyCent->Fill(fReconstructedE,fCentClass);
	      }
	      else{
		fHistAllGammasFoundRecoEnergyCent->Fill(fReconstructedE,fCentClass);
	      }
	      //cout<<"contains gamma"<<endl;
	    }
	    //else{cout<<"does not contain gamma"<<endl;}


	    fClusterMultGammas++;
	    written = kTRUE;
	    if(nottrackmatched){//Not removed and not supposed to be removed - signal
	      subtotalSignalEnergy += fReconstructedEt;
	      subtotalReconstructedSignalEnergy += fReconstructedEt;
	      if(fPrimaryCode==fgGammaCode){//if this comes from a primary electron - this is for reconstruction efficiency
		subtotalReconstructedGammaEnergy += fReconstructedEt;
	      }
// 	      else{
// 		if(fPrimaryCode!=fgEMinusCode && fPrimaryCode !=fgEPlusCode){
// 		  PrintFamilyTreeShort(iPart, stack);
// 		}
// 	      }
	      fEtNonRemovedGammas += clEt;
	      fMultNonRemovedGammas++;
	      fNeutralNotRemoved--;
	      fEnergyNeutralNotRemoved -= clEt;
	      fHistGammasAccepted->Fill(fDepositedEt);
	      if(fCalcTrackMatchVsMult){
		fHistGammasAcceptedMult->Fill(fDepositedEt,fClusterMult);
		if(fClusterMult<25){fHistGammasAcceptedPeripheral->Fill(fDepositedEt);}
	      }
	      if(containsGamma){
		if(!fSelector->CutGeometricalAcceptance(*primPart)){
		  fHistGammasFoundOutOfAccRecoEnergyCent->Fill(fReconstructedE,fCentClass);
		}
		else{
		  //cout<<"filling reference histogram"<<endl;
		  if(gammaEnergy>0) fHistGammaEnergyCrossCheckAlt->Fill(gammaEnergy,(fReconstructedEt-gammaEnergy)/fReconstructedEt);
		  fHistGammasFoundRecoEnergyCent->Fill(fReconstructedE,fCentClass);
		  fHistGammasFoundRecoEnergyTrueEnergyCent->Fill(fReconstructedE,gammaEnergy,fCentClass);
		  etGammaCrossCheckAlt += clEt;
		}
	      }
	      if(fDepositedCode==fgGammaCode){
		if(!fSelector->CutGeometricalAcceptance(*primPart)){
		  fHistGammasFoundOutOfAccCent->Fill(fReconstructedE,fCentClass);
		}
		else{
		  etGammaCrossCheck += clEt;
		  // cout<<"Found gamma et "<<fReconstructedEt<<" sim et "<< fPrimaryEt<<endl;
		  if(fPrimaryEt>0) fHistGammaEnergyCrossCheck->Fill(fPrimaryEt,(fReconstructedEt-fPrimaryEt)/fReconstructedEt);
		  if(fPrimaryEt>0) fHistGammaEnergyCrossCheckCent->Fill(fPrimaryEt,(fReconstructedEt-fPrimaryEt)/fReconstructedEt,fCentClass);
		  if((fReconstructedEt-fPrimaryEt)/fReconstructedEt>0.4 && fReconstructedEt>1.5){
		    Int_t n=caloCluster->GetNLabels() ;
		    cout<<"names ";
		    for(Int_t i=0;  i<n;  i++){
		      TParticle*  p=  stack->Particle(caloCluster->GetLabelAt(i)) ;
		      cout<<p->GetName()<<" ("<<p->Energy()<<") ";
		    }
		    cout<<endl;
		  }
		  fHistGammasFoundAltCent->Fill(fReconstructedE,fCentClass);
		}
	      }
	    }
	    else{//removed but shouldn't have been
	      Int_t trackindex = fSelector->GetLabel(caloCluster,stack);// (caloCluster->GetLabelsArray())->At(1);
	      if(caloCluster->GetLabel()!=trackindex){
		fHistBadTrackMatches->Fill(part->Pt(),fReconstructedE);
		fHistBadTrackMatchesdPhidEta->Fill(caloCluster->GetTrackDx(),caloCluster->GetTrackDz());
// 		cout<<"Track matched, label cluster "<<caloCluster->GetLabel()<<" track "<<trackindex<<endl;
// 		PrintFamilyTreeShort(trackindex, stack);
// 		cout<<"Cluster"<<endl;
	      }
	      fGammaRemoved++;
	      fGammaRemovedEt+=clEt; 
	      fHistGammasCut->Fill(fDepositedEt);
	      if(fCalcTrackMatchVsMult){
		fHistGammasCutMult->Fill(fDepositedEt,fClusterMult);
		if(fClusterMult<25){fHistGammasCutPeripheral->Fill(fDepositedEt);}
	      }
	      fHistMatchedTracksEvspTSignal->Fill(matchedTrackp,fReconstructedE);
	      if(fCalcTrackMatchVsMult){
		if(fClusterMult<25){fHistMatchedTracksEvspTSignalPeripheral->Fill(matchedTrackp,fReconstructedEt);}
		fHistMatchedTracksEvspTSignalvsCent->Fill(matchedTrackp,fReconstructedEt, fCentClass);
		fHistMatchedTracksEvspTSignalvsCentEffCorr->Fill(matchedTrackp,clEt, fCentClass);
	      }
	    }
	  }
	  //===================================END SIGNAL==================================================
	  //===================================BEGIN NEUTRONS==================================================
	  //all other cases - neutron, anti-neutron, not aware of other cases
	  if(!written){
// 	    cout<<"Neutron ";
// 	    PrintFamilyTreeShort(iPart, stack);
	    fNeutralNotRemoved++;
	    fEnergyNeutralNotRemoved += clEt;//this is the efficiency corrected energy
	    fHistRemovedOrNot->Fill(3.0, fCentClass);
	    if ((fDepositedCode == fgNeutronCode || fDepositedCode == fgAntiNeutronCode) && nottrackmatched){
	      etNeutronCrossCheckTrue += clEt;
	      subtotalNeutronEnergy += fReconstructedEt;
	      fHistDxDzNonRemovedNeutral->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());
	      fHistNeutronsEtVsCent->Fill(clEt,fCentClass);
	      nNeutrons++;
	      if(fDepositedCode == fgNeutronCode){
		etNeutronDeposited += clEt;
		etNeutronDeposited += fReconstructedEt;
	      }
	      if(fDepositedCode == fgAntiNeutronCode){
		etAntiNeutronDeposited += clEt;
		etAntiNeutronDeposited += fReconstructedEt;
	      }
	      //cout<<"I am a";
	      //if(fDepositedCode == fgAntiNeutronCode) cout<<"n anti-";
	      //else{cout<<" ";}
	      //cout<<"neutron!! pt "<<part->Pt()<<" eta "<<part->Eta()<<" phi "<<part->Phi()<<" Et "<<clEt<<endl;
	      //PrintFamilyTreeShort(iPart, stack);
	    }
	    else{
	      nNotNeutrons++;
// 	      cout<<"I am a";
// 	      if(fDepositedCode == fgAntiNeutronCode) cout<<"n anti-";
// 	      else{cout<<" ";}
// 	      cout<<"neutron!! pt "<<part->Pt()<<" eta "<<part->Eta()<<" phi "<<part->Phi()<<" Et "<<clEt<<endl;
// 	      PrintFamilyTreeShort(iPart, stack);
// 	      cout<<"I am not a neutron!!"<<endl;
// 	      PrintFamilyTreeShort(iPart, stack);
	    }
	  }
	  //===================================END NEUTRONS==================================================
	}
        if(fCuts->GetHistMakeTree()) fPrimaryTree->Fill();
    } // end of loop over clusters     


    fHistAllEnergy->Fill(fCentClass,subtotalAllEnergy);
    fHistSignalEnergy->Fill(fCentClass,subtotalSignalEnergy);
    fHistNeutronEnergy->Fill(fCentClass,subtotalNeutronEnergy);
    fHistKaonEnergy->Fill(fCentClass,subtotalKaonEnergy);
    fHistHadronEnergy->Fill(fCentClass,subtotalHadronEnergy);
    fHistSecondaryEnergy->Fill(fCentClass,subtotalSecondaryEnergy);
    fHistSecondaryChargedEnergy->Fill(fCentClass,subtotalSecondaryChargedEnergy);
    fHistSecondaryNeutronEnergy->Fill(fCentClass,subtotalSecondaryNeutronEnergy);
    fHistSecondaryGammaEnergy->Fill(fCentClass,subtotalSecondaryGammaEnergy);
    fHistSecondaryElectronEnergy->Fill(fCentClass,subtotalSecondaryElectronEnergy);
    fHistSecondaryOtherEnergy->Fill(fCentClass,subtotalSecondaryOtherEnergy);
    fHistReconstructedGammaEnergy->Fill(fCentClass,subtotalReconstructedGammaEnergy);
    fHistReconstructedSignalEnergy->Fill(fCentClass,subtotalReconstructedSignalEnergy);
    Float_t fracSignal =0;
    Float_t fracHadron = 0;
    Float_t fracKaon = 0;
    Float_t fracNeutron = 0;
    Float_t fracSecondary = 0;
    if(subtotalAllEnergy>0){
      fracSignal = subtotalReconstructedSignalEnergy/subtotalAllEnergy;
      fracHadron = subtotalHadronEnergy/subtotalAllEnergy;
      fracKaon = subtotalKaonEnergy/subtotalAllEnergy;
      fracNeutron = subtotalNeutronEnergy/subtotalAllEnergy;
      fracSecondary = subtotalSecondaryEnergy/subtotalAllEnergy;
    }
     fHistFracSignalVsNClusters->Fill(fClusterMult,fracSignal);
     fHistFracHadronsVsNClusters->Fill(fClusterMult,fracHadron);
     fHistFracKaonsVsNClusters->Fill(fClusterMult,fracKaon);
     fHistFracNeutronsVsNClusters->Fill(fClusterMult,fracNeutron);
     fHistFracSecondariesVsNClusters->Fill(fClusterMult,fracSecondary);
     fHistFracSignalVsNMultiplicity->Fill(multiplicity,fracSignal);
    fHistFracHadronsVsNMultiplicity->Fill(multiplicity,fracHadron);
    fHistFracKaonsVsNMultiplicity->Fill(multiplicity,fracKaon);
    fHistFracNeutronsVsNMultiplicity->Fill(multiplicity,fracNeutron);
    fHistFracSecondariesVsNMultiplicity->Fill(multiplicity,fracSecondary);
     fHistFracSignalVsNMatchedTracks->Fill(nChargedHadronsMeasured,fracSignal);
    fHistFracHadronsVsNMatchedTracks->Fill(nChargedHadronsMeasured,fracHadron);
    fHistFracKaonsVsNMatchedTracks->Fill(nChargedHadronsMeasured,fracKaon);
    fHistFracNeutronsVsNMatchedTracks->Fill(nChargedHadronsMeasured,fracNeutron);
    fHistFracSecondariesVsNMatchedTracks->Fill(nChargedHadronsMeasured,fracSecondary);
     fHistFracSignalVsNTotalTracks->Fill(nChargedHadronsTotal,fracSignal);
    fHistFracHadronsVsNTotalTracks->Fill(nChargedHadronsTotal,fracHadron);
    fHistFracKaonsVsNTotalTracks->Fill(nChargedHadronsTotal,fracKaon);
    fHistFracNeutronsVsNTotalTracks->Fill(nChargedHadronsTotal,fracNeutron);
    fHistFracSecondariesVsNTotalTracks->Fill(nChargedHadronsTotal,fracSecondary);



    fHistPiKPNotTrackMatchedDepositedVsNch->Fill(etPiKPDepositedNotTrackMatched,multiplicity);
    fHistPiKPDepositedVsNch->Fill(etPiKPDeposited,multiplicity);

    fHistNeutronsDepositedVsNch->Fill(etNeutronDeposited,multiplicity);
    fHistAntiNeutronsDepositedVsNch->Fill(etAntiNeutronDeposited,multiplicity);
    fHistProtonsDepositedVsNch->Fill(etProtonDeposited,multiplicity);
    fHistAntiProtonsDepositedVsNch->Fill(etAntiProtonDeposited,multiplicity);
    fHistProtonsNotTrackMatchedDepositedVsNch->Fill(etProtonDepositedNotTrackMatched,multiplicity);
    fHistAntiProtonsNotTrackMatchedDepositedVsNch->Fill(etAntiProtonDepositedNotTrackMatched,multiplicity);

    fHistNeutronsDepositedVsNcl->Fill(etNeutronDeposited,fClusterMult);
    fHistAntiNeutronsDepositedVsNcl->Fill(etAntiNeutronDeposited,fClusterMult);
    fHistProtonsDepositedVsNcl->Fill(etProtonDeposited,fClusterMult);
    fHistAntiProtonsDepositedVsNcl->Fill(etAntiProtonDeposited,fClusterMult);
    fHistProtonsNotTrackMatchedDepositedVsNcl->Fill(etProtonDepositedNotTrackMatched,fClusterMult);
    fHistAntiProtonsNotTrackMatchedDepositedVsNcl->Fill(etAntiProtonDepositedNotTrackMatched,fClusterMult);

    fHistSecondariesVsNch->Fill(etSecondaries,multiplicity);
    fHistSecondariesVsNcl->Fill(etSecondaries,fClusterMult);
    fHistSecondariesEffCorrVsNch->Fill(etSecondariesEffCorr,multiplicity);
    fHistSecondariesEffCorrVsNcl->Fill(etSecondariesEffCorr,fClusterMult);
    fHistSecondariesOutOfAccEffCorrVsNch->Fill(etSecondariesOutOfAccEffCorr,multiplicity);
    fHistSecondariesDetectorCoverEffCorrVsNch->Fill(etSecondariesDetectorCoverEffCorr,multiplicity);
    fHistCentVsNchVsNcl->Fill(fCentClass,multiplicity, fClusterMult);

    fHistNeutronsDepositedVsNchNoEffCorr->Fill(etNeutronDepositedNoEffCorr,multiplicity);
    fHistAntiNeutronsDepositedVsNchNoEffCorr->Fill(etAntiNeutronDepositedNoEffCorr,multiplicity);
    fHistNeutronsDepositedVsNclNoEffCorr->Fill(etNeutronDepositedNoEffCorr,fClusterMult);
    fHistAntiNeutronsDepositedVsNclNoEffCorr->Fill(etAntiNeutronDepositedNoEffCorr,fClusterMult);

    fHistProtonsDepositedVsNchNoEffCorr->Fill(etProtonDepositedNoEffCorr,multiplicity);
    fHistAntiProtonsDepositedVsNchNoEffCorr->Fill(etAntiProtonDepositedNoEffCorr,multiplicity);
    fHistProtonsNotTrackMatchedDepositedVsNchNoEffCorr->Fill(etProtonDepositedNotTrackMatchedNoEffCorr,multiplicity);
    fHistAntiProtonsNotTrackMatchedDepositedVsNchNoEffCorr->Fill(etAntiProtonDepositedNotTrackMatchedNoEffCorr,multiplicity);
    fHistProtonsDepositedVsNclNoEffCorr->Fill(etProtonDepositedNoEffCorr,fClusterMult);
    fHistAntiProtonsDepositedVsNclNoEffCorr->Fill(etAntiProtonDepositedNoEffCorr,fClusterMult);
    fHistProtonsNotTrackMatchedDepositedVsNclNoEffCorr->Fill(etProtonDepositedNotTrackMatchedNoEffCorr,fClusterMult);
    fHistAntiProtonsNotTrackMatchedDepositedVsNclNoEffCorr->Fill(etAntiProtonDepositedNotTrackMatchedNoEffCorr,fClusterMult);


    fHistNeutronsNumVsCent->Fill(nNeutrons,fCentClass);
    fHistNotNeutronsNumVsCent->Fill(nNotNeutrons,fCentClass);

    std::sort(foundGammas.begin(), foundGammas.end());
    for (Int_t iPart = 0; iPart < stack->GetNtrack(); iPart++)
    {

	if(!stack->IsPhysicalPrimary(iPart)) continue;
	//Some productions have signals added in.  This switch allows the explicit exclusion of these added in signals when running over the data.
	if(checkLabelForHIJING && !IsHIJINGLabel(iPart,mcEvent,stack) ) continue;
	
	TParticle *part = stack->Particle(iPart);

        if (!part)
        {
	  //Printf("ERROR: Could not get particle %d", iPart);
            continue;
        }
        TParticlePDG *pdg = part->GetPDG(0);

        if (!pdg)
        {
	  //Printf("ERROR: Could not get particle PDG %d", iPart);
            continue;
        }
        
        if(pdg->PdgCode()==fgGammaCode)// TMath::Abs(part->Eta()) < 0.12)
	{
	  if(fSelector->CutGeometricalAcceptance(*part)){
	    subtotalSimulatedGammaEnergy += part->Energy()*TMath::Sin(part->Theta());
	    if(fSelector->PassMinEnergyCut(*part)){
	      //cout<<"Sim gamma et "<< part->Energy()*TMath::Sin(part->Theta())<<endl;
	      etGammaCrossCheckTrue += part->Energy()*TMath::Sin(part->Theta());
	      subtotalSimulatedGammaEnergyAboveThreshold += part->Energy()*TMath::Sin(part->Theta());
	    }
	    fHistGammasGenerated->Fill(part->Energy());
	    fHistGammasGeneratedCent->Fill(part->Energy(),fCentClass);
	  }
	  if(std::binary_search(foundGammas.begin(),foundGammas.end(),iPart))
	  {
	    if(!fSelector->CutGeometricalAcceptance(*part)){
	      //cout<<"Gamma NOT in acceptance"<<endl;
	      fHistGammasFoundOutOfAccCent->Fill(part->Energy(),fCentClass);
	    }
	    else{
	      fHistGammasFound->Fill(part->Energy());
	      fHistGammasFoundCent->Fill(part->Energy(),fCentClass);
	      //cout<<"Gamma IN acceptance"<<endl;
	    }
	  }
	}
        if(pdg->PdgCode()==fgPiPlusCode || pdg->PdgCode()==fgPiMinusCode || pdg->PdgCode()==fgProtonCode || pdg->PdgCode()==fgAntiProtonCode){//section here for all hadrons generated
	  //if(fSelector->CutGeometricalAcceptance(*part)){
	    fHistHadronsAllCent->Fill(part->Pt(), fCentClass);
	    //}
	}
	
        
    }
    if(fCalcForKaonCorrection){
      Float_t etCuts[11] = {0.0,0.05,0.1,0.15,0.2,0.25, 0.3,0.35,0.4,0.45,0.5};
      Int_t nEtCuts = 11;
      //loop over simulated particles in order to find K0S
      for (Int_t iPart = 0; iPart < stack->GetNtrack(); iPart++){
	//Some productions have signals added in.  This switch allows the explicit exclusion of these added in signals when running over the data.
	if(checkLabelForHIJING && !IsHIJINGLabel(iPart,mcEvent,stack) ) continue;
	TParticle *part = stack->Particle(iPart);
	if (!part){
	  //Printf("ERROR: Could not get particle %d", iPart);
	  continue;
	}
	TParticlePDG *pdg = part->GetPDG(0);
	if (!pdg){
	  //Printf("ERROR: Could not get particle PDG %d", iPart);
	  continue;
	}
	//if(stack->IsPhysicalPrimary(iPart)){//if it is a K0 it might have decayed into four pions
	//fgK0SCode, fgGammaCode, fgPi0Code
	Int_t code = pdg->PdgCode();
	if(code == fgK0SCode || code==fgKPlusCode || code==fgKMinusCode ||code==fgK0LCode || code==fgK0Code){//this is a kaon
	  //cout<<"I am a kaon too! "<<stack->Particle(iPart)->GetName()<<" "<<code<<endl;
	  Float_t pTk = stack->Particle(iPart)->Pt();
	  if(TMath::Abs(stack->Particle(iPart)->Y())<0.5 && stack->IsPhysicalPrimary(iPart)){//these are particles which would be included in our spectra measurements
	    fHistSimKaonsInAcceptance->Fill(pTk);
	    if(code == fgK0SCode){fHistSimK0SInAcceptance->Fill(pTk);}
	    if(code == fgK0LCode){fHistSimK0LInAcceptance->Fill(pTk);}
	    if(code == fgKPlusCode){fHistSimKPlusInAcceptance->Fill(pTk);}
	    if(code == fgKMinusCode){fHistSimKMinusInAcceptance->Fill(pTk);}
	    if(code == fgK0Code){//Split K0's between the two
	      fHistSimK0SInAcceptance->Fill(pTk,0.5);
	      fHistSimK0LInAcceptance->Fill(pTk,0.5);
	    }
	  }
	  else{
	    fHistSimKaonsOutOfAcceptance->Fill(pTk);
	    // 	    if(!stack->IsPhysicalPrimary(iPart)){
	    // 	      PrintFamilyTree(iPart, stack);
	    // 	    }
	  }
	  Float_t totalGammaEts[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,  0.0};
	  Float_t totalClusterEts[11] = {0.0,0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0,  0.0};
	  for (int iCluster = 0; iCluster < nCluster; iCluster++ ){//if this cluster is from any of the decay daughters of any kaon...  but there is no easy way to look at this so we loop over clusters...
	    AliESDCaloCluster* caloCluster = ( AliESDCaloCluster* )caloClusters->At( iCluster );
	    if (!fSelector->CutGeometricalAcceptance(*caloCluster)) continue;
	    const Int_t myPart = fSelector->GetLabel(caloCluster,stack);
	    //const Int_t myPart = TMath::Abs(caloCluster->GetLabel());
	    //identify the primary particle which created this cluster
	    int primIdx = myPart;
	    if (!stack->IsPhysicalPrimary(myPart)){
	      primIdx = GetPrimMother(iPart, stack);
	    } // end of primary particle check
	    TParticle *hitPart = stack->Particle(myPart);
	    Bool_t hitsAsChargedKaon = kFALSE;
	    if(hitPart->GetPdgCode()== fgKPlusCode || hitPart->GetPdgCode()== fgKPlusCode){
	      if(myPart==primIdx){
		//The particle hits as a charged kaon and that kaon is a primary kaon - do not count because this is counted in the hadronic correction!
		hitsAsChargedKaon = kTRUE;
		//cout<<"Found primary charged kaon cluster!"<<endl;
	      }
	    }
	    if(primIdx==iPart && primIdx>0 && !hitsAsChargedKaon){//This cluster is from our primary particle and our primary particle is a kaon
	      //cout<<"I have a particle match! prim code"<<code<<" id "<<primIdx <<endl;
	      Float_t pos[3];
	      caloCluster->GetPosition(pos);
	      TVector3 cp(pos);
	      Double_t clEt = caloCluster->E()*TMath::Sin(cp.Theta());
	      Double_t clEtCorr = CorrectForReconstructionEfficiency(*caloCluster,fCentClass);
	      for(int l=0;l<nEtCuts;l++){//loop over cut values
		if(clEt>=etCuts[l]){
		  //cout<<", "<<clEt<<">="<<etCuts[l];
		  totalClusterEts[l] += clEtCorr;//if cluster et is above the cut off energy add it
		  totalGammaEts[l] += clEt;//if cluster et is above the cut off energy add it
		}
	      }
	    }
	  }
	  // 	  cout<<"Deposits:  pT: "<<pTk;
	  // 	  for(int l=0;l<nEtCuts;l++){//loop over cut values
	  // 	    cout<<" "<<totalClusterEts[l];
	  // 	  }
	  // 	  cout<<endl;
	  if(TMath::Abs(stack->Particle(iPart)->Y())<0.5 && stack->IsPhysicalPrimary(iPart)){//within the acceptance of our spectra and is a primary particle
	    if(totalClusterEts[0]>0.0){fHistSimKaonsInAcceptanceWithDepositsPrimaries->Fill(pTk);}
	    //cout<<"I have a particle match! prim code"<<code<<" id "<<iPart <<endl;
	    for(int l=0;l<nEtCuts;l++){
	      fHistK0EDepositsVsPtInAcceptance->Fill(pTk,totalClusterEts[l],etCuts[l]+0.001);
	      fHistK0EGammaVsPtInAcceptance->Fill(pTk,totalGammaEts[l],etCuts[l]+0.001);
	    }
	  }
	  else{//outside the acceptance of our spectra
	    if(totalClusterEts[0]>0.0){
	      if(stack->IsPhysicalPrimary(iPart)){fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries->Fill(pTk);}
	      else{fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries->Fill(pTk);}
	    }
	    for(int l=0;l<nEtCuts;l++){
	      fHistK0EDepositsVsPtOutOfAcceptance->Fill(pTk,totalClusterEts[l],etCuts[l]+0.001);
	      fHistK0EGammaVsPtOutOfAcceptance->Fill(pTk,totalGammaEts[l],etCuts[l]+0.001);
	    }
	  } 
	  
	}
      }
    }

    fHistSimulatedGammaEnergy->Fill(fCentClass,subtotalSimulatedGammaEnergy);
    fHistSimulatedGammaEnergyAboveThreshold->Fill(fCentClass,subtotalSimulatedGammaEnergyAboveThreshold);

    fHistMultChVsSignalVsMult->Fill(fClusterMultChargedTracks,fClusterMultGammas,fClusterMult);
    fHistNeutralRemovedSecondaryNumVsNCluster->Fill(nNeutralSecondariesRemoved,fClusterMult);
    fHistChargedRemovedSecondaryNumVsNCluster->Fill(nChargedSecondariesRemoved,fClusterMult);
    fHistNeutralNotRemovedSecondaryNumVsNCluster->Fill(nNeutralSecondariesNotRemoved,fClusterMult);
    fHistChargedNotRemovedSecondaryNumVsNCluster->Fill(nChargedSecondariesNotRemoved,fClusterMult);
    fHistNeutralRemovedSecondaryNumVsCent->Fill(nNeutralSecondariesRemoved,fCentClass);
    fHistChargedRemovedSecondaryNumVsCent->Fill(nChargedSecondariesRemoved,fCentClass);
    fHistNeutralNotRemovedSecondaryNumVsCent->Fill(nNeutralSecondariesNotRemoved,fCentClass);
    fHistChargedNotRemovedSecondaryNumVsCent->Fill(nChargedSecondariesNotRemoved,fCentClass);
    if(etGammaCrossCheckTrue>0)fHistGammaCrossCheck->Fill(fCentClass,(etGammaCrossCheck-etGammaCrossCheckTrue)/etGammaCrossCheckTrue);
    if(etGammaCrossCheckTrue>0)fHistGammaCrossCheckAlt->Fill(fCentClass,(etGammaCrossCheckAlt-etGammaCrossCheckTrue)/etGammaCrossCheckTrue);
    if(fTmCorrections->GetNeutronCorrection(fCentClass)>0)fHistNeutronCrossCheck->Fill(fCentClass,(fTmCorrections->GetNeutronCorrection(fCentClass)-etNeutronCrossCheckTrue)/fTmCorrections->GetNeutronCorrection(fCentClass));
    if(fTmCorrections->GetSecondaryCorrection(fCentClass)>0)fHistSecondaryCrossCheck->Fill(fCentClass,(fTmCorrections->GetSecondaryCorrection(fCentClass)-etSecondaryCrossCheckTrue)/fTmCorrections->GetSecondaryCorrection(fCentClass));
    if(fTmCorrections->GetHadronCorrection(fCentClass)>0)fHistHadronCrossCheck->Fill(fCentClass,(fTmCorrections->GetHadronCorrection(fCentClass)-etHadronCrossCheckTrue)/fTmCorrections->GetHadronCorrection(fCentClass));
    if(fTmCorrections->GetKaonCorrection(fCentClass)>0)fHistKaonCrossCheck->Fill(fCentClass,(fTmCorrections->GetKaonCorrection(fCentClass)-etKaonCrossCheckTrue)/fTmCorrections->GetKaonCorrection(fCentClass));

    fHistSimEmEtCent->Fill(fTotNeutralEt,fCentClass);
    fHistNeutronCorrection->Fill(fCentClass,etNeutronCrossCheckTrue);
    fHistSecondaryCorrection->Fill(fCentClass,etSecondaryCrossCheckTrue);
    fHistHadronCorrection->Fill(fCentClass,etHadronCrossCheckTrue);
    fHistKaonCorrection->Fill(fCentClass,etKaonCrossCheckTrue);
    //cout<<"Secondaries not removed: "<<fSecondaryNotRemoved<<" neutral secondaries not removed "<<nNeutralSecondariesNotRemoved<<" cluster mult "<<fClusterMult<<endl;
    FillHistograms();
    return 0;
}

void AliAnalysisEtMonteCarlo::Init()
{   // init
    AliAnalysisEt::Init();
}

void AliAnalysisEtMonteCarlo::ResetEventValues()
{   // reset event values
    AliAnalysisEt::ResetEventValues();
  if(!fIsMC) return;

    fTotEtSecondary = 0;
    fTotEtSecondaryFromEmEtPrimary = 0;
    fTotEtWithSecondaryRemoved = 0;

    // collision geometry defaults for p+p:
    fImpactParameter = 0;
    fNcoll = 1;
    fNpart = 2;

    fEtNonRemovedProtons = 0;
    fEtNonRemovedAntiProtons = 0;
    fEtNonRemovedPiPlus = 0;
    fEtNonRemovedPiMinus = 0;
    fEtNonRemovedKaonPlus = 0;
    fEtNonRemovedKaonMinus = 0;
    fEtNonRemovedK0S = 0;
    fEtNonRemovedK0L = 0;
    fEtNonRemovedLambdas = 0;
    fEtNonRemovedElectrons = 0;
    fEtNonRemovedPositrons = 0;
    fEtNonRemovedMuPlus = 0;
    fEtNonRemovedMuMinus = 0;
    fEtNonRemovedNeutrons = 0;
    fEtNonRemovedAntiNeutrons = 0;
    fEtNonRemovedGammas = 0;
    fEtNonRemovedGammasFromPi0 = 0;

    fEtRemovedProtons = 0;
    fEtRemovedAntiProtons = 0;
    fEtRemovedPiPlus = 0;
    fEtRemovedPiMinus = 0;
    fEtRemovedKaonPlus = 0;
    fEtRemovedKaonMinus = 0;
    fEtRemovedK0s = 0;
    fEtRemovedK0L = 0;
    fEtRemovedLambdas = 0;
    fEtRemovedElectrons = 0;
    fEtRemovedPositrons = 0;
    fEtRemovedMuPlus = 0;
    fEtRemovedMuMinus = 0;
    fEtRemovedNeutrons = 0;

    fEtRemovedGammasFromPi0 = 0;
    fEtRemovedGammas = 0;
    fEtRemovedNeutrons = 0;
    fEtRemovedAntiNeutrons = 0;

    fMultNonRemovedProtons = 0;
    fMultNonRemovedAntiProtons = 0;
    fMultNonRemovedPiPlus = 0;
    fMultNonRemovedPiMinus = 0;
    fMultNonRemovedKaonPlus = 0;
    fMultNonRemovedKaonMinus = 0;
    fMultNonRemovedK0s = 0;
    fMultNonRemovedK0L = 0;
    fMultNonRemovedLambdas = 0;
    fMultNonRemovedElectrons = 0;
    fMultNonRemovedPositrons = 0;
    fMultNonRemovedMuPlus = 0;
    fMultNonRemovedMuMinus = 0;
    fMultNonRemovedNeutrons = 0;
    fMultNonRemovedAntiNeutrons = 0;
    fMultNonRemovedGammas = 0;

    fMultRemovedProtons = 0;
    fMultRemovedAntiProtons = 0;
    fMultRemovedPiPlus = 0;
    fMultRemovedPiMinus = 0;
    fMultRemovedKaonPlus = 0;
    fMultRemovedKaonMinus = 0;
    fMultRemovedK0s = 0;
    fMultRemovedK0L = 0;
    fMultRemovedLambdas = 0;
    fMultRemovedElectrons = 0;
    fMultRemovedPositrons = 0;
    fMultRemovedMuPlus = 0;
    fMultRemovedMuMinus = 0;

    fMultRemovedGammas = 0;
    fMultRemovedNeutrons = 0;
    fMultRemovedAntiNeutrons = 0;

    fEnergyChargedNotRemoved = 0;
    fEnergyChargedRemoved = 0;
    fEnergyNeutralNotRemoved = 0;
    fEnergyNeutralRemoved = 0;

    fChargedNotRemoved = 0;
    fChargedRemoved = 0;
    fNeutralNotRemoved = 0;
    fNeutralRemoved = 0;
    fGammaRemoved = 0;
    fSecondaryNotRemoved = 0;

    fTrackMultInAcc = 0;

    fTotNeutralEtAfterMinEnergyCut = 0;
  
    fSecondaryNotRemoved = 0;
    
    fTotPx = 0;
    fTotPy = 0;
    fTotPz = 0;
    
    
}

void AliAnalysisEtMonteCarlo::CreateHistograms()
{   // histogram related additions
    AliAnalysisEt::CreateHistograms();
    if(!fIsMC) return;
    if (fEventSummaryTree) {
        fEventSummaryTree->Branch("fImpactParameter",&fImpactParameter,"fImpactParameter/D");
        fEventSummaryTree->Branch("fNcoll",&fNcoll,"fNcoll/I");
        fEventSummaryTree->Branch("fNpart",&fNpart,"fNpart/I");
        fEventSummaryTree->Branch("fTotEtWithSecondaryRemoved", &fTotEtWithSecondaryRemoved, "fTotEtWithSecondaryRemoved/D");
        fEventSummaryTree->Branch("fTotEtSecondaryFromEmEtPrimary", &fTotEtSecondaryFromEmEtPrimary, "fTotEtSecondaryFromEmEtPrimary/D");
        fEventSummaryTree->Branch("fTotEtSecondary", &fTotEtSecondary, "fTotEtSecondary/D");
        fEventSummaryTree->Branch("fTotNeutralEtAfterMinEnergyCut", &fTotNeutralEtAfterMinEnergyCut, "fTotNeutralEtAfterMinEnergyCut/D");
        fEventSummaryTree->Branch("fSecondaryNotRemoved", &fSecondaryNotRemoved, "fSecondaryNotRemoved/I");
	fEventSummaryTree->Branch("fChargedNotRemoved", &fChargedNotRemoved, "fChargedNotRemoved/I");
	fEventSummaryTree->Branch("fNeutralNotRemoved", &fNeutralNotRemoved, "fNeutralNotRemoved/I");
	fEventSummaryTree->Branch("fChargedRemoved", &fChargedRemoved, "fChargedRemoved/I");
	fEventSummaryTree->Branch("fNeutralRemoved", &fNeutralRemoved, "fNeutralRemoved/I");
	fEventSummaryTree->Branch("fGammaRemoved", &fGammaRemoved, "fGammaRemoved/I");
	fEventSummaryTree->Branch("fTotPx", &fTotPx, "fTotPx/D");
	fEventSummaryTree->Branch("fTotPy", &fTotPy, "fTotPy/D");
	fEventSummaryTree->Branch("fTotPz", &fTotPz, "fTotPz/D");
// 	fEventSummaryTree->Branch("f
    }
    Float_t scale = 1;//scale up histograms if EMCal 2011 so we have the right range
    if(fDataSet==2011   && !fHistogramNameSuffix.Contains("P")){
      scale = 3.0;
    }


    //fHistDecayVertexNonRemovedCharged = new TH3F("fHistDecayVertexNonRemovedCharged","fHistDecayVertexNonRemovedCharged", 500, -470, 30, 500, -300, 300, 40, -20, 20);
    //fHistDecayVertexRemovedCharged = new TH3F("fHistDecayVertexRemovedCharged","fHistDecayVertexRemovedCharged", 500, -470, 30, 500, -300, 300, 40, -20, 20);
    //fHistDecayVertexNonRemovedNeutral = new TH3F("fHistDecayVertexNonRemovedNeutral","fHistDecayVertexNonRemovedNeutral", 500, -470, 30, 500, -300, 300, 40, -20, 20);
    //fHistDecayVertexRemovedNeutral = new TH3F("fHistDecayVertexRemovedNeutral","fHistDecayVertexRemovedNeutral", 500, -470, 30, 500, -300, 300, 40, -20, 20);

    fHistRemovedOrNot = new TH2F("fHistRemovedOrNot", "fHistRemovedOrNot", 4, -0.5, 3.5, 10, -0.5, 9.5);

    //I don't think any of these are actually used...  flag for later removal
    Int_t nbins = 15;
    fHistEtNonRemovedProtons = new TH2F("fHistEtNonRemovedProtons", "fHistEtNonRemovedProtons", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedAntiProtons = new TH2F("fHistEtNonRemovedAntiProtons", "fHistEtNonRemovedAntiProtons", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedPiPlus = new TH2F("fHistEtNonRemovedPiPlus", "fHistEtNonRemovedPiPlus", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedPiMinus = new TH2F("fHistEtNonRemovedPiMinus", "fHistEtNonRemovedPiMinus", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedKaonPlus = new TH2F("fHistEtNonRemovedKaonPlus", "fHistEtNonRemovedKaonPlus", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedKaonMinus = new TH2F("fHistEtNonRemovedKaonMinus", "fHistEtNonRemovedKaonMinus", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedK0s = new TH2F("fHistEtNonRemovedK0s", "fHistEtNonRemovedK0s", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedK0L = new TH2F("fHistEtNonRemovedK0L", "fHistEtNonRemovedK0L", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedLambdas = new TH2F("fHistEtNonRemovedLambdas", "fHistEtNonRemovedLambdas", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedElectrons = new TH2F("fHistEtNonRemovedElectrons", "fHistEtNonRemovedElectrons", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedPositrons = new TH2F("fHistEtNonRemovedPositrons", "fHistEtNonRemovedPositrons", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedMuPlus = new TH2F("fHistEtNonRemovedMuPlus", "fHistEtNonRemovedMuPlus", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedMuMinus = new TH2F("fHistEtNonRemovedMuMinus", "fHistEtNonRemovedMuMinus", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedNeutrons = new TH2F("fHistEtNonRemovedNeutrons", "fHistEtNonRemovedNeutrons", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedAntiNeutrons = new TH2F("fHistEtNonRemovedAntiNeutrons", "fHistEtNonRemovedAntiNeutrons", nbins, 0, 30, 10, -0.5, 9.5);

    fHistEtNonRemovedGammas = new  TH2F("fHistEtNonRemovedGammas", "fHistEtNonRemovedGammas", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedGammasFromPi0 = new  TH2F("fHistEtNonRemovedGammasFromPi0", "fHistEtNonRemovedGammasFromPi0", nbins, 0, 30, 10, -0.5, 9.5);

    fHistEtRemovedGammas = new  TH2F("fHistEtRemovedGammas", "fHistEtRemovedGammas", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtRemovedNeutrons = new  TH2F("fHistEtRemovedNeutrons", "fHistEtRemovedNeutrons", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtRemovedAntiNeutrons = new  TH2F("fHistEtRemovedAntiNeutrons", "fHistEtRemovedAntiNeutrons", nbins, 0, 30, 10, -0.5, 9.5);

    fHistEtRemovedCharged = new  TH2F("fHistEtRemovedCharged", "fHistEtRemovedCharged", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtRemovedNeutrals = new  TH2F("fHistEtRemovedNeutrals", "fHistEtRemovedNeutrals", nbins, 0, 30, 10, -0.5, 9.5);

    fHistEtNonRemovedCharged = new  TH2F("fHistEtNonRemovedCharged", "fHistEtNonRemovedCharged", nbins, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedNeutrals = new  TH2F("fHistEtNonRemovedNeutrals", "fHistEtNonRemovedNeutrals", nbins, 0, 30, 10, -0.5, 9.5);

    Int_t nmultbins = 10;
    fHistMultNonRemovedProtons = new TH2F("fHistMultNonRemovedProtons", "fHistMultNonRemovedProtons", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedAntiProtons = new TH2F("fHistMultNonRemovedAntiProtons", "fHistMultNonRemovedAntiProtons", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedPiPlus = new TH2F("fHistMultNonRemovedPiPlus", "fHistMultNonRemovedPiPlus", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedPiMinus = new TH2F("fHistMultNonRemovedPiMinus", "fHistMultNonRemovedPiMinus", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedKaonPlus = new TH2F("fHistMultNonRemovedKaonPlus", "fHistMultNonRemovedKaonPlus", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedKaonMinus = new TH2F("fHistMultNonRemovedKaonMinus", "fHistMultNonRemovedKaonMinus", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedK0s = new TH2F("fHistMultNonRemovedK0s", "fHistMultNonRemovedK0s", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedK0L = new TH2F("fHistMultNonRemovedK0L", "fHistMultNonRemovedK0L", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedLambdas = new TH2F("fHistMultNonRemovedLambdas", "fHistMultNonRemovedLambdas", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedElectrons = new TH2F("fHistMultNonRemovedElectrons", "fHistMultNonRemovedElectrons", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedPositrons = new TH2F("fHistMultNonRemovedPositrons", "fHistMultNonRemovedPositrons", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedMuPlus = new TH2F("fHistMultNonRemovedMuPlus", "fHistMultNonRemovedMuPlus", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedMuMinus = new TH2F("fHistMultNonRemovedMuMinus", "fHistMultNonRemovedMuMinus", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedNeutrons = new TH2F("fHistMultNonRemovedNeutrons", "fHistMultNonRemovedNeutrons", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedAntiNeutrons = new TH2F("fHistMultNonRemovedAntiNeutrons", "fHistMultNonRemovedAntiNeutrons", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);

    fHistMultNonRemovedGammas = new  TH2F("fHistMultNonRemovedGammas", "fHistMultNonRemovedGammas", nmultbins, -0.5, 99.5, nmultbins, -0.5, 99.5);

    fHistMultRemovedGammas = new  TH2F("fHistMultRemovedGammas", "fHistMultRemovedGammas", nmultbins, -0.5, 99.5, nmultbins, -0.5, 99.5);
    fHistMultRemovedNeutrons = new  TH2F("fHistMultRemovedNeutrons", "fHistMultRemovedNeutrons", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultRemovedAntiNeutrons = new  TH2F("fHistMultRemovedAntiNeutrons", "fHistMultRemovedAntiNeutrons", nmultbins, -0.5, 99.5, 10, -0.5, 9.5);
    /*
      fHistMultRemovedCharged = new  TH2F("fHistMultRemovedCharged", "fHistMultRemovedCharged", nbins, 0, 30, 10, -0.5, 9.5);
      fHistMultRemovedNeutrals = new  TH2F("fHistMultRemovedNeutrals", "fHistMultRemovedNeutrals", nbins, 0, 30, 10, -0.5, 9.5);

      fHistMultNonRemovedCharged = new  TH2F("fHistMultNonRemovedCharged", "fHistMultNonRemovedCharged", nbins, 0, 30, 10, -0.5, 9.5);
      fHistMultNonRemovedNeutrals = new  TH2F("fHistMultNonRemovedNeutrals", "fHistMultNonRemovedNeutrals", nbins, 0, 30, 10, -0.5, 9.5);*/


    fHistMultRemovedCharged = new  TH2F("fHistMultRemovedCharged", "fHistMultRemovedCharged", nmultbins, -0.5, 99.5, nmultbins, -0.5, 99.5);
    fHistMultRemovedNeutrals = new  TH2F("fHistMultRemovedNeutrals", "fHistMultRemovedNeutrals", nmultbins, -0.5, 99.5, nmultbins, -0.5, 99.5);

    fHistMultNonRemovedCharged = new  TH2F("fHistMultNonRemovedCharged", "fHistMultNonRemovedCharged", nmultbins, -0.5, 99.5, nmultbins, -0.5, 99.5);
    fHistMultNonRemovedNeutrals = new  TH2F("fHistMultNonRemovedNeutrals", "fHistMultNonRemovedNeutrals", nmultbins, -0.5, 99.5, nmultbins, -0.5, 99.5);

    fHistTrackMultvsNonRemovedCharged = new TH2F("fHistTrackMultvsNonRemovedCharged", "fHistTrackMultvsNonRemovedCharged", nmultbins*10, -0.5, 999.5, nmultbins, -0.5, 99.5);
    fHistTrackMultvsNonRemovedNeutral = new TH2F("fHistTrackMultvsNonRemovedNeutral", "fHistTrackMultvsNonRemovedNeutral", nmultbins*10, -0.5, 999.5, nmultbins, -0.5, 99.5);
    fHistTrackMultvsRemovedGamma = new TH2F("fHistTrackMultvsRemovedGamma", "fHistTrackMultvsRemovedGamma", nmultbins*10, -0.5, 999.5, 100, -0.5, 99.5);

    fHistClusterMultvsNonRemovedCharged = new TH2F("fHistClusterMultvsNonRemovedCharged", "fHistClusterMultvsNonRemovedCharged", nmultbins*10, -0.5, 999.5, nmultbins, -0.5, 99.5);
    fHistClusterMultvsNonRemovedNeutral = new TH2F("fHistClusterMultvsNonRemovedNeutral", "fHistClusterMultvsNonRemovedNeutral", nmultbins*10, -0.5, 999.5, nmultbins, -0.5, 99.5);
    fHistClusterMultvsRemovedGamma = new TH2F("fHistClusterMultvsRemovedGamma", "fHistClusterMultvsRemovedGamma", nmultbins*10, -0.5, 999.5, nmultbins, -0.5, 99.5);

    Int_t ndnbins = 8;
    fHistDxDzNonRemovedCharged = new TH2F("fHistDxDzNonRemovedCharged", "fHistDxDzNonRemovedCharged", ndnbins, -200, 200, ndnbins, -200, 200);
    fHistDxDzRemovedCharged = new TH2F("fHistDxDzRemovedCharged", "fHistDxDzRemovedCharged", ndnbins, -200, 200, ndnbins, -200, 200);
    fHistDxDzNonRemovedNeutral = new TH2F("fHistDxDzNonRemovedNeutral", "fHistDxDzNonRemovedNeutral", ndnbins, -200, 200, ndnbins, -200, 200);
    fHistDxDzRemovedNeutral = new TH2F("fHistDxDzRemovedNeutral", "fHistDxDzRemovedNeutral", ndnbins, -200, 200, ndnbins, -200, 200);

    if(fCalcForKaonCorrection){
      Int_t nEtCut = 11;
      Float_t etCutAxis[12] = {0.00,0.05,0.10,0.15,0.20,  0.25,0.30,0.35,0.40,0.45, 0.50,0.55};
      fHistK0EDepositsVsPtInAcceptance = new TH3F("fHistK0EDepositsVsPtInAcceptance","Kaon deposits with corrections for kaons with y<0.5",fgNumOfPtBins,fgPtAxis,fgNumOfPtBins,fgPtAxis,nEtCut,etCutAxis);
      fHistK0EGammaVsPtInAcceptance = new TH3F("fHistK0EGammaVsPtInAcceptance","Kaon deposits without corrections for kaons with y<0.5",fgNumOfPtBins,fgPtAxis,fgNumOfPtBins,fgPtAxis,nEtCut,etCutAxis);
      fHistK0EDepositsVsPtOutOfAcceptance = new TH3F("fHistK0EDepositsVsPtOutOfAcceptance","Kaon deposits with corrections for kaons with y>0.5",fgNumOfPtBins,fgPtAxis,fgNumOfPtBins,fgPtAxis,nEtCut,etCutAxis);
      fHistK0EGammaVsPtOutOfAcceptance = new TH3F("fHistK0EGammaVsPtOutOfAcceptance","Kaon deposits without corrections for kaons with y>0.5",fgNumOfPtBins,fgPtAxis,fgNumOfPtBins,fgPtAxis,nEtCut,etCutAxis);
      fHistK0EDepositsVsPtInAcceptance->GetXaxis()->SetTitle("Kaon p_{T}");
      fHistK0EGammaVsPtInAcceptance->GetXaxis()->SetTitle("Kaon p_{T}");
      fHistK0EDepositsVsPtOutOfAcceptance->GetXaxis()->SetTitle("Kaon p_{T}");
      fHistK0EGammaVsPtOutOfAcceptance->GetXaxis()->SetTitle("Kaon p_{T}");
      fHistK0EDepositsVsPtInAcceptance->GetYaxis()->SetTitle("Deposited E_{T}");
      fHistK0EGammaVsPtInAcceptance->GetYaxis()->SetTitle("Deposited E_{T}");
      fHistK0EDepositsVsPtOutOfAcceptance->GetYaxis()->SetTitle("Deposited E_{T}");
      fHistK0EGammaVsPtOutOfAcceptance->GetYaxis()->SetTitle("Deposited E_{T}");
      fHistK0EDepositsVsPtInAcceptance->GetZaxis()->SetTitle("E_{T} cut");
      fHistK0EGammaVsPtInAcceptance->GetZaxis()->SetTitle("E_{T} cut");
      fHistK0EDepositsVsPtOutOfAcceptance->GetZaxis()->SetTitle("E_{T} cut");
      fHistK0EGammaVsPtOutOfAcceptance->GetZaxis()->SetTitle("E_{T} cut");

      fHistSimKaonsInAcceptance = new TH1F("fHistSimKaonsInAcceptance","Kaons with y<0.5",fgNumOfPtBins,fgPtAxis);
      fHistSimK0SInAcceptance = new TH1F("fHistSimK0SInAcceptance","Kaons with y<0.5",fgNumOfPtBins,fgPtAxis);
      fHistSimK0LInAcceptance = new TH1F("fHistSimK0LInAcceptance","Kaons with y<0.5",fgNumOfPtBins,fgPtAxis);
      fHistSimKPlusInAcceptance = new TH1F("fHistSimKPlusInAcceptance","Kaons with y<0.5",fgNumOfPtBins,fgPtAxis);
      fHistSimKMinusInAcceptance = new TH1F("fHistSimKMinusInAcceptance","Kaons with y<0.5",fgNumOfPtBins,fgPtAxis);
      fHistSimKaonsOutOfAcceptance = new TH1F("fHistSimKaonsOutOfAcceptance","Kaons with y>0.5",fgNumOfPtBins,fgPtAxis);
      fHistSimKaonsInAcceptanceWithDepositsPrimaries = new TH1F("fHistSimKaonsInAcceptanceWithDepositsPrimaries","Primary Kaons which deposited energy in calorimeter with y>0.5",fgNumOfPtBins,fgPtAxis);
      fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries = new TH1F("fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries","Secondary Kaons which deposited energy in calorimeter with y>0.5",fgNumOfPtBins,fgPtAxis);
      fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries = new TH1F("fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries","Primary Kaons which deposited energy in calorimeter with y>0.5",fgNumOfPtBins,fgPtAxis);
    }

    fHistPiPlusMult = new TH1F("fHistPiPlusMult", "fHistPiPlusMult", nmultbins, -0.5, 1999.5);
    fHistPiMinusMult = new TH1F("fHistPiMinusMult", "fHistPiMinusMult", nmultbins, -0.5, 1999.5);
    fHistPiZeroMult = new TH1F("fHistPiZeroMult", "fHistPiZeroMult", nmultbins, -0.5, 1999.5);

    fHistPiPlusMultAcc = new TH1F("fHistPiPlusMultAcc", "fHistPiPlusMultAcc", nmultbins, -0.5, 1999.5);
    fHistPiMinusMultAcc = new TH1F("fHistPiMinusMultAcc", "fHistPiMinusMultAcc", nmultbins, -0.5, 1999.5);
    fHistPiZeroMultAcc = new TH1F("fHistPiZeroMultAcc", "fHistPiZeroMultAcc", nmultbins, -0.5, 1999.5);

    if(fCuts->GetHistMakeTree())
    {
        TString treename = "fPrimaryTree" + fHistogramNameSuffix;
        fPrimaryTree = new TTree(treename, treename);

        fPrimaryTree->Branch("fTotEt",&fTotEt,"fTotEt/D");
        fPrimaryTree->Branch("fNeutralMultiplicity", &fNeutralMultiplicity, "fNeutralMultiplicity/I");
        fPrimaryTree->Branch("fCentClass",&fCentClass,"fCentClass/I");

        fPrimaryTree->Branch("fPrimaryCode", &fPrimaryCode, "fPrimaryCode/I");
        fPrimaryTree->Branch("fPrimaryCharge", &fPrimaryCharge, "fPrimaryCharge/I");

        fPrimaryTree->Branch("fPrimaryE", &fPrimaryE, "fPrimaryE/D");
        fPrimaryTree->Branch("fPrimaryEt", &fPrimaryEt, "fPrimaryEt/D");

        fPrimaryTree->Branch("fPrimaryPx", &fPrimaryPx, "fPrimaryPx/D");
        fPrimaryTree->Branch("fPrimaryPy", &fPrimaryPy, "fPrimaryPy/D");
        fPrimaryTree->Branch("fPrimaryPz", &fPrimaryPz, "fPrimaryPz/D");

        fPrimaryTree->Branch("fPrimaryVx", &fPrimaryVx, "fPrimaryVx/D");
        fPrimaryTree->Branch("fPrimaryVy", &fPrimaryVy, "fPrimaryVy/D");
        fPrimaryTree->Branch("fPrimaryVz", &fPrimaryVz, "fPrimaryVz/D");

        fPrimaryTree->Branch("fPrimaryAccepted", &fPrimaryAccepted, "fPrimaryAccepted/B");
        fPrimaryTree->Branch("fPrimaryMatched", &fPrimaryMatched, "fPrimaryMatched/B");


        fPrimaryTree->Branch("fDepositedCode", &fDepositedCode, "fDepositedCode/I");
        fPrimaryTree->Branch("fDepositedCharge", &fDepositedCharge, "fDepositedCharge/I");
	fPrimaryTree->Branch("fDepositedE", &fDepositedE, "fDepositedE/D");
        fPrimaryTree->Branch("fDepositedEt", &fDepositedEt, "fDepositedEt/D");

        fPrimaryTree->Branch("fDepositedVx", &fDepositedVx, "fDepositedVx/D");
        fPrimaryTree->Branch("fDepositedVy", &fDepositedVy, "fDepositedVy/D");
        fPrimaryTree->Branch("fDepositedVz", &fDepositedVz, "fDepositedVz/D");

	fPrimaryTree->Branch("fSecondary", &fSecondary, "fSecondary/I");

	
	fPrimaryTree->Branch("fReconstructedE", &fReconstructedE, "fReconstructedE/D");
        fPrimaryTree->Branch("fReconstructedEt", &fReconstructedEt, "fReconstructedEt/D");
	
	fPrimaryTree->Branch("fClusterMult", &fClusterMult,  "fClusterMult/I");
	
	
    }

    fHistGammasFound = new TH1F("fHistGammasFound", "fHistGammasFound",200, 0, 10);
    fHistGammasGenerated = new TH1F("fHistGammasGenerated", "fHistGammasGenerated",200, 0, 10);
    fHistGammasFoundOutOfAccCent = new TH2F("fHistGammasFoundOutOfAccCent", "fHistGammasFoundOutOfAccCent",200, 0, 10,20,-0.5,19.5);
    fHistGammasFoundCent = new TH2F("fHistGammasFoundCent", "fHistGammasFoundCent",200, 0, 10,20,-0.5,19.5);
    fHistSimEmEtCent = new TH2F("fHistSimEmEtCent", "fHistSimEmEtCent",200, 0, 200,20,-0.5,19.5);
    fHistGammasFoundOutOfAccAltCent = new TH2F("fHistGammasFoundOutOfAccAltCent", "fHistGammasFoundOutOfAccCent",200, 0, 10,20,-0.5,19.5);
    fHistGammasFoundRecoEnergyCent = new TH2F("fHistGammasFoundRecoEnergyCent", "fHistGammasFoundRecoEnergyCent",200, 0, 10,20,-0.5,19.5);
    fHistGammasFoundRecoEnergyTrueEnergyCent = new TH3F("fHistGammasFoundRecoEnergyTrueEnergyCent", "fHistGammasFoundRecoEnergyTrueEnergyCent",200, 0, 10,200, 0, 10,20,-0.5,19.5);
    fHistAllGammasFoundRecoEnergyCent = new TH2F("fHistAllGammasFoundRecoEnergyCent", "fHistAllGammasFoundRecoEnergyCent",200, 0, 10,20,-0.5,19.5);
    fHistGammasFoundOutOfAccRecoEnergyCent = new TH2F("fHistGammasFoundOutOfAccRecoEnergyCent", "fHistGammasFoundOutOfAccRecoEnergyCent",200, 0, 10,20,-0.5,19.5);
    fHistAllGammasFoundOutOfAccRecoEnergyCent = new TH2F("fHistAllGammasFoundOutOfAccRecoEnergyCent", "fHistAllGammasFoundOutOfAccRecoEnergyCent",200, 0, 10,20,-0.5,19.5);
    fHistGammasFoundAltCent = new TH2F("fHistGammasFoundAltCent", "fHistGammasFoundAltCent",200, 0, 10,20,-0.5,19.5);
    fHistGammasGeneratedCent = new TH2F("fHistGammasGeneratedCent", "fHistGammasGeneratedCent",200, 0, 10,20,-0.5,19.5);
    fHistChargedTracksCut = new TH1F("fHistChargedTracksCut", "fHistChargedTracksCut",100, 0, 5);
    fHistChargedTracksAccepted = new TH1F("fHistChargedTracksAccepted", "fHistChargedTracksAccepted",100, 0, 5);
    fHistGammasCut = new TH1F("fHistGammasTracksCut", "fHistGammasTracksCut",100, 0, 5);
    fHistGammasAccepted = new TH1F("fHistGammasTracksAccepted", "fHistGammasTracksAccepted",100, 0, 5);

    if(fCalcTrackMatchVsMult){
      fHistChargedTracksCutMult = new TH2F("fHistChargedTracksCutMult", "fHistChargedTracksCutMult",100, 0, 5,10*scale,0,100*scale);
      fHistChargedTracksAcceptedMult = new TH2F("fHistChargedTracksAcceptedMult", "fHistChargedTracksAcceptedMult",100, 0, 5,10*scale,0,100*scale);
      fHistChargedTrackDepositsAcceptedVsPt = new TH2F("fHistChargedTrackDepositsAcceptedVsPt", "fHistChargedTrackDepositsAcceptedVsPt",100, 0, 5,20,-0.5,19.5);
      fHistChargedTrackDepositsAllVsPt = new TH2F("fHistChargedTrackDepositsAllVsPt", "fHistChargedTrackDepositsAllVsPt",100, 0, 5,20,-0.5,19.5);
      fHistChargedTrackDepositsAcceptedVsPtEffCorr = new TH2F("fHistChargedTrackDepositsAcceptedVsPtEffCorr", "fHistChargedTrackDepositsAcceptedVsPtEffCorr",100, 0, 5,20,-0.5,19.5);
      fHistChargedTrackDepositsAllVsPtEffCorr = new TH2F("fHistChargedTrackDepositsAllVsPtEffCorr", "fHistChargedTrackDepositsAllVsPtEffCorr",100, 0, 5,20,-0.5,19.5);
      fHistChargedTracksAcceptedLowPtCentEffCorr = new TH2F("fHistChargedTracksAcceptedLowPtCentEffCorr", "fHistChargedTracksAcceptedLowPtCentEffCorr",100, 0, 5,20,-0.5,19.5);
      fHistChargedTracksAcceptedLowPtCent = new TH2F("fHistChargedTracksAcceptedLowPtCent", "fHistChargedTracksAcceptedLowPtCent",100, 0, 5,20,-0.5,19.5);
      fHistChargedTracksAcceptedLowPtCent500MeV = new TH2F("fHistChargedTracksAcceptedLowPtCent500MeV", "fHistChargedTracksAcceptedLowPtCent500MeV",100, 0, 5,20,-0.5,19.5);
      fHistChargedTracksAcceptedLowPtCentNoAntiProtons = new TH2F("fHistChargedTracksAcceptedLowPtCentNoAntiProtons", "fHistChargedTracksAcceptedLowPtCentNoAntiProtons",100, 0, 5,20,-0.5,19.5);
      fHistChargedTracksAcceptedLowPtCentAntiProtons = new TH2F("fHistChargedTracksAcceptedLowPtCentAntiProtons", "fHistChargedTracksAcceptedLowPtCentAntiProtons",100, 0, 5,20,-0.5,19.5);
      fHistGammasCutMult = new TH2F("fHistGammasTracksCutMult", "fHistGammasTracksCutMult",100, 0, 5,10*scale,0,100*scale);
      fHistGammasAcceptedMult = new TH2F("fHistGammasTracksAcceptedMult", "fHistGammasTracksAcceptedMult",100, 0, 5,10*scale,0,100*scale);
    }

    fHistBadTrackMatches = new TH1F("fHistBadTrackMatches", "fHistBadTrackMatches",100, 0, 5);
    fHistMatchedTracksEvspTBkgd = new TH2F("fHistMatchedTracksEvspTBkgd", "fHistMatchedTracksEvspTBkgd",100, 0, 3,100,0,3);
    fHistMatchedTracksEvspTSignal = new TH2F("fHistMatchedTracksEvspTSignal", "fHistMatchedTracksEvspTSignal",100, 0, 3,100,0,3);
    if(fCalcTrackMatchVsMult){
      fHistMatchedTracksEvspTBkgdPeripheral = new TH2F("fHistMatchedTracksEvspTBkgdPeripheral", "fHistMatchedTracksEvspTBkgd",100, 0, 3,100,0,3);
      fHistMatchedTracksEvspTSignalPeripheral = new TH2F("fHistMatchedTracksEvspTSignalPeripheral", "fHistMatchedTracksEvspTSignal",100, 0, 3,100,0,3);

      fHistMatchedTracksEvspTBkgdvsCent = new TH3F("fHistMatchedTracksEvspTBkgdvsCent", "fHistMatchedTracksEvspTBkgdvsCent",100, 0, 3,100,0,3,20,-0.5,19.5);
      fHistMatchedTracksEvspTSignalvsCent = new TH3F("fHistMatchedTracksEvspTSignalvsCent", "fHistMatchedTracksEvspTSignalvsCent",100, 0, 3,100,0,3,20,-0.5,19.5);
      fHistMatchedTracksEvspTBkgdvsCentEffCorr = new TH3F("fHistMatchedTracksEvspTBkgdvsCentEffCorr", "fHistMatchedTracksEvspTBkgdvsCent",100, 0, 3,100,0,3,20,-0.5,19.5);
      fHistMatchedTracksEvspTSignalvsCentEffCorr = new TH3F("fHistMatchedTracksEvspTSignalvsCentEffCorr", "fHistMatchedTracksEvspTSignalvsCent",100, 0, 3,100,0,3,20,-0.5,19.5);
    

      fHistChargedTracksCutPeripheral = new TH1F("fHistChargedTracksCutPeripheral", "fHistChargedTracksCut",100, 0, 5);
      fHistChargedTracksAcceptedPeripheral = new TH1F("fHistChargedTracksAcceptedPeripheral", "fHistChargedTracksAccepted",100, 0, 5);
      fHistGammasCutPeripheral = new TH1F("fHistGammasTracksCutPeripheral", "fHistGammasTracksCut",100, 0, 5);
      fHistGammasAcceptedPeripheral = new TH1F("fHistGammasTracksAcceptedPeripheral", "fHistGammasTracksAccepted",100, 0, 5);
    }
    fHistBadTrackMatchesdPhidEta = new TH2F("fHistBadTrackMatchesdPhidEta", "fHistBadTrackMatchesdPhidEta",20, -0.1, 0.1,20,-.1,0.1);
    fHistGoodTrackMatchesdPhidEta = new TH2F("fHistGoodTrackMatchesdPhidEta", "fHistGoodTrackMatchesdPhidEta",20, -0.1, 0.1,20,-.1,0.1);

    fHistHadronDepositsAll = new TH1F("fHistHadronDepositsAll","All Hadrons which deposited energy in calorimeter",fgNumOfPtBins,fgPtAxis);
    fHistHadronDepositsReco = new TH1F("fHistHadronDepositsReco","Reconstructed Hadrons which deposited energy in calorimeter",fgNumOfPtBins,fgPtAxis);
    //,10,0,100
      Int_t nMult = 50;
      Float_t nMultCuts[51] = { 0  , 5 ,10 ,15 ,20 , 25, 30, 35, 40, 45, 
				50 ,55 ,60 ,65 ,70 , 75, 80, 85, 90, 95,
				100,105,110,115,120,125,130,135,140,145,
				150,155,160,165,170,175,180,185,190,195,
				200,205,210,215,220,225,230,235,240,245,
				250};
      Int_t nCent = 20;
      Float_t nCentCuts[21] = { 0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};

      fHistHadronDepositsAllCent = new TH2F("fHistHadronDepositsAllCent","All Hadrons which deposited energy in calorimeter",fgNumOfPtBins,fgPtAxis,nCent,nCentCuts);
      fHistHadronDepositsAllCent500MeV = new TH2F("fHistHadronDepositsAllCent500MeV","All Hadrons which deposited energy in calorimeter pT>500MeV",fgNumOfPtBins,fgPtAxis,nCent,nCentCuts);
      fHistHadronDepositsRecoCent = new TH2F("fHistHadronDepositsRecoCent","Reconstructed Hadrons which deposited energy in calorimeter",fgNumOfPtBins,fgPtAxis,nCent,nCentCuts);

      fHistHadronDepositsAllvsECent = new TH2F("fHistHadronDepositsAllvsECent","All Hadrons which deposited energy in calorimeter",fgNumOfPtBins,fgPtAxis,nCent,nCentCuts);
      fHistHadronDepositsRecovsECent = new TH2F("fHistHadronDepositsRecovsECent","Reconstructed Hadrons which deposited energy in calorimeter",fgNumOfPtBins,fgPtAxis,nCent,nCentCuts);

      fHistHadronsAllCent = new TH2F("fHistHadronsAllCent","All Hadrons vs cluster mult",fgNumOfPtBins,fgPtAxis,nCent,nCentCuts);

      fHistMultChVsSignalVsMult = new TH3F("fHistMultChVsSignalVsMult","Charged particle Multiplicity vs Signal particle multiplicity vs Cluster Mult",nMult,nMultCuts,nMult,nMultCuts,nMult,nMultCuts);
      fHistNeutralRemovedSecondaryEtVsCent = new TH2F("fHistNeutralRemovedSecondaryEtVsCent","Neutral Removed Secondaries E_{T} vs centrality",100*scale,0.0,2.0*scale,20,-0.5,19.5);
      fHistChargedRemovedSecondaryEtVsCent = new TH2F("fHistChargedRemovedSecondaryEtVsCent","Charged Removed Secondaries E_{T} vs centrality",100*scale,0.0,2.0*scale,20,-0.5,19.5);
      fHistNeutralNotRemovedSecondaryEtVsCent = new TH2F("fHistNeutralNotRemovedSecondaryEtVsCent","Neutral NotRemoved Secondaries E_{T} vs centrality",100*scale,0.0,2.0*scale,20,-0.5,19.5);
      fHistChargedNotRemovedSecondaryEtVsCent = new TH2F("fHistChargedNotRemovedSecondaryEtVsCent","Charged NotRemoved Secondaries E_{T} vs centrality",100*scale,0.0,2.0*scale,20,-0.5,19.5);
      fHistNeutralRemovedSecondaryNumVsNCluster = new TH2F("fHistNeutralRemovedSecondaryNumVsNCluster","Neutral Removed Secondaries Number vs N_{cluster}",20,-0.5,19.5,250*scale,0,250*scale);
      fHistChargedRemovedSecondaryNumVsNCluster = new TH2F("fHistChargedRemovedSecondaryNumVsNCluster","Charged Removed Secondaries Number vs N_{cluster}",20,-0.5,19.5,250*scale,0,250*scale);
      fHistNeutralNotRemovedSecondaryNumVsNCluster = new TH2F("fHistNeutralNotRemovedSecondaryNumVsNCluster","Neutral NotRemoved Secondaries Number vs N_{cluster}",20,-0.5,19.5,250*scale,0,250*scale);
      fHistChargedNotRemovedSecondaryNumVsNCluster = new TH2F("fHistChargedNotRemovedSecondaryNumVsNCluster","Charged NotRemoved Secondaries Number vs N_{cluster}",20,-0.5,19.5,250*scale,0,250*scale);
      fHistNeutralRemovedSecondaryNumVsCent = new TH2F("fHistNeutralRemovedSecondaryNumVsCent","Neutral Removed Secondaries Number vs N_{cluster}",20,-0.5,19.5,20,-0.5,19.5);
      fHistChargedRemovedSecondaryNumVsCent = new TH2F("fHistChargedRemovedSecondaryNumVsCent","Charged Removed Secondaries Number vs N_{cluster}",20,-0.5,19.5,20,-0.5,19.5);
      fHistNeutralNotRemovedSecondaryNumVsCent = new TH2F("fHistNeutralNotRemovedSecondaryNumVsCent","Neutral NotRemoved Secondaries Number vs N_{cluster}",20,-0.5,19.5,20,-0.5,19.5);
      fHistChargedNotRemovedSecondaryNumVsCent = new TH2F("fHistChargedNotRemovedSecondaryNumVsCent","Charged NotRemoved Secondaries Number vs N_{cluster}",20,-0.5,19.5,20,-0.5,19.5);
      fHistNeutronsEtVsCent = new TH2F("fHistNeutronsEtVsCent","Neutrons and anti-neutrons - deposited ET vs Centrality bin",100*scale,0,4.0*scale,20,-0.5,19.5);
      fHistNeutronsNumVsCent = new TH2F("fHistNeutronsNumVsCent","Neutrons and anti-neutrons - number vs Centrality bin",20,-0.5,19.5,20,-0.5,19.5);
      fHistNotNeutronsNumVsCent = new TH2F("fHistNotNeutronsNumVsCent","Neutral particles not otherwise classified - number vs Centrality bin",20,-0.5,19.5,20,-0.5,19.5);
      Int_t nbinsEt = 150*scale;
      Float_t maxEtRange = 150*scale;
      Float_t maxEtRangeShort = 60*scale;
      Float_t minEtRange = 0;
      Int_t nbinsMult = 100*scale;
      Float_t maxMult = 3000*scale;
      Float_t minMult = 0;
      Int_t nbinsCl = 300*scale;
      Float_t maxCl = 600*scale;
      Float_t minCl = 0;
      fHistPiKPDepositedVsNch = new TH2F("fHistPiKPDepositedVsNch","#pi,K,p E_{T} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRange,nbinsMult,minMult,maxMult);
      fHistPiKPNotTrackMatchedDepositedVsNch = new TH2F("fHistPiKPNotTrackMatchedDepositedVsNch","#pi,K,p E_{T} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRange,nbinsMult,minMult,maxMult);

    fHistNeutronsDepositedVsNch = new TH2F("fHistNeutronsDepositedVsNch","n deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistAntiNeutronsDepositedVsNch = new TH2F("fHistAntiNeutronsDepositedVsNch","#bar{n} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistProtonsDepositedVsNch = new TH2F("fHistProtonsDepositedVsNch","p deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistAntiProtonsDepositedVsNch = new TH2F("fHistAntiProtonsDepositedVsNch","#bar{p} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistProtonsNotTrackMatchedDepositedVsNch = new TH2F("fHistProtonsNotTrackMatchedDepositedVsNch","p not track matched deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistAntiProtonsNotTrackMatchedDepositedVsNch = new TH2F("fHistAntiProtonsNotTrackMatchedDepositedVsNch","#bar{p} not track matched deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);



    fHistNeutronsDepositedVsNcl = new TH2F("fHistNeutronsDepositedVsNcl","n deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);
    fHistAntiNeutronsDepositedVsNcl = new TH2F("fHistAntiNeutronsDepositedVsNcl","#bar{n} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);
    fHistProtonsDepositedVsNcl = new TH2F("fHistProtonsDepositedVsNcl","p deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);
    fHistAntiProtonsDepositedVsNcl = new TH2F("fHistAntiProtonsDepositedVsNcl","#bar{p} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);
    fHistProtonsNotTrackMatchedDepositedVsNcl = new TH2F("fHistProtonsNotTrackMatchedDepositedVsNcl","p not track matched deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);
    fHistAntiProtonsNotTrackMatchedDepositedVsNcl = new TH2F("fHistAntiProtonsNotTrackMatchedDepositedVsNcl","#bar{p} not track matched deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);

    fHistSecondariesVsNch = new TH2F("fHistSecondariesVsNch","secondaries deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistSecondariesVsNcl = new TH2F("fHistSecondariesVsNcl","secondaries deposited in calorimeter vs number of clusters",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);
    fHistSecondariesEffCorrVsNch = new TH2F("fHistSecondariesEffCorrVsNch","efficiency corrected secondaries deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistSecondariesEffCorrVsNcl = new TH2F("fHistSecondariesEffCorrVsNcl","efficiency corrected secondaries deposited in calorimeter vs number of clusters",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);

    fHistSecondariesOutOfAccEffCorrVsNch = new TH2F("fHistSecondariesOutOfAccEffCorrVsNch","efficiency corrected secondaries deposited in calorimeter vs number of clusters for secondary particles out of detector acceptance",nbinsEt,minEtRange,maxEtRangeShort/10.0,nbinsMult,minMult,maxMult);
    fHistSecondariesDetectorCoverEffCorrVsNch = new TH2F("fHistSecondariesDetectorCoverEffCorrVsNch","efficiency corrected secondaries deposited in calorimeter vs number of clusters for secondaries from the detector cover",nbinsEt,minEtRange,maxEtRangeShort/10.0,nbinsMult,minMult,maxMult);
    //start
    fHistNeutronsDepositedVsNchNoEffCorr = new TH2F("fHistNeutronsDepositedVsNchNoEffCorr","n deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistAntiNeutronsDepositedVsNchNoEffCorr = new TH2F("fHistAntiNeutronsDepositedVsNchNoEffCorr","#bar{n} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistProtonsDepositedVsNchNoEffCorr = new TH2F("fHistProtonsDepositedVsNchNoEffCorr","p deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistAntiProtonsDepositedVsNchNoEffCorr = new TH2F("fHistAntiProtonsDepositedVsNchNoEffCorr","#bar{p} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistProtonsNotTrackMatchedDepositedVsNchNoEffCorr = new TH2F("fHistProtonsNotTrackMatchedDepositedVsNchNoEffCorr","p not track matched deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistAntiProtonsNotTrackMatchedDepositedVsNchNoEffCorr = new TH2F("fHistAntiProtonsNotTrackMatchedDepositedVsNchNoEffCorr","#bar{p} not track matched deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);



    fHistNeutronsDepositedVsNclNoEffCorr = new TH2F("fHistNeutronsDepositedVsNclNoEffCorr","n deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);
    fHistAntiNeutronsDepositedVsNclNoEffCorr = new TH2F("fHistAntiNeutronsDepositedVsNclNoEffCorr","#bar{n} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);
    fHistProtonsDepositedVsNclNoEffCorr = new TH2F("fHistProtonsDepositedVsNclNoEffCorr","p deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);
    fHistAntiProtonsDepositedVsNclNoEffCorr = new TH2F("fHistAntiProtonsDepositedVsNclNoEffCorr","#bar{p} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);
    fHistProtonsNotTrackMatchedDepositedVsNclNoEffCorr = new TH2F("fHistProtonsNotTrackMatchedDepositedVsNclNoEffCorr","p not track matched deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);
    fHistAntiProtonsNotTrackMatchedDepositedVsNclNoEffCorr = new TH2F("fHistAntiProtonsNotTrackMatchedDepositedVsNclNoEffCorr","#bar{p} not track matched deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);

    //end

    fHistCentVsNchVsNcl = new TH3F("fHistCentVsNchVsNcl","Cent bin vs Nch Vs NCl",20,-0.5,19.5,nbinsMult,minMult,maxMult,nbinsCl,minCl,maxCl);
    float maxpos = 500*scale;
     int nbinspos = 50*scale;
    fHistSecondaryPositionInDetector = new TH3F("fHistSecondaryPositionInDetector","Position of secondaries",nbinspos,-maxpos,maxpos,nbinspos,-maxpos,maxpos,nbinspos,-maxpos,maxpos);
    fHistSecondaryPositionInDetector->GetXaxis()->SetTitle("X");
    fHistSecondaryPositionInDetector->GetYaxis()->SetTitle("Y");
    fHistSecondaryPositionInDetector->GetZaxis()->SetTitle("Z");
//     fHistSecondaryPositionInDetectorMultiple = new TH3F("fHistSecondaryPositionInDetectorMultiple","Position of secondaries",nbinspos,-maxpos,maxpos,nbinspos,-maxpos,maxpos,nbinspos,-maxpos,maxpos);
//     fHistSecondaryPositionInDetectorMultiple->GetXaxis()->SetTitle("X");
//     fHistSecondaryPositionInDetectorMultiple->GetYaxis()->SetTitle("Y");
//     fHistSecondaryPositionInDetectorMultiple->GetZaxis()->SetTitle("Z");
    fClusterPositionWeird =  new TH2F("fClusterPositionWeird", "Position of weird secondary clusters",300, -TMath::Pi(),TMath::Pi(), 100, -0.7 , 0.7);

   fSecondaryClusterEnergy = new TH1F("fSecondaryClusterEnergy","fSecondaryClusterEnergy", 100, 0, 5);
   fHistGammaCrossCheck =  new TH2F("fHistGammaCrossCheck", "(E_{T}^{#gamma,rec}-E_{T}^{#gamma,sim})/E_{T}^{#gamma,rec}",20,-0.5,19.5,100,-1.5,4);
   fHistGammaCrossCheckAlt =  new TH2F("fHistGammaCrossCheckAlt", "(E_{T}^{#gamma,rec}-E_{T}^{#gamma,sim})/E_{T}^{#gamma,rec}",20,-0.5,19.5,100,-1.5,4);
   fHistGammaEnergyCrossCheck =  new TH2F("fHistGammaEnergyCrossCheck", "(E_{T}^{#gamma,rec}-E_{T}^{#gamma,sim})/E_{T}^{#gamma,rec}",100,0.0,10.0,100,-2,2);
   fHistGammaEnergyCrossCheckCent =  new TH3F("fHistGammaEnergyCrossCheckCent", "(E_{T}^{#gamma,rec}-E_{T}^{#gamma,sim})/E_{T}^{#gamma,rec} vs Cent",100,0.0,10.0,100,-2,2,20,-0.5,19.5);
   fHistGammaEnergyCrossCheckAlt =  new TH2F("fHistGammaEnergyCrossCheckAlt", "(E_{T}^{#gamma,rec}-E_{T}^{#gamma,sim})/E_{T}^{#gamma,rec}",100,0.0,10.0,100,-2,2);
   fHistNeutronCrossCheck =  new TH2F("fHistNeutronCrossCheck", "(E_{T}^{n,rec}-E_{T}^{n,sim})/E_{T}^{n,rec}",20,-0.5,19.5,100,-14.8,1.2);
   fHistSecondaryCrossCheck =  new TH2F("fHistSecondaryCrossCheck", "(E_{T}^{sec,rec}-E_{T}^{sec,sim})/E_{T}^{sec,rec}",20,-0.5,19.5,100,-14.8,1.2);
   fHistHadronCrossCheck =  new TH2F("fHistHadronCrossCheck", "(E_{T}^{h,rec}-E_{T}^{h,sim})/E_{T}^{h,rec}",20,-0.5,19.5,100,-7.8,1.2);
   fHistKaonCrossCheck =  new TH2F("fHistKaonCrossCheck", "(E_{T}^{K,rec}-E_{T}^{K,sim})/E_{T}^{K,rec}",20,-0.5,19.5,100,-7.8,1.2);
   fHistNeutronCorrection =  new TH2F("fHistNeutronCorrection", "E_{T}^{n,rec}",20,-0.5,19.5,200,0,50);
   fHistSecondaryCorrection =  new TH2F("fHistSecondaryCorrection", "E_{T}^{sec,rec}",20,-0.5,19.5,200,0,200);
   fHistHadronCorrection =  new TH2F("fHistHadronCorrection", "E_{T}^{h,rec}",20,-0.5,19.5,200,0,200);
   fHistKaonCorrection =  new TH2F("fHistKaonCorrection", "E_{T}^{K,rec}",20,-0.5,19.5,200,0,50);

   fHistAllEnergy = new TH1F("fHistAllEnergy","fHistAllEnergy",20,-0.5,19.5);
   fHistSignalEnergy = new TH1F("fHistSignalEnergy","fHistSignalEnergy",20,-0.5,19.5);
   fHistNeutronEnergy = new TH1F("fHistNeutronEnergy","fHistNeutronEnergy",20,-0.5,19.5);
   fHistKaonEnergy = new TH1F("fHistKaonEnergy","fHistKaonEnergy",20,-0.5,19.5);
   fHistHadronEnergy = new TH1F("fHistHadronEnergy","fHistHadronEnergy",20,-0.5,19.5);
   fHistSecondaryEnergy = new TH1F("fHistSecondaryEnergy","fHistSecondaryEnergy",20,-0.5,19.5);
   fHistSecondaryChargedEnergy = new TH1F("fHistSecondaryChargedEnergy","fHistSecondaryChargedEnergy",20,-0.5,19.5);
   fHistSecondaryNeutronEnergy = new TH1F("fHistSecondaryNeutronEnergy","fHistSecondaryNeutronEnergy",20,-0.5,19.5);
   fHistSecondaryGammaEnergy = new TH1F("fHistSecondaryGammaEnergy","fHistSecondaryGammaEnergy",20,-0.5,19.5);
   fHistSecondaryElectronEnergy = new TH1F("fHistSecondaryElectronEnergy","fHistSecondaryElectronEnergy",20,-0.5,19.5);
   fHistSecondaryOtherEnergy = new TH1F("fHistSecondaryOtherEnergy","fHistSecondaryOtherEnergy",20,-0.5,19.5);
   fHistSimulatedGammaEnergy = new TH1F("fHistSimulatedGammaEnergy","fHistSimulatedGammaEnergy",20,-0.5,19.5);
   fHistReconstructedGammaEnergy = new TH1F("fHistReconstructedGammaEnergy","fHistReconstructedGammaEnergy",20,-0.5,19.5);
   fHistSimulatedGammaEnergyAboveThreshold = new TH1F("fHistSimulatedGammaEnergyAboveThreshold","fHistSimulatedGammaEnergyAboveThreshold",20,-0.5,19.5);
   fHistReconstructedSignalEnergy = new TH1F("fHistReconstructedSignalEnergy","fHistReconstructedSignalEnergy",20,-0.5,19.5);

   Float_t fracMin = 0.0;
   Float_t fracMax = 1.0;
   Float_t fracMaxLow = 0.2;
   Float_t nBinsFrac = 50;
   fHistFracSignalVsNClusters = new TH2F("fHistFracSignalVsNClusters","fHistFracSignalVsNClusters",nbinsCl,minCl,maxCl,nBinsFrac,fracMin,fracMax);
   fHistFracHadronsVsNClusters = new TH2F("fHistFracHadronsVsNClusters","fHistFracHadronsVsNClusters",nbinsCl,minCl,maxCl,nBinsFrac,fracMin,fracMax);
   fHistFracNeutronsVsNClusters = new TH2F("fHistFracNeutronsVsNClusters","fHistFracNeutronsVsNClusters",nbinsCl,minCl,maxCl,nBinsFrac,fracMin,fracMaxLow);
   fHistFracKaonsVsNClusters = new TH2F("fHistFracKaonsVsNClusters","fHistFracKaonsVsNClusters",nbinsCl,minCl,maxCl,nBinsFrac,fracMin,fracMaxLow);
   fHistFracSecondariesVsNClusters = new TH2F("fHistFracSecondariesVsNClusters","fHistFracSecondariesVsNClusters",nbinsCl,minCl,maxCl,nBinsFrac,fracMin,fracMax);
   fHistFracSignalVsNMultiplicity = new TH2F("fHistFracSignalVsNMultiplicity","fHistFracSignalVsNMultiplicity",nbinsMult,minMult,maxMult,nBinsFrac,fracMin,fracMax);
   fHistFracHadronsVsNMultiplicity = new TH2F("fHistFracHadronsVsNMultiplicity","fHistFracHadronsVsNMultiplicity",nbinsMult,minMult,maxMult,nBinsFrac,fracMin,fracMax);
   fHistFracNeutronsVsNMultiplicity = new TH2F("fHistFracNeutronsVsNMultiplicity","fHistFracNeutronsVsNMultiplicity",nbinsMult,minMult,maxMult,nBinsFrac,fracMin,fracMaxLow);
   fHistFracKaonsVsNMultiplicity = new TH2F("fHistFracKaonsVsNMultiplicity","fHistFracKaonsVsNMultiplicity",nbinsMult,minMult,maxMult,nBinsFrac,fracMin,fracMaxLow);
   fHistFracSecondariesVsNMultiplicity = new TH2F("fHistFracSecondariesVsNMultiplicity","fHistFracSecondariesVsNMultiplicity",nbinsMult,minMult,maxMult,nBinsFrac,fracMin,fracMax);
   Int_t nBinsMatchedTracks = 100*scale;
   Float_t lowMatchedTracks = 0;
   Float_t highMatchedTracks = 100*scale;
   fHistFracSignalVsNMatchedTracks = new TH2F("fHistFracSignalVsNMatchedTracks","fHistFracSignalVsNMatchedTracks",nBinsMatchedTracks,lowMatchedTracks,highMatchedTracks,nBinsFrac,fracMin,fracMax);
   fHistFracHadronsVsNMatchedTracks = new TH2F("fHistFracHadronsVsNMatchedTracks","fHistFracHadronsVsNMatchedTracks",nBinsMatchedTracks,lowMatchedTracks,highMatchedTracks,nBinsFrac,fracMin,fracMax);
   fHistFracNeutronsVsNMatchedTracks = new TH2F("fHistFracNeutronsVsNMatchedTracks","fHistFracNeutronsVsNMatchedTracks",nBinsMatchedTracks,lowMatchedTracks,highMatchedTracks,nBinsFrac,fracMin,fracMaxLow);
   fHistFracKaonsVsNMatchedTracks = new TH2F("fHistFracKaonsVsNMatchedTracks","fHistFracKaonsVsNMatchedTracks",nBinsMatchedTracks,lowMatchedTracks,highMatchedTracks,nBinsFrac,fracMin,fracMaxLow);
   fHistFracSecondariesVsNMatchedTracks = new TH2F("fHistFracSecondariesVsNMatchedTracks","fHistFracSecondariesVsNMatchedTracks",nBinsMatchedTracks,lowMatchedTracks,highMatchedTracks,nBinsFrac,fracMin,fracMax);
   fHistFracSignalVsNTotalTracks = new TH2F("fHistFracSignalVsNTotalTracks","fHistFracSignalVsNTotalTracks",nBinsMatchedTracks,lowMatchedTracks,highMatchedTracks,nBinsFrac,fracMin,fracMax);
   fHistFracHadronsVsNTotalTracks = new TH2F("fHistFracHadronsVsNTotalTracks","fHistFracHadronsVsNTotalTracks",nBinsMatchedTracks,lowMatchedTracks,highMatchedTracks,nBinsFrac,fracMin,fracMax);
   fHistFracNeutronsVsNTotalTracks = new TH2F("fHistFracNeutronsVsNTotalTracks","fHistFracNeutronsVsNTotalTracks",nBinsMatchedTracks,lowMatchedTracks,highMatchedTracks,nBinsFrac,fracMin,fracMaxLow);
   fHistFracKaonsVsNTotalTracks = new TH2F("fHistFracKaonsVsNTotalTracks","fHistFracKaonsVsNTotalTracks",nBinsMatchedTracks,lowMatchedTracks,highMatchedTracks,nBinsFrac,fracMin,fracMaxLow);
   fHistFracSecondariesVsNTotalTracks = new TH2F("fHistFracSecondariesVsNTotalTracks","fHistFracSecondariesVsNTotalTracks",nBinsMatchedTracks,lowMatchedTracks,highMatchedTracks,nBinsFrac,fracMin,fracMax);
   fHistRCorrVsPtVsCent = new TH3F("fHistRCorrVsPtVsCent","fHistRCorrVsPtVsCent",72,0,2,50,0,10,20,-0.5,19.5);
   //NClusters
}

void AliAnalysisEtMonteCarlo::FillOutputList(TList *list)
{   //fill the output list
    AliAnalysisEt::FillOutputList(list);

  if(!fIsMC) return;
    if(fCuts->GetHistMakeTree())
    {
        list->Add(fPrimaryTree);
    }

    list->Add(fHistRemovedOrNot);

    list->Add(fHistEtNonRemovedProtons);
    list->Add(fHistEtNonRemovedAntiProtons);
    list->Add(fHistEtNonRemovedPiPlus);
    list->Add(fHistEtNonRemovedPiMinus);
    list->Add(fHistEtNonRemovedKaonPlus);
    list->Add(fHistEtNonRemovedKaonMinus);
    list->Add(fHistEtNonRemovedK0s);
    list->Add(fHistEtNonRemovedK0L);
    list->Add(fHistEtNonRemovedLambdas);
    list->Add(fHistEtNonRemovedElectrons);
    list->Add(fHistEtNonRemovedPositrons);
    list->Add(fHistEtNonRemovedMuPlus);
    list->Add(fHistEtNonRemovedMuMinus);
    list->Add(fHistEtNonRemovedNeutrons);
    list->Add(fHistEtNonRemovedAntiNeutrons);
    list->Add(fHistEtNonRemovedGammas);
    list->Add(fHistEtNonRemovedGammasFromPi0);

    list->Add(fHistEtRemovedGammas);
    list->Add(fHistEtRemovedNeutrons);
    list->Add(fHistEtRemovedAntiNeutrons);

    list->Add(fHistEtRemovedCharged);
    list->Add(fHistEtRemovedNeutrals);

    list->Add(fHistEtNonRemovedCharged);
    list->Add(fHistEtNonRemovedNeutrals);

    list->Add(fHistMultNonRemovedProtons);
    list->Add(fHistMultNonRemovedAntiProtons);
    list->Add(fHistMultNonRemovedPiPlus);
    list->Add(fHistMultNonRemovedPiMinus);
    list->Add(fHistMultNonRemovedKaonPlus);
    list->Add(fHistMultNonRemovedKaonMinus);
    list->Add(fHistMultNonRemovedK0s);
    list->Add(fHistMultNonRemovedK0L);
    list->Add(fHistMultNonRemovedLambdas);
    list->Add(fHistMultNonRemovedElectrons);
    list->Add(fHistMultNonRemovedPositrons);
    list->Add(fHistMultNonRemovedMuPlus);
    list->Add(fHistMultNonRemovedMuMinus);
    list->Add(fHistMultNonRemovedNeutrons);
    list->Add(fHistMultNonRemovedAntiNeutrons);
    list->Add(fHistMultNonRemovedGammas);

    list->Add(fHistMultRemovedGammas);
    list->Add(fHistMultRemovedNeutrons);
    list->Add(fHistMultRemovedAntiNeutrons);

    list->Add(fHistMultRemovedCharged);
    list->Add(fHistMultRemovedNeutrals);

    list->Add(fHistMultNonRemovedCharged);
    list->Add(fHistMultNonRemovedNeutrals);

    list->Add(fHistTrackMultvsNonRemovedCharged);
    list->Add(fHistTrackMultvsNonRemovedNeutral);
    list->Add(fHistTrackMultvsRemovedGamma);

    list->Add(fHistClusterMultvsNonRemovedCharged);
    list->Add(fHistClusterMultvsNonRemovedNeutral);
    list->Add(fHistClusterMultvsRemovedGamma);

    //list->Add(fHistDecayVertexNonRemovedCharged);
    //list->Add(fHistDecayVertexNonRemovedNeutral);
    //list->Add(fHistDecayVertexRemovedCharged);
    //list->Add(fHistDecayVertexRemovedNeutral);

    list->Add(fHistDxDzNonRemovedCharged);
    list->Add(fHistDxDzRemovedCharged);
    list->Add(fHistDxDzNonRemovedNeutral);
    list->Add(fHistDxDzRemovedNeutral);

    if(fCalcForKaonCorrection){
      list->Add(fHistK0EDepositsVsPtInAcceptance);
      list->Add(fHistK0EGammaVsPtInAcceptance);
      list->Add(fHistK0EDepositsVsPtOutOfAcceptance);
      list->Add(fHistK0EGammaVsPtOutOfAcceptance);
      list->Add(fHistSimKaonsInAcceptance);
      list->Add(fHistSimK0SInAcceptance);
      list->Add(fHistSimK0LInAcceptance);
      list->Add(fHistSimKPlusInAcceptance);
      list->Add(fHistSimKMinusInAcceptance);
      list->Add(fHistSimKaonsOutOfAcceptance);
      list->Add(fHistSimKaonsInAcceptanceWithDepositsPrimaries);
      list->Add(fHistSimKaonsOutOfAcceptanceWithDepositsSecondaries);
      list->Add(fHistSimKaonsOutOfAcceptanceWithDepositsPrimaries);
    }

    list->Add(fHistPiPlusMult);
    list->Add(fHistPiMinusMult);
    list->Add(fHistPiZeroMult);
    list->Add(fHistPiPlusMultAcc);
    list->Add(fHistPiMinusMultAcc);
    list->Add(fHistPiZeroMultAcc);
    
    list->Add(fHistSimEmEtCent);
    list->Add(fHistGammasFound);
    list->Add(fHistGammasGenerated);
    list->Add(fHistGammasFoundOutOfAccCent);
    list->Add(fHistGammasFoundCent);
    list->Add(fHistGammasFoundOutOfAccAltCent);
    list->Add(fHistGammasFoundAltCent);
    list->Add(fHistGammasGeneratedCent);
    list->Add(fHistGammasFoundRecoEnergyCent);
    list->Add(fHistGammasFoundRecoEnergyTrueEnergyCent);
    list->Add(fHistAllGammasFoundRecoEnergyCent);
    list->Add(fHistGammasFoundOutOfAccRecoEnergyCent);
    list->Add(fHistAllGammasFoundOutOfAccRecoEnergyCent);
    list->Add(fHistChargedTracksCut);
    list->Add(fHistChargedTracksAccepted);
    list->Add(fHistGammasCut);
    list->Add(fHistGammasAccepted);
    if(fCalcTrackMatchVsMult){
      list->Add(fHistChargedTracksCutMult);
      list->Add(fHistChargedTrackDepositsAcceptedVsPt);
      list->Add(fHistChargedTrackDepositsAllVsPt);
      list->Add(fHistChargedTrackDepositsAcceptedVsPtEffCorr);
      list->Add(fHistChargedTrackDepositsAllVsPtEffCorr);
      list->Add(fHistChargedTracksAcceptedMult);
      list->Add(fHistChargedTracksAcceptedLowPtCentEffCorr);
      list->Add(fHistChargedTracksAcceptedLowPtCent);
      list->Add(fHistChargedTracksAcceptedLowPtCent500MeV);
      list->Add(fHistChargedTracksAcceptedLowPtCentNoAntiProtons);
      list->Add(fHistChargedTracksAcceptedLowPtCentAntiProtons);
      list->Add(fHistGammasCutMult);
      list->Add(fHistGammasAcceptedMult);
    }
    list->Add(fHistBadTrackMatches);
    list->Add(fHistMatchedTracksEvspTBkgd);
    list->Add(fHistMatchedTracksEvspTSignal);
    if(fCalcTrackMatchVsMult){
      list->Add(fHistMatchedTracksEvspTBkgdPeripheral);
      list->Add(fHistMatchedTracksEvspTSignalPeripheral);
      list->Add(fHistMatchedTracksEvspTBkgdvsCent);
      list->Add(fHistMatchedTracksEvspTSignalvsCent);
      list->Add(fHistMatchedTracksEvspTBkgdvsCentEffCorr);
      list->Add(fHistMatchedTracksEvspTSignalvsCentEffCorr);
      list->Add(fHistChargedTracksCutPeripheral);
      list->Add(fHistChargedTracksAcceptedPeripheral);
      list->Add(fHistGammasCutPeripheral);
      list->Add(fHistGammasAcceptedPeripheral);
    }
    list->Add(fHistBadTrackMatchesdPhidEta);
    list->Add(fHistGoodTrackMatchesdPhidEta);
    list->Add(fHistHadronDepositsAll);
    list->Add(fHistHadronDepositsReco);
    list->Add(fHistHadronDepositsAllCent);
    list->Add(fHistHadronDepositsAllCent500MeV);
    list->Add(fHistHadronDepositsRecoCent);
    list->Add(fHistHadronDepositsAllvsECent);
    list->Add(fHistHadronDepositsRecovsECent);
    list->Add(fHistHadronsAllCent);
    list->Add(fHistMultChVsSignalVsMult);
    list->Add(fHistNeutralRemovedSecondaryEtVsCent);
    list->Add(fHistChargedRemovedSecondaryEtVsCent);
    list->Add(fHistNeutralNotRemovedSecondaryEtVsCent);
    list->Add(fHistChargedNotRemovedSecondaryEtVsCent);
    list->Add(fHistNeutralRemovedSecondaryNumVsNCluster);
    list->Add(fHistChargedRemovedSecondaryNumVsNCluster);
    list->Add(fHistNeutralNotRemovedSecondaryNumVsNCluster);
    list->Add(fHistChargedNotRemovedSecondaryNumVsNCluster);
    list->Add(fHistNeutralRemovedSecondaryNumVsCent);
    list->Add(fHistChargedRemovedSecondaryNumVsCent);
    list->Add(fHistNeutralNotRemovedSecondaryNumVsCent);
    list->Add(fHistChargedNotRemovedSecondaryNumVsCent);
    list->Add(fHistNeutronsEtVsCent);
    list->Add(fHistNeutronsNumVsCent);
    list->Add(fHistNotNeutronsNumVsCent);
    list->Add(fHistPiKPDepositedVsNch);
    list->Add(fHistPiKPNotTrackMatchedDepositedVsNch);

    list->Add(fHistNeutronsDepositedVsNch);
    list->Add(fHistAntiNeutronsDepositedVsNch);
    list->Add(fHistProtonsDepositedVsNch);
    list->Add(fHistAntiProtonsDepositedVsNch);
    list->Add(fHistProtonsNotTrackMatchedDepositedVsNch);
    list->Add(fHistAntiProtonsNotTrackMatchedDepositedVsNch);
    list->Add(fHistNeutronsDepositedVsNcl);
    list->Add(fHistAntiNeutronsDepositedVsNcl);
    list->Add(fHistProtonsDepositedVsNcl);
    list->Add(fHistAntiProtonsDepositedVsNcl);
    list->Add(fHistProtonsNotTrackMatchedDepositedVsNcl);
    list->Add(fHistAntiProtonsNotTrackMatchedDepositedVsNcl);
    list->Add(fHistSecondariesVsNch);
    list->Add(fHistSecondariesVsNcl);
    list->Add(fHistSecondariesEffCorrVsNch);
    list->Add(fHistSecondariesEffCorrVsNcl);
    list->Add(fHistSecondariesOutOfAccEffCorrVsNch);
    list->Add(fHistSecondariesDetectorCoverEffCorrVsNch);
    //start

    list->Add(fHistNeutronsDepositedVsNchNoEffCorr);
    list->Add(fHistAntiNeutronsDepositedVsNchNoEffCorr);
    list->Add(fHistProtonsDepositedVsNchNoEffCorr);
    list->Add(fHistAntiProtonsDepositedVsNchNoEffCorr);
    list->Add(fHistProtonsNotTrackMatchedDepositedVsNchNoEffCorr);
    list->Add(fHistAntiProtonsNotTrackMatchedDepositedVsNchNoEffCorr);
    list->Add(fHistNeutronsDepositedVsNclNoEffCorr);
    list->Add(fHistAntiNeutronsDepositedVsNclNoEffCorr);
    list->Add(fHistProtonsDepositedVsNclNoEffCorr);
    list->Add(fHistAntiProtonsDepositedVsNclNoEffCorr);
    list->Add(fHistProtonsNotTrackMatchedDepositedVsNclNoEffCorr);
    list->Add(fHistAntiProtonsNotTrackMatchedDepositedVsNclNoEffCorr);

    //end

    list->Add(fHistCentVsNchVsNcl);
    list->Add(fHistSecondaryPositionInDetector);
    list->Add(fClusterPositionWeird);
    //list->Add(fHistSecondaryPositionInDetectorMultiple);
    list->Add(fSecondaryClusterEnergy);
    list->Add(fHistGammaCrossCheck);
    list->Add(fHistGammaCrossCheckAlt);
    list->Add(fHistGammaEnergyCrossCheck);
    list->Add(fHistGammaEnergyCrossCheckCent);
    list->Add(fHistGammaEnergyCrossCheckAlt);
    list->Add(fHistNeutronCrossCheck);
    list->Add(fHistSecondaryCrossCheck);
    list->Add(fHistHadronCrossCheck);
    list->Add(fHistKaonCrossCheck);
    list->Add(fHistNeutronCorrection);
    list->Add(fHistSecondaryCorrection);
    list->Add(fHistHadronCorrection);
    list->Add(fHistKaonCorrection);


    list->Add(fHistAllEnergy);
    list->Add(fHistSignalEnergy);
    list->Add(fHistNeutronEnergy);
    list->Add(fHistKaonEnergy);
    list->Add(fHistHadronEnergy);
    list->Add(fHistSecondaryEnergy);
    list->Add(fHistSecondaryChargedEnergy);
    list->Add(fHistSecondaryNeutronEnergy);
    list->Add(fHistSecondaryGammaEnergy);
    list->Add(fHistSecondaryElectronEnergy);
    list->Add(fHistSecondaryOtherEnergy);
    list->Add(fHistSimulatedGammaEnergy);
    list->Add(fHistReconstructedGammaEnergy);
    list->Add(fHistSimulatedGammaEnergyAboveThreshold);
    list->Add(fHistReconstructedSignalEnergy);


    list->Add(fHistFracSignalVsNClusters);
    list->Add(fHistFracHadronsVsNClusters);
    list->Add(fHistFracNeutronsVsNClusters);
    list->Add(fHistFracKaonsVsNClusters);
    list->Add(fHistFracSecondariesVsNClusters);
    list->Add(fHistFracSignalVsNMultiplicity);
    list->Add(fHistFracHadronsVsNMultiplicity);
    list->Add(fHistFracNeutronsVsNMultiplicity);
    list->Add(fHistFracKaonsVsNMultiplicity);
    list->Add(fHistFracSecondariesVsNMultiplicity);
    list->Add(fHistFracSignalVsNMatchedTracks);
    list->Add(fHistFracHadronsVsNMatchedTracks);
    list->Add(fHistFracNeutronsVsNMatchedTracks);
    list->Add(fHistFracKaonsVsNMatchedTracks);
    list->Add(fHistFracSecondariesVsNMatchedTracks);
    list->Add(fHistFracSignalVsNTotalTracks);
    list->Add(fHistFracHadronsVsNTotalTracks);
    list->Add(fHistFracNeutronsVsNTotalTracks);
    list->Add(fHistFracKaonsVsNTotalTracks);
    list->Add(fHistFracSecondariesVsNTotalTracks);
    list->Add(fHistRCorrVsPtVsCent);

}


bool AliAnalysisEtMonteCarlo::TrackHitsCalorimeter(TParticle* part, Double_t magField)
{
    //  printf(" TrackHitsCalorimeter - magField %f\n", magField);
    AliESDtrack *esdTrack = new AliESDtrack(part);
    // Printf("MC Propagating track: eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());

    Bool_t prop = esdTrack->PropagateTo(fDetectorRadius, magField);

    // if(prop) Printf("Track propagated, eta: %f, phi: %f, pt: %f", esdTrack->Eta(), esdTrack->Phi(), esdTrack->Pt());

    bool status = prop && fSelector->CutGeometricalAcceptance(*esdTrack);
    delete esdTrack;

    return status;
}

void AliAnalysisEtMonteCarlo::FillHistograms()
{   // let base class fill its histograms, and us fill the local ones
    AliAnalysisEt::FillHistograms();
  if(!fIsMC) return;
    //std::cout << fEtNonRemovedPiPlus << " " << fCentClass << std::endl;

    fHistEtNonRemovedProtons->Fill(fEtNonRemovedProtons, fCentClass);
    fHistEtNonRemovedAntiProtons->Fill(fEtNonRemovedAntiProtons, fCentClass);
    fHistEtNonRemovedKaonPlus->Fill(fEtNonRemovedKaonPlus, fCentClass);
    fHistEtNonRemovedKaonMinus->Fill(fEtNonRemovedKaonMinus, fCentClass);
    fHistEtNonRemovedK0s->Fill(fEtNonRemovedK0S, fCentClass);
    fHistEtNonRemovedK0L->Fill(fEtNonRemovedK0L, fCentClass);
    fHistEtNonRemovedLambdas->Fill(fEtNonRemovedLambdas, fCentClass);
    fHistEtNonRemovedPiPlus->Fill(fEtNonRemovedPiPlus, fCentClass);
    fHistEtNonRemovedPiMinus->Fill(fEtNonRemovedPiMinus, fCentClass);
    fHistEtNonRemovedElectrons->Fill(fEtNonRemovedElectrons, fCentClass);
    fHistEtNonRemovedPositrons->Fill(fEtNonRemovedPositrons, fCentClass);
    fHistEtNonRemovedMuPlus->Fill(fEtNonRemovedMuPlus, fCentClass);
    fHistEtNonRemovedMuMinus->Fill(fEtNonRemovedMuMinus, fCentClass);
    fHistEtNonRemovedNeutrons->Fill(fEtNonRemovedNeutrons, fCentClass);
    fHistEtNonRemovedAntiNeutrons->Fill(fEtNonRemovedAntiNeutrons, fCentClass);
    fHistEtNonRemovedGammas->Fill(fEtNonRemovedGammas, fCentClass);
    fHistEtNonRemovedGammasFromPi0->Fill(fEtNonRemovedGammasFromPi0, fCentClass);

    fHistEtRemovedGammas->Fill(fEtRemovedGammas, fClusterMult);
    fHistEtRemovedNeutrons->Fill(fEtRemovedNeutrons, fCentClass);
    fHistEtRemovedAntiNeutrons->Fill(fEtRemovedAntiNeutrons, fCentClass);

//     fHistEtRemovedCharged->Fill(fEtRemovedAntiProtons+fEtRemovedElectrons+fEtRemovedKaonMinus+fEtRemovedKaonPlus
//                                             +fEtRemovedMuMinus+fEtRemovedMuPlus+fEtRemovedPiMinus+fEtRemovedPiPlus+fEtRemovedPositrons
//                                             +fEtRemovedProtons.
// 				fCentClass);
//     fHistEtRemovedNeutrals->Fill(fEtRemovedNeutrons+fEtRemovedAntiNeutrons, fCentClass);
//
//     fHistEtNonRemovedCharged->Fill(fEtNonRemovedAntiProtons+fEtNonRemovedElectrons+fEtNonRemovedKaonMinus+fEtNonRemovedKaonPlus
//                                             +fEtNonRemovedMuMinus+fEtNonRemovedMuPlus+fEtNonRemovedPiMinus+fEtNonRemovedPiPlus+fEtNonRemovedPositrons
//                                             +fEtNonRemovedProtons,
// 				fCentClass);
//     fHistEtRemovedNeutrals->Fill(fEtNonRemovedNeutrons+fEtNonRemovedAntiNeutrons, fCentClass);

    fHistEtRemovedCharged->Fill(fEnergyChargedRemoved, fClusterMult);
    fHistEtRemovedNeutrals->Fill(fEnergyNeutralRemoved, fClusterMult);
    fHistEtNonRemovedCharged->Fill(fEnergyChargedNotRemoved, fClusterMult);
    fHistEtNonRemovedNeutrals->Fill(fEnergyNeutralNotRemoved, fClusterMult);

    fHistMultRemovedCharged->Fill(fChargedRemoved, fClusterMult);
    fHistMultRemovedNeutrals->Fill(fNeutralRemoved, fClusterMult);
    fHistMultNonRemovedCharged->Fill(fChargedNotRemoved, fClusterMult);
    fHistMultNonRemovedNeutrals->Fill(fNeutralNotRemoved, fClusterMult);


    fHistMultNonRemovedProtons->Fill(fMultNonRemovedProtons, fCentClass);
    fHistMultNonRemovedAntiProtons->Fill(fMultNonRemovedAntiProtons, fCentClass);
    fHistMultNonRemovedKaonPlus->Fill(fMultNonRemovedKaonPlus, fCentClass);
    fHistMultNonRemovedKaonMinus->Fill(fMultNonRemovedKaonMinus, fCentClass);
    fHistMultNonRemovedK0s->Fill(fMultNonRemovedK0s, fCentClass);
    fHistMultNonRemovedK0L->Fill(fMultNonRemovedK0L, fCentClass);
    fHistMultNonRemovedLambdas->Fill(fMultNonRemovedLambdas, fCentClass);
    fHistMultNonRemovedPiPlus->Fill(fMultNonRemovedPiPlus, fCentClass);
    fHistMultNonRemovedPiMinus->Fill(fMultNonRemovedPiMinus, fCentClass);
    fHistMultNonRemovedElectrons->Fill(fMultNonRemovedElectrons, fCentClass);
    fHistMultNonRemovedPositrons->Fill(fMultNonRemovedPositrons, fCentClass);
    fHistMultNonRemovedMuPlus->Fill(fMultNonRemovedMuPlus, fCentClass);
    fHistMultNonRemovedMuMinus->Fill(fMultNonRemovedMuMinus, fCentClass);
    fHistMultNonRemovedNeutrons->Fill(fMultNonRemovedNeutrons, fCentClass);
    fHistMultNonRemovedAntiNeutrons->Fill(fMultNonRemovedAntiNeutrons, fCentClass);
    fHistMultNonRemovedGammas->Fill(fMultNonRemovedGammas, fCentClass);

    fHistMultRemovedGammas->Fill(fMultRemovedGammas, fCentClass);
    fHistMultRemovedNeutrons->Fill(fMultRemovedNeutrons, fCentClass);
    fHistMultRemovedAntiNeutrons->Fill(fMultRemovedAntiNeutrons, fCentClass);

    fHistTrackMultvsNonRemovedCharged->Fill(fTrackMultInAcc,
                                            fMultNonRemovedAntiProtons+fMultNonRemovedElectrons+fMultNonRemovedKaonMinus+fMultNonRemovedKaonPlus
                                            +fMultNonRemovedMuMinus+fMultNonRemovedMuPlus+fMultNonRemovedPiMinus+fMultNonRemovedPiPlus+fMultNonRemovedPositrons
                                            +fMultNonRemovedProtons);

    fHistTrackMultvsNonRemovedNeutral->Fill(fTrackMultInAcc,
                                            fMultNonRemovedNeutrons+fMultNonRemovedAntiNeutrons+fMultNonRemovedK0s+fMultNonRemovedK0L+fMultNonRemovedLambdas+fK0sMult);

    fHistTrackMultvsRemovedGamma->Fill(fTrackMultInAcc,
                                       fMultRemovedGammas);

    fHistClusterMultvsNonRemovedCharged->Fill(fClusterMult,
            fMultNonRemovedAntiProtons+fMultNonRemovedElectrons+fMultNonRemovedKaonMinus
            +fMultNonRemovedKaonPlus+fMultNonRemovedMuMinus+fMultNonRemovedMuPlus
            +fMultNonRemovedPiMinus+fMultNonRemovedPiPlus+fMultNonRemovedPositrons+fMultNonRemovedProtons);

    fHistClusterMultvsNonRemovedNeutral->Fill(fClusterMult,
            fMultNonRemovedNeutrons+fMultNonRemovedAntiNeutrons+fMultNonRemovedK0s+fMultNonRemovedK0L+fMultNonRemovedLambdas+fK0sMult);

    fHistClusterMultvsRemovedGamma->Fill(fClusterMult,
                                         fMultRemovedGammas);

}




Int_t AliAnalysisEtMonteCarlo::PrintFamilyTree(Int_t partIdx, AliStack* stack)
{ // print family tree
    TParticle *part = stack->Particle(partIdx);
//     if(part->GetPdgCode() == fgK0SCode)
    {
      std::cout << "This is index: " << partIdx << " (" << stack->Particle(partIdx)->GetName() <<") , is it primary: " << stack->IsPhysicalPrimary(partIdx)<<" is it from secondary interaction "<<fSelector->FromSecondaryInteraction(partIdx, *stack)<< std::endl;
        std::cout << "PID: " << part->GetPdgCode() << "/" << part->GetName() << std::endl;
        std::cout << "Energy: " << part->Energy() << std::endl;
	Float_t vtx = TMath::Sqrt( TMath::Power(part->Vx(),2) + TMath::Power(part->Vy(),2) + TMath::Power(part->Vz(),2) );
	Float_t vtxy = TMath::Sqrt( TMath::Power(part->Vx(),2) + TMath::Power(part->Vy(),2)  );
        std::cout << "Vertex: " << part->Vx() << ", " << part->Vy() << ", " << part->Vz() <<" |Vtx| "<<vtx <<" |Vtxy| "<<vtxy << std::endl;
    }
    return PrintMothers(partIdx, stack, 1);
}

Int_t AliAnalysisEtMonteCarlo::PrintMothers(Int_t partIdx, AliStack* stack, Int_t gen)
{ // print mothers
    char *tabs = new char[gen+1];
    for(Int_t i = 0; i < gen; ++i)
    {
        //std::cout << i << std::endl;
        tabs[i] = '\t';
    }
    tabs[gen] = '\0';
    Int_t mothIdx = stack->Particle(partIdx)->GetMother(0);
    if(mothIdx < 0)
    {
      delete [] tabs;
      return 0;
    }
    TParticle *mother = stack->Particle(mothIdx);
//     if(mother->GetPdgCode() == fgK0SCode)
    {
        //std::cout << tabs << "Mother of index: " << partIdx << " (" << stack->Particle(partIdx)->GetName() <<") is: " << mothIdx << ", is it primary: " << stack->IsPhysicalPrimary(mothIdx)<< std::endl;
        std::cout << tabs << "Index: " << mothIdx << std::endl;
        std::cout << tabs << "Primary: " << stack->IsPhysicalPrimary(mothIdx)<<" is it from secondary interaction "<<fSelector->FromSecondaryInteraction(partIdx, *stack) << std::endl;
        std::cout << tabs << "PID: " << mother->GetPdgCode() << "/" << mother->GetName() << std::endl;
        std::cout << tabs << "Energy: " << mother->Energy() << std::endl;
	if(mother->GetFirstMother() >= 0)
	{
	  std::cout << tabs << "Mother(s): " << stack->Particle(mother->GetFirstMother())->GetPdgCode();
	  if(mother->GetSecondMother() >= 0) std::cout << ", " << stack->Particle(mother->GetSecondMother())->GetPdgCode();
	  std::cout << std::endl;
	}
	Float_t vtx = TMath::Sqrt( TMath::Power(mother->Vx(),2) + TMath::Power(mother->Vy(),2) + TMath::Power(mother->Vz(),2) );
	Float_t vtxy = TMath::Sqrt( TMath::Power(mother->Vx(),2) + TMath::Power(mother->Vy(),2)  );
        std::cout<<tabs << "Vertex: " << mother->Vx() << ", " << mother->Vy() << ", " << mother->Vz() <<"|Vtx| "<<vtx<<" |Vtxy| "<<vtxy << std::endl;
    }
    if(mother->GetPdgCode() == fgK0SCode)
    {
//	std::cout << "K0S!!!!!!!!!!!!!11111!!!!!" << std::endl;
    }
//  std::cout << "Mother of index: " << partIdx << " (" << stack->Particle(partIdx)->GetName() <<") is: " << mothIdx << std::endl;
//  std::cout << "PID: " << mother->GetPdgCode() << "/" << mother->GetName() << std::endl;
//  std::cout << "Energy: " << mother->Energy() << std::endl;
//  std::cout << "Vertex: " << mother->Vx() << ", " << mother->Vy() << ", " << mother->Vz() << std::endl;

    delete [] tabs;
    return PrintMothers(mothIdx, stack, gen+1) + 1;
}

Int_t AliAnalysisEtMonteCarlo::GetPrimMother(Int_t partIdx, AliStack *stack)
{ // get primary mother
    if(partIdx >= 0)
    {
	//return stack->GetPrimary(partIdx);
	
        Int_t mothIdx = stack->Particle(partIdx)->GetMother(0);
        if(mothIdx < 0) return -1;
        TParticle *mother = stack->Particle(mothIdx);
        if(mother)
        {
            if(stack->IsPhysicalPrimary(mothIdx)) return mothIdx;
            else return GetPrimMother(mothIdx, stack);
        }
        else
        {
            return -1;
        }
    }
    return -1;
}

Int_t AliAnalysisEtMonteCarlo::GetK0InFamily(Int_t partIdx, AliStack* stack)
{ // get K0 in family
    if(partIdx >= 0)
    {
        if(stack->Particle(partIdx)->GetPdgCode() == fgK0SCode) return partIdx;
        Int_t mothIdx = stack->Particle(partIdx)->GetMother(0);
        if(mothIdx < 0) return -1;
        TParticle *mother = stack->Particle(mothIdx);
        if(mother)
        {
	  if(mother->GetPdgCode() == fgK0SCode)
            {
                return mothIdx;
            }
            return GetK0InFamily(mothIdx, stack);
        }
        else
        {
            return -1;
        }
    }
    return -1;
}
Int_t AliAnalysisEtMonteCarlo::PrintFamilyTreeShort(Int_t partIdx, AliStack* stack)
{ // print family tree
    TParticle *part = stack->Particle(partIdx);
    if(part){
	Float_t vtxy = TMath::Sqrt( TMath::Power(part->Vx(),2) + TMath::Power(part->Vy(),2)  );
	cout<<part->GetName()<<"( is scondary "<<fSelector->FromSecondaryInteraction(partIdx, *stack)<<", vtx "<<vtxy <<", index "<<partIdx<<")";
	//cout<<"<-"<<part->GetName()<<"("<<fSelector->FromSecondaryInteraction(partIdx, *stack)<<","<<vtxy <<")";
    }
    else{return 0;}
    Int_t value = PrintMothersShort(partIdx, stack, 1);
    cout<<endl;
    return value;
}

Int_t AliAnalysisEtMonteCarlo::PrintMothersShort(Int_t partIdx, AliStack* stack, Int_t gen)
{ // print mothers
    Int_t mothIdx = stack->Particle(partIdx)->GetMother(0);
    if(mothIdx < 0)
    {
      return 0;
    }
    TParticle *mother = stack->Particle(mothIdx);
    if(mother){
      //Float_t vtx = TMath::Sqrt( TMath::Power(mother->Vx(),2) + TMath::Power(mother->Vy(),2) + TMath::Power(mother->Vz(),2) );
      Float_t vtxy = TMath::Sqrt( TMath::Power(mother->Vx(),2) + TMath::Power(mother->Vy(),2)  );
      cout<<"<-"<<mother->GetName()<<"( is scondary "<<fSelector->FromSecondaryInteraction(mothIdx, *stack)<<", vtx "<<vtxy <<", index "<<mothIdx<<")";
      //std::cout<<tabs << "Vertex: " << mother->Vx() << ", " << mother->Vy() << ", " << mother->Vz() <<"|Vtx| "<<vtx<<" |Vtxy| "<<vtxy << std::endl;
    }
    else{return 0;}
    return PrintMothersShort(mothIdx, stack, gen+1) + 1;
}


void AliAnalysisEtMonteCarlo::SetGeneratorMinMaxParticles(AliMCEvent *eventMC){
  // In case of access only to hijing particles in cocktail
  // get the min and max labels
  // TODO: Check when generator is not the first one ...
  
  fNMCProducedMin = 0;
  fNMCProducedMax = 0;
  
  AliGenEventHeader * eventHeader = eventMC->GenEventHeader();
  
  AliGenCocktailEventHeader *cocktail = dynamic_cast<AliGenCocktailEventHeader *>(eventHeader);
  
  if(!cocktail) return ;
    
  TList *genHeaders = cocktail->GetHeaders();
  
  Int_t nGenerators = genHeaders->GetEntries();
  //printf("N generators %d \n", nGenerators);
  
  for(Int_t igen = 0; igen < nGenerators; igen++)
    {
      AliGenEventHeader * eventHeader2 = (AliGenEventHeader*)genHeaders->At(igen) ;
      TString name = eventHeader2->GetName();
      
      //printf("Generator %d: Class Name %s, Name %s, title %s \n",igen, eventHeader2->ClassName(), name.Data(), eventHeader2->GetTitle());
      
      fNMCProducedMin = fNMCProducedMax;
      fNMCProducedMax+= eventHeader2->NProduced();
      
      if(name.Contains("Hijing",TString::kIgnoreCase)){
	//cout<<"Found HIJING event and set range "<<fNMCProducedMin<<"-"<<fNMCProducedMax<<endl;
	return ;
      }
    }
        
}
AliGenEventHeader* AliAnalysisEtMonteCarlo::GetGenEventHeader(AliMCEvent *eventMC) const
{
  // Return pointer to Generated event header
  // If requested and cocktail, search for the hijing generator
  AliGenEventHeader * eventHeader = eventMC->GenEventHeader();
  AliGenCocktailEventHeader *cocktail = dynamic_cast<AliGenCocktailEventHeader *>(eventHeader);
  
  if(!cocktail) return 0x0 ;
  
  TList *genHeaders = cocktail->GetHeaders();
  
  Int_t nGenerators = genHeaders->GetEntries();
  //printf("N generators %d \n", nGenerators);
  
  for(Int_t igen = 0; igen < nGenerators; igen++)
    {
      AliGenEventHeader * eventHeader2 = (AliGenEventHeader*)genHeaders->At(igen) ;
      TString name = eventHeader2->GetName();
      
      //printf("Generator %d: Class Name %s, Name %s, title %s \n",igen, eventHeader2->ClassName(), name.Data(), eventHeader2->GetTitle());
      
      if(name.Contains("Hijing",TString::kIgnoreCase)) return eventHeader2 ;
    }
  
  return 0x0;
  
}
Bool_t AliAnalysisEtMonteCarlo::IsHIJINGLabel(Int_t label,AliMCEvent *eventMC,AliStack *stack)
{
 
  // Find if cluster/track was generated by HIJING
  
  AliGenHijingEventHeader*  hijingHeader =  dynamic_cast<AliGenHijingEventHeader *> (GetGenEventHeader(eventMC));
  
  //printf("header %p, label %d\n",hijingHeader,label);
  
  if(!hijingHeader || label < 0 ) return kFALSE;
  
  
  //printf("pass a), N produced %d\n",nproduced);
  
  if(label >= fNMCProducedMin && label < fNMCProducedMax)
  {
    //printf(" accept!, label is smaller than produced, N %d\n",nproduced);

    return kTRUE;
  }
  
  if(!stack) return kFALSE;
  
  Int_t nprimaries = stack->GetNtrack();
  
  if(label > nprimaries) return kFALSE;
    
  TParticle * mom = stack->Particle(label);
    
  Int_t iMom = label;
  Int_t iParent = mom->GetFirstMother();
  while(iParent!=-1){
    if(iParent >= fNMCProducedMin && iParent < fNMCProducedMax){
      return kTRUE;
    }
      
    iMom = iParent;
    mom = stack->Particle(iMom);
    iParent = mom->GetFirstMother();
  }
    
  return kFALSE ;
    
}




