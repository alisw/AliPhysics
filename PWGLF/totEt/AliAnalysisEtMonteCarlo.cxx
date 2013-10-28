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
using namespace std;

ClassImp(AliAnalysisEtMonteCarlo);


// ctor
AliAnalysisEtMonteCarlo::AliAnalysisEtMonteCarlo():AliAnalysisEt()
						  ,fIsMC(kTRUE)
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
						  ,fCalcTrackMatchVsMult(kFALSE)
						  ,fHistGammasFound(0)
						  ,fHistGammasGenerated(0)
						  ,fHistGammasFoundCent(0)
						  ,fHistGammasFoundOutOfAccCent(0)
						  ,fHistGammasGeneratedCent(0)
						  ,fHistChargedTracksCut(0)
						  ,fHistChargedTracksAccepted(0)
						  ,fHistGammasCut(0)
						  ,fHistGammasAccepted(0)
						  ,fHistChargedTracksCutMult(0)
						  ,fHistChargedTracksAcceptedMult(0)
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
						  ,fHistNeutronsDepositedVsNch(0)
						  ,fHistAntiNeutronsDepositedVsNch(0)
						  ,fHistProtonsDepositedVsNch(0)
						  ,fHistAntiProtonsDepositedVsNch(0)
						  ,fHistProtonsNotTrackMatchedDepositedVsNch(0)
						  ,fHistAntiProtonsNotTrackMatchedDepositedVsNch(0)
						  ,fHistSecondariesVsNch(0)
						  ,fHistSecondariesVsNcl(0)
						  ,fHistSecondariesEffCorrVsNch(0)
						  ,fHistSecondariesEffCorrVsNcl(0)
						  ,fHistSecondariesOutOfAccEffCorrVsNch(0)
						  ,fHistSecondariesDetectorCoverEffCorrVsNch(0)
						  ,fHistCentVsNchVsNcl(0)
						//,fHistSecondaryPositionInDetector(0)
						  ,fClusterPositionWeird(0)
						//,fHistSecondaryPositionInDetectorMultiple(0)
						  ,fSecondaryClusterEnergy(0)
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
    delete fHistGammasFound; // enter comment here
    delete fHistGammasGenerated; // enter comment here
    delete fHistGammasFoundOutOfAccCent; // enter comment here
    delete fHistGammasFoundCent; // enter comment here
    delete fHistGammasGeneratedCent; // enter comment here
    delete fHistChargedTracksCut;
    delete fHistChargedTracksAccepted;
    delete fHistGammasCut;
    delete fHistGammasAccepted;
    delete fHistChargedTracksCutMult;
    delete fHistChargedTracksAcceptedMult;
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
    delete fHistSecondariesVsNch;
    delete fHistSecondariesVsNcl;
    delete fHistSecondariesEffCorrVsNch;
    delete fHistSecondariesEffCorrVsNcl;
    delete fHistSecondariesOutOfAccEffCorrVsNch;
    delete fHistSecondariesDetectorCoverEffCorrVsNch;
    delete fHistCentVsNchVsNcl;
    //delete fHistSecondaryPositionInDetector;
    delete fClusterPositionWeird;
    //delete fHistSecondaryPositionInDetectorMultiple;
    delete fSecondaryClusterEnergy;
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

    Int_t partCount = 0;
    for (Int_t iPart = 0; iPart < nPrim; iPart++)
    {

        TParticle *part = stack->Particle(iPart);
	
	

        if (!part)
        {
            Printf("ERROR: Could not get particle %d", iPart);
            continue;
        }
        TParticlePDG *pdg = part->GetPDG(0);

        if (!pdg)
        {
            Printf("ERROR: Could not get particle PDG %d", iPart);
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
//             PrintFamilyTree(iPart, stack);
//
            if (TMath::Abs(part->Eta()) < fCuts->GetCommonEtaCut())
            {
                //Printf("Particle with eta: %f, pid: %d", part->Eta(), code);
                // calculate E_T
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


//                 // Fill up total E_T counters for each particle species
//                 if (code == fgProtonCode || code == fgAntiProtonCode)
//                 {
//                 }
//                 if (code == fgPiPlusCode || code == fgPiMinusCode)
//                 {
//                     if (code == fgPiPlusCode)
//                     {
//                     }
//                     else
//                     {
//                     }
//                 }
//                 if (code == fgGammaCode)
//                 {
//                 }
//                 if (code == fgKPlusCode)
//                 {
//                 }
//                 if(code == fgKMinusCode)
//                 {
//                 }
//                 if (code == fgMuPlusCode || code == fgMuMinusCode)
//                 {
//                 }
//                 if (code == fgEPlusCode || code == fgEMinusCode)
//                 {
//                 }
//                 // some neutrals also
//                 if (code == fgNeutronCode)
//                 {
//                 }
//                 if (code == fgAntiNeutronCode)
//                 {
//                 }
//                 if (code == fgGammaCode)
//                 {
//                 }

                // Neutral particles
                //if (TMath::Abs(pdg->Charge() - fCuts->GetMonteCarloNeutralParticle()) <1e-3 )

                if(code == fgGammaCode || code == fgPi0Code || code == fgEtaCode)
                {
		//  PrintFamilyTree(iPart,stack);
                    //Printf("Gamma, phi: %f, eta: %f, phi cut min: %f, phi cut max: %f, eta cut: %f", part->Phi(), part->Eta(), fPhiMinCutAcc, fPhiMaxCutAcc, fEtaCutAcc);
                    //if (et > fCuts->GetCommonClusterEnergyCut()) fTotNeutralEt += et;

                    // inside EMCal acceptance

                    //if (TMath::Abs(part->Eta()) < fEtaCutAcc && part->Phi() < fPhiMaxCutAcc && part->Phi() > fPhiMinCutAcc)

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
//Int_t AliAnalysisEtMonteCarlo::AnalyseEvent(AliMCEvent* mcEvent,AliESDEvent* realEvent)
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

    std::vector<Int_t> foundGammas;
    
    fSelector->SetEvent(realEvent);

    AnalyseEvent(ev);

    AliStack *stack = mcEvent->Stack();

    // get all detector clusters
    //  TRefArray* caloClusters = new TRefArray();

//    if (fDetector == fCuts->GetDetectorEmcal()) realEvent->GetEMCALClusters( caloClusters );
    //else if (fDetector == fCuts->GetDetectorPhos()) realEvent->GetPHOSClusters( caloClusters );
    //else {
    //AliFatal("Detector ID has not been specified");
    //return -1;
//    }

//Note that this only returns clusters for the selected detector.  fSelector actually calls the right GetClusters... for the detector
//It does not apply any cuts on these clusters
    TRefArray *caloClusters = fSelector->GetClusters();

    Int_t nCluster = caloClusters->GetEntries();
    fClusterMult = nCluster;
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
    Float_t multiplicity = fEsdtrackCutsTPC->GetReferenceMultiplicity(realEvent,kTRUE);
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
	Float_t matchedTrackp = 0.0;
	Float_t matchedTrackpt = 0.0;
	nottrackmatched = fSelector->PassTrackMatchingCut(*caloCluster);
	//by default ALL matched tracks are accepted, whether or not the match is good.  So we check to see if the track is good.
	if(!nottrackmatched){//if the track is trackmatched
	  Int_t trackMatchedIndex = caloCluster->GetTrackMatchedIndex();
	  if(trackMatchedIndex < 0) nottrackmatched=kTRUE;
	  AliESDtrack *track = realEvent->GetTrack(trackMatchedIndex);
	  matchedTrackp = track->P();
	  matchedTrackpt = track->Pt();
	  //if this is a good track, accept track will return true.  The track matched is good, so not track matched is false
	  nottrackmatched = !(fEsdtrackCutsTPC->AcceptTrack(track));//if the track is bad, this is track matched
	  //if(!nottrackmatched) cout<<"Matched track p: "<<matchedTrackp<<" sim "<<part->P()<<endl;
	}
     

	for(UInt_t i = 0; i < caloCluster->GetNLabels(); i++)
	{
	  Int_t pIdx = caloCluster->GetLabelAt(i);

	  //TParticle *p = stack->Particle(pIdx);
	  
	  if(!stack->IsPhysicalPrimary(pIdx))
	  {
// 	    PrintFamilyTree(pIdx, stack);
	    pIdx = GetPrimMother(pIdx, stack);
	  }
 	  if(fSelector->PassDistanceToBadChannelCut(*caloCluster))//&&fSelector->CutGeometricalAcceptance(*(stack->Particle(primIdx))))
	  {

//	    std::cout << "Gamma primary: " << primIdx << std::endl;
// 	    foundGammas.push_back(primIdx); 
	    if(nottrackmatched){
	      foundGammas.push_back(pIdx); 
	    }
	  }
	}
	fCutFlow->Fill(cf++);
        if(!fSelector->PassDistanceToBadChannelCut(*caloCluster)) continue;
        Double_t clEt = CorrectForReconstructionEfficiency(*caloCluster,fCentClass);
//	if(code == fgK0SCode) std::cout << "K0 energy: " << caloCluster->E() << std::endl;
        if(!fSelector->PassMinEnergyCut(*caloCluster)) continue;

	
        fCutFlow->Fill(cf++);
        Float_t pos[3];
	//PrintFamilyTree(
        caloCluster->GetPosition(pos);
        TVector3 cp(pos);

        TParticle *primPart = stack->Particle(primIdx);
        fPrimaryCode = primPart->GetPdgCode();
        fPrimaryCharge = (Int_t) primPart->GetPDG()->Charge();

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

        fDepositedCode = part->GetPdgCode();
	fDepositedE = part->Energy();
        fDepositedEt = part->Energy()*TMath::Sin(part->Theta());
        fDepositedCharge = (Int_t) part->GetPDG()->Charge();

        fDepositedVx = part->Vx();
        fDepositedVy = part->Vy();
        fDepositedVz = part->Vz();
	fReconstructedE = caloCluster->E();
	fReconstructedEt = caloCluster->E()*TMath::Sin(cp.Theta());

	//fSecondary = fSelector->FromSecondaryInteraction(*primPart, *stack);
	fSecondary =fSelector->FromSecondaryInteraction(*part, *stack);
	//if(fSecondary && fReconstructedEt<0.3) fSecondary = kFALSE;//patch to do quick cross check THIS SHOULD NOT GET COMMITTED!!!  IF IT DID YELL AT CHRISTINE ASAP
// 	if(fSecondary) 
// 	{
// 	  //std::cout << "Have secondary!" << std::endl;
// 	  //PrintFamilyTree(iPart, stack);
// 	}
	
	pdg = primPart->GetPDG(0);
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
	  written = kTRUE;
	  if(nottrackmatched){//secondaries not removed
// 	    Float_t vtx = TMath::Sqrt( TMath::Power(part->Vx(),2) + TMath::Power(part->Vy(),2) + TMath::Power(part->Vz(),2) );



	    if(!fSelector->CutGeometricalAcceptance(*part)){
	      if(TMath::Sqrt(part->Vx() * part->Vx() + part->Vy() * part->Vy() + part->Vz()* part->Vz())>100){
		// 	      fClusterPositionWeird->Fill(cp.Phi(), cp.PseudoRapidity());
		// 	      //Principal arc tangent of x, in the interval [-pi/2,+pi/2] radians
		// 	      float phiVertex = TMath::ATan(part->Vy()/part->Vx());
		// 	      //PHOS acceptance is -0.222 pi to -0.555 pi
		// 	      //but tan(phi) = tan(phi+pi)
		// 	      if(part->Vx()<0) phiVertex + TMath::Pi();
		// 	      //maximum angular distance from PHOS
		// 	      float dphi1 = phiVertex +40.0/180.0*TMath::Pi();
		// 	      float dphi2 = phiVertex +100.0/180.0*TMath::Pi();
		// 	      float dphi = dphi1;
		// 	      if(TMath::Abs(dphi1)>TMath::Abs(dphi2)) dphi= dphi2;
		// 	      cout<<"DPhi:  "<<dphi<<" phi "<<phiVertex<<" dphi1 "<<dphi1<<" dphi2 "<<dphi2<<endl;
		// 	      //cout<<"dphi "<<dphi1/TMath::Pi()<<"pi "<<dphi2/TMath::Pi()<<"pi phi "<<phiVertex/TMath::Pi()<<endl;
		// 	      PrintFamilyTree(iPart, stack);
		//out of acceptance clusters

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
	    //if(fReconstructedEt>0.3){//patch to do quick cross check THIS SHOULD NOT GET COMMITTED!!!  IF IT DID YELL AT CHRISTINE ASAP
// 	      fHistSecondaryPositionInDetector->Fill(part->Vx(),part->Vy(),part->Vz());
// 	      if(caloCluster->GetNLabels()>1){
// 		fHistSecondaryPositionInDetectorMultiple->Fill(part->Vx(),part->Vy(),part->Vz());
// 	      }
	      //}
// 	    if(vtx>300){
// 	      //cout<<"Vtx "<<vtx<<endl;
// 	      if(fPrimaryCode==fgProtonCode ||  fPrimaryCode==fgAntiProtonCode || fPrimaryCode==fgPiPlusCode || fPrimaryCode==fgPiMinusCode || fPrimaryCode==fgKPlusCode || fPrimaryCode==fgKMinusCode){
// 		cout<<"I think I am really a charged hadron!"<<endl;
// 	      }
// 	      //PrintFamilyTree(iPart, stack);
// 	      }
	    fSecondaryNotRemoved++;
	    etSecondaries += fReconstructedEt;
	    etSecondariesEffCorr += clEt;
	    if (fDepositedCharge != 0){//charged track not removed
	      fChargedNotRemoved++;
	      fEnergyChargedNotRemoved += clEt;
	      fHistRemovedOrNot->Fill(2.0, fCentClass);
	      fHistDxDzNonRemovedCharged->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());
	      fHistChargedNotRemovedSecondaryEtVsCent->Fill(fReconstructedEt,fCentClass);
	      nChargedSecondariesNotRemoved++;
	    }
	    else{
	      fHistNeutralNotRemovedSecondaryEtVsCent->Fill(fReconstructedEt,fCentClass);
	      nNeutralSecondariesNotRemoved++;
	    }
	  }
	  else{//secondaries removed
            if (fDepositedCharge != 0){
	      fChargedRemoved++;
	      fEnergyChargedRemoved += clEt;
	      fHistRemovedOrNot->Fill(0.0, fCentClass);
	      fHistDxDzRemovedCharged->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());
	      fHistChargedRemovedSecondaryEtVsCent->Fill(fReconstructedEt,fCentClass);
	      nChargedSecondariesRemoved++;
	      fHistChargedTracksCut->Fill(fDepositedEt);
	      if(fCalcTrackMatchVsMult){
		fHistChargedTracksCutMult->Fill(fDepositedEt,fClusterMult);
		if(fClusterMult<25){fHistChargedTracksCutPeripheral->Fill(fDepositedEt);}
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
	    else{//neutral energy removed
	      fNeutralRemoved++;
	      fEnergyNeutralRemoved += clEt;
	      fHistRemovedOrNot->Fill(1.0, fCentClass);
	      fHistDxDzRemovedNeutral->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());
	      fHistNeutralRemovedSecondaryEtVsCent->Fill(fReconstructedEt,fCentClass);

	      nNeutralSecondariesRemoved++;
	    }
	  }
	}
	else{//not a secondary

	  if (fDepositedCharge != 0 && fDepositedCode!=fgEMinusCode && fDepositedCode!=fgEPlusCode){//if the particle hitting the calorimeter is pi/k/p/mu
	    written = kTRUE;
	    fClusterMultChargedTracks++;
	    etPiKPDeposited += clEt;
	    if(fDepositedCode==fgProtonCode){
	      etProtonDeposited += clEt;
	    }
	    if(fDepositedCode==fgAntiProtonCode){
	      etAntiProtonDeposited += clEt;
	    }
	    if(nottrackmatched){//not removed but should be
	      etPiKPDepositedNotTrackMatched += clEt;
	      if(fDepositedCode==fgProtonCode){
		etProtonDepositedNotTrackMatched += clEt;
	      }
	      if(fDepositedCode==fgAntiProtonCode){
		etAntiProtonDepositedNotTrackMatched += clEt;
	      }
	      fHistHadronDepositsAll->Fill(part->Pt());
	      fHistHadronDepositsAllCent->Fill(part->Pt(), fCentClass);
	      if(fReconstructedEt>0.5) fHistHadronDepositsAllCent500MeV->Fill(part->Pt(), fCentClass);
	      fChargedNotRemoved++;
	      fEnergyChargedNotRemoved += clEt;
	      fHistRemovedOrNot->Fill(2.0, fCentClass);
	      fHistDxDzNonRemovedCharged->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());
	      fHistChargedTracksAccepted->Fill(fDepositedEt);
	      if(fCalcTrackMatchVsMult){
		if(matchedTrackpt<0.5){//if we could never have matched this because of its pt, how much energy did it deposit?
		  fHistChargedTracksAcceptedLowPtCent->Fill(fDepositedEt, fCentClass);
		  if(fDepositedEt>=0.5) fHistChargedTracksAcceptedLowPtCent500MeV->Fill(fDepositedEt, fCentClass);
		  if(pdg->PdgCode()!=fgAntiProtonCode){
		    fHistChargedTracksAcceptedLowPtCentNoAntiProtons->Fill(fDepositedEt, fCentClass);
		  }
		  else{
		    fHistChargedTracksAcceptedLowPtCentAntiProtons->Fill(fDepositedEt, fCentClass);
		  }
		}
		fHistChargedTracksAcceptedMult->Fill(fDepositedEt,fClusterMult);
		if(fClusterMult<25){fHistChargedTracksAcceptedPeripheral->Fill(fDepositedEt);}
	      }
	    }
	    else{//removed and should have been
	      Int_t trackindex =  fSelector->GetLabel(caloCluster,stack);// (caloCluster->GetLabelsArray())->At(0);
	      fHistHadronDepositsReco->Fill(part->Pt());
	      fHistHadronDepositsRecoCent->Fill(part->Pt(), fCentClass);
	      fHistHadronDepositsAll->Fill(part->Pt());
	      fHistHadronDepositsAllCent->Fill(part->Pt(), fCentClass);
	      if(fReconstructedEt>0.5) fHistHadronDepositsAllCent500MeV->Fill(part->Pt(), fCentClass);
	      if(caloCluster->GetLabel()!=trackindex){
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
	      fHistChargedTracksCut->Fill(fDepositedEt);
	      if(fCalcTrackMatchVsMult){
		fHistChargedTracksCutMult->Fill(fDepositedEt,fClusterMult);
		if(fClusterMult<25){fHistChargedTracksCutPeripheral->Fill(fDepositedEt);}
	      }
	      fHistMatchedTracksEvspTBkgd->Fill(matchedTrackp,fReconstructedE);
	      if(fCalcTrackMatchVsMult){
		if(fClusterMult<25){fHistMatchedTracksEvspTBkgdPeripheral->Fill(matchedTrackp,fReconstructedEt);}
		fHistMatchedTracksEvspTBkgdvsCent->Fill(matchedTrackp,fReconstructedEt, fCentClass);
		fHistMatchedTracksEvspTBkgdvsCentEffCorr->Fill(matchedTrackp,clEt, fCentClass);//fill with the efficiency corrected energy
	      }
	    }
	  }
	  //K0L and any neutral particles from the decay of K+/- or K0S
	  if(!written && (fPrimaryCode==fgKPlusCode || fPrimaryCode==fgKMinusCode || fPrimaryCode==fgK0SCode ||fPrimaryCode==fgK0LCode)){
	    written = kTRUE;//At this point we are not tracking them but we don't count them as neutrals accidentally removed
	  }
	  
	  if(!written && (fDepositedCode==fgGammaCode || fDepositedCode==fgEMinusCode || fDepositedCode ==fgEPlusCode)){//if the particle hitting the calorimeter is gamma, electron and not from a kaon
	    fClusterMultGammas++;
	    written = kTRUE;
	    if(nottrackmatched){//Not removed and not supposed to be removed - signal
	      fEtNonRemovedGammas += clEt;
	      fMultNonRemovedGammas++;
	      fNeutralNotRemoved--;
	      fEnergyNeutralNotRemoved -= clEt;
	      fHistGammasAccepted->Fill(fDepositedEt);
	      if(fCalcTrackMatchVsMult){
		fHistGammasAcceptedMult->Fill(fDepositedEt,fClusterMult);
		if(fClusterMult<25){fHistGammasAcceptedPeripheral->Fill(fDepositedEt);}
	      }
	    }
	    else{//removed but shouldn't have been
	      Int_t trackindex = fSelector->GetLabel(caloCluster,stack);// (caloCluster->GetLabelsArray())->At(1);
	      if(caloCluster->GetLabel()!=trackindex){
		fHistBadTrackMatches->Fill(part->Pt(),fReconstructedE);
		fHistBadTrackMatchesdPhidEta->Fill(caloCluster->GetTrackDx(),caloCluster->GetTrackDz());
// 		cout<<"Track matched, label cluster "<<caloCluster->GetLabel()<<" track "<<trackindex<<endl;
// 		PrintFamilyTree(trackindex, stack);
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
	  //all other cases - neutron, anti-neutron, not aware of other cases
	  if(!written){
	    fNeutralNotRemoved++;
	    fEnergyNeutralNotRemoved += clEt;//this is the efficiency corrected energy
	    fHistRemovedOrNot->Fill(3.0, fCentClass);
	    if (fDepositedCode == fgNeutronCode || fDepositedCode == fgAntiNeutronCode){
	      fHistDxDzNonRemovedNeutral->Fill(caloCluster->GetTrackDx(), caloCluster->GetTrackDz());
	      fHistNeutronsEtVsCent->Fill(clEt,fCentClass);
	      nNeutrons++;
	      if(fDepositedCode == fgNeutronCode){
		etNeutronDeposited += clEt;
	      }
	      if(fDepositedCode == fgAntiNeutronCode){
		etAntiNeutronDeposited += clEt;
	      }
	      //cout<<"I am a";
	      //if(fDepositedCode == fgAntiNeutronCode) cout<<"n anti-";
	      //else{cout<<" ";}
	      //cout<<"neutron!! pt "<<part->Pt()<<" eta "<<part->Eta()<<" phi "<<part->Phi()<<" Et "<<clEt<<endl;
	      //PrintFamilyTree(iPart, stack);
	    }
	    else{
	      nNotNeutrons++;
// 	      cout<<"I am a";
// 	      if(fDepositedCode == fgAntiNeutronCode) cout<<"n anti-";
// 	      else{cout<<" ";}
// 	      cout<<"neutron!! pt "<<part->Pt()<<" eta "<<part->Eta()<<" phi "<<part->Phi()<<" Et "<<clEt<<endl;
// 	      PrintFamilyTree(iPart, stack);
// 	      cout<<"I am not a neutron!!"<<endl;
// 	      PrintFamilyTree(iPart, stack);
	    }
	  }
	}
        fPrimaryTree->Fill();
    } // end of loop over clusters     
    fHistPiKPNotTrackMatchedDepositedVsNch->Fill(etPiKPDepositedNotTrackMatched,multiplicity);
    fHistPiKPDepositedVsNch->Fill(etPiKPDeposited,multiplicity);
    fHistNeutronsDepositedVsNch->Fill(etNeutronDeposited,multiplicity);
    fHistAntiNeutronsDepositedVsNch->Fill(etAntiNeutronDeposited,multiplicity);
    fHistProtonsDepositedVsNch->Fill(etProtonDeposited,multiplicity);
    fHistAntiProtonsDepositedVsNch->Fill(etAntiProtonDeposited,multiplicity);
    fHistProtonsNotTrackMatchedDepositedVsNch->Fill(etProtonDepositedNotTrackMatched,multiplicity);
    fHistAntiProtonsNotTrackMatchedDepositedVsNch->Fill(etAntiProtonDepositedNotTrackMatched,multiplicity);
    fHistSecondariesVsNch->Fill(etSecondaries,multiplicity);
    fHistSecondariesVsNcl->Fill(etSecondaries,fNClusters);
    fHistSecondariesEffCorrVsNch->Fill(etSecondariesEffCorr,multiplicity);
    fHistSecondariesEffCorrVsNcl->Fill(etSecondariesEffCorr,fNClusters);
    fHistSecondariesOutOfAccEffCorrVsNch->Fill(etSecondariesOutOfAccEffCorr,multiplicity);
    fHistSecondariesDetectorCoverEffCorrVsNch->Fill(etSecondariesDetectorCoverEffCorr,multiplicity);
    fHistCentVsNchVsNcl->Fill(fCentClass,multiplicity, fNClusters);

    fHistNeutronsNumVsCent->Fill(nNeutrons,fCentClass);
    fHistNotNeutronsNumVsCent->Fill(nNotNeutrons,fCentClass);

    std::sort(foundGammas.begin(), foundGammas.end());
    for (Int_t iPart = 0; iPart < stack->GetNtrack(); iPart++)
    {

	if(!stack->IsPhysicalPrimary(iPart)) continue;
	
	TParticle *part = stack->Particle(iPart);

        if (!part)
        {
            Printf("ERROR: Could not get particle %d", iPart);
            continue;
        }
        TParticlePDG *pdg = part->GetPDG(0);

        if (!pdg)
        {
            Printf("ERROR: Could not get particle PDG %d", iPart);
            continue;
        }
        
        if(pdg->PdgCode()==fgGammaCode)// TMath::Abs(part->Eta()) < 0.12)
	{
	  if(fSelector->CutGeometricalAcceptance(*part)){
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
	  fHistHadronsAllCent->Fill(part->Pt(), fCentClass);
	}
	
        
    }
    if(fCalcForKaonCorrection){
      Float_t etCuts[11] = {0.0,0.05,0.1,0.15,0.2,0.25, 0.3,0.35,0.4,0.45,0.5};
      Int_t nEtCuts = 11;
      //loop over simulated particles in order to find K0S
      for (Int_t iPart = 0; iPart < stack->GetNtrack(); iPart++){
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
    fHistMultChVsSignalVsMult->Fill(fClusterMultChargedTracks,fClusterMultGammas,fNClusters);
    fHistNeutralRemovedSecondaryNumVsNCluster->Fill(nNeutralSecondariesRemoved,fNClusters);
    fHistChargedRemovedSecondaryNumVsNCluster->Fill(nChargedSecondariesRemoved,fNClusters);
    fHistNeutralNotRemovedSecondaryNumVsNCluster->Fill(nNeutralSecondariesNotRemoved,fNClusters);
    fHistChargedNotRemovedSecondaryNumVsNCluster->Fill(nChargedSecondariesNotRemoved,fNClusters);
    fHistNeutralRemovedSecondaryNumVsCent->Fill(nNeutralSecondariesRemoved,fCentClass);
    fHistChargedRemovedSecondaryNumVsCent->Fill(nChargedSecondariesRemoved,fCentClass);
    fHistNeutralNotRemovedSecondaryNumVsCent->Fill(nNeutralSecondariesNotRemoved,fCentClass);
    fHistChargedNotRemovedSecondaryNumVsCent->Fill(nChargedSecondariesNotRemoved,fCentClass);
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

    //fHistDecayVertexNonRemovedCharged = new TH3F("fHistDecayVertexNonRemovedCharged","fHistDecayVertexNonRemovedCharged", 500, -470, 30, 500, -300, 300, 40, -20, 20);
    //fHistDecayVertexRemovedCharged = new TH3F("fHistDecayVertexRemovedCharged","fHistDecayVertexRemovedCharged", 500, -470, 30, 500, -300, 300, 40, -20, 20);
    //fHistDecayVertexNonRemovedNeutral = new TH3F("fHistDecayVertexNonRemovedNeutral","fHistDecayVertexNonRemovedNeutral", 500, -470, 30, 500, -300, 300, 40, -20, 20);
    //fHistDecayVertexRemovedNeutral = new TH3F("fHistDecayVertexRemovedNeutral","fHistDecayVertexRemovedNeutral", 500, -470, 30, 500, -300, 300, 40, -20, 20);

    fHistRemovedOrNot = new TH2F("fHistRemovedOrNot", "fHistRemovedOrNot", 4, -0.5, 3.5, 10, -0.5, 9.5);

    fHistEtNonRemovedProtons = new TH2F("fHistEtNonRemovedProtons", "fHistEtNonRemovedProtons", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedAntiProtons = new TH2F("fHistEtNonRemovedAntiProtons", "fHistEtNonRemovedAntiProtons", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedPiPlus = new TH2F("fHistEtNonRemovedPiPlus", "fHistEtNonRemovedPiPlus", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedPiMinus = new TH2F("fHistEtNonRemovedPiMinus", "fHistEtNonRemovedPiMinus", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedKaonPlus = new TH2F("fHistEtNonRemovedKaonPlus", "fHistEtNonRemovedKaonPlus", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedKaonMinus = new TH2F("fHistEtNonRemovedKaonMinus", "fHistEtNonRemovedKaonMinus", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedK0s = new TH2F("fHistEtNonRemovedK0s", "fHistEtNonRemovedK0s", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedK0L = new TH2F("fHistEtNonRemovedK0L", "fHistEtNonRemovedK0L", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedLambdas = new TH2F("fHistEtNonRemovedLambdas", "fHistEtNonRemovedLambdas", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedElectrons = new TH2F("fHistEtNonRemovedElectrons", "fHistEtNonRemovedElectrons", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedPositrons = new TH2F("fHistEtNonRemovedPositrons", "fHistEtNonRemovedPositrons", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedMuPlus = new TH2F("fHistEtNonRemovedMuPlus", "fHistEtNonRemovedMuPlus", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedMuMinus = new TH2F("fHistEtNonRemovedMuMinus", "fHistEtNonRemovedMuMinus", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedNeutrons = new TH2F("fHistEtNonRemovedNeutrons", "fHistEtNonRemovedNeutrons", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedAntiNeutrons = new TH2F("fHistEtNonRemovedAntiNeutrons", "fHistEtNonRemovedAntiNeutrons", 1500, 0, 30, 10, -0.5, 9.5);

    fHistEtNonRemovedGammas = new  TH2F("fHistEtNonRemovedGammas", "fHistEtNonRemovedGammas", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedGammasFromPi0 = new  TH2F("fHistEtNonRemovedGammasFromPi0", "fHistEtNonRemovedGammasFromPi0", 1500, 0, 30, 10, -0.5, 9.5);

    fHistEtRemovedGammas = new  TH2F("fHistEtRemovedGammas", "fHistEtRemovedGammas", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtRemovedNeutrons = new  TH2F("fHistEtRemovedNeutrons", "fHistEtRemovedNeutrons", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtRemovedAntiNeutrons = new  TH2F("fHistEtRemovedAntiNeutrons", "fHistEtRemovedAntiNeutrons", 1500, 0, 30, 10, -0.5, 9.5);

    fHistEtRemovedCharged = new  TH2F("fHistEtRemovedCharged", "fHistEtRemovedCharged", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtRemovedNeutrals = new  TH2F("fHistEtRemovedNeutrals", "fHistEtRemovedNeutrals", 1500, 0, 30, 10, -0.5, 9.5);

    fHistEtNonRemovedCharged = new  TH2F("fHistEtNonRemovedCharged", "fHistEtNonRemovedCharged", 1500, 0, 30, 10, -0.5, 9.5);
    fHistEtNonRemovedNeutrals = new  TH2F("fHistEtNonRemovedNeutrals", "fHistEtNonRemovedNeutrals", 1500, 0, 30, 10, -0.5, 9.5);

    fHistMultNonRemovedProtons = new TH2F("fHistMultNonRemovedProtons", "fHistMultNonRemovedProtons", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedAntiProtons = new TH2F("fHistMultNonRemovedAntiProtons", "fHistMultNonRemovedAntiProtons", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedPiPlus = new TH2F("fHistMultNonRemovedPiPlus", "fHistMultNonRemovedPiPlus", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedPiMinus = new TH2F("fHistMultNonRemovedPiMinus", "fHistMultNonRemovedPiMinus", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedKaonPlus = new TH2F("fHistMultNonRemovedKaonPlus", "fHistMultNonRemovedKaonPlus", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedKaonMinus = new TH2F("fHistMultNonRemovedKaonMinus", "fHistMultNonRemovedKaonMinus", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedK0s = new TH2F("fHistMultNonRemovedK0s", "fHistMultNonRemovedK0s", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedK0L = new TH2F("fHistMultNonRemovedK0L", "fHistMultNonRemovedK0L", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedLambdas = new TH2F("fHistMultNonRemovedLambdas", "fHistMultNonRemovedLambdas", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedElectrons = new TH2F("fHistMultNonRemovedElectrons", "fHistMultNonRemovedElectrons", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedPositrons = new TH2F("fHistMultNonRemovedPositrons", "fHistMultNonRemovedPositrons", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedMuPlus = new TH2F("fHistMultNonRemovedMuPlus", "fHistMultNonRemovedMuPlus", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedMuMinus = new TH2F("fHistMultNonRemovedMuMinus", "fHistMultNonRemovedMuMinus", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedNeutrons = new TH2F("fHistMultNonRemovedNeutrons", "fHistMultNonRemovedNeutrons", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultNonRemovedAntiNeutrons = new TH2F("fHistMultNonRemovedAntiNeutrons", "fHistMultNonRemovedAntiNeutrons", 100, -0.5, 99.5, 10, -0.5, 9.5);

    fHistMultNonRemovedGammas = new  TH2F("fHistMultNonRemovedGammas", "fHistMultNonRemovedGammas", 100, -0.5, 99.5, 100, -0.5, 99.5);

    fHistMultRemovedGammas = new  TH2F("fHistMultRemovedGammas", "fHistMultRemovedGammas", 100, -0.5, 99.5, 100, -0.5, 99.5);
    fHistMultRemovedNeutrons = new  TH2F("fHistMultRemovedNeutrons", "fHistMultRemovedNeutrons", 100, -0.5, 99.5, 10, -0.5, 9.5);
    fHistMultRemovedAntiNeutrons = new  TH2F("fHistMultRemovedAntiNeutrons", "fHistMultRemovedAntiNeutrons", 100, -0.5, 99.5, 10, -0.5, 9.5);
    /*
      fHistMultRemovedCharged = new  TH2F("fHistMultRemovedCharged", "fHistMultRemovedCharged", 1500, 0, 30, 10, -0.5, 9.5);
      fHistMultRemovedNeutrals = new  TH2F("fHistMultRemovedNeutrals", "fHistMultRemovedNeutrals", 1500, 0, 30, 10, -0.5, 9.5);

      fHistMultNonRemovedCharged = new  TH2F("fHistMultNonRemovedCharged", "fHistMultNonRemovedCharged", 1500, 0, 30, 10, -0.5, 9.5);
      fHistMultNonRemovedNeutrals = new  TH2F("fHistMultNonRemovedNeutrals", "fHistMultNonRemovedNeutrals", 1500, 0, 30, 10, -0.5, 9.5);*/


    fHistMultRemovedCharged = new  TH2F("fHistMultRemovedCharged", "fHistMultRemovedCharged", 100, -0.5, 99.5, 100, -0.5, 99.5);
    fHistMultRemovedNeutrals = new  TH2F("fHistMultRemovedNeutrals", "fHistMultRemovedNeutrals", 100, -0.5, 99.5, 100, -0.5, 99.5);

    fHistMultNonRemovedCharged = new  TH2F("fHistMultNonRemovedCharged", "fHistMultNonRemovedCharged", 100, -0.5, 99.5, 100, -0.5, 99.5);
    fHistMultNonRemovedNeutrals = new  TH2F("fHistMultNonRemovedNeutrals", "fHistMultNonRemovedNeutrals", 100, -0.5, 99.5, 100, -0.5, 99.5);

    fHistTrackMultvsNonRemovedCharged = new TH2F("fHistTrackMultvsNonRemovedCharged", "fHistTrackMultvsNonRemovedCharged", 1000, -0.5, 999.5, 100, -0.5, 99.5);
    fHistTrackMultvsNonRemovedNeutral = new TH2F("fHistTrackMultvsNonRemovedNeutral", "fHistTrackMultvsNonRemovedNeutral", 1000, -0.5, 999.5, 100, -0.5, 99.5);
    fHistTrackMultvsRemovedGamma = new TH2F("fHistTrackMultvsRemovedGamma", "fHistTrackMultvsRemovedGamma", 1000, -0.5, 999.5, 100, -0.5, 99.5);

    fHistClusterMultvsNonRemovedCharged = new TH2F("fHistClusterMultvsNonRemovedCharged", "fHistClusterMultvsNonRemovedCharged", 1000, -0.5, 999.5, 100, -0.5, 99.5);
    fHistClusterMultvsNonRemovedNeutral = new TH2F("fHistClusterMultvsNonRemovedNeutral", "fHistClusterMultvsNonRemovedNeutral", 1000, -0.5, 999.5, 100, -0.5, 99.5);
    fHistClusterMultvsRemovedGamma = new TH2F("fHistClusterMultvsRemovedGamma", "fHistClusterMultvsRemovedGamma", 1000, -0.5, 999.5, 100, -0.5, 99.5);

    fHistDxDzNonRemovedCharged = new TH2F("fHistDxDzNonRemovedCharged", "fHistDxDzNonRemovedCharged", 800, -200, 200, 800, -200, 200);
    fHistDxDzRemovedCharged = new TH2F("fHistDxDzRemovedCharged", "fHistDxDzRemovedCharged", 800, -200, 200, 800, -200, 200);
    fHistDxDzNonRemovedNeutral = new TH2F("fHistDxDzNonRemovedNeutral", "fHistDxDzNonRemovedNeutral", 800, -200, 200, 800, -200, 200);
    fHistDxDzRemovedNeutral = new TH2F("fHistDxDzRemovedNeutral", "fHistDxDzRemovedNeutral", 800, -200, 200, 800, -200, 200);

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

    fHistPiPlusMult = new TH1F("fHistPiPlusMult", "fHistPiPlusMult", 2000, -0.5, 1999.5);
    fHistPiMinusMult = new TH1F("fHistPiMinusMult", "fHistPiMinusMult", 2000, -0.5, 1999.5);
    fHistPiZeroMult = new TH1F("fHistPiZeroMult", "fHistPiZeroMult", 2000, -0.5, 1999.5);

    fHistPiPlusMultAcc = new TH1F("fHistPiPlusMultAcc", "fHistPiPlusMultAcc", 2000, -0.5, 1999.5);
    fHistPiMinusMultAcc = new TH1F("fHistPiMinusMultAcc", "fHistPiMinusMultAcc", 2000, -0.5, 1999.5);
    fHistPiZeroMultAcc = new TH1F("fHistPiZeroMultAcc", "fHistPiZeroMultAcc", 2000, -0.5, 1999.5);

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
    fHistGammasGeneratedCent = new TH2F("fHistGammasGeneratedCent", "fHistGammasGeneratedCent",200, 0, 10,20,-0.5,19.5);
    fHistChargedTracksCut = new TH1F("fHistChargedTracksCut", "fHistChargedTracksCut",100, 0, 5);
    fHistChargedTracksAccepted = new TH1F("fHistChargedTracksAccepted", "fHistChargedTracksAccepted",100, 0, 5);
    fHistGammasCut = new TH1F("fHistGammasTracksCut", "fHistGammasTracksCut",100, 0, 5);
    fHistGammasAccepted = new TH1F("fHistGammasTracksAccepted", "fHistGammasTracksAccepted",100, 0, 5);

    if(fCalcTrackMatchVsMult){
      fHistChargedTracksCutMult = new TH2F("fHistChargedTracksCutMult", "fHistChargedTracksCutMult",100, 0, 5,10,0,100);
      fHistChargedTracksAcceptedMult = new TH2F("fHistChargedTracksAcceptedMult", "fHistChargedTracksAcceptedMult",100, 0, 5,10,0,100);
      fHistChargedTracksAcceptedLowPtCent = new TH2F("fHistChargedTracksAcceptedLowPtCent", "fHistChargedTracksAcceptedLowPtCent",100, 0, 5,20,-0.5,19.5);
      fHistChargedTracksAcceptedLowPtCent500MeV = new TH2F("fHistChargedTracksAcceptedLowPtCent500MeV", "fHistChargedTracksAcceptedLowPtCent500MeV",100, 0, 5,20,-0.5,19.5);
      fHistChargedTracksAcceptedLowPtCentNoAntiProtons = new TH2F("fHistChargedTracksAcceptedLowPtCentNoAntiProtons", "fHistChargedTracksAcceptedLowPtCentNoAntiProtons",100, 0, 5,20,-0.5,19.5);
      fHistChargedTracksAcceptedLowPtCentAntiProtons = new TH2F("fHistChargedTracksAcceptedLowPtCentAntiProtons", "fHistChargedTracksAcceptedLowPtCentAntiProtons",100, 0, 5,20,-0.5,19.5);
      fHistGammasCutMult = new TH2F("fHistGammasTracksCutMult", "fHistGammasTracksCutMult",100, 0, 5,10,0,100);
      fHistGammasAcceptedMult = new TH2F("fHistGammasTracksAcceptedMult", "fHistGammasTracksAcceptedMult",100, 0, 5,10,0,100);
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
      Int_t nMult = 20;
      Float_t nMultCuts[21] = { 0, 5,10,15,20, 25,30,35,40,45, 
			       50,55,60,65,70, 75,80,85,90,95,
				100};
      Int_t nCent = 20;
      Float_t nCentCuts[21] = { 0, 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};

      fHistHadronDepositsAllCent = new TH2F("fHistHadronDepositsAllCent","All Hadrons which deposited energy in calorimeter",fgNumOfPtBins,fgPtAxis,nCent,nCentCuts);
      fHistHadronDepositsAllCent500MeV = new TH2F("fHistHadronDepositsAllCent500MeV","All Hadrons which deposited energy in calorimeter pT>500MeV",fgNumOfPtBins,fgPtAxis,nCent,nCentCuts);
      fHistHadronDepositsRecoCent = new TH2F("fHistHadronDepositsRecoCent","Reconstructed Hadrons which deposited energy in calorimeter",fgNumOfPtBins,fgPtAxis,nCent,nCentCuts);

      fHistHadronsAllCent = new TH2F("fHistHadronsAllCent","All Hadrons vs cluster mult",fgNumOfPtBins,fgPtAxis,nCent,nCentCuts);

      fHistMultChVsSignalVsMult = new TH3F("fHistMultChVsSignalVsMult","Charged particle Multiplicity vs Signal particle multiplicity vs Cluster Mult",nMult,nMultCuts,nMult,nMultCuts,nMult,nMultCuts);
      fHistNeutralRemovedSecondaryEtVsCent = new TH2F("fHistNeutralRemovedSecondaryEtVsCent","Neutral Removed Secondaries E_{T} vs centrality",100,0.0,2.0,20,-0.5,19.5);
      fHistChargedRemovedSecondaryEtVsCent = new TH2F("fHistChargedRemovedSecondaryEtVsCent","Charged Removed Secondaries E_{T} vs centrality",100,0.0,2.0,20,-0.5,19.5);
      fHistNeutralNotRemovedSecondaryEtVsCent = new TH2F("fHistNeutralNotRemovedSecondaryEtVsCent","Neutral NotRemoved Secondaries E_{T} vs centrality",100,0.0,2.0,20,-0.5,19.5);
      fHistChargedNotRemovedSecondaryEtVsCent = new TH2F("fHistChargedNotRemovedSecondaryEtVsCent","Charged NotRemoved Secondaries E_{T} vs centrality",100,0.0,2.0,20,-0.5,19.5);
      fHistNeutralRemovedSecondaryNumVsNCluster = new TH2F("fHistNeutralRemovedSecondaryNumVsNCluster","Neutral Removed Secondaries Number vs N_{cluster}",20,-0.5,19.5,250,0,250);
      fHistChargedRemovedSecondaryNumVsNCluster = new TH2F("fHistChargedRemovedSecondaryNumVsNCluster","Charged Removed Secondaries Number vs N_{cluster}",20,-0.5,19.5,250,0,250);
      fHistNeutralNotRemovedSecondaryNumVsNCluster = new TH2F("fHistNeutralNotRemovedSecondaryNumVsNCluster","Neutral NotRemoved Secondaries Number vs N_{cluster}",20,-0.5,19.5,250,0,250);
      fHistChargedNotRemovedSecondaryNumVsNCluster = new TH2F("fHistChargedNotRemovedSecondaryNumVsNCluster","Charged NotRemoved Secondaries Number vs N_{cluster}",20,-0.5,19.5,250,0,250);
      fHistNeutralRemovedSecondaryNumVsCent = new TH2F("fHistNeutralRemovedSecondaryNumVsCent","Neutral Removed Secondaries Number vs N_{cluster}",20,-0.5,19.5,20,-0.5,19.5);
      fHistChargedRemovedSecondaryNumVsCent = new TH2F("fHistChargedRemovedSecondaryNumVsCent","Charged Removed Secondaries Number vs N_{cluster}",20,-0.5,19.5,20,-0.5,19.5);
      fHistNeutralNotRemovedSecondaryNumVsCent = new TH2F("fHistNeutralNotRemovedSecondaryNumVsCent","Neutral NotRemoved Secondaries Number vs N_{cluster}",20,-0.5,19.5,20,-0.5,19.5);
      fHistChargedNotRemovedSecondaryNumVsCent = new TH2F("fHistChargedNotRemovedSecondaryNumVsCent","Charged NotRemoved Secondaries Number vs N_{cluster}",20,-0.5,19.5,20,-0.5,19.5);
      fHistNeutronsEtVsCent = new TH2F("fHistNeutronsEtVsCent","Neutrons and anti-neutrons - deposited ET vs Centrality bin",100,0,4.0,20,-0.5,19.5);
      fHistNeutronsNumVsCent = new TH2F("fHistNeutronsNumVsCent","Neutrons and anti-neutrons - number vs Centrality bin",20,-0.5,19.5,20,-0.5,19.5);
      fHistNotNeutronsNumVsCent = new TH2F("fHistNotNeutronsNumVsCent","Neutral particles not otherwise classified - number vs Centrality bin",20,-0.5,19.5,20,-0.5,19.5);
      Int_t nbinsEt = 125;
      Float_t maxEtRange = 125;
      Float_t maxEtRangeShort = 25;
      Float_t minEtRange = 0;
      Int_t nbinsMult = 100;
      Float_t maxMult = 3000;
      Float_t minMult = 0;
      Int_t nbinsCl = 175;
      Float_t maxCl = 350;
      Float_t minCl = 0;
      fHistPiKPDepositedVsNch = new TH2F("fHistPiKPDepositedVsNch","#pi,K,p E_{T} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRange,nbinsMult,minMult,maxMult);
      fHistPiKPNotTrackMatchedDepositedVsNch = new TH2F("fHistPiKPNotTrackMatchedDepositedVsNch","#pi,K,p E_{T} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRange,nbinsMult,minMult,maxMult);
    fHistNeutronsDepositedVsNch = new TH2F("fHistNeutronsDepositedVsNch","n deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistAntiNeutronsDepositedVsNch = new TH2F("fHistAntiNeutronsDepositedVsNch","#bar{n} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistProtonsDepositedVsNch = new TH2F("fHistProtonsDepositedVsNch","p deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistAntiProtonsDepositedVsNch = new TH2F("fHistAntiProtonsDepositedVsNch","#bar{p} deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistProtonsNotTrackMatchedDepositedVsNch = new TH2F("fHistProtonsNotTrackMatchedDepositedVsNch","p not track matched deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistAntiProtonsNotTrackMatchedDepositedVsNch = new TH2F("fHistAntiProtonsNotTrackMatchedDepositedVsNch","#bar{p} not track matched deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistSecondariesVsNch = new TH2F("fHistSecondariesVsNch","secondaries deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistSecondariesVsNcl = new TH2F("fHistSecondariesVsNcl","secondaries deposited in calorimeter vs number of clusters",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);
    fHistSecondariesEffCorrVsNch = new TH2F("fHistSecondariesEffCorrVsNch","efficiency corrected secondaries deposited in calorimeter vs multiplicity",nbinsEt,minEtRange,maxEtRangeShort,nbinsMult,minMult,maxMult);
    fHistSecondariesEffCorrVsNcl = new TH2F("fHistSecondariesEffCorrVsNcl","efficiency corrected secondaries deposited in calorimeter vs number of clusters",nbinsEt,minEtRange,maxEtRangeShort,nbinsCl,minCl,maxCl);

    fHistSecondariesOutOfAccEffCorrVsNch = new TH2F("fHistSecondariesOutOfAccEffCorrVsNch","efficiency corrected secondaries deposited in calorimeter vs number of clusters for secondary particles out of detector acceptance",nbinsEt,minEtRange,maxEtRangeShort/10.0,nbinsMult,minMult,maxMult);
    fHistSecondariesDetectorCoverEffCorrVsNch = new TH2F("fHistSecondariesDetectorCoverEffCorrVsNch","efficiency corrected secondaries deposited in calorimeter vs number of clusters for secondaries from the detector cover",nbinsEt,minEtRange,maxEtRangeShort/10.0,nbinsMult,minMult,maxMult);
    fHistCentVsNchVsNcl = new TH3F("fHistCentVsNchVsNcl","Cent bin vs Nch Vs NCl",20,-0.5,19.5,nbinsMult,minMult,maxMult,nbinsCl,minCl,maxCl);
    //float maxpos = 500;
    // int nbinspos = 200;
//     fHistSecondaryPositionInDetector = new TH3F("fHistSecondaryPositionInDetector","Position of secondaries",nbinspos,-maxpos,maxpos,nbinspos,-maxpos,maxpos,nbinspos,-maxpos,maxpos);
//     fHistSecondaryPositionInDetector->GetXaxis()->SetTitle("X");
//     fHistSecondaryPositionInDetector->GetYaxis()->SetTitle("Y");
//     fHistSecondaryPositionInDetector->GetZaxis()->SetTitle("Z");
//     fHistSecondaryPositionInDetectorMultiple = new TH3F("fHistSecondaryPositionInDetectorMultiple","Position of secondaries",nbinspos,-maxpos,maxpos,nbinspos,-maxpos,maxpos,nbinspos,-maxpos,maxpos);
//     fHistSecondaryPositionInDetectorMultiple->GetXaxis()->SetTitle("X");
//     fHistSecondaryPositionInDetectorMultiple->GetYaxis()->SetTitle("Y");
//     fHistSecondaryPositionInDetectorMultiple->GetZaxis()->SetTitle("Z");
    fClusterPositionWeird =  new TH2F("fClusterPositionWeird", "Position of weird secondary clusters",300, -TMath::Pi(),TMath::Pi(), 100, -0.7 , 0.7);

   fSecondaryClusterEnergy = new TH1F("fSecondaryClusterEnergy","fSecondaryClusterEnergy", 100, 0, 5);
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
    
    list->Add(fHistGammasFound);
    list->Add(fHistGammasGenerated);
    list->Add(fHistGammasFoundOutOfAccCent);
    list->Add(fHistGammasFoundCent);
    list->Add(fHistGammasGeneratedCent);
    list->Add(fHistChargedTracksCut);
    list->Add(fHistChargedTracksAccepted);
    list->Add(fHistGammasCut);
    list->Add(fHistGammasAccepted);
    if(fCalcTrackMatchVsMult){
      list->Add(fHistChargedTracksCutMult);
      list->Add(fHistChargedTracksAcceptedMult);
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
    list->Add(fHistSecondariesVsNch);
    list->Add(fHistSecondariesVsNcl);
    list->Add(fHistSecondariesEffCorrVsNch);
    list->Add(fHistSecondariesEffCorrVsNcl);
    list->Add(fHistSecondariesOutOfAccEffCorrVsNch);
    list->Add(fHistSecondariesDetectorCoverEffCorrVsNch);
    list->Add(fHistCentVsNchVsNcl);
    //list->Add(fHistSecondaryPositionInDetector);
    list->Add(fClusterPositionWeird);
    //list->Add(fHistSecondaryPositionInDetectorMultiple);
    list->Add(fSecondaryClusterEnergy);


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

    fHistEtRemovedGammas->Fill(fEtRemovedGammas, fNClusters);
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

    fHistEtRemovedCharged->Fill(fEnergyChargedRemoved, fNClusters);
    fHistEtRemovedNeutrals->Fill(fEnergyNeutralRemoved, fNClusters);
    fHistEtNonRemovedCharged->Fill(fEnergyChargedNotRemoved, fNClusters);
    fHistEtNonRemovedNeutrals->Fill(fEnergyNeutralNotRemoved, fNClusters);

    fHistMultRemovedCharged->Fill(fChargedRemoved, fNClusters);
    fHistMultRemovedNeutrals->Fill(fNeutralRemoved, fNClusters);
    fHistMultNonRemovedCharged->Fill(fChargedNotRemoved, fNClusters);
    fHistMultNonRemovedNeutrals->Fill(fNeutralNotRemoved, fNClusters);


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

    fHistClusterMultvsNonRemovedCharged->Fill(fNClusters,
            fMultNonRemovedAntiProtons+fMultNonRemovedElectrons+fMultNonRemovedKaonMinus
            +fMultNonRemovedKaonPlus+fMultNonRemovedMuMinus+fMultNonRemovedMuPlus
            +fMultNonRemovedPiMinus+fMultNonRemovedPiPlus+fMultNonRemovedPositrons+fMultNonRemovedProtons);

    fHistClusterMultvsNonRemovedNeutral->Fill(fNClusters,
            fMultNonRemovedNeutrons+fMultNonRemovedAntiNeutrons+fMultNonRemovedK0s+fMultNonRemovedK0L+fMultNonRemovedLambdas+fK0sMult);

    fHistClusterMultvsRemovedGamma->Fill(fNClusters,
                                         fMultRemovedGammas);

}




Int_t AliAnalysisEtMonteCarlo::PrintFamilyTree(Int_t partIdx, AliStack* stack)
{ // print family tree
    TParticle *part = stack->Particle(partIdx);
//     if(part->GetPdgCode() == fgK0SCode)
    {
        std::cout << "This is index: " << partIdx << " (" << stack->Particle(partIdx)->GetName() <<") , is it primary: " << stack->IsPhysicalPrimary(partIdx)<< std::endl;
        std::cout << "PID: " << part->GetPdgCode() << "/" << part->GetName() << std::endl;
        std::cout << "Energy: " << part->Energy() << std::endl;
	Float_t vtx = TMath::Sqrt( TMath::Power(part->Vx(),2) + TMath::Power(part->Vy(),2) + TMath::Power(part->Vz(),2) );
        std::cout << "Vertex: " << part->Vx() << ", " << part->Vy() << ", " << part->Vz() <<"|Vtx| "<<vtx << std::endl;
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
        std::cout << tabs << "Primary: " << stack->IsPhysicalPrimary(mothIdx) << std::endl;
        std::cout << tabs << "PID: " << mother->GetPdgCode() << "/" << mother->GetName() << std::endl;
        std::cout << tabs << "Energy: " << mother->Energy() << std::endl;
	if(mother->GetFirstMother() >= 0)
	{
	  std::cout << tabs << "Mother(s): " << stack->Particle(mother->GetFirstMother())->GetPdgCode();
	  if(mother->GetSecondMother() >= 0) std::cout << ", " << stack->Particle(mother->GetSecondMother())->GetPdgCode();
	  std::cout << std::endl;
	}
	Float_t vtx = TMath::Sqrt( TMath::Power(mother->Vx(),2) + TMath::Power(mother->Vy(),2) + TMath::Power(mother->Vz(),2) );
        std::cout<<tabs << "Vertex: " << mother->Vx() << ", " << mother->Vy() << ", " << mother->Vz() <<"|Vtx| "<<vtx << std::endl;
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

