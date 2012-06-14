//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD analysis, for PHOS
//  - reconstruction output
// implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________
#include "AliAnalysisEtReconstructedPhos.h"
#include "AliAnalysisEtSelectorPhos.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include <fstream>
#include <iostream>
#include <AliPHOSGeoUtils.h>
#include "AliPHOSGeometry.h"
#include <vector>
#include "TH2I.h"

using namespace std;

ClassImp(AliAnalysisEtReconstructedPhos);

/*// Worst case (protons and neutrons):
const Double_t kMEANCHARGED = 0.335;
const Double_t kMEANNEUTRAL = 0.434;
const Double_t kMEANGAMMA = 0.374;

// Best case (pions and K0s):
const Double_t kMEANCHARGED = 0.304;
const Double_t kMEANNEUTRAL = 0.3356;
const Double_t kMEANGAMMA = 0.374;
*/
// Simulated case:
//corrEnergy =cluster->E()/(0.51 + 0.02*cluster->E());
const Double_t kMEANCHARGED = 0.307/(0.51 + 0.02*0.307);
const Double_t kMEANNEUTRAL = 0.407/(0.51 + 0.02*0.407);
const Double_t kMEANGAMMA = 0.374/(0.51 + 0.02*0.374);


AliAnalysisEtReconstructedPhos::AliAnalysisEtReconstructedPhos() :
        AliAnalysisEtReconstructed()
{
    fHistogramNameSuffix = TString("PhosRec");

    fChargedContributionCorrectionParameters[0] = -0.017;
    fChargedContributionCorrectionParameters[1] = 0.065;

    fNeutralContributionCorrectionParameters[0] = -0.002;
    fNeutralContributionCorrectionParameters[1] = 0.017;
    fNeutralContributionCorrectionParameters[2] = -3.6e-5;

    fRemovedGammaContributionCorrectionParameters[0] = 0.001;
    fRemovedGammaContributionCorrectionParameters[1] = 0.37e-5;
    fRemovedGammaContributionCorrectionParameters[2] = 0.0003;

}

AliAnalysisEtReconstructedPhos::~AliAnalysisEtReconstructedPhos()
{
}

void AliAnalysisEtReconstructedPhos::Init()
{ // Init
    AliAnalysisEtReconstructed::Init();

    fSelector = new AliAnalysisEtSelectorPhos(fCuts);
    
    fSelector->Init(137366);
    
    fDetectorRadius = fCuts->GetGeometryPhosDetectorRadius();
    fEtaCutAcc = fCuts->GetGeometryPhosEtaAccCut();
    fPhiCutAccMax = fCuts->GetGeometryPhosPhiAccMaxCut() * TMath::Pi()/180.;
    fPhiCutAccMin = fCuts->GetGeometryPhosPhiAccMinCut() * TMath::Pi()/180.;
    fClusterEnergyCut = fCuts->GetReconstructedPhosClusterEnergyCut();
    fSingleCellEnergyCut = fCuts->GetReconstructedPhosSingleCellEnergyCut();

    fClusterType = fCuts->GetPhosClusterType();
    fTrackDistanceCut = fCuts->GetPhosTrackDistanceCut();
    fTrackDxCut = fCuts->GetPhosTrackDxCut();
    fTrackDzCut = fCuts->GetPhosTrackDzCut();

    fDetector = fCuts->GetDetectorPhos();
  

}

bool AliAnalysisEtReconstructedPhos::TrackHitsCalorimeter(AliVParticle* track, Double_t magField)
{
    return  AliAnalysisEtReconstructed::TrackHitsCalorimeter(track, magField);
}

Double_t AliAnalysisEtReconstructedPhos::GetChargedContribution(Int_t clusterMult)
{ // Charged contrib
    if (clusterMult > 0)
    {
        Double_t nPart = fChargedContributionCorrectionParameters[0] + fChargedContributionCorrectionParameters[1]*clusterMult;

        Double_t contr = nPart*kMEANCHARGED;

        return contr;
    }
    return 0;

}

Double_t AliAnalysisEtReconstructedPhos::GetNeutralContribution(Int_t clusterMult)
{ // Neutral contrib
    if (clusterMult > 0)
    {
        Double_t nPart = fNeutralContributionCorrectionParameters[0] + fNeutralContributionCorrectionParameters[1]*clusterMult + fNeutralContributionCorrectionParameters[2]*clusterMult*clusterMult;

        Double_t contr = nPart*kMEANNEUTRAL;

        return contr;
    }
    return 0;
}

Double_t AliAnalysisEtReconstructedPhos::GetGammaContribution(Int_t clusterMult)
{ // Gamma contrib
    if (clusterMult > 0)
    {
        Double_t nPart = fRemovedGammaContributionCorrectionParameters[0] + fRemovedGammaContributionCorrectionParameters[1]*clusterMult + fRemovedGammaContributionCorrectionParameters[2]*clusterMult*clusterMult;

        Double_t contr = nPart*kMEANGAMMA;

        return contr;
    }
    return 0;
}




