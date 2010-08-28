#include "AliAnalysisEtMonteCarloEmcal.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"

AliAnalysisEtMonteCarloEmcal::AliAnalysisEtMonteCarloEmcal()
{
   fHistogramNameSuffix = TString("EmcalMC");
}


void AliAnalysisEtMonteCarloEmcal::Init()
{
    AliAnalysisEtMonteCarlo::Init();
    
    fDetectorRadius = EtGeometryCutsEmcal::kDetectorRadius;
    fEtaCutAcc = EtGeometryCutsEmcal::kEtaAccCut;
    fPhiCutAccMax = EtGeometryCutsEmcal::kPhiAccMaxCut*TMath::Pi()/180.;
    fPhiCutAccMin = EtGeometryCutsEmcal::kPhiAccMinCut*TMath::Pi()/180.;
    fDetectorRadius = EtGeometryCutsEmcal::kDetectorRadius;
    fClusterEnergyCut = EtReconstructedCutsEmcal::kClusterEnergyCut;
    fSingleCellEnergyCut = EtReconstructedCutsEmcal::kSingleCellEnergyCut;
    
}
