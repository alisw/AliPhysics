#include "AliAnalysisEtMonteCarloPhos.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"

AliAnalysisEtMonteCarloPhos::AliAnalysisEtMonteCarloPhos()
{
   fHistogramNameSuffix = TString("PhosMC");
}


void AliAnalysisEtMonteCarloPhos::Init()
{
    AliAnalysisEtMonteCarlo::Init();
    
    fDetectorRadius = EtGeometryCutsPhos::kDetectorRadius;
    fEtaCutAcc = EtGeometryCutsPhos::kEtaAccCut;
    fPhiCutAccMax = EtGeometryCutsPhos::kPhiAccMaxCut*TMath::Pi()/180.;
    fPhiCutAccMin = EtGeometryCutsPhos::kPhiAccMinCut*TMath::Pi()/180.;
    fDetectorRadius = EtGeometryCutsPhos::kDetectorRadius;
    fClusterEnergyCut = EtReconstructedCutsPhos::kClusterEnergyCut;
    fSingleCellEnergyCut = EtReconstructedCutsPhos::kSingleCellEnergyCut;
    
}
