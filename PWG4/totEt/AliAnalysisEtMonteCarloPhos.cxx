//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis, for PHOS
//  - MC output
//  implementation file 
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________
#include "AliAnalysisEtMonteCarloPhos.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"

AliAnalysisEtMonteCarloPhos::AliAnalysisEtMonteCarloPhos()
{
   fHistogramNameSuffix = TString("PhosMC");
}

AliAnalysisEtMonteCarloPhos::~AliAnalysisEtMonteCarloPhos()
{
}


void AliAnalysisEtMonteCarloPhos::Init()
{ // Init
    AliAnalysisEtMonteCarlo::Init();
    
    fDetectorRadius = EtGeometryCutsPhos::kDetectorRadius;
    fEtaCutAcc = EtGeometryCutsPhos::kEtaAccCut;
    fPhiCutAccMax = EtGeometryCutsPhos::kPhiAccMaxCut*TMath::Pi()/180.;
    fPhiCutAccMin = EtGeometryCutsPhos::kPhiAccMinCut*TMath::Pi()/180.;
    fDetectorRadius = EtGeometryCutsPhos::kDetectorRadius;
    fClusterEnergyCut = EtReconstructedCutsPhos::kClusterEnergyCut;
    fSingleCellEnergyCut = EtReconstructedCutsPhos::kSingleCellEnergyCut;
    
}
