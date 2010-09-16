//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for MC analysis, for EMCAL
//  - MC output
//  implementation file 
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________
#include "AliAnalysisEtMonteCarloEmcal.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"

AliAnalysisEtMonteCarloEmcal::AliAnalysisEtMonteCarloEmcal()
{
   fHistogramNameSuffix = TString("EmcalMC");
}

AliAnalysisEtMonteCarloEmcal::~AliAnalysisEtMonteCarloEmcal()
{
}


void AliAnalysisEtMonteCarloEmcal::Init()
{ // Init
    AliAnalysisEtMonteCarlo::Init();
    
    fDetectorRadius = EtGeometryCutsEmcal::kDetectorRadius;
    fEtaCutAcc = EtGeometryCutsEmcal::kEtaAccCut;
    fPhiCutAccMax = EtGeometryCutsEmcal::kPhiAccMaxCut*TMath::Pi()/180.;
    fPhiCutAccMin = EtGeometryCutsEmcal::kPhiAccMinCut*TMath::Pi()/180.;
    fDetectorRadius = EtGeometryCutsEmcal::kDetectorRadius;
    fClusterEnergyCut = EtReconstructedCutsEmcal::kClusterEnergyCut;
    fSingleCellEnergyCut = EtReconstructedCutsEmcal::kSingleCellEnergyCut;
    
}
