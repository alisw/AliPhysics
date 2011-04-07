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

using namespace std;

ClassImp(AliAnalysisEtMonteCarloEmcal);


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
  
  fDetectorRadius = fCuts->GetGeometryEmcalDetectorRadius();
  fEtaCutAcc = fCuts->GetGeometryEmcalEtaAccCut();
  fPhiCutAccMax = fCuts->GetGeometryEmcalPhiAccMaxCut() * TMath::Pi()/180.;
  fPhiCutAccMin = fCuts->GetGeometryEmcalPhiAccMinCut() * TMath::Pi()/180.;
  fClusterEnergyCut = fCuts->GetReconstructedEmcalClusterEnergyCut();
  fSingleCellEnergyCut = fCuts->GetReconstructedEmcalSingleCellEnergyCut();
  
  fDetector = fCuts->GetDetectorEmcal();

}
