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
#include "AliAnalysisEtSelectorEmcal.h"
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
  fSelector = new AliAnalysisEtSelectorEmcal(fCuts);
  fDetectorRadius = fCuts->GetGeometryEmcalDetectorRadius();
  fSingleCellEnergyCut = fCuts->GetReconstructedEmcalSingleCellEnergyCut();
}
