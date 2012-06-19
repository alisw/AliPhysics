//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD analysis, for EMCAL
//  - reconstruction output
//  implementation file 
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________
#include "AliAnalysisEtReconstructedEmcal.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"

using namespace std;

ClassImp(AliAnalysisEtReconstructedEmcal);


AliAnalysisEtReconstructedEmcal::AliAnalysisEtReconstructedEmcal() :
AliAnalysisEtReconstructed()
{
   fHistogramNameSuffix = TString("EmcalRec");    
}

AliAnalysisEtReconstructedEmcal::~AliAnalysisEtReconstructedEmcal() 
{
}


void AliAnalysisEtReconstructedEmcal::Init()
{ // Init
  AliAnalysisEtReconstructed::Init();
    
  fDetectorRadius = fCuts->GetGeometryEmcalDetectorRadius();
  fSingleCellEnergyCut = fCuts->GetReconstructedEmcalSingleCellEnergyCut();

}

bool AliAnalysisEtReconstructedEmcal::TrackHitsCalorimeter(AliVParticle* track, Double_t magField)
{
  return  AliAnalysisEtReconstructed::TrackHitsCalorimeter(track, magField);
}
