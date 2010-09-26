//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD analysis, for PHOS
//  - reconstruction output
// implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________
#include "AliAnalysisEtReconstructedPhos.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"

using namespace std;

ClassImp(AliAnalysisEtReconstructedPhos);


AliAnalysisEtReconstructedPhos::AliAnalysisEtReconstructedPhos() :
AliAnalysisEtReconstructed()
{
   fHistogramNameSuffix = TString("PhosRec");    
}

AliAnalysisEtReconstructedPhos::~AliAnalysisEtReconstructedPhos() 
{
}

void AliAnalysisEtReconstructedPhos::Init()
{ // Init
  AliAnalysisEtReconstructed::Init();
    
  fDetectorRadius = fCuts->GetGeometryPhosDetectorRadius();
  fEtaCutAcc = fCuts->GetGeometryPhosEtaAccCut();
  fPhiCutAccMax = fCuts->GetGeometryPhosPhiAccMaxCut() * TMath::Pi()/180.;
  fPhiCutAccMin = fCuts->GetGeometryPhosPhiAccMinCut() * TMath::Pi()/180.;
  fClusterEnergyCut = fCuts->GetReconstructedPhosClusterEnergyCut();
  fSingleCellEnergyCut = fCuts->GetReconstructedPhosSingleCellEnergyCut();
  
  fClusterType = fCuts->GetReconstructedPhosClusterType();
  fTrackDistanceCut = fCuts->GetReconstructedPhosTrackDistanceCut();

}

bool AliAnalysisEtReconstructedPhos::TrackHitsCalorimeter(AliVParticle* track, Double_t magField)
{
  return  AliAnalysisEtReconstructed::TrackHitsCalorimeter(track, magField);
}

