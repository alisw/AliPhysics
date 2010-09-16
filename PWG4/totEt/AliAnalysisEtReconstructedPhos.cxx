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
    
    fClusterType = EtReconstructedCutsPhos::kClusterType;
    fDetectorRadius = EtGeometryCutsPhos::kDetectorRadius;
    fEtaCutAcc = EtGeometryCutsPhos::kEtaAccCut;
    fPhiCutAccMax = EtGeometryCutsPhos::kPhiAccMaxCut*TMath::Pi()/180.;
    fPhiCutAccMin = EtGeometryCutsPhos::kPhiAccMinCut*TMath::Pi()/180.;
    fClusterEnergyCut = EtReconstructedCutsPhos::kClusterEnergyCut;
    fSingleCellEnergyCut = EtReconstructedCutsPhos::kSingleCellEnergyCut;
    fTrackDistanceCut = EtReconstructedCutsPhos::kTrackDistanceCut;
	 
}

bool AliAnalysisEtReconstructedPhos::TrackHitsCalorimeter(AliVParticle* track, Double_t magField)
{
  return  AliAnalysisEtReconstructed::TrackHitsCalorimeter(track, magField);
}

