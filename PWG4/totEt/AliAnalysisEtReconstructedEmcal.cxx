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
    
    fClusterType = EtReconstructedCutsEmcal::kClusterType;
    fDetectorRadius = EtGeometryCutsEmcal::kDetectorRadius;    
    fEtaCutAcc = EtGeometryCutsEmcal::kEtaAccCut;
    fPhiCutAccMax = EtGeometryCutsEmcal::kPhiAccMaxCut*TMath::Pi()/180.;
    fPhiCutAccMin = EtGeometryCutsEmcal::kPhiAccMinCut*TMath::Pi()/180.;
    fClusterEnergyCut = EtReconstructedCutsEmcal::kClusterEnergyCut;
    fSingleCellEnergyCut = EtReconstructedCutsEmcal::kSingleCellEnergyCut;
    fTrackDistanceCut = EtReconstructedCutsEmcal::kTrackDistanceCut;
	 
}

bool AliAnalysisEtReconstructedEmcal::TrackHitsCalorimeter(AliVParticle* track, Double_t magField)
{
  return  AliAnalysisEtReconstructed::TrackHitsCalorimeter(track, magField);
}
