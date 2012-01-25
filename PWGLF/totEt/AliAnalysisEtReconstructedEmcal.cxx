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
  fEtaCutAcc = fCuts->GetGeometryEmcalEtaAccCut();
  fPhiCutAccMax = fCuts->GetGeometryEmcalPhiAccMaxCut() * TMath::Pi()/180.;
  fPhiCutAccMin = fCuts->GetGeometryEmcalPhiAccMinCut() * TMath::Pi()/180.;
  fClusterEnergyCut = fCuts->GetReconstructedEmcalClusterEnergyCut();
  fSingleCellEnergyCut = fCuts->GetReconstructedEmcalSingleCellEnergyCut();

  fClusterType = fCuts->GetReconstructedEmcalClusterType();
  fTrackDistanceCut = fCuts->GetEmcalTrackDistanceCut();
  
  fDetector = fCuts->GetDetectorEmcal();
	 
}

bool AliAnalysisEtReconstructedEmcal::TrackHitsCalorimeter(AliVParticle* track, Double_t magField)
{
  return  AliAnalysisEtReconstructed::TrackHitsCalorimeter(track, magField);
}
