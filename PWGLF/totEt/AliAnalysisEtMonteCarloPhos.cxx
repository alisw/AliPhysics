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

using namespace std;

ClassImp(AliAnalysisEtMonteCarloPhos);


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
    
  fDetectorRadius = fCuts->GetGeometryPhosDetectorRadius();
  fEtaCutAcc = fCuts->GetGeometryPhosEtaAccCut();
  fPhiCutAccMax = fCuts->GetGeometryPhosPhiAccMaxCut() * TMath::Pi()/180.;
  fPhiCutAccMin = fCuts->GetGeometryPhosPhiAccMinCut() * TMath::Pi()/180.;
  fClusterEnergyCut = fCuts->GetReconstructedPhosClusterEnergyCut();
  fSingleCellEnergyCut = fCuts->GetReconstructedPhosSingleCellEnergyCut();
  
  fTrackDistanceCut = fCuts->GetPhosTrackDistanceCut();
  fTrackDxCut = fCuts->GetPhosTrackDxCut();
  fTrackDzCut = fCuts->GetPhosTrackDzCut();
    
  fDetector = fCuts->GetDetectorPhos();
    
}
