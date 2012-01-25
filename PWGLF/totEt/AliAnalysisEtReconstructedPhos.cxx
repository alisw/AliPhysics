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

/*// Worst case (protons and neutrons):
const Double_t kMEANCHARGED = 0.335;
const Double_t kMEANNEUTRAL = 0.434;
const Double_t kMEANGAMMA = 0.374;

// Best case (pions and K0s): 
const Double_t kMEANCHARGED = 0.304;
const Double_t kMEANNEUTRAL = 0.3356;
const Double_t kMEANGAMMA = 0.374;
*/
// Simulated case:
const Double_t kMEANCHARGED = 0.307;
const Double_t kMEANNEUTRAL = 0.407;
const Double_t kMEANGAMMA = 0.374;


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
  fTrackDistanceCut = fCuts->GetPhosTrackDistanceCut();
  fTrackDxCut = fCuts->GetPhosTrackDxCut();
  fTrackDzCut = fCuts->GetPhosTrackDzCut();
  
  fDetector = fCuts->GetDetectorPhos();
  
  fGeomCorrection = 1.0/0.036;
  
  fEMinCorrection = 1.0;

}

bool AliAnalysisEtReconstructedPhos::TrackHitsCalorimeter(AliVParticle* track, Double_t magField)
{
  return  AliAnalysisEtReconstructed::TrackHitsCalorimeter(track, magField);
}

Double_t AliAnalysisEtReconstructedPhos::GetChargedContribution(Int_t clusterMult)
{ // Charged contrib
  if(clusterMult > 0)
  {
    Double_t nPart = 0.067 + 0.137*clusterMult;
  
    Double_t contr = nPart*kMEANCHARGED;
  
    return contr;
  }
  return 0;
  
}

Double_t AliAnalysisEtReconstructedPhos::GetNeutralContribution(Int_t clusterMult)
{ // Neutral contrib
  if(clusterMult > 0)
  {
    Double_t nPart = 0.012 + 0.024*clusterMult - 0.00006*clusterMult*clusterMult;
  
    Double_t contr = nPart*kMEANNEUTRAL;
  
    return contr;
  }
  return 0;
}

Double_t AliAnalysisEtReconstructedPhos::GetGammaContribution(Int_t clusterMult)
{ // Gamma contrib
  if(clusterMult > 0)
  {
    Double_t nPart = -0.008 + 0.0057*clusterMult + 0.0002*clusterMult*clusterMult;
  
    Double_t contr = nPart*kMEANGAMMA;
  
    return contr;
  }
  return 0;
}






