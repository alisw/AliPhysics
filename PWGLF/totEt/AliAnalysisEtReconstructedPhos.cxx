//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD analysis, for PHOS
//  - reconstruction output
// implementation file
//
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________
#include "AliAnalysisEtReconstructedPhos.h"
#include "AliAnalysisEtSelectorPhos.h"
#include "AliAnalysisEtCuts.h"
#include "AliESDtrack.h"
#include <fstream>
#include <iostream>
#include <AliPHOSGeoUtils.h>
#include "AliPHOSGeometry.h"
#include <vector>
#include "TH2I.h"

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
//corrEnergy =cluster->E()/(0.31 + 0.02*cluster->E());
const Double_t kMEANCHARGED = 0.48/(0.31 + 0.00*0.48);
const Double_t kMEANNEUTRAL = 0.53/(0.31 + 0.00*0.53);
const Double_t kMEANGAMMA = 0.51/(0.31 + 0.00*0.31);


AliAnalysisEtReconstructedPhos::AliAnalysisEtReconstructedPhos() :
        AliAnalysisEtReconstructed()
{ // ctor
    fHistogramNameSuffix = TString("PhosRec");

//     fChargedContributionCorrectionParameters[0] = -0.017;
//     fChargedContributionCorrectionParameters[1] = 0.065;
// 
//     fNeutralContributionCorrectionParameters[0] = -0.002;
//     fNeutralContributionCorrectionParameters[1] = 0.017;
//     fNeutralContributionCorrectionParameters[2] = -3.6e-5;
// 
//     fRemovedGammaContributionCorrectionParameters[0] = 0.001;
//     fRemovedGammaContributionCorrectionParameters[1] = 0.37e-5;
//     fRemovedGammaContributionCorrectionParameters[2] = 0.0003;

    fChargedContributionCorrectionParameters[0] = -0.06;
    fChargedContributionCorrectionParameters[1] = 0.316;
    fChargedContributionCorrectionParameters[2] = 0.0022;

    fNeutralContributionCorrectionParameters[0] = -0.003;
    fNeutralContributionCorrectionParameters[1] = 0.232;
    fNeutralContributionCorrectionParameters[2] = 0.002;

    fRemovedGammaContributionCorrectionParameters[0] = 0.001;
    fRemovedGammaContributionCorrectionParameters[1] = 0.009;
    fRemovedGammaContributionCorrectionParameters[2] = 0.0;

    fSecondaryContributionCorrectionParameters[0] = -0.03;
    fSecondaryContributionCorrectionParameters[1] = 0.221;
    fSecondaryContributionCorrectionParameters[2] = 0.002;

}

AliAnalysisEtReconstructedPhos::~AliAnalysisEtReconstructedPhos()
{
}

void AliAnalysisEtReconstructedPhos::Init()
{ // Init
    AliAnalysisEtReconstructed::Init();

    
    fDetectorRadius = fCuts->GetGeometryPhosDetectorRadius();
    fSingleCellEnergyCut = fCuts->GetReconstructedPhosSingleCellEnergyCut();
    fSelector = new AliAnalysisEtSelectorPhos(fCuts);

}

bool AliAnalysisEtReconstructedPhos::TrackHitsCalorimeter(AliVParticle* track, Double_t magField)
{ 
    return  AliAnalysisEtReconstructed::TrackHitsCalorimeter(track, magField);
}



void AliAnalysisEtReconstructedPhos::CreateHistograms()
{ // add some extra histograms & objects to the ones from base class
  if(!fSelector){
    cout<<__FILE__<<" "<<"Creating new fSelector"<<endl;
    fSelector = new AliAnalysisEtSelectorPhos(fCuts);
  }
  AliAnalysisEtReconstructed::CreateHistograms();
}

Double_t AliAnalysisEtReconstructedPhos::GetCorrectionModification(const AliESDCaloCluster& cluster,Int_t nonLinCorr, Int_t effCorr){//nonLinCorr 0 = nominal 1 = high -1 = low, effCorr  0 = nominal 1 = high -1 = low
  Double_t factor = 1.0;
  //Double_t E = fReCorrections->CorrectedEnergy(cluster.E());
  if(nonLinCorr!=0){
    if(nonLinCorr==1){//high bound on nonlinearity
      factor *=1.05;
    }
    else{//nonLinCorr==-1
      factor *=0.95;
    }
  }
  if(effCorr!=0){
    if(effCorr==1){//high bound
      factor *=1.05;
    }
    else{//low bound
      factor *=0.95;
    }
  }
  return factor;
}
