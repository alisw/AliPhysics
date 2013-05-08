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
#include "AliAnalysisEtSelectorEmcal.h"
#include "AliESDtrack.h"
#include "AliAnalysisEtRecEffCorrection.h"

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

  fSelector = new AliAnalysisEtSelectorEmcal(fCuts);
}

bool AliAnalysisEtReconstructedEmcal::TrackHitsCalorimeter(AliVParticle* track, Double_t magField)
{
  return  AliAnalysisEtReconstructed::TrackHitsCalorimeter(track, magField);
}

void AliAnalysisEtReconstructedEmcal::CreateHistograms()
{ // add some extra histograms & objects to the ones from base class
  if(!fSelector){
    cout<<__FILE__<<" "<<"Creating new fSelector"<<endl;
    fSelector = new AliAnalysisEtSelectorEmcal(fCuts);
  }
  AliAnalysisEtReconstructed::CreateHistograms();
}


Double_t AliAnalysisEtReconstructedEmcal::GetCorrectionModification(const AliESDCaloCluster& cluster,Int_t nonLinCorr, Int_t effCorr){//nonLinCorr 0 = nominal 1 = high -1 = low, effCorr  0 = nominal 1 = high -1 = low
  Double_t factor = 1.0;
  Double_t E = fReCorrections->CorrectedEnergy(cluster.E());
  if(nonLinCorr!=0){
    Double_t p0 = 9.90780e-01;
    Double_t p1 = 1.61503e-01;
    Double_t p2 = 6.55150e-01;
    Double_t p3 = 1.34100e-01;
    Double_t p4 = 1.63282e+02;
    Double_t p5 = 2.36902e+01;
    Double_t nominal = p0*(1./(1.+p1*TMath::Exp(-E/p2))*1./(1.+p3*TMath::Exp((E-p4)/p5)));
    Double_t alt = 1;
    if(nonLinCorr==1){//high bound on nonlinearity
      p0 = .984;
      p1 = 0.42;
      p2 = 0.35;
      alt = (p0)*1./(1.+p1*TMath::Exp(-E/p2));
    }
    else{//nonLinCorr==-1
      p0 = .992;
      p1 =0.115;
      p2 = 0.68;
      alt = (p0)*1./(1.+p1*TMath::Exp(-E/p2));
    }
    factor *=alt/nominal;
  }
  if(effCorr!=0){
    if(effCorr==1){//high bound
      factor *=1.02;
    }
    else{//low bound
      factor *=0.98;
    }
  }
  //cout<<"Factor:  "<<factor<<endl;
  return factor;
}
