// AliAnalysisMultPbCentralitySelector 
// Interface class to centrality estimators for the PbPb
// track-multiplicity analysis
// Michele Floris, CERN

#include "AliAnalysisMultPbCentralitySelector.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliESDVZERO.h"
#include <iostream>
#include "AliMultiplicity.h"

using namespace std;



ClassImp(AliAnalysisMultPbCentralitySelector)

Bool_t AliAnalysisMultPbCentralitySelector::IsCentralityBinSelected(AliESDEvent* aEsd, AliESDtrackCuts * trackCuts) {

  // Centrality selection
  // On MC cuts on the number of good tracks,
  // On data cuts using AliESDCentrality and the cut requested in ntracks

  //  cout << "Tracks " << trackCuts->CountAcceptedTracks(aEsd) << endl;
  ///  cout << "CENTRALITY " << fUseV0CutRange << " " << fUseMultRange << " " << fMultMin << " " << fMultMax << endl;
  
  if (fUseMultRange && fUseV0CutRange && fUseSPDOuterRange) {
    AliFatal(Form("Cannot use multiple estimators at once: fUseMultRange [%d], fUseV0CutRange[%d], fUseSPDOuterRange[%d]!!",
		  fUseMultRange , fUseV0CutRange , fUseSPDOuterRange)); 
  }

  if (fUseV0CutRange) {

    Float_t multV0=0;
    Float_t dummy=0;
    if (fUseCorrV0) {
      multV0 = GetCorrV0(aEsd,dummy);      
    }
    else { 
      AliESDVZERO* esdV0 = aEsd->GetVZEROData();
      Float_t multV0A=esdV0->GetMTotV0A();
      Float_t multV0C=esdV0->GetMTotV0C();
      multV0 = multV0A+multV0C;
    }

    if (fIsMC) multV0 = 0.85871 * multV0;
    if (multV0 < fMultMin) return kFALSE;
    if (multV0 > fMultMax) return kFALSE;
    //    cout << "ok" << endl;

  } 
  else if (fUseSPDOuterRange) {

    const AliMultiplicity * mult = aEsd->GetMultiplicity();
    Float_t outerLayerSPD = mult->GetNumberOfITSClusters(1);  
    
    if (outerLayerSPD < fMultMin) return kFALSE;
    if (outerLayerSPD > fMultMax) return kFALSE;
    //    cout << "ok" << endl;

  }
 else if(fUseMultRange) {
    if(!trackCuts){
      AliFatal("Track cuts object is invalid");
    }
    Float_t ntracks = trackCuts->CountAcceptedTracks(aEsd);
    //    cout << "Hey! " << fCentrBin << " " << ntracks << " " << fMultMin <<" - " << fMultMax << endl;
    
    if (fCentrBin == -1 && !fUseMultRange) return kTRUE;
    if (ntracks < fMultMin) return kFALSE;
    if (ntracks > fMultMax) return kFALSE;						       
  } else {

    AliCentrality *centrality = aEsd->GetCentrality();
    if(!centrality && !fUseMultRange) {
      AliFatal("Centrality object not available"); 
    }
    else {
      Int_t centrBin = centrality->GetCentralityClass5(fCentrEstimator.Data()) ;    
      if (centrBin != fCentrBin && fCentrBin != -1 && !fUseMultRange) return kFALSE;
    }
  }

  //  cout << "Selected" << endl;
  

  return kTRUE;

}

void AliAnalysisMultPbCentralitySelector::Print(Option_t* option ) const {
  // Print some information

  Printf("AliAnalysisMultPbCentralitySelector [%s]", option);
  Printf(" - Centrality estimator [%s]",fCentrEstimator.Data());
  Printf(" - Centrality bin       [%d]",fCentrBin);
  Printf(" - File1 used for centrality estimate: [%s]", fFile1.Data());
  Printf(" - File2 used for centrality estimate: [%s]", fFile2.Data());
  if ( fUseMultRange ) {
    Printf ("Using multiplicity range [%1.1f - %1.1f]",fMultMin,fMultMax);
  }
  if ( fUseV0CutRange ) {
    Printf ("Using V0 range [%1.1f - %1.1f]",fMultMin,fMultMax);
  }
  if ( fIsMC ) {    
    Printf("Running on Monte Carlo, actual cut was on tracks multiplicity [%1.1f - %1.1f]",fMultMin,fMultMax);    
  } 
  
}
//________________________________________________________________________
Float_t AliAnalysisMultPbCentralitySelector::GetCorrV0(const AliESDEvent* esd, float &v0CorrResc) const 
{
  // correct V0 non-linearity, prepare a version rescaled to SPD2 corr
  Double_t *par0;
  Double_t *par1;
  Double_t *par2;
  
  Double_t par0_137161[64] = { 6.71e-02 , 6.86e-02 , 7.06e-02 , 6.32e-02 , 
			5.91e-02 , 6.07e-02 , 5.78e-02 , 5.73e-02 , 5.91e-02 , 6.22e-02 , 
			5.90e-02 , 6.11e-02 , 5.55e-02 , 5.29e-02 , 5.19e-02 , 5.56e-02 , 
			6.25e-02 , 7.03e-02 , 5.64e-02 , 5.81e-02 , 4.57e-02 , 5.30e-02 , 
			5.13e-02 , 6.43e-02 , 6.27e-02 , 6.48e-02 , 6.07e-02 , 1.01e-01 , 
			6.68e-02 , 7.16e-02 , 6.36e-02 , 5.95e-02 , 2.52e-02 , 2.82e-02 , 
			2.56e-02 , 2.86e-02 , 2.82e-02 , 2.10e-02 , 2.13e-02 , 2.32e-02 , 
			2.75e-02 , 4.34e-02 , 3.78e-02 , 4.52e-02 , 4.11e-02 , 3.89e-02 , 
			4.10e-02 , 3.73e-02 , 4.51e-02 , 5.07e-02 , 5.42e-02 , 4.74e-02 , 
			4.33e-02 , 4.44e-02 , 4.64e-02 , 3.01e-02 , 6.38e-02 , 5.26e-02 , 
			4.99e-02 , 5.26e-02 , 5.47e-02 , 3.84e-02 , 5.00e-02 , 5.20e-02 };
  Double_t par1_137161[64] = { -6.68e-05 , -7.78e-05 , -6.88e-05 , -5.92e-05 , 
			-2.43e-05 , -3.54e-05 , -2.91e-05 , -1.99e-05 , -1.40e-05 , -4.01e-05 , 
			-2.29e-05 , -3.68e-05 , -2.53e-05 , -2.44e-06 , -9.22e-06 , -1.51e-05 , 
			-2.80e-05 , -2.34e-05 , -1.72e-05 , -1.81e-05 , -1.29e-05 , -2.65e-05 , 
			-1.61e-05 , -2.86e-05 , -1.74e-05 , -4.23e-05 , -3.41e-05 , -1.05e-04 , 
			-2.76e-05 , -4.71e-05 , -3.06e-05 , -2.32e-05 , -1.55e-06 , 2.15e-05 , 
			1.40e-05 , 2.16e-05 , 1.21e-05 , 3.05e-06 , 1.67e-05 , -3.84e-06 , 
			3.09e-06 , 1.50e-05 , 3.47e-06 , 4.87e-06 , -3.71e-07 , -1.75e-06 , 
			-1.80e-06 , 9.99e-06 , -6.46e-06 , -4.91e-06 , 1.33e-05 , -2.52e-07 , 
			-3.85e-06 , 4.94e-06 , -2.48e-07 , -1.20e-05 , 2.07e-06 , 6.12e-06 , 
			-1.18e-06 , 4.54e-06 , -1.54e-05 , -1.25e-05 , 1.46e-06 , -6.67e-06 };
  Double_t par2_137161[64] = { 1.29e-08 , 1.51e-08 , 1.43e-08 , 1.11e-08 , 
			5.04e-09 , 6.99e-09 , 5.58e-09 , 4.15e-09 , 4.00e-09 , 8.22e-09 , 
			4.97e-09 , 7.66e-09 , 4.91e-09 , 1.10e-09 , 2.64e-09 , 3.64e-09 , 
			5.76e-09 , 5.46e-09 , 3.38e-09 , 3.47e-09 , 2.43e-09 , 4.13e-09 , 
			2.80e-09 , 5.80e-09 , 3.86e-09 , 7.46e-09 , 5.98e-09 , 2.58e-08 , 
			5.50e-09 , 8.72e-09 , 5.23e-09 , 4.37e-09 , 2.33e-09 , -6.01e-10 , 
			3.99e-11 , -2.02e-10 , 7.67e-10 , 2.03e-09 , 1.17e-10 , 2.56e-09 , 
			1.16e-09 , -4.75e-10 , 1.28e-09 , 1.23e-09 , 1.62e-09 , 1.61e-09 , 
			1.93e-09 , 2.97e-10 , 2.21e-09 , 2.16e-09 , 5.22e-10 , 1.03e-09 , 
			1.56e-09 , 5.00e-10 , 1.01e-09 , 2.93e-09 , 1.05e-09 , 9.96e-11 , 
			1.21e-09 , 7.45e-10 , 3.07e-09 , 2.31e-09 , 6.70e-10 , 1.89e-09 };

  Double_t par0_137366[64] = { 7.12e-02 , 7.34e-02 , 7.39e-02 , 6.54e-02 , 6.11e-02 , 6.31e-02 , 6.15e-02 , 
			       6.00e-02 , 6.10e-02 , 6.49e-02 , 6.17e-02 , 6.33e-02 , 6.00e-02 , 5.48e-02 , 
			       5.44e-02 , 5.81e-02 , 6.49e-02 , 7.07e-02 , 5.91e-02 , 6.18e-02 , 4.82e-02 , 
			       5.67e-02 , 5.36e-02 , 6.60e-02 , 6.37e-02 , 6.78e-02 , 6.31e-02 , 1.04e-01 , 
			       6.91e-02 , 7.32e-02 , 6.61e-02 , 6.16e-02 , 2.64e-02 , 2.81e-02 , 2.64e-02 , 
			       2.85e-02 , 2.87e-02 , 2.18e-02 , 2.19e-02 , 2.43e-02 , 2.81e-02 , 4.37e-02 , 
			       3.90e-02 , 4.66e-02 , 4.24e-02 , 4.09e-02 , 4.21e-02 , 3.88e-02 , 4.83e-02 , 
			       5.23e-02 , 5.44e-02 , 4.85e-02 , 4.42e-02 , 4.58e-02 , 4.74e-02 , 3.14e-02 , 
			       6.31e-02 , 5.30e-02 , 5.01e-02 , 5.33e-02 , 5.70e-02 , 3.95e-02 , 4.98e-02 , 5.31e-02 };
  Double_t par1_137366[64] = { -6.99e-05 , -6.99e-05 , -6.94e-05 , -6.55e-05 , -3.55e-05 , -4.50e-05 , 
			       -3.10e-05 , -2.81e-05 , -2.29e-05 , -3.89e-05 , -2.53e-05 , -4.25e-05 ,
			       -1.87e-05 , -2.01e-05 , -1.53e-05 , -2.14e-05 , -2.86e-05 , -4.70e-05 ,
			       -2.23e-05 , -3.30e-05 ,-9.74e-06 , -2.62e-05 , -1.76e-05 , -2.38e-05 , 
			       -2.40e-05 , -3.43e-05 , -2.75e-05 , -6.86e-05 ,-2.35e-05 , -4.45e-05 , 
			       -2.51e-05 , -2.20e-05 , -1.25e-16 , -2.04e-17 , -2.06e-17 , -3.74e-19 ,
			       -1.18e-18 , -2.02e-15 , -3.78e-06 , -1.26e-06 , -2.71e-06 , -6.23e-17 , 
			       -7.39e-08 , -1.76e-16 , -8.98e-06 , -4.10e-18 , -1.34e-05 , -1.06e-16 , 
			       -3.34e-06 , -1.04e-05 , -5.28e-06 , -7.34e-06 , -1.05e-05 , -7.68e-06 ,
			       -1.78e-05 , -1.19e-05 , -1.78e-05 , -1.34e-06 , -9.23e-06 , -3.34e-06 ,
			       -8.02e-06 , -1.39e-05 , -1.38e-05 , -1.40e-05 };
  Double_t par2_137366[64] = { 1.41e-08 , 1.47e-08 , 1.48e-08 , 1.24e-08 , 6.82e-09 , 8.73e-09 , 6.26e-09 , 
			       5.53e-09 , 5.40e-09 , 7.93e-09 , 5.49e-09 , 8.77e-09 , 4.21e-09 , 3.93e-09 , 
			       3.60e-09 , 4.67e-09 , 5.59e-09 , 8.81e-09 , 3.89e-09 , 6.19e-09 , 1.97e-09 , 
			       4.38e-09 , 3.26e-09 , 5.00e-09 , 4.58e-09 , 6.39e-09 , 5.03e-09 , 1.30e-08 , 
			       4.95e-09 , 8.26e-09 , 4.57e-09 , 4.10e-09 , 2.35e-09 , 2.30e-09 , 2.15e-09 , 
			       2.27e-09 , 2.17e-09 , 2.27e-09 , 2.97e-09 , 2.25e-09 , 1.69e-09 , 1.44e-09 , 
			       1.66e-09 , 1.75e-09 , 2.88e-09 , 1.82e-09 , 3.64e-09 , 1.80e-09 , 1.71e-09 , 
			       2.66e-09 , 3.01e-09 , 1.95e-09 , 2.64e-09 , 2.42e-09 , 3.68e-09 , 2.66e-09 , 
			       3.92e-09 , 1.18e-09 , 2.26e-09 , 1.57e-09 , 2.02e-09 , 2.71e-09 , 2.99e-09 , 3.04e-09 }; 
  
  
  if (esd->GetRunNumber() <= 137165) {
    par0=par0_137161;
    par1=par1_137161;
    par2=par2_137161;
  }  else  {
    par0=par0_137366;
    par1=par1_137366;
    par2=par2_137366;
 }
  //
  Float_t multCorr = 0;
  Float_t multCorr2 = 0;
  Float_t multChCorr[64];
  AliESDVZERO* esdV0 = esd->GetVZEROData();
  for(Int_t i = 0; i < 64; ++i) {
    Double_t b = (esdV0->GetMultiplicity(i)*par1[i]-par0[i]);
    Double_t s = (b*b-4.*par2[i]*esdV0->GetMultiplicity(i)*esdV0->GetMultiplicity(i));
    Double_t n = (s<0) ? -b : (-b + TMath::Sqrt(s));
    multChCorr[i] = 2.*esdV0->GetMultiplicity(i)/n*par0[i];
    multCorr += multChCorr[i];
    multCorr2 += (multChCorr[i]/par0[i]/64.);
  }
  v0CorrResc =  multCorr2;
  return multCorr;
}
