// AliAnalysisMultPbCentralitySelector 
// Interface class to centrality estimators for the PbPb
// track-multiplicity analysis
// Michele Floris, CERN

#include "AliAnalysisCentralitySelector.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliESDVZERO.h"
#include <iostream>
#include "AliMultiplicity.h"

using namespace std;



ClassImp(AliAnalysisCentralitySelector)

Bool_t AliAnalysisCentralitySelector::IsCentralityBinSelected(AliESDEvent* aEsd, AliESDtrackCuts * trackCuts) {

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
    AliESDVZERO* esdV0 = aEsd->GetVZEROData();
    Float_t multV0A=esdV0->GetMTotV0A();
    Float_t multV0C=esdV0->GetMTotV0C();
    multV0 = multV0A+multV0C;
    
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
  } 
  else if(fUsePercentile) {
    AliCentrality *centrality = (AliCentrality*) aEsd->GetCentrality(); 
    return centrality->IsEventInCentralityClass(fMultMin, fMultMax, fCentrEstimator.Data()) ;

  }
  
  else {

   AliCentrality *centrality = (AliCentrality*) aEsd->GetCentrality();
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

void AliAnalysisCentralitySelector::Print(Option_t* option ) const {
  // Print some information

  Printf("AliAnalysisCentralitySelector [%s]", option);
  Printf(" - Centrality estimator [%s]",fCentrEstimator.Data());
  Printf(" - Centrality bin       [%d]",fCentrBin);
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
