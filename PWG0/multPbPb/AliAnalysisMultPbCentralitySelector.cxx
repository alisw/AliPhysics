#include "AliAnalysisMultPbCentralitySelector.h"
#include "AliESDtrackCuts.h"
#include "AliESDCentrality.h"
#include "AliESDEvent.h"
#include "AliLog.h"

ClassImp(AliAnalysisMultPbCentralitySelector)

Bool_t AliAnalysisMultPbCentralitySelector::IsCentralityBinSelected(AliESDEvent* aEsd, AliESDtrackCuts * trackCuts) {

  // Centrality selection
  // On MC cuts on the number of good tracks,
  // On data cuts using AliESDCentrality and the cut requested in ntracks
  if(fIsMC) {
    if(!trackCuts){
      AliFatal("Track cuts object is invalid");
    }
    if (trackCuts->CountAcceptedTracks(aEsd) < fMultMin) return kFALSE;
    if (trackCuts->CountAcceptedTracks(aEsd) > fMultMax) return kFALSE;						       
  }

  AliESDCentrality *centrality = aEsd->GetCentrality();
  if(!centrality) {
    AliFatal("Centrality object not available"); 
  }
  else {
    Int_t centrBin = centrality->GetCentralityClass5(fCentrEstimator.Data()) ;    
    if (centrBin != fCentrBin && fCentrBin != -1) return kFALSE;
  }

  return kTRUE;

}

void AliAnalysisMultPbCentralitySelector::Print(Option_t* option ) const {
  // Print some information

  Printf("AliAnalysisMultPbCentralitySelector");
  Printf(" - Centrality estimator [%s]",fCentrEstimator.Data());

  if ( fIsMC ) {    
    Printf("Running on Monte Carlo, actual cut was on tracks multiplicity [%d - %d]",fMultMin,fMultMax);    
  } 
  
}
