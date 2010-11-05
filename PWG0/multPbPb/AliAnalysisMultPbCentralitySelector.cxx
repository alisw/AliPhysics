// AliAnalysisMultPbCentralitySelector 
// Interface class to centrality estimators for the PbPb
// track-multiplicity analysis
// Michele Floris, CERN

#include "AliAnalysisMultPbCentralitySelector.h"
#include "AliESDtrackCuts.h"
#include "AliESDCentrality.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include <iostream>

using namespace std;


// FIXME: bookkeep here all parameters of centrality estimate (files, estimator, selected bin...)

ClassImp(AliAnalysisMultPbCentralitySelector)

Bool_t AliAnalysisMultPbCentralitySelector::IsCentralityBinSelected(AliESDEvent* aEsd, AliESDtrackCuts * trackCuts) {

  // Centrality selection
  // On MC cuts on the number of good tracks,
  // On data cuts using AliESDCentrality and the cut requested in ntracks

  //  cout << "Tracks " << trackCuts->CountAcceptedTracks(aEsd) << endl;
  

  if(fIsMC || fUseMultRange) {
    if(!trackCuts){
      AliFatal("Track cuts object is invalid");
    }
    //    cout << "Hey!" << endl;
    
    if (fCentrBin == -1) return kTRUE;
    if (trackCuts->CountAcceptedTracks(aEsd) < fMultMin) return kFALSE;
    if (trackCuts->CountAcceptedTracks(aEsd) > fMultMax) return kFALSE;						       
  }

  AliESDCentrality *centrality = aEsd->GetCentrality();
  if(!centrality && !fUseMultRange) {
    AliFatal("Centrality object not available"); 
  }
  else {
    Int_t centrBin = centrality->GetCentralityClass5(fCentrEstimator.Data()) ;    
    if (centrBin != fCentrBin && fCentrBin != -1 && !fUseMultRange) return kFALSE;
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
    Printf ("Using multiplicity range [%d - %d]",fMultMin,fMultMax);
  }
  if ( fIsMC ) {    
    Printf("Running on Monte Carlo, actual cut was on tracks multiplicity [%d - %d]",fMultMin,fMultMax);    
  } 
  
}
