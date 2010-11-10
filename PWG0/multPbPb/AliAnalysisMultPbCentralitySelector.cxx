// AliAnalysisMultPbCentralitySelector 
// Interface class to centrality estimators for the PbPb
// track-multiplicity analysis
// Michele Floris, CERN

#include "AliAnalysisMultPbCentralitySelector.h"
#include "AliESDtrackCuts.h"
#include "AliESDCentrality.h"
#include "AliESDEvent.h"
#include "AliLog.h"
#include "AliESDVZERO.h"
#include <iostream>

using namespace std;



ClassImp(AliAnalysisMultPbCentralitySelector)

Bool_t AliAnalysisMultPbCentralitySelector::IsCentralityBinSelected(AliESDEvent* aEsd, AliESDtrackCuts * trackCuts) {

  // Centrality selection
  // On MC cuts on the number of good tracks,
  // On data cuts using AliESDCentrality and the cut requested in ntracks

  //  cout << "Tracks " << trackCuts->CountAcceptedTracks(aEsd) << endl;
  ///  cout << "CENTRALITY " << fUseV0CutRange << " " << fUseMultRange << " " << fMultMin << " " << fMultMax << endl;
  

  if (fUseV0CutRange) {

    AliESDVZERO* esdV0 = aEsd->GetVZEROData();
    Float_t multV0A=esdV0->GetMTotV0A();
    Float_t multV0C=esdV0->GetMTotV0C();
    Float_t multV0 = multV0A+multV0C;
    //    cout << "V0 Mult: " << multV0 << " " << fMultMin << " " << fMultMax << endl;
    
    if (multV0 < fMultMin) return kFALSE;
    if (multV0 > fMultMax) return kFALSE;
    //    cout << "ok" << endl;

  }
  else if(fIsMC || fUseMultRange) {
    if(!trackCuts){
      AliFatal("Track cuts object is invalid");
    }
    //    cout << "Hey! " << fCentrBin << " " << fMultMin <<" - " << fMultMax << endl;
    
    if (fCentrBin == -1) return kTRUE;
    if (trackCuts->CountAcceptedTracks(aEsd) < fMultMin) return kFALSE;
    if (trackCuts->CountAcceptedTracks(aEsd) > fMultMax) return kFALSE;						       
  } else {

    AliESDCentrality *centrality = aEsd->GetCentrality();
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
