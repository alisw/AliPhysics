#include "AliPHOSTracker.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSPIDv1.h"
#include "AliRunLoader.h"
#include "AliESD.h"

//-------------------------------------------------------------------------
//                          PHOS tracker.
// Matches ESD tracks with the PHOS and makes the PID.  
// Currently, has only one function implemented : PropagateBack(AliESD*)
//-------------------------------------------------------------------------

ClassImp(AliPHOSTracker)

Bool_t AliPHOSTracker::fgDebug = kFALSE ; 

Int_t AliPHOSTracker::PropagateBack(AliESD *esd) {
  // Called by AliReconstruction 
  // Creates the tracksegments and Recparticles
  // Makes the PID
  
  Int_t eventNumber = fRunLoader->GetEventNumber() ;

  TString headerFile(fRunLoader->GetFileName()) ; 
  TString branchName(fRunLoader->GetEventFolder()->GetName()) ;  

  AliPHOSTrackSegmentMakerv1 tsm(headerFile, branchName);
  tsm.SetESD(esd) ; 
  AliPHOSPIDv1 pid(headerFile, branchName);

  // do current event; the loop over events is done by AliReconstruction::Run()
  tsm.SetEventRange(eventNumber, eventNumber) ; 
  pid.SetEventRange(eventNumber, eventNumber) ; 
  if ( Debug() ) {
   tsm.ExecuteTask("deb all") ;
   pid.ExecuteTask("deb all") ;
  }
  else {
    tsm.ExecuteTask("") ;
    pid.ExecuteTask("") ;
  }
  
  return 0;
}
