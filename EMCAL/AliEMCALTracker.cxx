#include "AliEMCALTracker.h"
#include "AliEMCALPIDv1.h"
#include "AliRunLoader.h"
#include "AliESD.h"

//-------------------------------------------------------------------------
//                          EMCAL tracker.
// Matches ESD tracks with the EMCAL and makes the PID.  
// Currently, has only one function implemented : PropagateBack(AliESD*)
//-------------------------------------------------------------------------

ClassImp(AliEMCALTracker)

Bool_t AliEMCALTracker::fgDebug = kFALSE ; 

Int_t AliEMCALTracker::PropagateBack(AliESD *esd) {
  // Makes the Particle Identification

  esd=0; // This is to avoid a compilation warning. 
         // This pointer is reserved for future needs

  Int_t eventNumber = fRunLoader->GetEventNumber() ;

  TString headerFile(fRunLoader->GetFileName()) ; 
  TString branchName(fRunLoader->GetEventFolder()->GetName()) ;  

  AliEMCALPIDv1 pid(headerFile, branchName);

  // do current event; the loop over events is done by AliReconstruction::Run()
  pid.SetEventRange(eventNumber, eventNumber) ; 
  if ( Debug() ) 
    pid.ExecuteTask("deb all") ;
  else 
    pid.ExecuteTask("") ;
 
  return 0;
}
