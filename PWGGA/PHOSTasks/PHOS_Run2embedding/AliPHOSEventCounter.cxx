#include "TChain.h"
#include "AliLog.h"
#include "AliPHOSEventCounter.h"

ClassImp(AliPHOSEventCounter)

//________________________________________________________________________
AliPHOSEventCounter::AliPHOSEventCounter(const char *name) 
  : AliAnalysisTaskSE(name),
		fEventCounter(0)
{
  // Constructor

}

//________________________________________________________________________
void AliPHOSEventCounter::UserCreateOutputObjects()
{

}
//________________________________________________________________________
void AliPHOSEventCounter::UserExec(Option_t *)
{
  // Main loop, called for each event
  
//  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());//for general information before embedding
//  if (!event) {
//    Printf("ERROR: Could not retrieve event");
//    PostData(0, fTreeOut);
//    return;
//  }


  fEventCounter++;


}
//________________________________________________________________________
void AliPHOSEventCounter::Terminate(Option_t *)
{
  AliInfo(Form("The number of analyszed events is %d",fEventCounter));

}
//________________________________________________________________________
