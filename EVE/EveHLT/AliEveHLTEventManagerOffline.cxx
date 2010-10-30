// Author: 2010 Svein Lindal <slindal@fys.uio.no>                        *
//         for The ALICE HLT Project.                                    *

#include "AliHLTEveHLT.h"
#include "AliHLTEvePhos.h"
#include "AliHLTEveEmcal.h"
#include "TEveManager.h"

#include "AliESDEvent.h"
#include "AliEveHLTEventManager.h"
#include "AliEveEventBufferOffline.h"
#include "AliEveHLTEventManagerOffline.h"


ClassImp(AliEveHLTEventManagerOffline)

AliEveHLTEventManagerOffline::AliEveHLTEventManagerOffline() : 
  AliEveHLTEventManager(),
  fEvent(NULL),
  fEventBuffer(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  
}

AliEveHLTEventManagerOffline::AliEveHLTEventManagerOffline(TString filename) : 
  AliEveHLTEventManager(),
  fEvent(NULL),
  fEventBuffer(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

  fEventBuffer = new AliEveEventBufferOffline(filename);
}
 
AliEveHLTEventManagerOffline::~AliEveHLTEventManagerOffline() {

  //DestroyElements();
  //DestroyDetectorElements();  

  if(fEvent)
    delete fEvent;
  fEvent = NULL;

  if(fEventBuffer)
    delete fEventBuffer;
  fEventBuffer = NULL;
  
}

void AliEveHLTEventManagerOffline::NextEvent() {
  //See header file for documentation
  fEvent = dynamic_cast<AliESDEvent*>(fEventBuffer->NextEvent());
  if(fEvent) {
    //Int_t eventId = fBuffer->GetEventId();
    ProcessEvent(fEvent);
  } else {
    cout << "couldn't get the event"<<endl;
  }
}

void AliEveHLTEventManagerOffline::NavigateFwd() {
  //See header file for documentation
  fEvent = dynamic_cast<AliESDEvent*>(fEventBuffer->Fwd());
  if(fEvent) {
    ProcessEvent(fEvent);
  } else {
    cout << "couldn't get the fwd event"<<endl;
  }
}

void AliEveHLTEventManagerOffline::NavigateBack() {
  //See header file for documentation
  fEvent = dynamic_cast<AliESDEvent*>(fEventBuffer->Back());
  if(fEvent) {
    ProcessEvent(fEvent);
  } else {
    cout << "couldn't get the back event"<<endl;
  }
}

