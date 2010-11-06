// Author: 2010 Svein Lindal <slindal@fys.uio.no>                        *
//         for The ALICE HLT Project.                                    *

#include "AliHLTEveHLT.h"
#include "AliHLTEvePhos.h"
#include "AliHLTEveEmcal.h"


#include "AliESDEvent.h"
#include "AliEveHLTEventManager.h"
#include "AliEveEventBufferOffline.h"
#include "AliEveHLTEventManagerHomer.h"
#include "TList.h"
#include "AliEveHOMERManager.h"
#include "TEveManager.h"

ClassImp(AliEveHLTEventManagerHomer)

AliEveHLTEventManagerHomer::AliEveHLTEventManagerHomer() : 
  AliEveHLTEventManager(),
  fEventBuffer(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  fEventBuffer = new AliEveEventBufferHomer();
  fEventBuffer->StartBufferMonitor();
  //AliEveHOMERManager * hm = fEventBuffer->GetHomerManager();
  //if(hm) GetEveManager()->AddToListTree(hm, kTRUE);
  

}

 
AliEveHLTEventManagerHomer::~AliEveHLTEventManagerHomer() {

  //DestroyElements();
  //DestroyDetectorElements();  
  if(fEventBuffer)
    delete fEventBuffer;
  fEventBuffer = NULL;

}

///________________________________________________________________________________
void AliEveHLTEventManagerHomer::ProcessList(TList * blockList) {

  ProcessEvent(blockList);
  UpdateDisplay();

}



///________________________________________________________________________________
void AliEveHLTEventManagerHomer::NextEvent() {
  //See header file for documentation
  if(fEventBuffer->GetBusy() ) {
    cout << "event buffer already busy"<<endl;
    return;
  }else {
    fEventBuffer->SetBusy(kTRUE);
  }

  TList * aSyncEvent = fEventBuffer->GetASyncEvent();
  TList * event = static_cast<TList*>(fEventBuffer->NextEvent());
  if(event) {
    cout << "Got the event, reset the display " <<endl;
    ResetDisplay();
    cout << "Process event"<<endl;
    ProcessEvent(event);
    if(aSyncEvent) {
      cout  << "Process asynchroneous event" << endl;
      ProcessEvent(aSyncEvent);
    }  else {
      cout << "Could not get async event"<<endl;
    }
    
    cout << "Upate the display"<<endl;
    UpdateDisplay();
  
  } else {
    cout << "couldn't get the sync event"<<endl;
  }
  
  //  cout  << "doint async block"<<endl;
  //  TList * async = static_cast<TList*>(fEventBuffer->GetAList());
  //  if(async) {
  // 	ProcessEvent(async);
  //   }  else {
  // 	 cout << "No async bloc"<<endl;
  // }

  fEventBuffer->SetBusy(kFALSE);
}


///____________________________________________________________________________________
void AliEveHLTEventManagerHomer::NavigateFwd() {
  //See header file for documentation
  TList * fEvent = dynamic_cast<TList*>(fEventBuffer->Fwd());
  if(fEvent) {
    ResetDisplay();
    ProcessEvent(fEvent);
    UpdateDisplay();
  } else {
    cout << "couldn't get the fwd event"<<endl;
  }
}


///______________________________________________________________________________________
void AliEveHLTEventManagerHomer::NavigateBack() {
  //See header file for documentation
  TList * fEvent = dynamic_cast<TList*>(fEventBuffer->Back());
  if(fEvent) {
    ResetDisplay();
    ProcessEvent(fEvent);
    UpdateDisplay();
  } else {
    cout << "couldn't get the back event"<<endl;
  }
}

