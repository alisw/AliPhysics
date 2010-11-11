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
#include "TTimer.h"
#include "TGLOverlayButton.h"

ClassImp(AliEveHLTEventManagerHomer)

AliEveHLTEventManagerHomer::AliEveHLTEventManagerHomer() : 
  AliEveHLTEventManager(),
  fEventBuffer(NULL), 
  fNextEventTimer(NULL),
  fInfoButton(NULL)
{
  // see header file for class documentation

  fEventBuffer = new AliEveEventBufferHomer();
  fEventBuffer->StartBufferMonitor();
    
  fNextEventTimer = new TTimer();
  fNextEventTimer->Connect("Timeout()", "AliEveHLTEventManagerHomer", this, "TryNextEvent()" );
   
  

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
void AliEveHLTEventManagerHomer::NextEvent() {
  fNextEventTimer->Start(1000);
}

///________________________________________________________________________________
void AliEveHLTEventManagerHomer::TryNextEvent() {
  //See header file for documentation
    
  if ( fEventBuffer->LockMutex() ) {
    cout << "try again in 1 sec"<<endl;
    return;
  }
  
  fNextEventTimer->Stop();
    
  cout << "Mutex is freeee!!"<<endl;
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
  
  
  fEventBuffer->UnLockMutex();
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

