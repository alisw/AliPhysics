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
#include "TGLViewer.h"

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
   
  fInfoButton = new TGLOverlayButton(dynamic_cast<TGLViewerBase*>(gEve->GetDefaultGLViewer()),  "", 0, 540, 210, 25);
    fInfoButton->SetAlphaValues(0.0, 0.8);
  
  

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
    fInfoButton->SetAlphaValues(0.8, 0.8);
    fInfoButton->SetText("Waiting for buffer...");
    gEve->Redraw3D(kFALSE);
    cout << "try again in 1 sec"<<endl;
    return;
  }
  fInfoButton->SetAlphaValues(0.8, 0.8);

  fNextEventTimer->Stop();
    
  cout << "Mutex is freeee!!"<<endl;
  TList * aSyncEvent = fEventBuffer->GetASyncEvent();
  TList * event = static_cast<TList*>(fEventBuffer->NextEvent());
  if(event) {
    cout << "Got the event, reset the display " <<endl;
    fInfoButton->SetText("Reset display..");
    ResetDisplay();
    cout << "Process event"<<endl;
    fInfoButton->SetText("Processing event..");
    ProcessEvent(event);
    if(aSyncEvent) {
      cout  << "Process asynchroneous event" << endl;
      ProcessEvent(aSyncEvent);
    }  else {
      cout << "Could not get async event"<<endl;
    }
    
    cout << "Upate the display"<<endl;
    fInfoButton->SetText("Updating display...");
    UpdateDisplay();
  
  } else {
    cout << "couldn't get the sync event"<<endl;
    fInfoButton->SetAlphaValues(0.8, 0.8);
    fInfoButton->SetText("Waiting for buffer...");
    fEventBuffer->UnLockMutex();
    fEventBuffer->CreateBufferThread();
    fNextEventTimer->Start(3000);
    return;
  }
  
  fInfoButton->SetAlphaValues(0.0, 0.0);
  fInfoButton->SetText("Done..");
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

