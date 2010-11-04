#include <iostream>

#include "TObjArray.h"
#include "TTimer.h"
#include "TThread.h"
#include "AliEveEventBuffer.h"


//Not needed, only for debug
#include "AliESDEvent.h"

using namespace std;

ClassImp(AliEveEventBuffer)

///_______________________________________________________________________
AliEveEventBuffer::AliEveEventBuffer() :
  fBufferSize(10),
  fPreBuffer(4),
  fBusy(kFALSE),
  fEventBuffer(NULL),
  fCurrentEvent(NULL),
  fBIndex(),
  fTimer(NULL),
  fEventId(),
  fBufferMonStarted(kFALSE)
 {
  // see header file for class documentation
  fEventBuffer = new TObjArray(fBufferSize, 0);
  fEventBuffer->SetOwner(kFALSE);
  
  for(int id = 0; id < kSize; id++) {
    fBIndex[id] = -1;
  }
  
  fTimer = new TTimer();
  fTimer->Connect("Timeout()", "AliEveEventBuffer", this, "CreateBufferThread()");

  fEventId = new Int_t[fBufferSize];


}



///_______________________________________________________________________
AliEveEventBuffer::~AliEveEventBuffer() {
  // see header file for class documentation
  
  if ( fEventBuffer ) {
    fEventBuffer->Clear();
    delete fEventBuffer;
  }
  fEventBuffer = NULL;

  if(fCurrentEvent)
    delete fCurrentEvent;
  fCurrentEvent = NULL;

}

void AliEveEventBuffer::CreateBufferThread() {
  //  cout << "hereherehere"<<endl;
  TThread * fThread = new TThread(AliEveEventBuffer::BufferThread, (void*) this);
  fThread->Run();

}

///___________________________________________________________________________
void * AliEveEventBuffer::BufferThread(void * buffer) {
  cout <<"BufferThread : " <<endl;
  if(buffer) {
    if (!reinterpret_cast<AliEveEventBuffer*>(buffer)->GetBusy()) {
      reinterpret_cast<AliEveEventBuffer*>(buffer)->MonitorBuffer();
    } else {
      cout << "busy"<<endl;
    }
    
  } else {
    cout << "no buffer"<<endl;
  }
  return (void*)0;
}

///_____________________________________________________________________________
void AliEveEventBuffer::MonitorBuffer() {
  cout << "Monitorbuffer: " << endl;
  if ( (CalculateDifference(fBIndex[kTop],fBIndex[kLast]) < fPreBuffer) ) {
    if(GetBusy()) {
      cout << "Already called FetchEvent()" << endl;
      return;
    } else {
      fBusy = kTRUE;
      FetchEvent();
      fBusy = kFALSE;
    }
  } else {
    //StopBufferMonitor();
    
    fBusy = kFALSE;
  }
}


///_______________________________________________________________________________
TObject * AliEveEventBuffer::NextEvent() {
  //See header file for documentation
  cout << "In enxtevent"<<endl;


  // if(fBusy) {
  //   cout << "Event Buffer busy"<<endl;
  //   return NULL;
  // }
  //  else SetBusy(kTRUE);
  TObject * nextEvent = GetNextUnSeen();
  //SetBusy(kFALSE);
  return nextEvent;
}

///______________________________________________________________________________
TObject * AliEveEventBuffer::Back() {
  cout << "go back"<<endl;
  PrintIndeces();
  Int_t prevId = CalculatePrevious(fBIndex[kCurrent]);
  if(prevId == fBIndex[kTop]) {
    cout << "returning NULL" << endl;
    return NULL;
  } else {
    fBIndex[kCurrent] = prevId;
    PrintIndeces();
    cout <<"returning: "<< fBIndex[kCurrent] << " " << fEventBuffer->At(fBIndex[kCurrent]);
    return fEventBuffer->At(fBIndex[kCurrent]);
  }
}

///______________________________________________________________________________
TObject * AliEveEventBuffer::Fwd() {
  PrintIndeces();
  if (fBIndex[kCurrent] == fBIndex[kLast]) {
    cout<<  "returning NULL"<<endl;
    return NULL;
  }
  
  fBIndex[kCurrent] = CalculateNext(fBIndex[kCurrent]);
  TObject * event = fEventBuffer->At(fBIndex[kCurrent]);
  return event;
}



///________________________________________________________________________________
TObject * AliEveEventBuffer::GetNextUnSeen() {
  //See header file for documentation
  cout << "GetNextUnSeen"<<endl;
  PrintIndeces();
  if(CalculateDifference(fBIndex[kTop], fBIndex[kLast])) {
    fBIndex[kLast] = CalculateNext(fBIndex[kLast]);
    fBIndex[kCurrent] = fBIndex[kLast];
    PrintIndeces();
    return fEventBuffer->At(fBIndex[kCurrent]);      
  } else {
    cout << "No new event available, only events in buffer available!"<<endl;
    return NULL;
  } 
}
///_________________________________________________________________________________
void AliEveEventBuffer::PrintIndeces() {
  for(Int_t i = 0; i < kSize; i++) {
    cout << i << ": " << fBIndex[i] << endl;
  }
}
///_________________________________________________________________________________
void AliEveEventBuffer::PrintBuffer() {
  for(Int_t i = 0; i < 10; i++) {
    AliESDEvent * event = dynamic_cast<AliESDEvent*>(fEventBuffer->At(i));
    if(event) {
      cout << i << ": " <<event << " " << event->GetEventNumberInFile() << endl;;
    }
  }
}

///____________________________________________________________________________________
void AliEveEventBuffer::FetchEvent() {
  cout << "FetchEvent " << endl;
  TObject * event = GetEventFromSource();
  if(event) AddToBuffer(event);
  PrintIndeces();
  cout << "FetchedEvent " << endl;
  
}

///_________________________________________________________________________________
void AliEveEventBuffer::AddToBuffer(TObject * event) {
  cout << "Add to buffer"<<endl;
  if(!event) return;

  fBIndex[kTop] = CalculateNext(fBIndex[kTop]);
  //Delete the event already there (ok to delete as object, not aliesdevent, TList?)
  //TObject * object = fEventBuffer->At(fBIndex[kTop]);
  fEventBuffer->RemoveAt(fBIndex[kTop]);
  //if (object) delete object;
  fEventBuffer->AddAt(event, fBIndex[kTop]);
}



///_____________________________________________________________________________________
Int_t AliEveEventBuffer::CalculateNext(Int_t current) {
  //See header file for documentation
  current++;
  if(current == fBufferSize) current = 0;
  return current;
}


///_____________________________________________________________________________________
Int_t AliEveEventBuffer::CalculatePrevious(Int_t current) {
  //See header file for documentation
  cout << "CalculatePrev:  " << current; 
  current--;
  if(current == -1) current += fBufferSize;
  cout << "... " << current << endl;
  return current;
}

///__________________________________________________________________________________
Int_t AliEveEventBuffer::CalculateDifference(Int_t top, Int_t low) {
  //See header file for documentation
  if (top > low) {
    //    cout << "top > low"<<endl;
    return (top - low);
  } else if (top < low) {
    // cout << "low < top"<<endl;
    return (fBufferSize - low + top);
  } else {
    //cout << "calculated to 0"<<endl;
    return 0;
  }
}

///___________________________________________________________________________________
void AliEveEventBuffer::StartBufferMonitor() {
  //cout << "NOT !!! starting buffer mon"<<endl;
  cout << "starting buffer mon"<<endl;
  CreateBufferThread();
  SetBufferMonStarted(kTRUE);
  fTimer->Start(10000);
}
///___________________________________________________________________________________
void AliEveEventBuffer::StopBufferMonitor() {
  cout << "Stopping buffer mon"<<endl;
  SetBufferMonStarted(kFALSE);
  fTimer->Stop();
}


// //_________________________________________________________________________________
// Int_t AliEveEventBuffer::NavigateEventBufferBack() { 
//   // see header file for class documentation

//   // -- reached the end of the buffer
//   if ( fNavigateBufferIdx == fBufferLowIdx )
//     return -1;

//   Int_t newIdx = fNavigateBufferIdx - 1;
//   if ( newIdx == -1 )
//     newIdx = BUFFERSIZE-1;

//   fCurrentBufferIdx = fNavigateBufferIdx = newIdx;

//   return newIdx;
// }

// //_______________________________________________________________
// Int_t AliEveEventBuffer::NavigateEventBufferFwd() {
//   // see header file for class documentation

//   // -- reached the top of the buffer
//   if ( fNavigateBufferIdx == fBufferTopIdx )
//     return -1;

//   Int_t newIdx = fNavigateBufferIdx + 1;
//   if ( newIdx == BUFFERSIZE )
//     newIdx = 0;
  
//   fCurrentBufferIdx = fNavigateBufferIdx = newIdx;

//   return newIdx;
// }

// void AliEveEventBuffer::MonitorBuffer() {
//   //See header file for documentation
//   if( GetNAvailableEvents() < 10) {
//     StopBufferChecker();
//     StartLoop();
//   }
// }

// void AliEveEventBuffer::StartLoop() {
//   //See header file for documentation
//   fTimer->Start(2000);
// }
// void AliEveEventBuffer::StopLoop() {
//   //See header file for documentation
//   fTimer->Stop();
// }

// void AliEveEventBuffer::StartBufferChecker() {
//   //See header file for documentation
//   fBufferTimer->Start(2000);
// }
// void AliEveEventBuffer::StopBufferChecker() {
//   //See header file for documentation
//   fBufferTimer->Stop();
// }

// AliESDEvent * GetNextEvent() {
  
//   tree->GetEntry(fEvent++);

//   AliESDEvent * event = new AliESDEvent();
//   event->ReadFromTree(fTree);
//   if (event) {
//     return event;
//   } else {
//     cout << "error getting event" << endl;
//     return NULL;
//   }
// }
