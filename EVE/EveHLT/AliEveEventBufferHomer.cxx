#if __GNUC__>= 3
   using namespace std;
#endif

#include "AliEveEventBufferHomer.h"
#include "AliEveHOMERManager.h"

#include "TList.h"
#include "TFile.h"

#include <iostream>

ClassImp(AliEveEventBufferHomer);

///_______________________________________________________________________
AliEveEventBufferHomer::AliEveEventBufferHomer() :
  fHomer(NULL),
  fEventNo(0),
  fAsyncList(NULL)
{
  // see header file for class documentation
  fHomer = new AliEveHOMERManager();
  Initialize();

}


///____________________________________________________________________
AliEveEventBufferHomer::~AliEveEventBufferHomer() {
  
  if(fHomer)
    delete fHomer;
  fHomer = NULL;

}


///____________________________________________________________
ULong64_t AliEveEventBufferHomer::GetEventIdFromSource() {
  ///see header file for documentation
  return fHomer->GetEventID();
}

///______________________________________________________________________
TObject * AliEveEventBufferHomer::GetEventFromSource() {
  //see header file for documentation
  if(fHomer) {
    cout << "calling nexthomerevent"<<endl;
    TList * blockList = fHomer->NextHOMEREvent();
    cout  << "getting async block list"<<endl;
    TList * aList = fHomer->GetAsyncBlockList();
    fAsyncList = aList;
    if(blockList)  return dynamic_cast<TObject*>(blockList);
    else return NULL;
  } 

  cout << "AliEveEventBufferHomer::GetEventFromSource() : fHomer is null pointer!!"<<endl;
  return NULL;
}

///___________________________________________________________________________
void AliEveEventBufferHomer::AddToBuffer(TObject * event) {
  cout << "AliEveEventBufferHomer::Add to buffer"<<endl;
  if(!event) {
    cout << "event is NULL pointer, return" << endl;
    return;
  }


  TList * listIn = dynamic_cast<TList*>(event);
  if(!listIn || !(listIn->GetSize() > 0)) {
    cout  << "AliEveEventBufferHomer::AddToBuffer(): List Is empty, not added to buffer!"<<endl;
    return;
  }
  
  fBIndex[kTop] = CalculateNext(fBIndex[kTop]);
  TList * list = dynamic_cast<TList*>(fEventBuffer->At(fBIndex[kTop]));
  if(list) {
    list->Delete();
    list->AddAll(dynamic_cast<TList*>(event));
  }
}



///______________________________________________________________________
void AliEveEventBufferHomer::Initialize() {
  //Add TList pointers to the buffer array

  //Create TLists in all of buffer
  for(Int_t i = 0; i < fBufferSize; i++) {
    TList * list = new TList();
    list->SetOwner(kTRUE);
    fEventBuffer->AddAt(list, i);
  }

  //Create the async list
  fAsyncList = new TList();
  fAsyncList->SetOwner(kTRUE);

  Int_t iResult = fHomer->Initialize();
  if(iResult) {
    cout << "Error Initializing HOMER!!!"<<endl;
    return;
  }


  fHomer->SetRetryCount(1,15);
  fHomer->SetBlockOwner(kFALSE);
  //fHomer->StartEveSourceListLoop();
}


///______________________________________________________________________
void AliEveEventBufferHomer::ConnectToSource () {

  fHomer->SetRetryCount(200, 15);
  fHomer->ReConnectHOMER();
}

///_____________________________________________________________________
void AliEveEventBufferHomer::WriteToFile(Int_t runnumber){

  TFile * file = TFile::Open(Form("Run_%d_0x%016llu_ESD.root", runnumber, GetEventId()), "RECREATE"); 
  fEventBuffer->At(fBIndex[kCurrent])->Write("blockList", TObject::kSingleKey);
  file->Close();
  
  if(fAsyncList) {
    TFile * afile = TFile::Open(Form("Run_%d_0x%016llu_Async.root",  runnumber, GetEventId()), "RECREATE"); 
    fAsyncList->Write("blockList", TObject::kSingleKey);
    //aFile-
    afile->Close();
  }
}	     
