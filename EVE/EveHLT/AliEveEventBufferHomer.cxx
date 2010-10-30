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




///______________________________________________________________________
TObject * AliEveEventBufferHomer::GetEventFromSource() {
  //see header file for documentation
  if(fHomer) {
    cout << "calling nexthomerevent"<<endl;
    TList * blockList = fHomer->NextHOMEREvent();
    if(blockList)  return dynamic_cast<TObject*>(blockList);
    else return NULL;
  } 

  cout << "AliEveEventBufferHomer::GetEventFromSource() : fHomer is null pointer!!"<<endl;
  return NULL;
}
///___________________________________________________________________________
void AliEveEventBufferHomer::AddToBuffer(TObject * event) {
  cout << "Add to buffer"<<endl;
  if(!event) return;
  fBIndex[kTop] = CalculateNext(fBIndex[kTop]);
  TList * list = dynamic_cast<TList*>(fEventBuffer->At(fBIndex[kTop]));
  if(list) {
    list->Clear();
    list->AddAll(dynamic_cast<TList*>(event));
  }
}


///______________________________________________________________________
void AliEveEventBufferHomer::Initialize() {
  //Add TList pointers to the buffer array
  for(Int_t i = 0; i < fBufferSize; i++) {
    TList * list = new TList();
    list->SetOwner(kTRUE);
    fEventBuffer->AddAt(list, i);
  }


  Int_t iResult = fHomer->Initialize();
  if(iResult) {
    cout << "Error Initializing HOMER!!!"<<endl;
    return;
  }


  fHomer->SetRetryCount(1,15);
  fHomer->SetBlockOwner(kFALSE);
  fHomer->StartEveSourceListLoop();
}


///______________________________________________________________________
void AliEveEventBufferHomer::ConnectToSource () {

  fHomer->SetRetryCount(200, 15);
  fHomer->ReConnectHOMER();
}

///_____________________________________________________________________
void AliEveEventBufferHomer::WriteToFile(){

  TFile * file = TFile::Open(Form("Event_0x%016X_ITS.root", 100), "RECREATE"); 
  fEventBuffer->At(fBIndex[kTop])->Write("blockList", TObject::kSingleKey);
  file->Close();
  
  if(fAsyncList) {
    TFile * afile = TFile::Open(Form("Event_0x%016X_Async.root", 100), "RECREATE"); 
    fAsyncList->Write("blockList", TObject::kSingleKey);
    afile->Close();
  }
}	     
