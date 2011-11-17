#if __GNUC__>= 3
   using namespace std;
#endif

#include "AliEveEventBufferOffline.h"

#include "AliESDEvent.h"
#include "AliEveEventBufferOffline.h"
#include "AliEveEventBuffer.h"

#include "TTimer.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

ClassImp(AliEveEventBufferOffline);

///_______________________________________________________________________
AliEveEventBufferOffline::AliEveEventBufferOffline() :
  fFile(NULL),
  fNEntries(0),
  fEventNo(0),
  fEvent(NULL),
  fTree(NULL)
{
  // see header file for class documentation
  //Not Allowed
}

///_______________________________________________________________________
AliEveEventBufferOffline::AliEveEventBufferOffline(TString filename)  : 
  fFile(NULL),
  fNEntries(0),
  fEventNo(0),
  fEvent(NULL),
  fTree(NULL)
{
  
  fEvent = new AliESDEvent();

  cout <<"opening file " << filename << " what?" <<endl;
  fFile = TFile::Open(filename, "READ");
  if(!fFile) {
    cout << "Couldn't open file, crashing hard! Soon?"<<endl;
    return;
  }

  fTree = dynamic_cast<TTree *>(fFile->Get("HLTesdTree"));
  if(fTree) {
    fNEntries = fTree->GetEntries();
    cout << "File has " << fNEntries << "events" << endl;
    fEvent->ReadFromTree(fTree);
  }
}




///____________________________________________________________________
AliEveEventBufferOffline::~AliEveEventBufferOffline() {
  
  if(fFile)
    delete fFile;
  fFile = NULL;

  if(fTree)
    delete fTree;
  fTree = NULL;
  
  if(fEvent)
    delete fEvent;
  fEvent = NULL;

}


///______________________________________________________________________
TObject * AliEveEventBufferOffline::GetEventFromSource() {
  //see header file for documentation
  if(fTree) {
    fTree->GetEntry(fEventNo++);
    if(fEventNo == fNEntries)
      fEventNo = 0;
  }  else {
    cout << "Tree not found, probably bad file!!"<<endl;
    return NULL;
  }

  //Copy event into new event (must be deleted when no longer needed in list!!
  //AliESDEvent * event = new AliESDEvent();
  //fEvent->Copy(*event);
  //cout << event->GetNumberOfCaloClusters() << endl;

  if (fEvent) {
    return dynamic_cast<TObject*>(fEvent);
  } else {
    cout << "error getting event" << endl;
    return NULL;
  }
}

///___________________________________________________________________
void AliEveEventBufferOffline::ConnectToSource() {
  //Needed for homer version
  // see header file for class documentation
  return;
}


///_____________________________________________________________________
void AliEveEventBufferOffline::WriteToFile(Int_t runnumber){
  //Inherited from AliEveEventBuffer
  TFile * file = TFile::Open(Form("%d_0x%016llu_ESD.root", runnumber, GetEventId()), "RECREATE"); 
  fEventBuffer->At(fBIndex[kCurrent])->Write("blockList", TObject::kSingleKey);
  file->Close();
}	     
