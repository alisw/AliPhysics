#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TObjArray.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliT0CalibSeasonTimeShift.h"
#include "AliT0CalibOffsetChannelsTask.h"

#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"

// Task should calculate channels offset 
// Authors: Alla 

ClassImp(AliT0CalibOffsetChannelsTask)
//________________________________________________________________________
AliT0CalibOffsetChannelsTask::AliT0CalibOffsetChannelsTask() 
  : AliAnalysisTaskSE(),  fESD(0), fTzeroObject(0),fRunNumber(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(1, TObjArray::Class());
  fTzeroObject = new TObjArray(0);
  fTzeroObject->SetOwner(kTRUE);
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TList::Class());
}


//________________________________________________________________________
AliT0CalibOffsetChannelsTask::AliT0CalibOffsetChannelsTask(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fTzeroObject(0),fRunNumber(0)
{
  // Constructor
  
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(1, TObjArray::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliT0CalibOffsetChannelsTask::~AliT0CalibOffsetChannelsTask() 
{
  // Destructor
 printf("AliT0CalibOffsetChannels~AliT0CalibOffsetChannels() ");
 if( fTzeroObject )fTzeroObject->Delete();
}

//________________________________________________________________________
void AliT0CalibOffsetChannelsTask::ConnectInputData(Option_t *) {
  //
  //
  //
  TTree* tree=dynamic_cast<TTree*>(GetInputData(0));
  if (!tree) {
    printf("ERROR: Could not read chain from input slot 0");
  } 
  else {
    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if (!esdH) {
      printf ("ERROR: Could not get ESDInputHandler");
    } 
    else {
      fESD = esdH->GetEvent();
      printf ("*** CONNECTED NEW EVENT ****");
    }
  }
}

//________________________________________________________________________
void AliT0CalibOffsetChannelsTask::UserCreateOutputObjects()
{
  // Create histograms
  for (Int_t i=0; i<24; i++) {
    fTimeDiff[i]   = new TH1F (Form("CFD1minCFD%d",i+1),"fTimeDiff",300, -300, 300);
    fCFD[i]   = new TH1F("CFD","CFD",500, 6000, 7000);
  }
  fTzeroObject = new TObjArray(0);
  fTzeroObject->SetOwner(kTRUE);
 
  PostData(1, fTzeroObject);
 
  // Called once
}

//________________________________________________________________________
void AliT0CalibOffsetChannelsTask::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // Post output data.
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }

  const Double32_t* time = fESD->GetT0time();
  for (Int_t i=0; i<12; i++) {
    if( time[i]>1 ){
      fCFD[i]->Fill( time[i]);
      if(  time[0]>1 ) 
	fTimeDiff[i]->Fill( time[i]-time[0]);
    }
  }
  for (Int_t i=12; i<24; i++) {
    if( time[i]>1) {
      fCFD[i]->Fill( time[i]);
      if( time[12]>1 ) 
	fTimeDiff[i]->Fill( time[i]-time[12]);
    }
  }
  fRunNumber =  fESD->GetRunNumber() ; 
  
  // printf("%lf   %lf  %lf\n",orA,orC,time);
  PostData(1, fTzeroObject);
}      
 //________________________________________________________________________
void AliT0CalibOffsetChannelsTask::Terminate(Option_t *) 
{
  
  // Called once at the end of the query
  for (Int_t i=0; i<24; i++)
    fTzeroObject->AddAtAndExpand(fTimeDiff[i],i);

  for (Int_t i=24; i<48; i++)
    fTzeroObject->AddAtAndExpand(fCFD[i],i);

}

