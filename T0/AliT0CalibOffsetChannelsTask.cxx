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

#include "AliT0CalibOffsetChannelsTask.h"

//#include "AliCDBMetaData.h"
//#include "AliCDBId.h"
//#include "AliCDBEntry.h"
//#include "AliCDBManager.h"
//#include "AliCDBStorage.h"

// Task should calculate channels offset 
// Authors: Alla 

ClassImp(AliT0CalibOffsetChannelsTask)
//________________________________________________________________________
AliT0CalibOffsetChannelsTask::AliT0CalibOffsetChannelsTask() 
  : AliAnalysisTaskSE(),  fESD(0x0), fTzeroObject(0x0),
  fTzeroORA(0x0), fTzeroORC(0x0), fResolution(0x0), fTzeroORAplusORC(0x0),
  fRunNumber(0)
{
  // Constructor

  for( int ip=0; ip < 24; ip++){
    fTimeDiff[ip] = 0;
    fCFD[ip]      = 0;
    fCDBdelays[ip]= 0;
    fCDBcfds[ip]= 0;
    fCDBT0s[ip]= 0;
  }

  // Define input and output slots here
  // Input slot #0 works with a TChain
  //  DefineInput(0,  TChain::Class());
  //  DefineOutput(1, TObjArray::Class());
}


//________________________________________________________________________
AliT0CalibOffsetChannelsTask::AliT0CalibOffsetChannelsTask(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fTzeroObject(0),
  fTzeroORA(0x0), fTzeroORC(0x0), fResolution(0x0), fTzeroORAplusORC(0x0),
    fRunNumber(0)
{
  // Constructor
 
  for( int ip=0; ip < 24; ip++){
    fTimeDiff[ip] = 0;
    fCFD[ip]      = 0;
    fCDBdelays[ip]= 0;
    fCDBcfds[ip]= 0;
    fCDBT0s[ip]= 0;

  }
 
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  DefineOutput(1, TObjArray::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
}

//________________________________________________________________________
AliT0CalibOffsetChannelsTask::~AliT0CalibOffsetChannelsTask() 
{
  // Destructor
  // printf("AliT0CalibOffsetChannels~AliT0CalibOffsetChannels() ");
  delete fTzeroORA;
  delete fTzeroORC;
  delete fResolution;
  delete fTzeroORAplusORC;
  for( Int_t  ip=0; ip < 24; ip++){
    delete fTimeDiff[ip];
    delete fCFD[ip];
  }
  
  delete fTzeroObject;
}

//________________________________________________________________________
/*void AliT0CalibOffsetChannelsTaskX::ConnectInputData(Option_t *) {
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
*/
//________________________________________________________________________
void AliT0CalibOffsetChannelsTask::UserCreateOutputObjects()
{
  // Create histograms
  for (Int_t i=0; i<24; i++) {
    fTimeDiff[i]   = new TH1F (Form("CFD1minCFD%d",i+1),"fTimeDiff",150, -300, 300);
    fCFD[i]        = new TH1F(Form("CFD%d",i+1),"CFD",250, 2000, 3000);//6000, 7000);
    //    fCFD[i]        = new TH1F(Form("CFD%d",i+1),"CFD",250, -1000, 1000);//6000, 7000);
  }

  fTzeroORAplusORC = new TH1F("fTzeroORAplusORC","ORA+ORC /2",200,-2500,2500);   //or A plus or C 
  fResolution      = new TH1F("fResolution","fResolution",400,-2500,2500);// or A minus or C spectrum
  fTzeroORA        = new TH1F("fTzeroORA","fTzeroORA",200,-2500,2500);// or A spectrum
  fTzeroORC        = new TH1F("fTzeroORC","fTzeroORC",200,-2500,2500);// or C spectrum

  
  fTzeroObject     = new TObjArray(0);
  fTzeroObject->SetOwner(kTRUE);
  
  for (Int_t i=0; i<24; i++)
    fTzeroObject->AddAtAndExpand(fTimeDiff[i],i);

  for (Int_t i=0; i<24; i++)
    fTzeroObject->AddAtAndExpand(fCFD[i],i+24); //24 - 48

  fTzeroObject->AddAtAndExpand(fTzeroORAplusORC, 48);
  fTzeroObject->AddAtAndExpand(fResolution, 49);
  fTzeroObject->AddAtAndExpand(fTzeroORA, 50);
  fTzeroObject->AddAtAndExpand(fTzeroORC, 51);

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

  fRunNumber =  fESD->GetRunNumber() ; 

  const Double32_t* time = fESD->GetT0time();
  const Double32_t* amp = fESD->GetT0amplitude();
  Double32_t diff;
  for (Int_t i=0; i<12; i++) {
    if( time[i] !=0  && amp[i]>0.1 ){
      fCFD[i]->Fill( time[i] );
      //  printf(" time %f ocdb %f \n", time[i],fCDBcfds[i]); 
      
      if(  time[0] != 0 ) {
	diff =  time[i]-time[0] + fCDBdelays[i];
	fTimeDiff[i]->Fill( diff);
      }
    }
  }
  for (Int_t i=12; i<24; i++) {
    if( time[i] != 0 && amp[i]>0.1) {
      fCFD[i]->Fill( time[i]);
      //  printf(" time %f ocdb %f \n", time[i],fCDBcfds[i]); 
       if( time[12] != 0 ) {
	diff =  time[i]-time[12] + fCDBdelays[i];
	fTimeDiff[i]->Fill( diff);
      }
    }
  }
  const Double32_t* mean = fESD->GetT0TOF();
  Double32_t meanTOF = mean[0] +  fCDBT0s[0];
  Double32_t orA = mean[1]  +  fCDBT0s[1];
  Double32_t orC = mean[2] +  fCDBT0s[2];
 
  if(orA<9999) fTzeroORA->Fill(orA);
  if(orC<9999) fTzeroORC->Fill(orC);
  if(orA<9999 && orC<9999) fResolution->Fill((orA-orC)/2.);
  if(orA<9999 && orC<9999) fTzeroORAplusORC->Fill(meanTOF); 

  //  printf("%f   %f  %f\n",orA,orC,meanTOF);
  PostData(1, fTzeroObject);
}      
 //________________________________________________________________________
void AliT0CalibOffsetChannelsTask::Terminate(Option_t *) 
{
  
   // Called once at the end of the query
}
 
