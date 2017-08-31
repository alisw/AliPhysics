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
  fRunNumber(0),fRefPMTA(12), fRefPMTC(0),
  fEvent(0),fStartTime(0), fEndTime(0)
{
  // Constructor

  for( int ip=0; ip < 24; ip++){
    fTimeDiff[ip] = 0;
    fCFD[ip]      = 0;
    fCDBdelays[ip]= 0;
    fCDBcfds[ip]= 0;
    fCFDvsTimestamp[ip] = NULL;
    if (ip<4 ) {
      fCDBT0s[ip]= 0;
      fT0s[ip] =NULL;
    }
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
    fRunNumber(0),fRefPMTA(12), fRefPMTC(0), fEvent(0),
    fStartTime(0), fEndTime(0)
{
  // Constructor
 
  for( int ip=0; ip < 24; ip++){
    fTimeDiff[ip] = 0;
    fCFD[ip]      = 0;
    fCDBdelays[ip]= 0;
    fCDBcfds[ip]= 0;
    fCFDvsTimestamp[ip] = NULL;
		    
    if (ip<4 ) {
      fCDBT0s[ip]= 0;
      fT0s[ip] =NULL;
    }
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
    delete fCFDvsTimestamp[ip];
  }
  for( Int_t  ip=0; ip < 4; ip++) delete fT0s[ip];
  
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
  Float_t low = fCDBcfds[fRefPMTC] - 400;
  Float_t high = fCDBcfds[fRefPMTA] + 400;
  printf (" AliT0CalibOffsetChannelsTask %f %f \n",low,high);
  
  Float_t timestart = Float_t (fStartTime - 300);
  Float_t timeend = Float_t (fEndTime +300);
  Int_t nTimeBins = (fEndTime - fStartTime+600)/600;
  for (Int_t i=0; i<24; i++) {
    fTimeDiff[i]   = new TH1F (Form("CFD1minCFD%d",i+1),"fTimeDiff",150, -300, 300);
    fCFD[i]        = new TH1F(Form("CFD%d",i+1),"CFD",250,low, high);//6000, 7000);
    // fCFDvsTimestamp[i] = new TH2F( Form("hCFDvsTimestamp%i", i+1), Form("CFD vs timestamp%i", i+1), nTimeBins, timestart,  timeend,250,low, high );
    fCFDvsTimestamp[i] = new TH2F( Form("hCFDvsTimestamp%i", i+1), Form("CFD vs timestamp%i", i+1), nTimeBins,fStartTime - 300,  fEndTime +300 ,250,low, high );
    // fCFDvsTimestamp[i]->SetName(Form("hCFDvsTimestamp%i", i+1));
    // fCFDvsTimestamp[i]->SetTitle(Form("CFD vs timestamp%i", i+1));
  }

  fTzeroORAplusORC = new TH1F("fTzeroORAplusORC","ORA+ORC /2",200,-4000,4000);   //or A plus or C 
  fResolution      = new TH1F("fResolution","fResolution",200,-2000,2000);// or A minus or C spectrum
  fTzeroORA        = new TH1F("fTzeroORA","fTzeroORA",200,-4000,4000);// or A spectrum
  fTzeroORC        = new TH1F("fTzeroORC","fTzeroORC",200,-4000,4000);// or C spectrum
  TString histname[4] = {"hT0AC","hT0A","hT0C","hResolution"};
  for (int icase=0; icase<4; icase++) 
    fT0s[icase] = new TH2F(histname[icase].Data(), histname[icase].Data(), 100, 0, 200, 200, -1000, 1000);

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
  for (int icase=0; icase<4; icase++) 
    fTzeroObject->AddAtAndExpand(fT0s[icase], 52+icase);
  for (int icase=0; icase<24; icase++) 
    fTzeroObject->AddAtAndExpand(fCFDvsTimestamp[icase], 56+icase);

  PostData(1, fTzeroObject);
  fEvent=0;
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
  UInt_t timestamp=fESD->GetTimeStamp();
  /*  if (fEvent==0 ) 
    for (int iii=0; iii<24; iii++)
    fCFDvsTimestamp[iii]->SetBins(121, timestamp-60, timestamp+72000, 250,fCDBcfds[fRefPMTC]-500, fCDBcfds[fRefPMTA] + 500); */
   fEvent++;
 
            
  AliESDTZERO* tz= (AliESDTZERO*) fESD->GetESDTZERO();
  Int_t trigT0 = fESD->GetT0Trig();
  Float_t tvdctr = tz->GetTVDC(0);
  Bool_t eq = kTRUE;
  fRunNumber =  fESD->GetRunNumber() ;
  if( fRunNumber<165747) eq = kFALSE;
    
  const Double32_t* time = fESD->GetT0time();
  const Double32_t* amp = fESD->GetT0amplitude();
  
  if(tvdctr>-5 && tvdctr<5 && tvdctr!=0) { //event selection
    //    cout<<" tvdc "<<tvdctr<<endl;
    Double32_t diff;
    for (Int_t i=0; i<24; i++) {
      if( time[i] > 0   ){
	if (eq)	{
	  fCFD[i]->Fill( time[i] );//////!!!!!
	  fCFDvsTimestamp[i]->Fill(timestamp,time[i]);
	  if(  time[fRefPMTC] > 0 && i<12)   {
	    diff =  time[i]-time[fRefPMTC];
	    fTimeDiff[i]->Fill( diff);
	  }
	  if(  time[fRefPMTA] >0  && i>11)  {
	    diff =  time[i]-time[fRefPMTA] ;
	    fTimeDiff[i]->Fill( diff);
	  }
	} //eq=1
	else  {
	  fCFD[i]->Fill( time[i] + fCDBdelays[i] );
	  if(  time[fRefPMTC] > 0 && i<12) {
	    diff =  time[i]-time[fRefPMTC] + fCDBdelays[i];
	    fTimeDiff[i]->Fill( diff);
	  } //C
	  if(  time[fRefPMTA] >0  && i>11) {
	    diff =  time[i]-time[fRefPMTA] + fCDBdelays[i];
	    fTimeDiff[i]->Fill( diff);
	  } //A
	} //eq=0
      }
      
    }
    const Double32_t* mean = fESD->GetT0TOF();
    Double32_t meanTOF = mean[0]  +  fCDBT0s[0] ;
    Double32_t orA = mean[1]  +  fCDBT0s[1] ;
    Double32_t orC = mean[2] + fCDBT0s[2] ;
    Int_t ncont = fESD->GetPrimaryVertexSPD()->GetNContributors();
   
    if(orA<99999) {
      fTzeroORA->Fill(orA);
       if (ncont>0) fT0s[1]->Fill(ncont, orA);
    }
    if(orC<99999) {
      fTzeroORC->Fill(orC);
       if (ncont>0) fT0s[2]->Fill(ncont, orC);
    }
    if(orA<99999 && orC<99999) {
      fResolution->Fill((orA-orC)/2.);
      fTzeroORAplusORC->Fill(meanTOF); 
      if (ncont>0)   {
	fT0s[0]->Fill(ncont,meanTOF );
	fT0s[3]->Fill(ncont,(orA-orC)/2. );
      }
    }
  } //if TVDC on
  PostData(1, fTzeroObject);
}      
//________________________________________________________________________
  void AliT0CalibOffsetChannelsTask::Terminate(Option_t *) 
{
  
   // Called once at the end of the query
}
/* 
@@@@ start 1495913485 end 1495916316
@@  start 1495913344.000000 end 1495916672.000000 bins 5 
*/
