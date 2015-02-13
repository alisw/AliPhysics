/// \file AliTPCDDL.C

#if !defined(__CINT__)
#include <TFile.h>
#include <TTree.h>
#include "AliTPCParamSR.h"
#include "AliTPCDigitsArray.h"
#include "AliSimDigits.h"
#include "AliTPCBuffer.h"
#endif


void AliTPCDDL(Int_t eventNumber=0, Int_t eth=0){
  /// eth is a threshold.
  /// Digits stored into a file have an amplitude value greater than "eth"

  const char * inFile_new = "galice.root";
  AliRunLoader *rl = AliRunLoader::Open(inFile_new,"Event","read");

  Int_t nevents=rl->GetNumberOfEvents();
  cout<<"Number of Events:"<<nevents<<endl;
  while (eventNumber<=0 || eventNumber>nevents){
    cout<<"Insert the event number:";
    cin>>eventNumber;
    cout<<endl;
  }
  rl->GetEvent(eventNumber-1);
  AliLoader *tpcloader=rl->GetLoader("TPCLoader");
  tpcloader->LoadDigits();
  TTree *digitsTree=tpcloader->TreeD();

  AliSimDigits digrows, *dummy=&digrows;
  digitsTree->GetBranch("Segment")->SetAddress(&dummy);
  Stat_t nrows = digitsTree->GetEntries();
  cout<<"Number of entries (rows):"<<nrows<<endl;
  // get the TPC parameters
  rl->CdGAFile();
  AliTPCParamSR* param = AliTPC::LoadTPCParam(gFile);
  if (!param)
    cout<<"No TPC parameter"<<endl;
  AliTPCDigitsArray *digarr=new AliTPCDigitsArray;
  digarr->Setup(param);
  digarr->ConnectTree(digitsTree);

  AliTPCBuffer *b=new AliTPCBuffer("AliTPCDDL.dat");

  //Verbose level
  // 0: Silent
  // 1: cout messages
  // 2: txt files with digits 
  //BE CAREFUL, verbose level 2 MUST be used only for debugging and
  //it is highly suggested to use this mode only for debugging digits files
  //reasonably small, because otherwise the size of the txt files can reach
  //quickly several MB wasting time and disk space.
  b->SetVerbose(0);


  nrows=Int_t(digarr->GetTree()->GetEntries());
  cout<<"Number of entries "<<nrows<<endl;
  Int_t PSector=-1;
  Int_t SubSec=0;
  for (Int_t n=0; n<nrows; n++) {
    AliSimDigits *digrow=(AliSimDigits*)digarr->LoadEntry(n);

    Int_t sec,row; // sector and row number (in the TPC)
    param->AdjustSectorRow(digrow->GetID(),sec,row);   
    // cout<<sec<<" row "<<row<<endl;
    if(PSector!=sec){
      SubSec=0;
      PSector=sec;
    }//end if

    if(sec<36){
      //inner sector [0;35]
      if(row!=30)
	//the whole row is written into the output file
	b->WriteRowBinary(eth,digrow,0,0,0,sec,SubSec,row);
      else{
	//only the pads in the range [37;48] are written into the output file
	b->WriteRowBinary(eth,digrow,37,48,1,sec,SubSec,row);
	SubSec=1;
	//only the pads outside the range [37;48] are written into the output file
	b->WriteRowBinary(eth,digrow,37,48,2,sec,SubSec,row);
      }//end else
    }//end if
    else{
      //outer sector [36;71]
      if(row==54)SubSec=2;
      if((row!=27)&&(row!=76))
	b->WriteRowBinary(eth,digrow,0,0,0,sec,SubSec,row);
      else{
	if(row==27){
	  //only the pads outside the range [43;46] are written into the output file
	  b->WriteRowBinary(eth,digrow,43,46,2,sec,SubSec,row);
	  SubSec=1;
	  //only the pads in the range [43;46] are written into the output file
	  b->WriteRowBinary(eth,digrow,43,46,1,sec,SubSec,row);
	}
	if(row==76){
	  //only the pads outside the range [33;88] are written into the output file
	  b->WriteRowBinary(eth,digrow,33,88,2,sec,SubSec,row);
	  SubSec=3;
	  //only the pads in the range [33;88] are written into the output file
	  b->WriteRowBinary(eth,digrow,33,88,1,sec,SubSec,row);
	}
      }
    }//end else
  }//end for
  cout<<"File created !"<<endl;
  cout<<"Total number of digits: "<<b->GetDigNumber()<<endl;
  delete b;
  return;
}//end AliTPCDataChallenge
