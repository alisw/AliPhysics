#if !defined(__CINT__)
#include <Riostream.h>
#include "AliITSDDLRawData.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "AliITS.h"
#endif

/*
Before running this macro it is necessary to comment the following line of the method
AddDigit in the class AliITSsimulationSDD
//if( fResponse->Do10to8() ) signal = Convert8to10( signal ); 
In this way the amplitude value for signal coming from SDD takes only 8 bits and not 10.
*/
//DigitsFile is the input file that contains digits

void AliITSDDLRawData(Int_t eventNumber=0){

  Int_t spdLDCs=2;
  Int_t sddLDCs=4;
  Int_t ssdLDCs=2;
  Int_t eventn=0;
  const char * inFile_new = "galice.root";
  AliRunLoader *rl = AliRunLoader::Open(inFile_new,"Event","read");
  rl->LoadgAlice();
  gAlice=rl->GetAliRun();
  Int_t nevents=rl->GetNumberOfEvents();
  cout<<"Number of Events:"<<nevents<<endl;
  while (eventNumber<=0 || eventNumber>nevents){
    cout<<"Insert the event number:";
    cin>>eventNumber;
    cout<<endl;
  }
  rl->GetEvent(eventNumber-1);
  AliLoader *itsloader=rl->GetLoader("ITSLoader");
  itsloader->LoadDigits();
  TTree *TD=itsloader->TreeD();
  gAlice=rl->GetAliRun(); 
  if(!gAlice){
    cout<<"gAlice is null"<<endl;
    return;
  } 
  AliITS *ITS  = (AliITS*)gAlice->GetDetector("ITS");

  Int_t nmodules;
  ITS->InitModules(-1,nmodules);
  ITS->GetDetTypeSim()->SetTreeAddressD(TD,"ITS");

    AliITSDDLRawData *util=new AliITSDDLRawData();
    //Verbose level
    // 0: Silent
    // 1: cout messages
    // 2: txt files with digits 
    //BE CAREFUL, verbose level 2 MUST be used only for debugging and
    //it is highly suggested to use this mode only for debugging digits files
    //reasonably small, because otherwise the size of the txt files can reach
    //quickly several MB wasting time and disk space.
    util->SetVerbose(0);
    TStopwatch timer;
    
    //SILICON PIXEL DETECTOR
    cout<<"Formatting data for SPD"<<endl;
    timer.Start();
    util->RawDataSPD(ITS,TD,spdLDCs,eventNumber);
    timer.Stop();
    timer.Print();
    //ONLY FOR DEBUGGING 
    //    util->TestFormat(eventNumber);
    
    //SILICON DRIFT DETECTOR
    cout<<"Formatting data for SDD"<<endl;
    timer.Start();
    util->RawDataSDD(ITS,TD,sddLDCs,eventNumber);
    timer.Stop();
    timer.Print();
    
    //SILICON STRIP DETECTOR
    cout<<"Formatting data for SSD"<<endl;
    timer.Start();
    util->RawDataSSD(ITS,TD,ssdLDCs,eventNumber);
    timer.Stop();
    timer.Print();
    
    delete util;

    return;
}
