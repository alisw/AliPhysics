#if !defined(__CINT__)
#include <Riostream.h>
#include "AliITSDDLRawData.h"
#endif

/*
Before running this macro it is necessary to comment the following line of the method
AddDigit in the class AliITSsimulationSDD
//if( fResponse->Do10to8() ) signal = Convert8to10( signal ); 
In this way the amplitude value for signal coming from SDD takes only 8 bits and not 10.
*/
//DigitsFile is the input file that contains digits

void AliITSDDLRawData(char* DigitsFile="galiceD.root"){
#ifdef __NOCOMPILED__
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
  else {
#endif
    if(gAlice){
      delete gAlice;
      gAlice=0;
    }
#ifdef __NOCOMPILED__
  }
#endif
  Int_t eventNumber=0;
  Int_t spdLDCs=2;
  Int_t sddLDCs=4;
  Int_t ssdLDCs=2;
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(DigitsFile);
  if (!file){
    file = new TFile(DigitsFile);
  }//end if
  file->ls();

  // Get AliRun object from file
    if (!gAlice){
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice)cout<<"AliRun object found on file "<<DigitsFile<<endl;
      if(!gAlice){
	cout<<"Can't access AliRun object on file "<<DigitsFile<<endl;
	cout<<"Macro execution stopped!!!"<<endl;
	exit(1);
      }//end if
    }//end if 
    gAlice->SetTreeDFileName(DigitsFile);
    //  Long_t nparticles = gAlice->GetEvent(0);
  
    //Int_t nparticles = gAlice->GetEvent(0);
    // 
    // ITS
    AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
    Int_t nmodules;
    ITS->InitModules(-1,nmodules);
    cout<<"Number of ITS modules= "<<nmodules<<endl;
    //cout<<"Filling modules... It takes a while, now. Please be patient"<<endl;
    //ITS->FillModules(0,0,nmodules," "," ");
    //cout<<"ITS modules .... DONE!"<<endl;
    // DIGITS
    
    
    TTree* TD = (TTree*)file->Get("TreeD0");
    if (TD == 0x0){
      ::Error("DDLRawData","Can not find tree with ITS digits");
      return;
    }//end if
    ITS->SetTreeAddressD(TD);
    
    
    //TTree *TD = gAlice->TreeD();
    cout<<"Insert the event number:";
    cin>>eventNumber;
    cout<<endl;

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
