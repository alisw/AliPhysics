#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iostream.h>
#include "AliITSDDLRawData.h"
#endif

/*
Before running this macro it is necessary comment the following line of the method
AddDigit in the class AliITSsimulationSDD
//if( fResponse->Do10to8() ) signal = Convert8to10( signal ); 
In this way the amplitude value for signal coming from SDD takes only 8 bits and not 10.
*/

void AliITSDDLRawData(char* DigitsFile="galice.root"){
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
  // Connect the Root input  file containing Geometry, Kine and Hits
  // galice.root file by default
  char* filename="galice.root";
  // TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(DigitsFile);
  if (!file){
    //    file = new TFile(filename);
    file = new TFile(DigitsFile);
    cout<<"NEW FILE CREATED !!!"<<endl;
  }//end if
  file->ls();
  // Get AliRun object from file
  if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice)cout<<"AliRun object found on file "<<filename<<endl;
    if(!gAlice){
      cout<<"Can't access AliRun object on file "<<filename<<endl;
      cout<<"Macro execution stopped!!!"<<endl;
    }
  }
  //  gAlice->SetTreeDFileName("digits.root");
  gAlice->SetTreeDFileName(DigitsFile);
  Int_t nparticles = gAlice->GetEvent(0);
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
  TTree *TD = gAlice->TreeD();

  AliITSDDLRawData *util=new AliITSDDLRawData();
  TStopwatch timer;
  
  //SILICON PIXEL DETECTOR
  cout<<"Formatting data for SPD"<<endl;
  timer.Start();
  util->RawDataSPD(ITS,TD);
  timer.Stop();
  timer.Print();
    
  //SILICON DRIFT DETECTOR
  cout<<"Formatting data for SDD"<<endl;
  timer.Start();
  util->RawDataSDD(ITS,TD);
  timer.Stop();
  timer.Print();
  
  //SILICON STRIP DETECTOR
  cout<<"Formatting data for SSD"<<endl;
  timer.Start();
  util->RawDataSSD(ITS,TD);
  timer.Stop();
  timer.Print();
  
  delete util;
  return;
}
