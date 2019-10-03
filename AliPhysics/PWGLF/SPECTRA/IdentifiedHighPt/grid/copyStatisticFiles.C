//using namespace std;

#include <iostream>

#include "ListDirectories.h"



const char * GetLocalFileName1(Int_t run,  const char * path);
const char * GetLocalFileName2(Int_t run,  const char * path);
//---------------------------------------------------------


void copyStatisticFiles() {
  
  
  
  
  //loading libraries
  loadlibs();
 
  
  // connect to grid
  TGrid::Connect("alien://");  
  
  // do not use scientific notation for run number
  TGaxis::SetMaxDigits(7)  ;
  
  
  // loop over all files
  Int_t ifile =-1;
  Int_t ifileGood = 0;
  Int_t ifileNotEmpty  = 0;
  while (runs[++ifile] > 0) {
    
    
    //loop over two root files
    for(Int_t i=0;i<2;++i){
      
      Long_t *id,*size,*flags,*mt;
      
      TString file;
      TFile *fr=0;
      TString file2 ;
      TFile *fr2=0;
      
      TFile *fc=0; // centrality, only in local mode for the time being
      
      cout<<"location.Data()="<<location.Data()<<endl;
      cout<<"runs[ifile]="<<runs[ifile]<<endl;
      //cout<<" output.Data()="<<output.Data()<<endl;
      
      
      switch(i){
      case 0:{
	file.Form("alien://%s/000%d/HighPtDeDx_Tree.root",location.Data(),runs[ifile] );
	
	Printf("\nBegin of reading: %s", file.Data());    
	
	gSystem->Exec(Form("alien_cp %s %s",file.Data(), GetLocalFileName1(runs[ifile],  localPath)));
	cout << Form("alien_cp %s %s",file.Data(), GetLocalFileName1(runs[ifile], localPath)) <<endl;
      }break;
      case 1:{
	file.Form("alien://%s/000%d/HighPtDeDxV0_Tree.root",location.Data(),runs[ifile] );
	
	Printf("\nBegin of reading: %s", file.Data());    
	
	gSystem->Exec(Form("alien_cp %s %s",file.Data(), GetLocalFileName2(runs[ifile],  localPath)));
	cout << Form("alien_cp %s %s",file.Data(), GetLocalFileName2(runs[ifile], localPath)) <<endl;
      }break;
	
	
	
      }
      
      
    }
    //gSystem->Exec(Form("alien_cp %s %s",file2.Data(), GetLocalFileName2(runs[ifile], localSuffix, localPath)));
    //cout << Form("alien_cp %s %s",file2.Data(), GetLocalFileName2(runs[ifile], localSuffix, localPath)) <<endl;   
    
    
  }
}

const char * GetLocalFileName1(Int_t run, const char * path) {
  // returns the filename of the local copy of the event_stat file
  static TString name;
  //  name.Form("%s/event_stat_%s_%d.root", path, suffix, run);
  name.Form("%s/HighPtDeDx_Tree_%d.root", path, run);
  return name.Data();

}
const char * GetLocalFileName2(Int_t run, const char * path) {
  // returns the filename of the local copy of the event_stat file
  static TString name;
  //  name.Form("%s/event_stat_%s_%d.root", path, suffix, run);
  name.Form("%s/HighPtDeDxV0_Tree_%d.root", path, run);
  return name.Data();

}

void loadlibs()
{
  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libMinuit");
  gSystem->Load("libPWG2spectra");
  gSystem->Load("libPWG0base"); 
}
