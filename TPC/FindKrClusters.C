/****************************************************************************
 *           Origin: A.Matyja amatyja@cern.ch                               *
 ****************************************************************************/

/*

  macro to create array of clusters from TPC digits
  input files - galice.root 
                digits.root - file with digits - usualy use link to galice.root
		            - in splitted mode - neccesary to create link to proper file
			    
   output file - TPC.RecPoints.root


//  Warning - if cluster file AliTPCclusters.root already exist - macro exit and don't produce anything
	       
 
*/


#ifndef __CINT__
#include <iostream.h>
#include "AliRun.h"
#include "AliTPCv4.h"
#include "AliTPCParam.h"
#include "AliTPCclusterKr.h"
#include "AliTPCclustererKr.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TTree.h"
#endif

Int_t FindKrClusters(){

  AliRunLoader* rl = AliRunLoader::Open("galice.root");
  if (rl == 0x0) {
    cerr<<"Can not open session"<<endl;
    return 1;
  }
  
  AliTPCLoader *tpcl = (AliTPCLoader*)rl->GetLoader("TPCLoader");
  if (tpcl == 0x0) {
    cerr<<"Can not get TPC Loader"<<endl;
    return 1;
  }

  if (tpcl->LoadDigits()) {
    cerr<<"Error occured while loading digits"<<endl;
    return 1;
  }

  if (rl->LoadgAlice()) {
    cerr<<"Error occured while LoadgAlice"<<endl;
    return 1;
  }
  
  gAlice=rl->GetAliRun();
  if (!gAlice) {
    cerr<<"Can't get gAlice !\n";
    return 1;
  }

  TDirectory *cwd = gDirectory;

  AliTPCv4 *tpc = (AliTPCv4*)gAlice->GetDetector("TPC");
  Int_t ver = tpc->IsVersion(); 
  cerr<<"TPC version "<<ver<<" has been found !\n";

  rl->CdGAFile();
  
  AliTPCParam *param=(AliTPCParamSR *)gDirectory->Get("75x40_100x60_150x60");
  if (!param) {cerr<<"TPC parameters have not been found !\n"; return 4;}
  
  AliTPCDigitsArray *digarr=new AliTPCDigitsArray;
  digarr->Setup(param);

  cerr<<"It has begun"<<endl;  
  TStopwatch timer;
  timer.Start();
  cwd->cd();
  
  Int_t nevmax=rl->GetNumberOfEvents();//number of events in run
  for(Int_t nev=0;nev<nevmax /*&& nev<1*/ ;nev++){
    rl->GetEvent(nev);
    
    TTree* input_tree= tpcl->TreeD();//tree with digits
    if (input_tree == 0x0){
      cerr << "Can not get TreeD for event " <<nev<<endl;
      continue;
    }
    
    digarr->ConnectTree(input_tree);
    
    TTree *output_tree =tpcl->TreeR();
    if(output_tree==0x0){
      tpcl->MakeTree("R");
      output_tree = tpcl->TreeR();
      if (output_tree == 0x0){
	cerr << "Problems with output tree (TreeR) for event "<<nev<<endl;
	continue;
      }
    }
    
    cout<<"Processing event "<<nev<<endl;

    //test
    //cout<<"nentries"<<input_tree->GetEntries();

    AliTPCclustererKr *clusters = new AliTPCclustererKr();
    clusters->SetParam(param);
    clusters->SetInput(input_tree);
    clusters->SetOutput(output_tree);
    clusters->SetDigArr(digarr);
    clusters->finderIO();

    tpcl->WriteRecPoints("OVERWRITE");
  }
  

  timer.Stop(); timer.Print();
  
  delete rl;//cleans everything
  
  return 0;
}
