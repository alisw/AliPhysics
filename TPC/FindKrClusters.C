/// \file FindKrClusters.C
/// \brief Macro to create array of clusters from TPC digits
///
/// Input files:
///
/// * galice.root
/// * digits.root: file with digits - usualy use link to galice.root
///   in splitted mode - neccesary to create link to proper file
///
/// Output file:
///
/// * TPC.RecPoints.root
///
/// Warning - if cluster file AliTPCclusters.root already exist - macro exit and don't produce anything
///
/// \author A.Matyja amatyja@cern.ch

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

  // Load DataBase
  //char *ocdbpath ="local:///afs/cern.ch/alice/tpctest/OCDB";
  //char *ocdbpath ="local:///home/matyja/baza/OCDB";
  char *ocdbpath ="local:///data/baza/OCDB";
  if (ocdbpath==0){
    ocdbpath="alien://folder=/alice/data/2007/LHC07w/OCDB/";
  }
  printf("OCDB PATH = %s\n",ocdbpath); 
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(ocdbpath);
  man->SetRun(0);

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

  TTree *output_tree;
  
  AliTPCclustererKr *clusters = new AliTPCclustererKr();
  clusters->SetParam(param);
  clusters->SetOutput(output_tree);

  clusters->SetMinAdc(3);//signal threshold (everything below is treated as 0)
  clusters->SetMinTimeBins(2);//number of neighbouring timebins
  clusters->SetMaxPadRangeCm(5.);//distance of the cluster center to the center of a pad (in cm)
  clusters->SetMaxRowRangeCm(5.);//distance of the cluster center to the center of a padrow (in cm)
  clusters->SetMaxTimeRange(7.);//distance of the cluster center to the max time bin on a pad (in tackts)
  //ie. fabs(centerT - time)<7

  clusters->SetIsolCut(3);//set isolation cut threshold
  clusters->SetValueToSize(3.1);//cut reduce peak at 0

  Int_t nevmax=rl->GetNumberOfEvents();//number of events in run
  for(Int_t nev=0;nev<nevmax /*&& nev<1*/ ;nev++){
    rl->GetEvent(nev);
    
    TTree* input_tree= tpcl->TreeD();//tree with digits
    if (input_tree == 0x0){
      cerr << "Can not get TreeD for event " <<nev<<endl;
      continue;
    }
    digarr->ConnectTree(input_tree);
    clusters->SetInput(input_tree);
    clusters->SetDigArr(digarr);
    cout<<"Processing event "<<nev<<endl;
    clusters->FinderIO();

  }
  delete clusters;
  timer.Stop(); timer.Print();
  
  delete rl;//cleans everything
  
  return 0;
}
