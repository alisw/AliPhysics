//
// Make event list
// Purpose: prepare  event list for the raw data filtering from the filtered trees
//  
// Input: 
//       root file or list of root files  with  filtered data trees (highPt,V0s,Laser, CosmicPairs, high dEdx )
// Output:
//       event list  to the terminal; input for the  makeEventList.sh script
// Usage:
// 
// Authors:
//    modifications:
//       mesut.arslandok@cern.ch
//       marian.ivanov@cern.ch
//

/*
   Example usage:

   //1.) define env varaibles - sourcing config file e.g:
   source $ALICE_PHYSICS/../src/PWGPP/rawmerge/makeEventList.config
   //2.) run aliroot
   .L $ALICE_PHYSICS/../src/PWGPP/rawmerge/makeEventList.C+
   makeEventList("alien:///alice/data/2012/LHC12h/000189406/pass1/PWGPP/AdHocCalibration/239_20150413-1059/FilterEvents_Trees.root",8.0,4.0)
   // makeEventList("/hera/alice/local/filtered/alice/data/2012/LHC12h/000189406/pass1/155/root_archive.zip#FilterEvents_Trees.root",8.0,4.0)
   
*/

#include <TFile.h>
#include "TSystem.h"
#include "TBranch.h"
#include "TGrid.h"
#include "TH1.h"
#include "TFrame.h"
#include "TROOT.h"
#include "TKey.h"
#include "TClass.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include "TCut.h"
#include "AliXRDPROOFtoolkit.h"
#include "AliSysInfo.h"
#include "TChain.h"
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;
using std::cout;
using std::setw;

// variables related to config file
TString rawDataPrefix;
TString commonPrefix;
TString alienPrefix="alien://";

Long_t  nEvents=0;
TFile  *fInputFile;
TString fInputFileName;
TString period;
Int_t   run=0;
Int_t   year=0;
Int_t   runType=0;

TTree *treeHighPt=NULL;
TTree *treeV0s=NULL;
TTree *treeLaser=NULL;
TTree *treeCosmicPairs=NULL;
TTree *treedEdx=NULL;

TString treeNameHighPt      = "highPt";
TString treeNameV0s         = "V0s";
TString treeNameLaser       = "Laser";
TString treeNameCosmicPairs = "CosmicPairs";
TString treeNamedEdx        = "dEdx";


TString scanOutputHighPt;      
TString scanOutputV0s;         
TString scanOutputLaser;      
TString scanOutputCosmicPairs;
TString scanOutputdEdx;

Int_t counterHighPtAll;   // number of all hightpt events
Int_t counterV0sAll;      // number of all V0s events
Int_t counterCosmicAll;   // number of all Cosmic events
Int_t counterLaserAll;    // number of all Laser events
Int_t counterdEdxAll;     // number of all dEdx events

Int_t counterHighPtSelected;  // number of selected hightpt events
Int_t counterV0sSelected;     // number of selected V0s events
Int_t counterCosmicSelected;  // number of selected Cosmic events
Int_t counterLaserSelected;   // number of selected Laser events
Int_t counterdEdxSelected;    // number of selected dEdx events

TCut cutHighPt="";
TCut cutV0s="";
TCut cutCosmic="";
TCut cutLaser="";
TCut cutdEdx;



// Helper Functions
void    ReadTrees();
Int_t   GuessRunNumber();
Int_t   FindNEventsPerTrigger(TTree *t, TCut selection);
void    SetDumpOutputFile(Int_t run);
void    DumpCounters();
void    MakeRawList(const char *scanOutput);


// Main Function
void makeEventList(TString fFileName="FilteredESDs.root", Double_t ptMinHighPt = 8., Double_t ptMinV0s = 3.)
{
  // 
  // Reads the filtered trees and creates the raw data list 
  // 

  rawDataPrefix = gSystem->GetFromPipe("echo $rawDataPrefix");
  commonPrefix  = gSystem->GetFromPipe("echo $commonPrefix");
  runType       = atoi(gSystem->Getenv("isCosmic"));

  // check alien connection
  if (fFileName.Contains("alien://")!=0) {
    TGrid * grid = TGrid::Connect("alien://");
    if (!grid){
      cout << " Error:: makeEventList.C  --> Alien not available " << endl;
      return;
    }
  }
  
  // Read input:  Filtered tree file and get the trees to be used for "Scan()"
  fInputFileName = fFileName; 
  ReadTrees();
  run = GuessRunNumber();
  SetDumpOutputFile(run);
  // Base cut for the high pt tracks:
  // 1.) rough pointing to the primary vertex 
  // 2.) reasonable q/pt resolution
  // 3.) at minimum one of the outer detectors - to clean pile-up
    
  // Analyse Hight pt tree
  if (treeHighPt)
  {
    if (treeHighPt->GetEntries()>0)
    { 
      // Dump the tree scan in to file
      printf("offlineTrigger: highPt\n");
      cutHighPt="abs(esdTrack.fdTPC)<4&&abs(esdTrack.fzTPC)<4&&esdTrack.fzTPC!=0&&sqrt(esdTrack.fC[14])<0.03&&((esdTrack.fFlags&0x4401)>0)";
      treeHighPt->GetUserInfo()->AddFirst(new TNamed("highPt","highPt"));
      treeHighPt->SetScanField(nEvents);
      //      tree->Scan("fileName.GetString():evtNumberInFile:triggerName:gid:evtTimeStamp:esdTrack.Pt()",Form("esdTrack.Pt()>%lf",ptMinHighPt),"col=.2f:8.d:8.d:130.s:15.lu:12.d");
      AliSysInfo::AddStamp("highPtFilterScanBegin",1,1);
      treeHighPt->Scan("fileName.GetString():evtNumberInFile:This->GetUserInfo()->At(0)->GetName():gid:evtTimeStamp:esdTrack.Pt()",Form("esdTrack.Pt()>%lf",ptMinHighPt)+cutHighPt,"col=130.s:14.d:30.s:20.lu:20.lu");
      AliSysInfo::AddStamp("highPtFilterScanEnd",1,2);
      
      // Process the Scan output further to prepare the raw data list
      MakeRawList(scanOutputHighPt.Data());
      counterHighPtAll      = FindNEventsPerTrigger(treeHighPt,"");
      counterHighPtSelected = FindNEventsPerTrigger(treeHighPt,cutHighPt);    
    } else {     
      counterHighPtAll      = 0;
      counterHighPtSelected = 0;
    }
  } 
  AliSysInfo::AddStamp("highPtFilterListEnd",1,3);
  
  // Analyse V0s tree
  if (treeV0s)
  {
    if (treeV0s->GetEntries()>0)
    { 
      // Dump the tree scan in to file
      printf("offlineTrigger: V0s\n");
      cutV0s="abs(track0.fzTPC)<20&&abs(track1.fzTPC)<20&&sqrt(track0.fC[14])<0.1&&sqrt(track1.fC[14])<0.03&&((track0.fFlags&0x4401)+(track1.fFlags&0x4401)>0)";
      treeV0s->GetUserInfo()->AddFirst(new TNamed("V0s","V0s"));
      treeV0s->SetScanField(nEvents);
      AliSysInfo::AddStamp("highV0sFilterScanBegin",2,1);
      treeV0s->Scan("fileName.GetString():evtNumberInFile:This->GetUserInfo()->At(0)->GetName():gid:evtTimeStamp",Form("v0.Pt()>%lf",ptMinV0s)+cutV0s,"col=130.s:14.d:30.s:20.lu:20.lu"); 
      AliSysInfo::AddStamp("highV0sFilterScanEnd",2,2);
      
      // Process the Scan output further to prepare the raw data list
      MakeRawList(scanOutputV0s.Data()); 
      counterV0sAll      = FindNEventsPerTrigger(treeV0s,"");
      counterV0sSelected = FindNEventsPerTrigger(treeV0s,cutV0s);              
    } else {     
      counterV0sAll      = 0;
      counterV0sSelected = 0;
    }
  }
  
  // Analyse Laser tree
  if (treeLaser)
  {
    if (treeLaser->GetEntries()>0)
    { 
      // Dump the tree scan in to file
      printf("offlineTrigger: Laser\n");
      treeLaser->GetUserInfo()->AddFirst(new TNamed("Laser","Laser"));
      treeLaser->SetScanField(nEvents);
      AliSysInfo::AddStamp("laserFilterBegin",3,1);
      treeLaser->Scan("fileName.GetString():evtNumberInFile:This->GetUserInfo()->At(0)->GetName():gid:evtTimeStamp","","col=130.s:14.d:30.s:20.lu:20.lu"); 
      AliSysInfo::AddStamp("laserFilterEnd",3,2);
      
      // Process the Scan output further to prepare the raw data list
      MakeRawList(scanOutputLaser.Data());    
      counterLaserAll      = FindNEventsPerTrigger(treeLaser,"");
      counterLaserSelected = FindNEventsPerTrigger(treeLaser,cutLaser);          
    } else {     
      counterLaserAll      = 0;
      counterLaserSelected = 0;
    }
  }
  
  // Analyse Laser tree
  if (treeCosmicPairs)
  {
    if (treeCosmicPairs->GetEntries()>0)
    { 
      cutCosmic="1";
      // Additional cuts for cosmics
      if (gSystem->Getenv("isCosmic")!=NULL &&strstr(gSystem->Getenv("isCosmic"),"0")!=0){
	treeCosmicPairs->SetAlias("alphaPrime","abs(t0.fAlpha)-pi/2");
	treeCosmicPairs->SetAlias("alphaPrimeFit","TMath::Gaus(alphaPrime,0,0.573+0)");
	treeCosmicPairs->SetAlias("alphaPrimeDownscale","alphaPrimeFit*rndm<0.05");  
	treeCosmicPairs->SetAlias("itsFiducial","min(abs(t0.fP[0]),abs(t1.fP[0]))<16&&min(abs(t0.fP[1]),abs(t1.fP[1]))<30");
	cutCosmic="itsFiducial || alphaPrimeDownscale";
      }
      else{
	TCut ptCut="abs(t0.fP[4])<0.33"; //cut on 1/pt < 0.33
	TCut cutDCA="abs(0.5*(t0.fD-t1.fD))>5&&abs(0.5*(t0.fD-t1.fD))<80"; //tracks crossing the inner field cage (80cm)
	TCut cutCross="t0.fOp.fP[1]*t1.fOp.fP[1]<0"; //tracks crossing central electrode
	cutCosmic= ptCut && cutDCA && cutCross;
      }
      // Dump the tree scan in to file
      printf("offlineTrigger: CosmicPairs\n");
      treeCosmicPairs->GetUserInfo()->AddFirst(new TNamed("CosmicPairs","CosmicPairs"));
      treeCosmicPairs->SetScanField(nEvents);
      AliSysInfo::AddStamp("cosmicFilterBegin",4,1);
      treeCosmicPairs->Scan("fileName.GetString():evtNumberInFile:This->GetUserInfo()->At(0)->GetName():gid:evtTimeStamp", cutCosmic, "col=130.s:14.d:30.s:20.lu:20.lu"); 
      AliSysInfo::AddStamp("cosmicFilterEnd",4,2);
      
      // Process the Scan output further to prepare the raw data list
      MakeRawList(scanOutputCosmicPairs.Data()); 
      counterCosmicAll      = FindNEventsPerTrigger(treeCosmicPairs,"");
      counterCosmicSelected = FindNEventsPerTrigger(treeCosmicPairs,cutCosmic);     
    } else {     
      counterCosmicAll      = 0;
      counterCosmicSelected = 0;
    }
  }
  
  // Analyse Laser tree
  if (treedEdx)
  {
    if (treedEdx->GetEntries()>0)
    { 
      // Dump the tree scan in to file
      printf("offlineTrigger: dEdx\n");
      treedEdx->GetUserInfo()->AddFirst(new TNamed("dEdx","dEdx"));
      treedEdx->SetScanField(nEvents);
      AliSysInfo::AddStamp("highdEdxFilterBegin",5,1);
      treedEdx->Scan("fileName.GetString():evtNumberInFile:This->GetUserInfo()->At(0)->GetName():gid:evtTimeStamp","","col=130.s:14.d:30.s:20.lu:20.lu"); 
      AliSysInfo::AddStamp("highdEdxFilterBegin",5,2);
      
      // Process the Scan output further to prepare the raw data list
      MakeRawList(scanOutputdEdx.Data());  
      counterdEdxAll      = FindNEventsPerTrigger(treedEdx,"");
      counterdEdxSelected = FindNEventsPerTrigger(treedEdx,cutdEdx);
    }
    else {     
      counterdEdxAll      = 0;
      counterdEdxSelected = 0;
    }
  }
 
 gSystem->Exec("cat rawEvents*.list >> event.list");
 DumpCounters();
 cout << "------------------------- Event list creation is successful -------------------------" << endl;

}
// -----------------------------------------------------------------------------------------------------
void ReadTrees()
{
  
  //
  // Read the trees from the .root file if the file is a list of root files then make the chain of files 
  //
  
  if (fInputFileName.Contains(".root")){
    fInputFile = TFile::Open(fInputFileName);
    cout << " make the tree " << endl;
    treeHighPt      = (TTree*)fInputFile->Get(treeNameHighPt);
    treeV0s         = (TTree*)fInputFile->Get(treeNameV0s);
    treeLaser       = (TTree*)fInputFile->Get(treeNameLaser);
    treeCosmicPairs = (TTree*)fInputFile->Get(treeNameCosmicPairs);  
    treedEdx        = (TTree*)fInputFile->Get(treeNamedEdx);      
  } else {
    cout << " make the chain " << endl;
    treeHighPt      = (TTree*)AliXRDPROOFtoolkit::MakeChainRandom(fInputFileName,treeNameHighPt     ,0,500000,0);  
    treeV0s         = (TTree*)AliXRDPROOFtoolkit::MakeChainRandom(fInputFileName,treeNameV0s        ,0,500000,0);  
    treeLaser       = (TTree*)AliXRDPROOFtoolkit::MakeChainRandom(fInputFileName,treeNameLaser      ,0,500000,0);  
    treeCosmicPairs = (TTree*)AliXRDPROOFtoolkit::MakeChainRandom(fInputFileName,treeNameCosmicPairs,0,500000,0); 
    treedEdx        = (TTree*)AliXRDPROOFtoolkit::MakeChainRandom(fInputFileName,treeNamedEdx       ,0,500000,0);   
  }
  
   treeHighPt->SetCacheSize(100000000);
   treeV0s->SetCacheSize(100000000);
   treeLaser->SetCacheSize(100000000);
   treeCosmicPairs->SetCacheSize(100000000);
   treedEdx->SetCacheSize(100000000);
}
// -----------------------------------------------------------------------------------------------------
Int_t GuessRunNumber()
{

  //
  // read run number from file
  //
  
  TTree * tree;
  if (treeHighPt && treeHighPt->GetEntries()>0) tree=treeHighPt;
  else if (treeV0s && treeV0s->GetEntries()>0) tree=treeV0s;
  else if (treeCosmicPairs && treeCosmicPairs->GetEntries()>0) tree=treeCosmicPairs;
  else if (treeLaser && treeLaser->GetEntries()>0) tree=treeLaser;
  else if (treedEdx && treedEdx->GetEntries()>0) tree=treedEdx;
  else return 0;
  
  // get run number 
  tree->Draw("runNumber>>htmp","","goff",100);
  TH1D * htmp = (TH1D*)tree->GetHistogram();
  return (Int_t)htmp->GetMean();
  
}
// -----------------------------------------------------------------------------------------------------
Int_t FindNEventsPerTrigger(TTree *t, TCut selection)
{

  //
  // guess number events of a ttree from the global id information
  //
  
  Int_t nEntries = t->Draw("gid>>his(100,0,1)",selection,"goff");
  std::map<int,int> gids;
  for (Int_t i=0; i<nEntries; i++) gids[long(t->GetV1()[i])]+=1;
  return gids.size();
}
// -----------------------------------------------------------------------------------------------------
void SetDumpOutputFile(Int_t run)
{

  //
  // Prepare the dump of output file
  //
  
  scanOutputHighPt      ="HighPt.list";
  scanOutputV0s         ="V0s.list";
  scanOutputLaser       ="Laser.list";
  scanOutputCosmicPairs ="CosmicPairs.list";
  scanOutputdEdx        ="dEdx.list";

  if (treeHighPt && treeHighPt->GetEntries()>0) {
    ((TTreePlayer*)(treeHighPt     ->GetPlayer()))->SetScanRedirect(true);
    ((TTreePlayer*)(treeHighPt     ->GetPlayer()))->SetScanFileName(scanOutputHighPt);
  }
  
  if (treeV0s && treeV0s->GetEntries()>0) {
    ((TTreePlayer*)(treeV0s        ->GetPlayer()))->SetScanRedirect(true);
    ((TTreePlayer*)(treeV0s        ->GetPlayer()))->SetScanFileName(scanOutputV0s);
  }
  
  if (treeLaser && treeLaser->GetEntries()>0) { 
    ((TTreePlayer*)(treeLaser      ->GetPlayer()))->SetScanRedirect(true);
    ((TTreePlayer*)(treeLaser      ->GetPlayer()))->SetScanFileName(scanOutputLaser);
  }
  
  if (treeCosmicPairs && treeCosmicPairs->GetEntries()>0) {
    ((TTreePlayer*)(treeCosmicPairs->GetPlayer()))->SetScanRedirect(true);
    ((TTreePlayer*)(treeCosmicPairs->GetPlayer()))->SetScanFileName(scanOutputCosmicPairs);
  }
  
  if (treedEdx && treedEdx->GetEntries()>0) { 
    ((TTreePlayer*)(treedEdx       ->GetPlayer()))->SetScanRedirect(true);
    ((TTreePlayer*)(treedEdx       ->GetPlayer()))->SetScanFileName(scanOutputdEdx);
  }

}
// -----------------------------------------------------------------------------------------------------
void MakeRawList(const char *scanOutput)
{
  
  //
  // Create a new list replacing the esd path with raw path of the chunk
  //

  // Create the rawlist output file
  ofstream rawlist;
  TString rawListName = Form("rawEvents_%s",scanOutput);
  rawlist.open(rawListName);   // rawEvents_<TriggerName>.list

  // open Tree::Scan output  
  std::ifstream file(scanOutput);
  std::string line;
  
  // Loop over all lines of the Tree::Scan output
  Int_t ncount=0;
  while (std::getline(file, line)) {
    
    TString sline = (line);
    if (!sline.Contains(".root")) continue;
    ncount++;
    
    // tokenize each line wrt "*" to get read of "*"
    TObjArray *objArr1  = sline.Tokenize("*");
    TString esdFilePath = (objArr1->At(1)->GetName());
        
    // Tokenize the esd file path wrt "/" to read run number and root file name
    TObjArray *objArr2  = esdFilePath.Tokenize("/");
    Int_t nElements = objArr2->GetLast();
    TString rootFile;
    Int_t rootFileIndex = 0;
    Int_t runNumberIndex = 0;
    for (Int_t i=0;i<(Int_t)objArr2->GetLast();i++){
      TString tmp = TString((objArr2->At(i))->GetName());
      if ( tmp.Contains("000") && tmp.Contains(".") ) {
        rootFile = tmp.Append(".root");
        rootFileIndex = i;
      }
      if ( tmp.Contains("000") && !(tmp.Contains(".")) ) {
        runNumberIndex = i;
      }
    }
         
    // Retrieve runNumber, period, and year from the esdfilePath
    TString runStr  = TString((objArr2->At(runNumberIndex))->GetName());
    TString yearStr = TString((objArr2->At(runNumberIndex-2))->GetName());
    period    = TString((objArr2->At(runNumberIndex-1))->GetName());
    year      = atoi((objArr2->At(runNumberIndex-2))->GetName());
    
    // prepare the raw file path
    rootFile.Prepend(Form("%s%s/%s/%s/%s/raw/",rawDataPrefix.Data(),commonPrefix.Data(),yearStr.Data(),period.Data(),runStr.Data()));
    
    // Other variables needed for raw list
    Int_t   eventNumber = atoi(objArr1->At(2)->GetName());
    TString trigName    = (TString)(objArr1->At(3)->GetName());
    TString gid         = (objArr1->At(4)->GetName());
    TString timeStamp   = (objArr1->At(5)->GetName());
     
    // Dump info into final raw list
    rawlist  << rootFile << setw(15) << eventNumber << setw(1) << trigName << setw(1) << gid << setw(1) << timeStamp << endl;
  } 
  
  rawlist.close();
  
  // Avoid duplicated event numbers in file
  gSystem->Exec(Form("{ rm %s && uniq > %s; } < %s ",rawListName.Data(),rawListName.Data(),rawListName.Data()));
}
// -----------------------------------------------------------------------------------------------------
void DumpCounters(){
  
  //
  //  Dump event counters and necessary info
  //
  cout << "-----------------------------------------------" << endl; 
  cout << "KeyValue_Index0="  << run << endl;   
  cout << "KeyValue_Index1="  << fInputFileName << endl;     
  cout << "KeyValue_year="    << year << endl;    
  cout << "KeyValue_period="  << period << endl;      
  cout << "KeyValue_runType=" << runType << endl;      
  cout << "-----------------------------------------------" << endl; 
  cout << "KeyValue_CounterHighPtAll="  << counterHighPtAll << endl; 
  cout << "KeyValue_CounterHighPtSelected="  << counterHighPtSelected << endl;  
  if (counterHighPtAll) cout << "KeyValue_HighPtRatio=" << counterHighPtSelected/Double_t(counterHighPtAll) << endl; 
  else  cout << "KeyValue_HighPtRatio=0" << endl;
  cout << "-----------------------------------------------" << endl; 
  cout << "KeyValue_CounterV0sAll="     << counterV0sAll << endl;  
  cout << "KeyValue_CounterV0sSelected="     << counterV0sSelected << endl;     
  if (counterV0sAll) cout << "KeyValue_V0sRatio=" << counterV0sSelected/Double_t(counterV0sAll) << endl;     
  else  cout << "KeyValue_V0sRatio=0" << endl;
  cout << "-----------------------------------------------" << endl; 
  cout << "KeyValue_CounterCosmicAll="  << counterCosmicAll << endl; 
  cout << "KeyValue_CounterCosmicSelected="  << counterCosmicSelected << endl;
  if (counterCosmicAll) cout << "KeyValue_CosmicRatio=" << counterCosmicSelected/Double_t(counterCosmicAll) << endl;  
  else  cout << "KeyValue_CosmicRatio=0" << endl;
  cout << "-----------------------------------------------" << endl; 
  cout << "KeyValue_CounterLaserAll="   << counterLaserAll << endl;
  cout << "KeyValue_CounterLaserSelected="   << counterLaserSelected << endl;   
  if (counterLaserAll) cout << "KeyValue_LaserRatio=" << counterLaserSelected/Double_t(counterLaserAll) << endl;   
  else  cout << "KeyValue_LaserRatio=0" << endl;
  cout << "-----------------------------------------------" << endl; 
  cout << "KeyValue_CounterdEdxAll="    << counterdEdxAll << endl;     
  cout << "KeyValue_CounterdEdxSelected="    << counterdEdxSelected << endl;  
  if (counterdEdxAll) cout << "KeyValue_dEdxRatio=" << counterdEdxSelected/Double_t(counterdEdxAll) << endl;    
  else  cout << "KeyValue_dEdxRatio=0" << endl;
}
