//--------------------------------------------------------------------------
// Macro for QA monitoring.
//
// In case it is not run with full aliroot, it needs the following libraries to compile:
//  - libSTEERBase.so
//  - libESD.so
//  - libAOD.so
//  - libANALYSIS.so
//  - libANALYSISalice.so
//  - libCORRFW.so
//  - libPWG3muon.so
//
//  TString includePath = "-I${ALICE_ROOT}/PWG3/base/ ";  gSystem->SetIncludePath(includePath.Data());


// The macro reads the PWG1 QA train output, produces a merged root files for the full period
// for event and track counters as well as separate root files run per run with all MUON_TRK related histograms.
// The results is stored under the directory "results". Then use PlotMUONQA.C, to draw QA histograms.
//
// Author: Cynthia Hadjidakis - IPN Orsay
//--------------------------------------------------------------------------

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>

// ROOT includes
#include "TROOT.h"
#include "TMath.h"
#include "TGrid.h"
#include "TGridCollection.h"
#include "TGridResult.h"
#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"

// ALIROOT includes
#include "AliCounterCollection.h"

#endif

TObjArray * GetListOfFiles(const char* baseDir, const char * trainName, const char* inFile);
TObjArray * GetListOfRuns(const char* runList, TObjArray *&listoffiles);

// .x TransferMuonQATrain.C("alien:///alice/data/2010/LHC10e","QA50","/Users/cynthia/Documents/alice/data/MuonQA/LHC10e/pass2/runlist_period3_test3.txt")
// .x TransferMuonQATrain.C("alien:///alice/cern.ch/user/s/suire/LHC10hV2116/output","","mylist.txt","AnalysisResults.root")
//--------------------------------------------------------------------------
Bool_t TransferMuonQATrain_v4(const char* baseDir, const char * trainName, const char* runList,char* inputFile = "QAresults.root")
{
#if defined(__CINT__) && !defined(__MAKECINT__)
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWG3base");
  gSystem->Load("libPWG3muon");
#endif
  
  TString sbaseDir = baseDir;
  if (sbaseDir.Contains("alien:") && !TGrid::Connect("alien://")) {
    Error("MuonQATrain","cannot connect to grid");
    return 0;
  }	
  
  //----------------------------------------------------------- //
  //          Build the list of files, the list of runs         //
  //					to be processed
  //----------------------------------------------------------- //
  
  TObjArray *listoffiles = (TObjArray*) GetListOfFiles(baseDir,trainName,inputFile);
  if(!listoffiles) return kFALSE;
  TObjArray *runs = (TObjArray*) GetListOfRuns(runList,listoffiles);
  if(!runs||!listoffiles){
    Error("TransferMuonQATrain","cannot get a list of selected runs or files");
    return kFALSE;
  }
  printf("TransferMuonQATrain: Files found to be processed = %d\n",runs->GetEntriesFast());
  if(runs->GetEntriesFast()==0)	return kFALSE;
  
  //-------------------------------------//
  //    Loop over the files on grid      //
  //-------------------------------------//
  
  TFile *file = 0;
  TIter next0(listoffiles);
  TObject *nextfile;
  TString snextfile;
  Int_t nFiles = -1;
  AliCounterCollection *trackCounters = 0;
  AliCounterCollection *eventCounters = 0;
  AliCounterCollection *mergedTrackCounters = 0;
  AliCounterCollection *mergedEventCounters = 0;	
  TObjArray*  outputList; 
  TObjArray*  outputListExpert;
  TObjArray*  outputListNorm;  
  TFile	outputHistoFile;	
  TString fileName;
  
  TString command = "mkdir results";
  cout<<"Shell command = "<<command<<endl;
  gSystem->Exec(command.Data());
  
  while ((nextfile=next0())) {
    
    nFiles++;
    snextfile = nextfile->GetName();
    //Open the file
    file = TFile::Open(snextfile.Data());
    if(!file) continue;
    
    outputList = (TObjArray *)file->Get("MUON_QA/general1");
    outputListExpert = (TObjArray *)file->Get("MUON_QA/expert");
    outputListNorm = (TObjArray *)file->Get("MUON_QA/general2");
    
    trackCounters = (AliCounterCollection *) file->Get("MUON_QA/trackCounters");
    eventCounters = (AliCounterCollection *) file->Get("MUON_QA/eventCounters");
    
    if(!trackCounters || !eventCounters){
      Error("TransferMuonQATrain","Objects not found for that file");
      continue;
    }
    
    //-------------------------------------//
    //    Merge the AliCounterCollection
    //-------------------------------------//
    if(nFiles==0){
      mergedTrackCounters = (AliCounterCollection*) trackCounters->Clone();
      mergedEventCounters = (AliCounterCollection*) eventCounters->Clone();
    }		
    else{
      mergedTrackCounters->Add(trackCounters);
      mergedEventCounters->Add(eventCounters);
    }
    
    
    //-------------------------------------//
    //     Save the MUONQA histos in AnalysisResults.root for each run number
    //-------------------------------------//
    fileName = "AnalysisResults.root";
    outputHistoFile.Open(fileName,"recreate");
    new TDirectoryFile("MUON_QA","MUON_QA");
    outputHistoFile.Cd("MUON_QA");
    
    TDirectory *dir0 = (TDirectory*) file->GetDirectory("MUON_QA");
    TObjArray* general1 = static_cast<TObjArray*>(dir0->FindObjectAny("general1"));
    TObjArray* expert = static_cast<TObjArray*>(dir0->FindObjectAny("expert"));
    TObjArray* general2 = static_cast<TObjArray*>(dir0->FindObjectAny("general2"));
    
    if(general1) general1->Write("general1",TObject::kSingleKey);
    if(general2) general2->Write("general2",TObject::kSingleKey);		
    if(expert) expert->Write("expert",TObject::kSingleKey);
    outputHistoFile.Close();
    
    command = "mkdir -p results/";
    command+=((TObjString*)runs->UncheckedAt(nFiles))->GetString();
    cout<<"Shell command = "<<command<<endl;
    gSystem->Exec(command.Data());
    command=" mv AnalysisResults.root results/";
    command+=((TObjString*)runs->UncheckedAt(nFiles))->GetString();
    command+= "/.";
    cout<<"Shell command = "<<command<<endl;
    gSystem->Exec(command.Data());
  } //end of loop over files
  
  //-------------------------------------//
  //      Save the AliCounterCollection in MergedAnalysisResults.root
  //-------------------------------------//
  
  TFile outputFile("MergedAnalysisResults.root","recreate");
  new TDirectoryFile("MUON_QA","MUON_QA");
  outputFile.Cd("MUON_QA");
  
  mergedTrackCounters->Write();
  mergedEventCounters->Write();
  outputFile.Close();
  
  command = "mv MergedAnalysisResults.root results/.";
  cout<<"Shell command = "<<command<<endl;
  gSystem->Exec(command.Data());
  
  return kTRUE;
}

TObjArray * GetListOfRuns(const char* runList, TObjArray *&listoffiles)
{

  TObjArray * runs = new TObjArray();
  runs->SetOwner();
  
  if(!runList){
    Error("GetListOfruns","runList is not defined... exit");
    return 0;
  }
  else {
    // only the ones in the runList
    ifstream inFile(runList);
    if (!inFile.is_open()) {
      Error("GetListOfRuns",Form("unable to open file %s", runList));
      return 0;
    }
    
    TString currRun;
    while (!inFile.eof()) {
      currRun.ReadLine(inFile, kTRUE);
      if (currRun.IsNull()) continue;
      if (!currRun.IsDigit()) {
	Error("GetListOfRuns","invalid run number: %s", currRun.Data());
	return 0;
      }
      runs->AddLast(new TObjString(Form("%09d", currRun.Atoi())));
    }
    
    inFile.close();
  }
	
  printf("GetListOfRuns: Nr of runs in the runlist = %d and in the list of files = %d\n",runs->GetEntriesFast(),listoffiles->GetEntriesFast());
  
  if(runList && listoffiles){	
    //Filter the selected runs and modify listoffiles
    TObjArray*  runsFound = new TObjArray();
    runsFound->SetOwner();	
    
    //filter the selected runs	
    TIter next0(listoffiles);
    TObject *nextfile;
    TString snextfile;
		Int_t isFound = 0;
		
    TObjArray *listoffilestmp = new TObjArray();	
    listoffilestmp->SetOwner();
    while ((nextfile=next0())) {//loop over files found on alien
      snextfile = nextfile->GetName();
			isFound = 0;
      for ( Int_t irun=0; irun<runs->GetEntriesFast(); irun++ ) { //loop over selected runs
	TString run = ((TObjString*)runs->UncheckedAt(irun))->GetString();
	if(snextfile.Contains(run)){
	  listoffilestmp->Add(nextfile);
	  runsFound->AddLast(new TObjString(Form("%09d", run.Atoi())));
		isFound=1;
	}
      }
	if(isFound==0) printf("GetListOfRuns: run = %s not found in the list of files.... continue....\n",snextfile.Data());	

    }
    runs = runsFound;
    listoffiles->Clear();
    listoffiles = (TObjArray*) listoffilestmp->Clone();
    
    printf("GetListOfRuns Nr of selected runs corresponding to the list of files = %d \n",runs->GetEntriesFast());		
  }
  
  return runs;
  
}

TObjArray* GetListOfFiles(const char* baseDir, const char * trainName, const char* inFile)
{
  
  TString sbaseDir = baseDir;
  TString strainName = trainName;
  TString inputFile = inFile;
  TString command;
  
  if(!sbaseDir.Contains("alien://")){
    Error("GetListOfFiles","Not implemented for files not on alien-->exit");
    return 0;
  }
  
  sbaseDir.ReplaceAll("alien://", "");
  
  TObjArray *listoffiles = new TObjArray();
  
  if (sbaseDir.Contains(".xml")) {
    // Read files pointed by the xml 
    TGridCollection *coll = (TGridCollection*)gROOT->ProcessLine(Form("TAlienCollection::Open(\"%s\");", sbaseDir.Data()));
    if (!coll) {
      ::Error("GetListOfFiles", "Input XML collection empty.");
      return 0;
      }
    // Iterate grid collection
    while (coll->Next()) {
      TString fname = gSystem->DirName(coll->GetTURL());
      fname += "/";
      fname += inputFile;      
      listoffiles->Add(new TNamed(fname.Data(),""));
    }   
  }
  else {   
    command = Form("find %s/ *%s/%s", sbaseDir.Data(), strainName.Data(), inputFile.Data());
    printf("command: %s\n", command.Data());
    TGridResult *res = gGrid->Command(command);
    if (!res) {
      ::Error("GetListOfFiles","No result for the find command\n");
			cout << sbaseDir.Data() << "   " <<  strainName.Data() << endl ; 
      delete listoffiles;
      return 0;
    }     
    TIter nextmap(res);
    TMap *map = 0;
    while ((map=(TMap*)nextmap())) {
      TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
      if (!objs || !objs->GetString().Length()) {
	// Nothing found - skip this output
	delete res;
	delete listoffiles;
	return 0;
      }
      listoffiles->Add(new TNamed(objs->GetName(),""));
    }
    delete res;
  }
  if (!listoffiles->GetEntries()) {
    ::Error("GetListOfFiles","No files from the find command=%s\n",command.Data());
      delete listoffiles;
      return 0;
  }     
  else printf("GetListOfFiles: Number of files found %d\n",listoffiles->GetEntries());
  
  return listoffiles;
  
}
