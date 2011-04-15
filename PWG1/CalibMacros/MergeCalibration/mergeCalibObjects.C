/*
  Merge calibration entry macro:
  
  Example usage:
  
  .L $ALICE_ROOT/ANALYSIS/CalibMacros/mergeCalibObjects.C
  mergeCalibObjects()
  
*/

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <fstream>
#include "TSystem.h"
#include "TFile.h"
#include "TObjArray.h"
#include "AliSysInfo.h"
#include "AliTPCcalibBase.h"
#include "TH1F.h"
#include "TMethodCall.h"
#include "TGrid.h"
#include "TGridResult.h"

#endif

void IterTXT( const char * fileList="calib.list",Bool_t separate=kFALSE);
void IterAlien(const char* outputDir, const char* outputFileName="AliESDfriends_v1.root");
void Merge(TFile* fileIn, TObjArray * array);
void StoreResults(TObjArray * array);
void StoreSeparateResults(TObjArray * array);

void mergeCalibObjects( const char * fileList="calib.list",Bool_t separate=kFALSE){

  // main function

  IterTXT(fileList,separate );
}

void LoadLib(){

  // Loading the necessary libraries

  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib"); 
  TH1::AddDirectory(0);
}


void IterXML(){ 

  // iterating over the files coming from an XML collection 
  // to be implemented

  LoadLib();
}


void IterAlien(const char* outputDir, const char* outputFileName){

	//
	// Merge the files coming out of the calibration job
	// 

	LoadLib();
	TObjArray * mergeArray= new TObjArray;

	TString outputFile(outputFileName);
	TString command;
	// looking for files to be merged in the output directory
	command = Form("find %s/ *%s", outputDir, outputFile.Data());
	printf("command: %s\n", command.Data());
	TGrid::Connect("alien://");
	TGridResult *res = gGrid->Command(command);
	if (!res) return;
	TIter nextmap(res);
	TMap *map = 0;
	// loop over the results
	while((map=(TMap*)nextmap())) {
		// getting the turl
		TObjString *objs = dynamic_cast<TObjString*>(map->GetValue("turl"));
		if (!objs || !objs->GetString().Length()) {
			// Nothing found - skip this output
			delete res;
			break;
		} 
		printf("looking for file %s\n",(objs->GetString()).Data());
		TFile* currentFile=TFile::Open((objs->GetString()).Data());
		Merge(currentFile, mergeArray);
	}
	Bool_t separate = kFALSE;
	if (separate) 
		StoreSeparateResults(mergeArray);
	else
		StoreResults(mergeArray);
	delete mergeArray;
}
		


void IterTXT( const char * fileList, Bool_t separate){

  // Merge the files indicated in the list - fileList
  // ASCII file opition example: 
  // find `pwd`/ | grep AliESDfriends_v1.root > calib.list

  LoadLib();
  TObjArray * mergeArray= new TObjArray;
  
  // Open the input stream
  
  ifstream in;
  in.open(fileList);
  // Read the input list of files 
  TString objfile;
  Int_t counter=0;
  while(in.good()) {
    in >> objfile;
    if (!objfile.Contains("root")) continue; // protection    
    printf("Open file:Counter\t%d\tMerging file %s\n",counter++,objfile.Data());
    TFile currentFile(objfile.Data());
    Merge(&currentFile, mergeArray);
  }
  if (separate) 
    StoreSeparateResults(mergeArray);
  else
    StoreResults(mergeArray);
  delete mergeArray;
}

void StoreResults(TObjArray * array){

  // Storing the results in one single file

  TFile *f = new TFile("CalibObjects.root","recreate");
  for (Int_t i=0; i<array->GetEntries(); i++){
    TObject *object0 = array->At(i);
    if (!object0) continue;
    object0->Write();
  }
  f->Close();
  delete f;
}


void StoreSeparateResults(TObjArray * array){

  // Store the results in separate files (one per object)

  for (Int_t i=0; i<array->GetEntries(); i++){
    TObject *object0 = array->At(i);
    if (!object0) continue;
    TFile *f = new TFile(Form("CalibObjects_%s.root",object0->GetName()),"recreate");
    object0->Write();
    f->Close();
    delete f;
  }
}



void Merge(TFile* fileIn, TObjArray * array){
  
  // Merging procedure
  
  TObjArray *carray = new TObjArray;   //array of the objects inside current file
  carray->SetOwner(kTRUE);

  // load all objects to  memory
  
  TList *farr = fileIn->GetListOfKeys();
  if (!farr) return;
  for (Int_t ical=0; ical<farr->GetEntries(); ical++){
    if (!farr->At(ical)) continue;
    TObject *obj = fileIn->Get(farr->At(ical)->GetName());
    if (obj) carray->AddLast(obj);
    AliSysInfo::AddStamp(farr->At(ical)->GetName(),1,ical);  
  }

  if (carray->GetEntries()==0) return;
  TMethodCall callEnv;

  for (Int_t i=0; i<carray->GetEntries(); i++){

    TObjArray *templist = new TObjArray(1);
    templist->SetOwner(kFALSE);
    TObject *currentObject = carray->At(i);
    if (!currentObject) continue;
    printf("%s\n",currentObject->GetName());
    callEnv.InitWithPrototype(currentObject->IsA(), "Merge", "TCollection*");
    if (!callEnv.IsValid()) {continue;}
    TString oname=currentObject->GetName();
    if (oname.Contains("esdFriendTree")) continue;
    TObject *mergedObject = array->FindObject(currentObject->GetName());
    if (!mergedObject) {
      array->AddLast(currentObject);
      carray->RemoveAt(i);
      continue;
    }
    templist->AddLast(currentObject);
    callEnv.SetParam((Long_t) templist);
    callEnv.Execute(mergedObject);
    delete templist;
  }
  delete carray;
}
