/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//  marian.ivanov@cern.ch
//  Utilities for file merging.
//  Additional functionality on top of the standard TFileMerger:
//
//  1. Possibility to Set the reject/accept list.
//     1.a)  Only entries selected in accept list are merged. By default all entries are selected
//           use AddAccept 0 to specify your desired entry
//     1.b)  Entries selected in reject list are not merged. By default the reject list is empty.
//
//  2. syswatch.log is created diring mergin procedure. 
//     Memeory consumption - for reading and for merging can be monitored


/*
  Usage:
  // Libraries for all classes to be merged should be loaded before using the class
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSIScalib");
  gSystem->Load("libTPCcalib"); 
  TH1::AddDirectory(0);

  //Example usage starting from the input data list in text file:
  //
  AliFileMerger merger;
  merger.AddReject("esdFriend");
  merger.IterTXT("calib.list","CalibObjects.root",kFALSE);
  //

*/
//////////////////////////////////////////////////////////////////////////
 

#include <fstream>
#include "TSystem.h"
#include "TFile.h"
#include "TGrid.h"
#include "TGridResult.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TMethodCall.h"

#include "AliSysInfo.h"
#include "AliFileMerger.h"

ClassImp(AliFileMerger)

////////////////////////////////////////////////////////////////////////

AliFileMerger::AliFileMerger():
  TNamed(),
  fRejectMask(0),
  fAcceptMask(0)
{
  //
  // Default constructor
  //
}

//______________________________________________________________________

AliFileMerger::AliFileMerger(const char* name):
  TNamed(name,name),
  fRejectMask(0),
  fAcceptMask(0)
{
  //
  // 
  //
}


void AliFileMerger::IterAlien(const char* outputDir, const char* outputFileName, const char* pattern){

  //
  // Merge the files coming out of the calibration job
  // 
  
  TObjArray * mergeArray= new TObjArray;
  
  TString outputFile(outputFileName);
  TString command;
  // looking for files to be merged in the output directory
  command = Form("find %s/ *%s", outputDir, pattern);
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
    StoreSeparateResults(mergeArray,outputFileName);
  else
    StoreResults(mergeArray,outputFileName);
  delete mergeArray;
}



void AliFileMerger::IterTXT( const char * fileList,  const char* outputFileName, Bool_t separate){
  
  // Merge the files indicated in the list - fileList
  // ASCII file option example: 
  // find `pwd`/ | grep AliESDfriends_v1.root > calib.list
  
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
    StoreSeparateResults(mergeArray, outputFileName);
  else
    StoreResults(mergeArray, outputFileName);
  delete mergeArray;
}

void AliFileMerger::StoreResults(TObjArray * array, const char* outputFileName){
  //
  // Storing the results in one single file
  //
  TFile *f = new TFile(outputFileName,"recreate");
  for (Int_t i=0; i<array->GetEntries(); i++){
    TObject *object0 = array->At(i);
    if (!object0) continue;
    object0->Write();
  }
  f->Close();
  delete f;
}


void AliFileMerger::StoreSeparateResults(TObjArray * array, const char* outputFileName){
  //
  // Store the results in separate files (one per object)
  //
  for (Int_t i=0; i<array->GetEntries(); i++){
    TObject *object0 = array->At(i);
    if (!object0) continue;
    TFile *f = new TFile(Form("%s_%s.root",outputFileName,object0->GetName()),"recreate");
    object0->Write();
    f->Close();
    delete f;
  }
}



void AliFileMerger::Merge(TFile* fileIn, TObjArray * array){
  //
  // Merging procedure
  //
  static Int_t counter=-1;
  counter++;
  TObjArray *carray = new TObjArray;   //array of the objects inside current file
  carray->SetOwner(kTRUE);
  
  // load all objects to  memory
  
  TList *farr = fileIn->GetListOfKeys();
  if (!farr) return;
  for (Int_t ical=0; ical<farr->GetEntries(); ical++){
    if (!farr->At(ical)) continue;
    TString name(farr->At(ical)->GetName());
    if (!IsAccepted(name)) continue;                        // skip not accepted entries
    TObject *obj = fileIn->Get(name.Data());
    if (obj) carray->AddLast(obj);
    AliSysInfo::AddStamp(name.Data(),1,ical,counter);  
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
    TObject *mergedObject = array->FindObject(currentObject->GetName());
    if (!mergedObject) {
      array->AddLast(currentObject);
      carray->RemoveAt(i);
      continue;
    }
    templist->AddLast(currentObject);
    callEnv.SetParam((Long_t) templist);
    callEnv.Execute(mergedObject);
    AliSysInfo::AddStamp(currentObject->GetName(),2,i,counter);  
    delete templist;
  }
  delete carray;
}

Bool_t AliFileMerger::IsAccepted(TString name){
  //
  // Accept/reject logic
  // name - name of the entry
  //
  //  if fAcceptMask specified   - entry has to be in list of selected
  //  if fRejectMask speciefied  - entry with name speciief in the list are rejected 
  //
  Bool_t accept=kTRUE;
  if (fAcceptMask){
    //
    accept=kFALSE;
    for (Int_t iaccept=0; iaccept<fAcceptMask->GetEntries(); iaccept++){
      if (name.Contains(fAcceptMask->At(iaccept)->GetName())) accept=kTRUE;   // entry was selected
    }
  }
  if (!accept) return kFALSE;

  if (fRejectMask){
    //
    for (Int_t ireject=0; ireject<fRejectMask->GetEntries(); ireject++){
      if (name.Contains(fRejectMask->At(ireject)->GetName())) accept=kFALSE;   // entry was rejected
    }
  }
  return accept;
}




void AliFileMerger::AddReject(const char *reject){
  //
  // add reject string to the list of entries to be rejected for merging
  //
  if (!fRejectMask) fRejectMask = new TObjArray;
  fRejectMask->AddLast(new TObjString(reject));
}
void AliFileMerger::AddAccept(const char *accept){
  //
  // add reject string to the list of entries to be rejected for merging
  //
  if (!fAcceptMask) fAcceptMask = new TObjArray;
  fAcceptMask->AddLast(new TObjString(accept));


}

