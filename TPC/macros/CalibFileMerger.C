/// \file CalibFileMerger.C
/// Macro for merging of the calibration classes:
/// \author marain.ivanov@cern.ch
///
///   Example usage:
///   1. Make a list of the calibration files
///   2. Load libraries
///   3. Merge - the output is stored in the current directory
///
/// ~~~{.cpp}
/// gSystem->AddIncludePath("-I$ALICE_ROOT/TPC")
/// gSystem->Load("libANALYSIS");
/// gSystem->Load("libTPCcalib");
/// .L $ALICE_ROOT/TPC/macros/CalibFileMerger.C+g 
///  CalibFileMerger();
/// ~~~

#include <fstream>
#include "TSystem.h"
#include "TFile.h"
#include "TObjArray.h"
#include "AliSysInfo.h"
#include "AliTPCcalibBase.h"

void CalibFileMerger(const char *oname, Int_t nmax=1000,const char* filename = "mergelist.txt") {

  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  TObjArray * mergeArray= 0;
  //
  // Open the input stream
  ifstream in;
  in.open(filename);
  
  // Read the input list of files 
  TString objfile;
  Int_t counter=0;
  while(in.good()) {
    in >> objfile;
    if (!objfile.Contains("root")) continue; // protection    
    TFile currentFile(objfile.Data());
    printf("Open file:Counter\t%d\tMerging file %s\n",counter,objfile.Data());
    TObjArray * carray =  (TObjArray*)currentFile.Get("TPCCalib");
    if (!carray) {
      carray = new TObjArray;
      TList *farr = currentFile.GetListOfKeys();
      if (!farr) continue;
      for (Int_t ical=0; ical<farr->GetEntries(); ical++){
	if (!farr->At(ical)) continue;
	TObject *obj = currentFile.Get(farr->At(ical)->GetName());
	if (obj) carray->AddLast(obj);
	AliSysInfo::AddStamp(farr->At(ical)->GetName(),1,ical,counter);
	// reading entries
      }
      if (carray->GetEntries()==0) continue;
    }
   
    if (!mergeArray) { mergeArray = carray; continue;}
    printf("Counter\t%d\tMerging file %s\n",counter,objfile.Data());
    for (Int_t icalib=0; icalib<mergeArray->GetEntries(); icalib++){
      AliTPCcalibBase *calib = ( AliTPCcalibBase *)mergeArray->At(icalib);
      //disable calib align
      if (!calib->InheritsFrom("AliTPCcalibBase")) continue;
      if (calib->InheritsFrom("AliTPCcalibTrigger")) continue;  // too big
      //
    //   if (counter>50&&counter%2>0){
// 	if (calib->InheritsFrom("AliTPCcalibCosmic")) continue;  // too big
//       }
//       if (counter>100&&counter%4>0){
// 	if (calib->InheritsFrom("AliTPCcalibCosmic")) continue;  // too big
//       }
//       if (counter>200&&counter%8>0){
// 	if (calib->InheritsFrom("AliTPCcalibCosmic")) continue;  // too big
//       }
      if (!calib) continue;
      printf("Merging\t%s\n",calib->GetName());
      TObjArray toMerge(1);
      TObject *objectMerge = carray->FindObject(calib->GetName());
      if (!objectMerge) {
	printf("Object not found\n");
	continue;
      }
      toMerge.AddAt(objectMerge,0);
      toMerge.SetOwner(kFALSE);
      calib->Merge(&toMerge);
      AliSysInfo::AddStamp(calib->GetName(),2,icalib,counter);
    }
    AliSysInfo::AddStamp(objfile.Data(),counter,counter,counter);
    if (carray){
      carray->SetOwner(kTRUE);
      carray->Delete();
      // for (Int_t i=0; i<carray->GetEntriesFast(); i++){
// 	TObject *o = carray->At(i);
// 	if (!o) continue;
// 	printf("Delete %s\n",o->GetName());
// 	AliSysInfo::AddStamp(o->GetName(),3,i,counter);
// 	delete o;
//       }
      delete carray;
    }
    counter++;
    currentFile.Close();
    if (counter>nmax) break;
  }

  in.close();
  if (mergeArray){
    for (Int_t icalib=0; icalib<mergeArray->GetEntries(); icalib++){
      TFile * output = new TFile(oname, "UPDATE");
      TObject * obj = mergeArray->At(icalib); 
      if (!obj) continue;
      obj->Write();
      delete obj;
      output->Close();
    }
  }
  //
  // problem to stop root
  // Make suicide
  //
  printf("Kill ourself\n");
  printf("Kill ourself\n");
  printf("Kill ourself\n");
  printf("Kill ourself\n");
  gSystem->Exec(Form("kill -9 %d",gSystem->GetPid()));
}


