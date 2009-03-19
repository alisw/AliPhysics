/*
  Macro for merging of the calibration classes:
  marain.ivanov@cern.ch
  //
  //
  Example usage:
  1. Make a list of the calibration files
  2. Load libraries
  3. Merge - the output is stored in the current directory
  //
  .x ~/rootlogon.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  .L $ALICE_ROOT/TPC/macros/CalibFileMerger.C+
  CalibFileMerger();
*/

#include <fstream>
#include "TFile.h"
#include "TObjArray.h"
#include "AliSysInfo.h"
#include "AliTPCcalibBase.h"

void CalibFileMerger(Int_t nmax=200,const char* filename = "mergelist.txt") {


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
    if (!carray) continue;
    if (!mergeArray) { mergeArray = carray; continue;}
    printf("Counter\t%d\tMerging file %s\n",counter,objfile.Data());
    for (Int_t icalib=0; icalib<mergeArray->GetEntries(); icalib++){
      AliTPCcalibBase *calib = ( AliTPCcalibBase *)mergeArray->At(icalib);
      //disable calib align
      if (calib->InheritsFrom("AliTPCcalibUnlinearity")) continue;
      //if (calib->InheritsFrom("AliTPCcalibCosmic")) continue;
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
    }
    AliSysInfo::AddStamp(objfile.Data(),counter,counter,counter);
    carray->SetOwner(kTRUE);
    carray->Delete();
    delete carray;
    counter++;
    currentFile.Close();
    if (counter>nmax) break;
  }

  in.close();
  TFile * output = new TFile("CalibObjects.root", "RECREATE");
  mergeArray->Write("TPCCalib",TObject::kSingleKey);
  output->Close();

  //

}


