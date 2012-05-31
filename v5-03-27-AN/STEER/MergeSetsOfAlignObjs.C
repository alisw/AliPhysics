#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliAlignObjParams.h"
#include <TString.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TError.h>
#endif


void MergeSetsOfAlignObjs(const char* filename1, const char* filename2, const char* det="ITS")
{
  // example macro: building an array by merging the non-SSD entries
  // from one file (or OCDB entry) with the remaining SSD entries taken
  // from another file (or OCDB entry); the first two arguments can be local filenames
  // or URLs of the OCDB folders 
  //  
  const char* macroname = "MergeSetsOfAlignObjs";
  
  TClonesArray* array1 = 0;
  TClonesArray* array2 = 0;

  TString arName(det);
  arName+="AlignObjs";
  TString path(det);
  path+="/Align/Data";
  
  TString f1(filename1);
  TString f2(filename2);
  
  AliCDBStorage* stor1 = 0;
  AliCDBStorage* stor2 = 0;
  
  Bool_t fromOcdb1=kFALSE;
  Bool_t fromOcdb2=kFALSE;

  if(f1.Contains("alien://folder=") || f1.Contains("local://")) fromOcdb1=kTRUE;
  if(f2.Contains("alien://folder=") || f2.Contains("local://")) fromOcdb2=kTRUE;

  
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  
  if(fromOcdb1){
    stor1 = cdb->GetStorage(f1.Data());
    AliCDBEntry* entry = stor1->Get(path.Data(),0);
    array1 = (TClonesArray*) entry->GetObject();
  }else{
    TFile* filein1 = TFile::Open(f1.Data(),"READ");
    if(!filein1)
    {
      Info(macroname,Form("Unable to open file %s! Exiting ...", f1.Data()));
      return;
    }
    array1 = (TClonesArray*) filein1->Get(arName.Data());
  }
  if(array1){
    Info(macroname,Form("First array has %d entries", array1->GetEntriesFast()));
  }else{
    Info(macroname,"Unable to get first array! Exiting ...");
    return;
  }

  if(fromOcdb2){
    stor2 = cdb->GetStorage(f2.Data());
    AliCDBEntry* entry = stor2->Get(path.Data(),0);
    array2 = (TClonesArray*) entry->GetObject();
  }else{
    TFile* filein2 = TFile::Open(f2.Data(),"READ");
    if(!filein2)
    {
      Info(macroname,Form("Unable to open file %s! Exiting ...", f2.Data()));
      return;
    }
    array2 = (TClonesArray*) filein2->Get(arName.Data());
  }
  if(array2){
    Info(macroname,Form("Second array has %d entries", array2->GetEntriesFast()));
  }else{
    Info(macroname,"Unable to get second array! Exiting ...");
    return;
  }


  TClonesArray *mergedArr = new TClonesArray("AliAlignObjParams",3000);

  Info(macroname,"Merging objects for SPD and SDD from the first array ...");
  // SSD starts from 500
  for(Int_t i=0; i<500; i++)
  {
    (*mergedArr)[i] = (AliAlignObjParams*) array1->UncheckedAt(i);
  }
  
  Info(macroname,"Merging objects for SSD from the second array ...");
  for(Int_t i=500; i<array2->GetEntriesFast(); i++)
  {
    (*mergedArr)[i] = (AliAlignObjParams*) array2->UncheckedAt(i);
  }
  
  TString foutName("merged");
  foutName+=det;
  foutName+="Alignment.root";
  Info(macroname,Form("... in a single array into the file %s", foutName.Data()));
  TFile* fileout = TFile::Open(foutName.Data(),"RECREATE");
  fileout->cd();
  fileout->WriteObject(mergedArr,arName.Data(),"kSingleKey");
  fileout->Close();
  
  mergedArr->Delete();
  
}
