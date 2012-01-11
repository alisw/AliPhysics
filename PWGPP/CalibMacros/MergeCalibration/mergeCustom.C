//
// macro to merge the files 
//
void mergeCustom(const char * inputList, const char *outputName, const char *accept=0, const char *reject=0, Bool_t singleKey=kTRUE){

  //
  //
  // 
  gROOT->Macro("$ALICE_ROOT/ANALYSIS/CalibMacros/Pass0/LoadLibraries.C");    
  TString hasAlien = gSystem->GetFromPipe("cat  calib.list | grep -c  alien");
  if (hasAlien.Length()>0){
    TGrid::Connect("alien://");
  }
  AliFileMerger merger;
  if (reject) {
    TString rstr(reject);
    TObjArray * arrayReject=rstr.Tokenize("++");
    for (Int_t i=0; i<arrayReject->GetEntries(); i++){
      merger.AddReject(arrayReject->At(i)->GetName());
    }    
  }
  if (accept){
    TString astr(accept);
    TObjArray * arrayAccept=astr.Tokenize("++");
    for (Int_t i=0; i<arrayAccept->GetEntries(); i++){
      merger.AddAccept(arrayAccept->At(i)->GetName());
    }
  }
  merger.IterTXT(inputList,outputName,kFALSE,singleKey);
}

