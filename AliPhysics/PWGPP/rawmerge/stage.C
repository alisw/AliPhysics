/*
  author: marian.ivanov@cern.ch
   
  This macro touch the files.
  In case files are on castor file staging is invoked
  
  Usage: 
  .L $ALICE_PHYSICS/../src/PWGPP/rawmerge/stage.C
  stage("uniqule.list",2,1);
  stage files from the unique.list
  // timeOut=2s for touching of  files

*/

void stage(const char *input="unique.list", Int_t timeOut=2, Int_t skip=1){
  //
  // Stage input files
  // Write sys info
  //
  gSystem->Setenv("XRDCLIENTMAXWAIT",Form("%d",timeOut));
  gEnv->SetValue("XNet.RequestTimeout", timeOut);
  gEnv->SetValue("XNet.ConnectTimeout", timeOut);
  gEnv->SetValue("XNet.TransactionTimeout", timeOut);
  TFile::SetOpenTimeout(timeOut);    
  TGrid::Connect("alien");  
  TString inputExpr=TString::Format("cat %s",input).Data();
  TObjArray * array = gSystem->GetFromPipe(inputExpr.Data()).Tokenize("\n");
  Int_t nfiles= array->GetEntries();
  TFile::SetOpenTimeout(timeOut);
  //
  for (Int_t jfile=0; jfile<nfiles; jfile+=skip){
    if (TString(array->At(jfile)->GetName()).Length()>10){
      TFile * f = TFile::Open(array->At(jfile)->GetName());
      AliSysInfo::AddStamp(array->At(jfile)->GetName(),jfile);
      delete f;
    }
  }						
  delete array;
}

/* 
  example cration of the file list:

  for a in ` alien_find  /alice/data/2012/LHC12h/000189406/raw root`; do echo alien://$a; done  > unique.list
  for a in ` alien_find  /alice/data/2012/LHC12h/000189231/raw root`; do echo alien://$a; done  > unique.list
  aliroot -b -q $ALICE_PHYSICS/../src/PWGPP/rawmerge/stage.C\(\"unique.list\",5\); 2>&1 | tee staging.log  

*/
