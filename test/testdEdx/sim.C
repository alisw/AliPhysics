//
// parameter to take from config file
//
const char * recoStorage="local:///cvmfs/alice.gsi.de/alice/data/2010/OCDB";
Int_t run=0;


void ModifyRecoParam(TObjArray* recoArray, Bool_t useIonTail, Double_t crossTalkCorrection){
  //
  // Modify reco param - and store it in the OCDB in local directory
  //
  AliCDBManager * man  =  AliCDBManager::Instance();
  for (Int_t i=0; i<4; i++){
    AliTPCRecoParam* p = ( AliTPCRecoParam*)recoArray->At(i);
    p->SetUseIonTailCorrection(useIonTail);
    p->SetCrosstalkCorrection(crossTalkCorrection);
  }
  TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDBsim";
  AliCDBStorage*pocdbStorage = AliCDBManager::Instance()->GetStorage(localStorage.Data());  
  AliCDBMetaData *metaData= new AliCDBMetaData();
  AliCDBId*   id1=new AliCDBId("TPC/Calib/RecoParam/", man->GetRun(), AliCDBRunRange::Infinity());
  pocdbStorage->Put(recoArray, (*id1), metaData);
  AliCDBManager::Instance()->SetSpecificStorage("TPC/Calib/RecoParam/",localStorage.Data());
}


void sim(Int_t nev, Bool_t useIonTail, Double_t crossTalkCorrection) {
  gSystem->Load("liblhapdf");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libpythia6");
  gSystem->Load("libAliPythia6");
  gSystem->Load("libHIJING");
  gSystem->Load("libTHijing");
  gSystem->Load("libgeant321");

  if (nev<0){
    AliCDBManager * man = AliCDBManager::Instance();
    man->SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
    man->SetSpecificStorage("TPC/Calib/RecoParam/",recoStorage);
    man->SetRun(run);
    AliCDBEntry* e = man->Get("TPC/Calib/RecoParam/",run); // get default
    // modify content
    TObjArray* recoArray = (TObjArray*)e->GetObject();
    ModifyRecoParam(recoArray, useIonTail, crossTalkCorrection);
    return;
  }

  if (gSystem->Getenv("EVENT")) nev = atoi(gSystem->Getenv("EVENT")) ;   
  
  AliSimulation simulator;
  simulator.SetMakeSDigits("ITS TPC TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0");
  //simulator.SetMakeDigitsFromHits( "ITS TPC");
  simulator.SetWriteRawData("ALL","raw.root",kTRUE);

  simulator.SetDefaultStorage("local://$ALIROOT_OCDB_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
  TString localStorage = "local://"+gSystem->GetFromPipe("pwd")+"/OCDBsim";
  simulator.SetSpecificStorage("TPC/Calib/RecoParam/",localStorage.Data());
  
  simulator.SetRunQA(":") ; 
  
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
  //
  // Print the OCDB setup which we used
  //
  AliCDBManager * man = AliCDBManager::Instance();
  AliCDBEntry* ocdbEntry = man->Get("TPC/Calib/RecoParam/",run);
  TObjArray* recoArray = (TObjArray*)ocdbEntry->GetObject();
  for (Int_t i=0; i<4; i++){
    AliTPCRecoParam* recoParam = ( AliTPCRecoParam*)recoArray->At(i);
    printf("ipar=%d\t%d\t%f\n",i,recoParam->GetUseIonTailCorrection(), recoParam->GetCrosstalkCorrection());
  } 
  printf("End of the simulation\n");

   
}
