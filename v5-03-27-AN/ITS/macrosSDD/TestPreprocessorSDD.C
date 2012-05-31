
void TestPreprocessorSDD(Char_t *optRunType="PULSER"){
  // macro to Test SDD preprocessor
  // needs:
  // - 4 tar files with simulated output of PULSER DA
  // - 4 tar files with simulated output of INJECTOR DA
  // - 1 root file with simulated output of DCS
  // all these files can be found on 
  gSystem->Load("$ALICE_ROOT/SHUTTLE/TestShuttle/libTestShuttle.so");

  // The parameters are run, startTime, endTime
  AliTestShuttle* shuttle = new AliTestShuttle(7, 0, 1);


  AliTestShuttle::SetMainCDB("local:///home/prino/alice/SDD/Calibration/preprocessor/OCDB");
  AliTestShuttle::SetMainRefStorage("local:///home/prino/alice/SDD/Calibration/preprocessor/OCDB");

  printf("Test OCDB storage Uri: %s\n", AliShuttleInterface::GetMainCDB().Data());
  printf("Test Reference storage Uri: %s\n", AliShuttleInterface::GetMainRefStorage().Data());

  // DCS input file
  TFile *fil=new TFile("DCSAliasMap_3h_SIM.root");
  TMap* dcsAliasMap =(TMap*)fil->Get("DCSAliasMap");
  shuttle->SetDCSInput(dcsAliasMap);

  // DA input files
  gSystem->Exec("rm -v OCDB/ITS/DCS/DataSDD/Run*.root");
  if(optRunType=="PULSER"){
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "SDD","SDD_Calib","LDC1","SDDbase_LDC1.tar");
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "SDD","SDD_Calib","LDC2","SDDbase_LDC2.tar");
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "SDD","SDD_Calib","LDC3","SDDbase_LDC3.tar");
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "SDD","SDD_Calib","LDC4","SDDbase_LDC4.tar");
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "SDD","SDD_Calib","LDC5","SDDbase_LDC5.tar");
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "SDD","SDD_Calib","LDC6","SDDbase_LDC6.tar");
    gSystem->Exec("rm -v OCDB/ITS/Calib/CalibSDD/Run*.root");
  }else if(optRunType=="INJECTOR"){ 
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "SDD","SDD_Injec","LDC1","SDDinj_LDC1.tar");
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "SDD","SDD_Injec","LDC2","SDDinj_LDC2.tar");
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "SDD","SDD_Injec","LDC3","SDDinj_LDC3.tar");
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "SDD","SDD_Injec","LDC4","SDDinj_LDC4.tar");
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "SDD","SDD_Injec","LDC5","SDDinj_LDC5.tar");
    shuttle->AddInputFile(AliShuttleInterface::kDAQ, "SDD","SDD_Injec","LDC6","SDDinj_LDC6.tar");
    gSystem->Exec("rm -v OCDB/ITS/Calib/DriftSpeedSDD/Run*.root");
  }

  shuttle->SetInputRunType(optRunType);

  shuttle->AddInputRunParameter("totalEvents", "30000");
  shuttle->AddInputRunParameter("NumberOfGDCs", "15");
  cout<<"Input run parameters added"<<endl;


  Bool_t hltStatus=kFALSE;
  //  shuttle->SetInputHLTStatus(hltStatus);


  // Call preprocessor
  AliPreprocessor* test = new AliITSPreprocessorSDD(shuttle);
  printf("Call SDD Preprocessor\n");
  shuttle->Process();
  printf("Preprocessor OK\n");

  // Check the file which should have been created
  Char_t theDir[100];
  Bool_t doCheck=kFALSE;
  if(optRunType=="PULSER"){
    gSystem->Exec("rm SDDbase_ddl*.data");
    gSystem->Exec("rm fee.conf");
    sprintf(theDir,"ITS/Calib/CalibSDD");
    doCheck=kTRUE;
  }else if(optRunType=="INJECTOR"){ 
    gSystem->Exec("rm SDDinj_ddl*.data");
    sprintf(theDir,"ITS/Calib/DriftSpeedSDD");
    doCheck=kTRUE;
  } 

  if(doCheck){
    AliCDBEntry* chkEntry = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())->Get(theDir, 7);
    if (!chkEntry){
      printf("The Calib file is not there. Something went wrong.\n");
    }else{
      chkEntry->PrintMetaData();
      TObjArray* arr=(TObjArray*)chkEntry->GetObject();
      arr->Inspect();
    }
  }

  sprintf(theDir,"ITS/DCS/DataSDD");
  AliCDBEntry* chkEntryDCS = AliCDBManager::Instance()->GetStorage(AliShuttleInterface::GetMainCDB())
  			->Get(theDir, 7);
  if (!chkEntryDCS){
    printf("The DCS data points file is not there. Something went wrong.\n");
  }else{
    chkEntryDCS->PrintMetaData();
    TObjArray* arrdcs=(TObjArray*)chkEntryDCS->GetObject();
    arrdcs->Inspect();
  }
}


