void ConfigLegoTrainPWGJE(int iFlag = 0){

  // Adds the Global Variables depending on the flag
  // 

  if(iFlag==1008){ // 10h
    AliAnalysisManager::SetGlobalStr("kJetRunPeriod","LHC10h");
  }
  else if(iFlag==1108){ // 11h
    AliAnalysisManager::SetGlobalStr("kJetRunPeriod","LHC11h");
  }

}
