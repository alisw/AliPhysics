void
MakeRecoParam(AliRecoParam::EventSpecie_t thedefault=AliRecoParam::kLowMult)
{
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);
  
  TObjArray *recoParamArray = new TObjArray();

  {
    // Default reconstruction parameters
    AliFMDRecoParam * fmdRecoParam = new AliFMDRecoParam();
    fmdRecoParam->SetEventSpecie(AliRecoParam::kDefault);
    fmdRecoParam->SetName("Default");
    recoParamArray->AddLast(fmdRecoParam);
  }   
 {
    // Default reconstruction parameters
    AliFMDRecoParam * fmdRecoParam = new AliFMDRecoParam();
    fmdRecoParam->SetEventSpecie(AliRecoParam::kCosmic);
    fmdRecoParam->SetName("Cosmic");
    recoParamArray->AddLast(fmdRecoParam);
  }   
  {
    // Default reconstruction parameters
    AliFMDRecoParam * fmdRecoParam = new AliFMDRecoParam();
    fmdRecoParam->SetEventSpecie(AliRecoParam::kLowMult);
    fmdRecoParam->SetName("LowMult");
    recoParamArray->AddLast(fmdRecoParam);
  }   
  {
    // Default reconstruction parameters
    AliFMDRecoParam * fmdRecoParam = new AliFMDRecoParam();
    fmdRecoParam->SetEventSpecie(AliRecoParam::kHighMult);
    fmdRecoParam->SetName("HighMult");
    recoParamArray->AddLast(fmdRecoParam);
  }   
  
  // Set the default
  Bool_t defaultIsSet = kFALSE;
  for(Int_t i =0; i < recoParamArray->GetEntriesFast(); i++) {
    AliDetectorRecoParam *param = (AliDetectorRecoParam *)recoParamArray->UncheckedAt(i);
    if (!param) continue;
    if (thedefault == param->GetEventSpecie()) {
      param->SetAsDefault();
      defaultIsSet = kTRUE;
    }
  }

  if (!defaultIsSet) {
    Error(macroname,"The default reconstruction parameters are not set! Exiting...");
    return;
  }

  //AliFMDRecoParam param;
  //param.SetEventSpecie(AliRecoParam::kLowMult) ;
  AliCDBId        id("FMD/Calib/RecoParam",0,999999999);
  AliCDBMetaData  meta;
  
  meta = new AliCDBMetaData;					    
  meta.SetResponsible(gSystem->GetUserInfo()->fRealName.Data());	
  meta.SetAliRootVersion(gROOT->GetVersion()); 
  meta.SetBeamPeriod(1);						
  meta.SetComment("Data for testing"); 
  meta.SetProperty("key1", recoParamArray);
  
  cdb->Put(recoParamArray, id, &meta);
}
  
 
