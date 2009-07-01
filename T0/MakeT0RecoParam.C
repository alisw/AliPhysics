void MakeT0RecoParam(Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity(), AliRecoParam::EventSpecie_t defaultParam = AliRecoParam::kHighMult)
{
  // Create T0 Calibration Object for Ideal calibration and
  // write it on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(startRun);

  TObjArray *recoParamArray = new TObjArray();
  AliT0RecoParam* t0RecoParam;

  t0RecoParam = AliT0RecoParam::GetHighFluxParam();
  t0RecoParam->SetEventSpecie(AliRecoParam::kHighMult);
  // t0RecoParam->Dump();
  cout<<" t0RecoParam->GetEventSpecie "<< t0RecoParam->GetEventSpecie()<<endl;
  //  t0RecoParam->Dump();
  //      t0RecoParam->PrintParameters();
  recoParamArray->AddLast(t0RecoParam);



  t0RecoParam = AliT0RecoParam::GetLowFluxParam();
  t0RecoParam->SetEventSpecie(AliRecoParam::kLowMult);
  cout<<" t0RecoParam->GetEventSpecie "<< t0RecoParam->GetEventSpecie()<<endl;
  // t0RecoParam->Dump();
 //  t0RecoParam->PrintParameters();
  recoParamArray->AddLast(t0RecoParam);

  t0RecoParam = AliT0RecoParam::GetLaserTestParam();
  t0RecoParam->SetEventSpecie(AliRecoParam::kCalib);
  // t0RecoParam->Dump();
  cout<<" t0RecoParam->GetEventSpecie "<< t0RecoParam->GetEventSpecie()<<endl;
  //  t0RecoParam->Dump();
  //      t0RecoParam->PrintParameters();
  recoParamArray->AddLast(t0RecoParam);

  // Set the default
 Bool_t defaultIsSet = kFALSE;
  cout<<"recoParamArray->GetEntriesFast() "<<recoParamArray.GetEntriesFast()<<endl;
  TIter next(recoParamArray.MakeIterator());
  while ( (param = static_cast<AliT0RecoParam*>(next())) ) {
    if (!param) continue;
      if (defaultParam ==  param->GetEventSpecie()) {
      cout<<" Specie "<<param->GetEventSpecie()<<endl;
          param->SetEventSpecie(param->GetEventSpecie());
        param->SetAsDefault(); 
       defaultIsSet = kTRUE;
      }
        param->Print("FULL");
  }


  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Alla");
  md->SetComment("Reconstruction parameters T0");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("T0/Calib/RecoParam",startRun,endRun);
  man->GetDefaultStorage()->Put(recoParamArray,id, md);

}


