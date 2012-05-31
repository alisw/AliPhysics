void MakeT0RecoParam(Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity(), AliRecoParam::EventSpecie_t defaultParam = AliRecoParam::kHighMult)
//void MakeT0RecoParam(Int_t startRun = 122000, Int_t endRun = 126437, AliRecoParam::EventSpecie_t defaultParam = AliRecoParam::kHighMult)
{
  // Create T0 Calibration Object for Ideal calibration and
  // write it on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  //  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetDefaultStorage("local:///scratch/alla/alice/Oct11/TestCDBf");
  man->SetRun(startRun);
  
  TObjArray *recoParamArray = new TObjArray();
  AliT0RecoParam* t0RecoParam;

  t0RecoParam = AliT0RecoParam::GetHighFluxParam();
  t0RecoParam->SetEventSpecie(AliRecoParam::kHighMult);
   Float_t low[500] , high[500] ;
   Float_t pos[24] = {2472, 2492, 2501, 2540, 2526, 2525, 2501, 2461, 2556, 2544, 2533, 2487, 
		      2485, 2502, 2453, 2470, 2520, 2486, 2512, 2523, 2505, 2562, 2533, 2481}; 
    for (Int_t i=0; i<107; i++)
      {
	low[i] = 500;
	high[i]=50000;
 	if(i>0  && i<13) { low[i] = pos[i-1] - 200; high[i] = pos[i-1] + 200;}
	if(i>56 && i<69) { low[i] = pos[i-57] - 200.; high[i]=pos[i-57] + 200;}
	//	if(i<13) { low[i]=pos[i]-500.;high[i]=pos[i]+500.;}
	//	if(i>56 && i<69) { low[i]=pos[i+13-57]-500.;high[i]=pos[i-57+13]+500.;}
  }

  // bad PMT
  t0RecoParam->SetRefPoint(-1);
  // amplitude range for recontruction
  low[200]=0.8;
  high[200]=80;
  for (Int_t i=0; i<24; i++) t0RecoParam-> SetBadChannels(i,0);
  
  for (Int_t i=0; i<500; i++) {
    t0RecoParam->SetLow(i,low[i]);
    t0RecoParam->SetHigh(i,high[i]);
  }
  t0RecoParam->SetEq(1);
  t0RecoParam->SetSatelliteThresholds(-15, -1.7);

  t0RecoParam->PrintParameters();
  recoParamArray->AddLast(t0RecoParam);
  
  //-----------------------------------------

  t0RecoParam = AliT0RecoParam::GetLowFluxParam();
  t0RecoParam->SetEventSpecie(AliRecoParam::kLowMult);
    Float_t low[500], high[500] ;
    for (Int_t i=0; i<500; i++)
      {
	low[i] = 500;
	high[i]=50000;
  	if(i>0  && i<13) { low[i] = pos[i-1] - 200; high[i] = pos[i-1] + 200;}
	if(i>56 && i<69) { low[i] = pos[i-57] - 200.; high[i]=pos[i-57] + 200;}
	//	if(i<13) { low[i]=pos[i]-500.;high[i]=pos[i]+500.;}
	//	if(i>56 && i<69) { low[i]=pos[i+13-57]-500.;high[i]=pos[i-57+13]+500.;}
    }

  // bad PMT
    //  t0RecoParam->SetRefPoint(-1);
    for (Int_t i=0; i<24; i++) t0RecoParam->SetBadChannels(i,0);
  t0RecoParam->SetRefPoint(-1);
  // amplitude range for recontruction
   t0RecoParam->SetEq(1);
 low[200]=0.8;
  high[200]=80;
  for (Int_t i=0; i<500; i++) {
    t0RecoParam->SetLow(i,low[i]);
    t0RecoParam->SetHigh(i,high[i]);
  }
  //  t0RecoParam->PrintParameters();
  recoParamArray->AddLast(t0RecoParam);
  //--------------------------------------
  t0RecoParam = AliT0RecoParam::GetLaserTestParam();
  t0RecoParam->SetEventSpecie(AliRecoParam::kCalib);
  // t0RecoParam->Dump();
  cout<<" t0RecoParam->GetEventSpecie "<< t0RecoParam->GetEventSpecie()<<endl;
  for (Int_t i=0; i<107; i++) {
    t0RecoParam->SetLow(i, 0);
    t0RecoParam->SetHigh(i, 50000);
  }
  t0RecoParam->SetEq(1);
  t0RecoParam->SetSatelliteThresholds(-15, -1.7);
  //  t0RecoParam->Dump();
  t0RecoParam->SetRefPoint(1);
  //  t0RecoParam->PrintParameters();
  recoParamArray->AddLast(t0RecoParam);


  // Set the default
 Bool_t defaultIsSet = kFALSE;
 //  cout<<"recoParamArray->GetEntriesFast() "<<recoParamArray.GetEntriesFast()<<endl;
  TIter next(recoParamArray->MakeIterator());
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


