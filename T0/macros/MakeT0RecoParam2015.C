void MakeT0RecoParam2015(Int_t startRun = 0, Int_t endRun = AliCDBRunRange::Infinity(), AliRecoParam::EventSpecie_t defaultParam = AliRecoParam::kLowMult)
//void MakeT0RecoParam(Int_t startRun = 122000, Int_t endRun = 126437, AliRecoParam::EventSpecie_t defaultParam = AliRecoParam::kHighMult)
{
  // Create T0 Calibration Object for Ideal calibration and
  // write it on CDB
  AliCDBManager *man = AliCDBManager::Instance();
  //   man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
   man->SetDefaultStorage("local:///home/alla/alice/Nov15/OCDB");
  man->SetRun(startRun);
  
  TObjArray *recoParamArray = new TObjArray();
  AliT0RecoParam* t0RecoParam;
  //coefficients for new QTC 00->10
  Float_t a[26] ={0.842277,  0.9588,   0.957138, 0.985494,0.983393, 0.960508,   
		  0.929485,  0.941351,  0.987975,  0.994604, 0.95595, 0.964502,   
		  0.964731,    0.97356,  0.996221,  0.981255, 0.943558, 1.00606,   
		  0.934549,   0.969384, 0.97097, 0.989021, 0.95537, 0.987018,
		  1,    1};   
  Float_t b[26] = { 35.431,  14.6214,  3.83161, 48.624, 18.9751, 31.1847 ,  
		    40.1204, 24.1932 ,20.8494,  17.2447, 23.0073, 20.1561,   
		    22.3656, 27.0031, 19.4473,  18.5315, 20.9247, 15.5292,   
		    42.4446, -3.91989, 11.0756, 16.9568,  25.8498, 21.2587, 0,0};   

  
  Float_t low[500], high[500] ;
  for (Int_t i=0; i<500; i++)
    {
      low[i] = 0;
      high[i]=50000;
    }
  
  t0RecoParam = AliT0RecoParam::GetHighFluxParam();
  t0RecoParam->SetEventSpecie(AliRecoParam::kHighMult);
  // bad PMT
  t0RecoParam->SetRefPoint(-1);
   // amplitude range for recontruction
  //C side
  low[200]=0.55;
  //A side
  low[201] = 0.55;
  //  amplitude cutted 10% spectra  
   // for runs 2015
  for (Int_t i=0; i<26; i++) {
    low[i+130]=a[i];
    low[i+156]=b[i];
  }
  for (Int_t i=0; i<24; i++) t0RecoParam-> SetBadChannels(i,0);
  //corridor : mean position +-corridor
  low[300]=70; 
  low[310]=0;
  low[311]=16;
  low[312]=20; //according LHC15k1
  
  t0RecoParam->SetEq(1);
  t0RecoParam->SetSatelliteThresholds(-15, -1.7);
  
   for (Int_t i=0; i<500; i++) {
    t0RecoParam->SetLow(i,low[i]);
    t0RecoParam->SetHigh(i,high[i]);
  }
  recoParamArray->AddLast(t0RecoParam);
 
  //-----------------------------------------
  
  t0RecoParam = AliT0RecoParam::GetLowFluxParam();
  t0RecoParam->SetEventSpecie(AliRecoParam::kLowMult);
 // bad PMT
  //  t0RecoParam->SetRefPoint(-1);
  for (Int_t i=0; i<24; i++) t0RecoParam->SetBadChannels(i,0);
  t0RecoParam->SetRefPoint(-1);
  // for runs 2015
   t0RecoParam->SetEq(1); //after  LHC11f
  //C side
  low[200]=0.55;
  //A side
  low[201]=0.55;
  //corridor : mean position +-corridor
  low[300]=70; // for 2015
  //shift simulated  T0A, T0C, (T0A+T0C)/2
  low[310]=40;
  low[311]=30;
  low[312]=45;
   //new QTC start position
  for (Int_t i=0; i<26; i++) {
    low[i+130]=a[i];
    low[i+156]=b[i];
  }
  t0RecoParam->SetSatelliteThresholds(-15, -1.5);
   for (Int_t i=0; i<500; i++) {
    t0RecoParam->SetLow(i,low[i]);
    t0RecoParam->SetHigh(i,high[i]);
  }
  recoParamArray->AddLast(t0RecoParam);
  //--------------------------------------
  t0RecoParam = AliT0RecoParam::GetLaserTestParam();
  t0RecoParam->SetEventSpecie(AliRecoParam::kCalib);
  // t0RecoParam->Dump();
  cout<<" t0RecoParam->GetEventSpecie "<< t0RecoParam->GetEventSpecie()<<endl;
  Float_t low[500], high[500] ;
  for (Int_t i=0; i<500; i++)
    {
      low[i] = 500;
      high[i]=50000;
    }
  t0RecoParam->SetEq(1);
  t0RecoParam->SetSatelliteThresholds(-15, -1.7);
  //  t0RecoParam->Dump();
  t0RecoParam->SetRefPoint(1);
  //corridor : mean position +-corridor
  low[300]=50;
  low[200]=0.8; //// for runs 2009-2010
  high[200]=60;
  //A side
  //low[201]=0.55;
  low[201]=0.8;// for runs 2009-2010
  high[201]=60;
  for (Int_t i=0; i<500; i++) {
    t0RecoParam->SetLow(i,low[i]);
    t0RecoParam->SetHigh(i,high[i]);
  }
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


