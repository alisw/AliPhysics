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
  Float_t a[26] = { 0.862486, 0.982631, 0.969829, 0.987521, 0.992425, 
		    0.991012, 0.955364, 0.973607, 0.994484, 1.000748,
		    0.965197, 0.978035, 0.981249, 0.993485, 1.013696, 
		    0.980984, 0.952803, 0.989909, 0.946670, 0.972652, 
		    0.984168, 0.986350, 0.952677, 0.986929,
		    1.1,    0.9};
  Float_t  b[26] = {33.057700, 11.560650, 3.381142, 48.665796, 
		    18.108813, 27.025644, 37.289124, 20.609302, 
		    19.673466, 16.640889, 20.773348, 18.378632, 
		    19.674898, 24.323490, 17.358610, 18.515412,
		    20.156329, 15.466606, 40.983598, -5.113665, 
		    10.714077, 17.500848, 26.717135, 21.718537,
		    29.8, 33.3};

  
  
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
  low[200]=200;
  //A side
  low[201] = 200;
  //  amplitude cutted 10% spectra  
   // for runs 2015
  for (Int_t i=0; i<26; i++) {
    low[i+130]=a[i];
    low[i+156]=b[i];
  }
  for (Int_t i=0; i<24; i++) t0RecoParam-> SetBadChannels(i,0);
  //corridor : mean position +-corridor
  low[300]=70; 
  low[310]=80;
  low[311]=56;
  low[312]=70; //according LHC15k1
  
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
  low[200]=100;
  //A side
  low[201]=100;
  //corridor : mean position +-corridor
  low[300]=70; // for 2015
  //shift simulated  T0A, T0C, (T0A+T0C)/2 for 13TeV
  // low[310]=40;
  // low[311]=30;
  // low[312]=45;
  //shift simulated  T0A, T0C, (T0A+T0C)/2 for 5 TeV
  //  low[310]=55;
  // low[311]=60;
  // low[312]=70;
  //shift simulated  T0A, T0C, (T0A+T0C)/2 for 13 TeV
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


