void MakeVZEROThresholdsEntryRun2()
{
  Double_t alpha[64][2] = {
 0.849114 , 1.182377 ,
 1.336399 , 0.888964 ,
 1.318963 , 0.970869 ,
 1.251016 , 1.043030 ,
 1.442306 , 1.020564 ,
 1.053460 , 1.201740 ,
 0.829986 , 1.022696 ,
 1.036325 , 1.165465 ,
 1.542461 , 0.872956 ,
 1.515863 , 0.839843 ,
 1.269698 , 1.215328 ,
 1.020312 , 0.937311 ,
 0.327139 , 1.469078 ,
 1.073537 , 1.235343 ,
 1.119917 , 1.188781 ,
 1.285190 , 1.238887 ,
 1.579386 , 1.005736 ,
 1.242708 , 0.888605 ,
 1.002295 , 0.894210 ,
 0.0,       3.0,  // 8.814347 , 7.313963 ,
 0.839430 , 1.180685 ,
 1.123069 , 0.847190 ,
 1.038351 , 1.079055 ,
 1.378613 , 1.213129 ,
 1.471212 , 1.103290 ,
 1.445550 , 1.089916 ,
 0.690587 , 1.271476 ,
 1.660245 , 1.221649 ,
 0.458279 , 1.331632 ,
 0.826517 , 1.277696 ,
 1.639883 , 1.265328 ,
 1.422140 , 1.045043 ,
 0.769852 , 1.007222 ,
 0.728513 , 1.027257 ,
 0.709045 , 0.966391 ,
 0.892321 , 1.455283 ,
 0.215711 , 1.273256 ,
 0.334897 , 1.231276 ,
 0.952247 , 1.055319 ,
 1.166976 , 1.057428 ,
 1.558643 , 1.068893 ,
 1.245368 , 1.195704 ,
 0.774615 , 1.338693 ,
 1.097145 , 1.084111 ,
 0.949590 , 1.063873 ,
 1.370406 , 0.913476 ,
 0.766243 , 0.678602 ,
 1.114483 , 1.171308 ,
 0.642971 , 1.093758 ,
 1.286589 , 1.162979 ,
 0.359324 , 1.166903 ,
 0.994052 , 1.088550 ,
 1.165003 , 1.248492 ,
 0.812264 , 1.023748 ,
 1.154415 , 1.213281 ,
 0.894247 , 1.220942 ,
 1.427910 , 1.133300 ,
 0.777129 , 1.358940 ,
 1.255835 , 1.162724 ,
 1.596315 , 0.999355 ,
 1.336408 , 1.209334 ,
 1.507150 , 1.186404 ,
 0.635214 , 1.599607 ,
 1.102436 , 1.141243
  };

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://./OCDB");

  // Creation of the functions which are used
  // in order to transform the threshold values in FEE to
  // actual thresholds in units of charge (ADC)
  TObjArray *arr = new TObjArray(64);
  arr->SetOwner(1);
  for(Int_t i = 0; i < 64; ++i) {
    TF1 *func = NULL;
    if (i != 19)
      func = new TF1(Form("thrFunc_%d",i),"[0]+[1]*x",0,10);
    else
      func = new TF1(Form("thrFunc_%d",i),"(x>=1./[1])?[0]+[1]*x:1.",0,10);

    func->SetParameter(0,alpha[i][0]);
    func->SetParameter(1,alpha[i][1]);
    arr->AddAt(func,i);
  }
  TObjString str("Discriminator thresholds calibration object");

  AliCDBMetaData *md= new AliCDBMetaData(); // metaData describing the object
  md->SetResponsible("Brigitte Cheynis");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("Discriminator thresholds calibration for Run2 and used in reconstruction and MC simulation");
  md->PrintMetaData();

  AliCDBStorage *storLoc = man->GetDefaultStorage();
  AliCDBId id("VZERO/Calib/Thresholds",243800,AliCDBRunRange::Infinity());

  storLoc->Put(arr, id, md);

  storLoc->Delete();
  delete md;

}
