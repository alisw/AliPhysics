void CreateOfflineCalibPars(){
  // Create TOF Offline Calibration Object for reconstruction
  // and write it on CDB
  AliTOFcalib *tofcalib = new AliTOFcalib();
  tofcalib->CreateCalArrays();
  TObjArray *tofCalOffline = (TObjArray*) tofcalib->GetTOFCalArrayOffline(); 

  TFile f("$ALICE_ROOT/TOF/data/spectrumScaled.root","READ");  

  TH1F *hTimeToTFit=  (TH1F*)f.Get("hTimeToTScaled");
  TF1  *fit=hTimeToTFit->GetFunction("pol5");
  
  // Slewing parameters (same for all channels)

  Float_t delay=0.;
  Float_t meanDelay=0.3;
  Float_t sigmaDelay=0.08;
  TRandom *rnd   = new TRandom(4357);

  Float_t par[6] = {0.,0.,0.,0.,0.,0.};
  for(Int_t i =0;i<6;i++){
    par[i]=fit->GetParameter(i);
    cout << " Slewing parameter " <<i<<" =" << par[i] << endl;
  }


  Int_t nChannels = AliTOFGeometry::NSectors()*(2*(AliTOFGeometry::NStripC()+AliTOFGeometry::NStripB())+AliTOFGeometry::NStripA())*AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX();

  for (Int_t ipad = 0 ; ipad<nChannels; ipad++){
    AliTOFChannelOffline *calChannelOffline = (AliTOFChannelOffline*)tofCalOffline->At(ipad);
    delay = rnd->Gaus(meanDelay,sigmaDelay);
    par[0]+=delay.; // adding time delay
    calChannelOffline->SetSlewPar(par);
  }
  // Write the valid offline calibration object on CDB

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE");
  tofcalib->WriteParOfflineOnCDB("TOF/Calib","valid",0,0);
  f.Close();
  return;

}
