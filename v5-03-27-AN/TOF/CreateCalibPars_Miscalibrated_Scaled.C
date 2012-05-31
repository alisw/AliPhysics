void CreateCalibPars_Miscalibrated_Scaled(){
  // Create TOF Calibration Object for miscalibrated detector
  // and write it on CDB
  AliTOFcalib *tofcalib = new AliTOFcalib();
  tofcalib->CreateCalArrays();
  TObjArray *tofCalOffline = (TObjArray*) tofcalib->GetTOFCalArrayOffline(); 
  // Input data for decalibration

  TFile f("$ALICE_ROOT/TOF/data/spectrumScaled.root","READ");  

  TH1F *hTimeToTFit=  (TH1F*)f.Get("hTimeToTScaled");
  TF1  *fit=hTimeToTFit->GetFunction("pol5");
  
  // Slewing parameters (same for all channels)

  Float_t par[6] = {0.,0.,0.,0.,0.,0.};
  Float_t kpar=0;
  for(Int_t i =0;i<6;i++){
    par[i]=fit->GetParameter(i);
    kpar=par[0];
    cout << " Slewing parameter " <<i<<" =" << par[i] << endl;
  }

  // Global time offset (randomly gen, gaussian with mean = 0.3, sig=0.08 ns)

  Float_t delay=0.;
  Float_t meanDelay=0.3;
  Float_t sigmaDelay=0.08;

  // ToT spectrum 

  TH1F *hToT=  (TH1F*)f.Get("hToTScaled");

  // Fill the Sim calibration object

  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  TRandom *rnd   = new TRandom(4357);
  Int_t nChannels = AliTOFGeometry::NSectors()*(2*(AliTOFGeometry::NStripC()+AliTOFGeometry::NStripB())+AliTOFGeometry::NStripA())*AliTOFGeometry::NpadZ()*AliTOFGeometry::NpadX();
  for (Int_t ipad = 0 ; ipad<nChannels; ipad++){
    AliTOFChannelOffline *calChannelOffline = (AliTOFChannelOffline*)tofCalOffline->At(ipad);
    delay=rnd->Gaus(meanDelay,sigmaDelay);
    par[0]=kpar+delay;
    calChannelOffline->SetSlewPar(par);
  }
  tofcalib->WriteParOfflineOnCDB("TOF/Calib","valid",0,AliCDBRunRange::Infinity());
  tofcalib->WriteSimHistoOnCDB("TOF/Calib",0,AliCDBRunRange::Infinity(),hToT);
  f.Close();
}


