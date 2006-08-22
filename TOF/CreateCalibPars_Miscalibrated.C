void CreateCalibPars_Miscalibrated(){
  // Create TOF Calibration Object for miscalibrated detector
  // and write it on CDB
  AliTOFGeometry *geom = new AliTOFGeometryV5(); 
  AliTOFcalib *tofcalib = new AliTOFcalib(geom);
  AliTOFCal *tofCal= new AliTOFCal(geom);
  tofCal->CreateArray();

  // Input data for decalibration

  TFile f("$ALICE_ROOT/TOF/data/spectrum.root","READ");  

  TH1F *hTimeToTFit=  (TH1F*)f.Get("TimeToTFit");
  TF1  *fit=hTimeToTFit->GetFunction("poly5");
  
  // Slewing parameters (same for all channels)

  Float_t par[6] = {0.,0.,0.,0.,0.,0.};
  for(Int_t i =0;i<6;i++){
    par[i]=fit->GetParameter(i);
    //    cout << " Slewing parameters=" << par[i] << endl;
  }

  // Global time offset (randomly gen, gaussian with mean = 0.3, sig=0.08 ns)

  Float_t delay=0.;
  Float_t meanDelay=0.3;
  Float_t sigmaDelay=0.08;

  // ToT spectrum 

  TH1F *hToT=  (TH1F*)f.Get("ToT");

  // Fill the Sim calibration object

  TRandom *rnd   = new TRandom(4357);
  for (Int_t ipad = 0 ; ipad<tofCal->NPads(); ipad++){
    AliTOFChannel *calChannel = tofCal->GetChannel(ipad);
    calChannel->SetSlewPar(par);
    delay=rnd->Gaus(meanDelay,sigmaDelay);
    // cout << " delay=" << delay << endl;
    calChannel->SetDelay(delay);
  }
  tofcalib->WriteSimParOnCDB("TOF/Calib",0,0,tofCal,hToT);
  f.Close();
}


