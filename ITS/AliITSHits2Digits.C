Int_t AliITSHits2Digits()
{

  // Connect the Root Galice file containing Geometry, Kine and Hits

  const char * inFile = "galice.root";  
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(inFile);
  if (file) {file->Close(); delete file;}
  printf("Hits2Digits\n");
  file = new TFile(inFile,"UPDATE");
  if (!file->IsOpen()) {
    cerr<<"Can't open "<<inFile<<" !\n";
    return 1;
  }
  file->ls();

  // Get AliRun object from file or return if not on file
  if (gAlice) delete gAlice;
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice) {
    cerr<<"ITSHits2Digits.C : AliRun object not found on file\n";
    return 2;
  }

  gAlice->GetEvent(0);
  AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");      
  if (!ITS) {
    cerr<<"ITSHits2Digits.C : AliITS object not found on file\n";
    return 3;
  }

// Set the simulation models for the three detector types
  AliITSgeom *geom = ITS->GetITSgeom();

  // SPD
  AliITSDetType *iDetType=ITS->DetType(0);
  AliITSsegmentationSPD *seg0=(AliITSsegmentationSPD*)iDetType->GetSegmentationModel();
  AliITSresponseSPD *res0 = (AliITSresponseSPD*)iDetType->GetResponseModel();
  AliITSsimulationSPD *sim0=new AliITSsimulationSPD(seg0,res0);
  ITS->SetSimulationModel(0,sim0);
  // test
  printf("SPD dimensions %f %f \n",seg0->Dx(),seg0->Dz());
  printf("SPD npixels %d %d \n",seg0->Npz(),seg0->Npx());
  printf("SPD pitches %d %d \n",seg0->Dpz(0),seg0->Dpx(0));
  // end test

  // SDD
  //Set response functions
  Float_t baseline = 10.;
  Float_t noise = 1.75;

  // SDD compression param: 2 fDecrease, 2fTmin, 2fTmax or disable, 2 fTolerance
  //Float_t fCutAmp = baseline + 2.*noise;

  Float_t maxadc = res1->MaxAdc();    
  Float_t topValue = res1->MagicValue();
  Float_t norm = maxadc/topValue;

  Float_t fCutAmp = baseline + 2.*noise;
  fCutAmp *= norm;  
  
  Int_t cp[8]={0,0,fCutAmp,fCutAmp,0,0,0,0};

  AliITSDetType *iDetType=ITS->DetType(1);
  AliITSresponseSDD *res1 = (AliITSresponseSDD*)iDetType->GetResponseModel();
  if (!res1) {
    res1=new AliITSresponseSDD();
    ITS->SetResponseModel(1,res1);
  }
  //res1->SetZeroSupp("2D");
  res1->SetZeroSupp("1D");
  res1->SetNoiseParam(noise,baseline);
  res1->SetDo10to8(kTRUE);
  res1->SetCompressParam(cp);
  res1->SetMinVal(4);
  res1->SetDiffCoeff(3.6,40.);
  res1->SetMagicValue(96.95);

  AliITSsegmentationSDD *seg1=(AliITSsegmentationSDD*)iDetType->GetSegmentationModel();
  if (!seg1) {
    seg1 = new AliITSsegmentationSDD(geom,res1);
    ITS->SetSegmentationModel(1,seg1);
  }
  AliITSsimulationSDD *sim1=new AliITSsimulationSDD(seg1,res1);
  sim1->SetDoFFT(1);
  sim1->SetCheckNoise(kFALSE);
  ITS->SetSimulationModel(1,sim1);

  // SSD
  AliITSDetType *iDetType=ITS->DetType(2);
  AliITSsegmentationSSD *seg2=(AliITSsegmentationSSD*)iDetType->GetSegmentationModel();
  AliITSresponseSSD *res2 = (AliITSresponseSSD*)iDetType->GetResponseModel();
  res2->SetSigmaSpread(3.,2.);
  AliITSsimulationSSD *sim2=new AliITSsimulationSSD(seg2,res2);
  ITS->SetSimulationModel(2,sim2);


  cerr<<"Digitizing ITS...\n";
  
  TStopwatch timer;
  timer.Start();
  ITS->HitsToDigits(0,0,-1," ","All"," ");
  timer.Stop(); timer.Print();

  delete sim0;
  delete sim1;
  delete sim2;


  delete gAlice;   gAlice=0;
  file->Close(); 
  delete file;
  return 0;
};

