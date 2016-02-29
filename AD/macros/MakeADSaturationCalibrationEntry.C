// -*- C++ -*-
// $Id:$

TTree* MakeTTree();

void MakeADSaturationCalibrationEntry(const char *outputCDB = "local://$ALICE_ROOT/OCDB")
{
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("C. Mayer");
  md->SetBeamPeriod(0);
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetComment("AD Saturation");
  AliCDBId id("AD/Calib/Saturation", 0, AliCDBRunRange::Infinity());
  
  AliCDBManager *man = AliCDBManager::Instance();  
  man->SetDefaultStorage(outputCDB);

  TTree *t = MakeTTree();
  man->Put(t, id, md);
}

const Int_t gOffline2Online[16] = {
  15,13,11,9,14,12,10,8,
  7,5,3,1,6,4,2,0
};

TTree* MakeTTree() {
  TTree *t = new TTree;

  Int_t chOffline=0, chOnline=0;
  t->Branch("chOffline", &chOffline);
  t->Branch("chOnline",  &chOnline);

  Float_t extrapolationThresholds[21];
  Bool_t  doExtrapolation[21];
  for (Int_t bc=0; bc<21; ++bc) {
    extrapolationThresholds[bc] = -999.9f;
    doExtrapolation[bc]         = kFALSE;
  }
  t->Branch("extrapolationThresholds",  &extrapolationThresholds, "val[21]/F");
  t->Branch("doExtrapolation",          &doExtrapolation,         "val[21]/O");

  Float_t chargeEqualizationFactor = 1.0f;
  t->Branch("chargeEqualizationFactor", &chargeEqualizationFactor);

  for (Int_t bc=0; bc<21; ++bc) {
    extrapolationThresholds[bc] = -1.0f;
    doExtrapolation[bc] = kFALSE;
  }
  TClonesArray f_Int0("TF1", 21);
  TClonesArray f_Int1("TF1", 21);
  t->Branch("f_Int0", &f_Int0, 32000, 0);
  t->Branch("f_Int1", &f_Int0, 32000, 0);
  f_Int0.BypassStreamer();
  f_Int1.BypassStreamer();
  for (chOffline=0; chOffline<16; ++chOffline) {
    chOnline = gOffline2Online[chOffline];    
    f_Int0.Clear("C");
    f_Int1.Clear("C");
    for (Int_t bc=0; bc<21; ++bc) {
      new (f_Int0[bc]) TF1(Form("f_Ch%02d_BC%02d_int0", chOffline, bc), "x");
      new (f_Int1[bc]) TF1(Form("f_Ch%02d_BC%02d_int1", chOffline, bc), "x");
    }    
    t->Fill();
  }  
  return t;
}

