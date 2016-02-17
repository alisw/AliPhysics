TH2* h=0;
TH1F *hdcar=0,*hdcaz=0;

void test()
{
  gSystem->Load("libITSUpgradeBase.so");
  gSystem->Load("libITSUpgradeSim.so");
  gSystem->Load("libITSUpgradeRec.so");
  gROOT->ProcessLine(".L FT2.cxx+");
  //
  TParticle prt;
  //
  FT2* det = new FT2();
  det->InitEnvLocal();
  det->InitDetector();
  det->SetSimMat(kTRUE); // simulate mat.effects in probe preparation
  //
  det->SetMinTPCHits(60);
  det->SetMinITSLrHit(4);
  //
  AliESDVertex* vtx = new AliESDVertex();

  h = new TH2F("h","h",100,0,TMath::Pi()*2,100,1,-1);
  hdcar = new TH1F("hdcar","dcaXY",100,-100e-4,100e-4);
  hdcaz = new TH1F("hdcaz","dcaZ",100,-100e-4,100e-4);

  double pt = 1;//0.2;
  for (int i=0;i<50000;i++) {
    double phi = gRandom->Rndm()*TMath::Pi()*2;
    double eta = gRandom->Rndm()*0.8;
    double theta = 2*TMath::ATan(TMath::Exp(-eta));
    double pz = pt/TMath::Tan(theta);
    double px = pt*TMath::Cos(phi);
    double py = pt*TMath::Sin(phi);
    //    double pxyz[3]={0.2,0.7,0.4};
    double en = TMath::Sqrt(pt*pt+pz*pz+0.14*0.14);
    prt.SetPdgCode(211);
    prt.SetMomentum(px,py,pz,en);
    if (det->ProcessTrack(&prt,vtx)) {
      AliExternalTrackParam& prb = det->GetProbe();
      h->Fill(prb.Phi(),prb.GetSigmaY2());
      const double* dca = det->GetDCA();
      hdcar->Fill(dca[0]);
      hdcaz->Fill(dca[1]);
    }
  }

}
