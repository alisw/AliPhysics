void AliTPCCreatePRFGEM()
{ 
  //create prf
  //with given parameters
  TFile *f=TFile::Open("AliTPCprf2dGEM.root","recreate");
  AliTPCPRF2D prf;
  prf.SetPad(0.4,0.75);
  prf.SetChevron(0.75,0,0);
  prf.SetGauss(0.025,0.0251,1 );
  prf.SetY(-0.36,0.35,50);
  prf.Update();
  prf.DrawPRF(-1,1,-1,1,50,50);
  prf.Write("prf0");

  prf.SetPad(0.6,1.);
  prf.SetChevron(1.,0,0);
  prf.SetGauss(0.025,0.0251,1 );
  prf.SetY(-0.5,0.5,50);
  prf.Update();
  prf.Write("prf1");

  prf.SetPad(0.6,1.5);
  prf.SetChevron(1.5,0,0);
  prf.SetGauss(0.025,0.025,1 );
  prf.SetY(-0.75,0.75,50);
  prf.Update();
  prf.Write("prf2");


  f->Close();
}




