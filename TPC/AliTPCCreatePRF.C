void AliTPCCreatePRF()
{ 
  //create prf
  //with given parameters
  TFile *f=TFile::Open("AliTPCprf2d.root","recreate");
  AliTPCPRF2D prf;
  prf.SetPad(0.4,0.75);
  prf.SetChevron(0.75,0,0);
  prf.SetGati(0.56,0.68,0.2);
  prf.SetY(-0.5,0.5,5);
  prf.Update();
  prf->Write("prf_07504_Gati_056068_d02");

  prf.SetPad(0.6,1.);
  prf.SetChevron(1.,0,0);
  prf.SetGati(0.47,0.51,0.3);
  prf.SetY(-0.875,0.875,8);
  prf.Update();
  prf.Write("prf_10006_Gati_047051_d03");

  prf.SetPad(0.6,1.5);
  prf.SetChevron(1.,0,0);
  prf.SetGati(0.47,0.51,0.3);
  prf.SetY(-1.125,1.125,10);
  prf.Update();
  prf.Write("prf_15006_Gati_047051_d03");


  f->Close();
}




