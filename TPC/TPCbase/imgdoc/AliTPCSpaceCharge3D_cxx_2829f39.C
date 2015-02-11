  {
  gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
  TCanvas *c2 = new TCanvas("cAliTPCSpaceCharge3D","cAliTPCSpaceCharge3D",500,400);
  AliTPCSpaceCharge3D sc;
  sc.WriteChargeDistributionToFile("SC_zr2_GGleaks.root");
  sc.SetSCDataFileName("SC_zr2_GGleaks.root");
  sc.SetOmegaTauT1T2(0,1,1); // B=0
  sc.InitSpaceCharge3DDistortion();
  sc.CreateHistoDRinXY(15,300,300)->Draw("colz");
  return c2;
  }
