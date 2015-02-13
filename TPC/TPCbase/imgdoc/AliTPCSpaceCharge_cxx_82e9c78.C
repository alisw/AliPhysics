   {
   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
   TCanvas *c2 = new TCanvas("cAliTPCSpaceCharge","cAliTPCSpaceCharge",500,300);
   AliTPCSpaceCharge sc;
   sc.SetOmegaTauT1T2(-0.32,1,1); // B=0.5 Tesla
   sc.SetCorrectionFactor(0.0015);
   sc.CreateHistoDRinZR(0.)->Draw("surf2");
   return c2;
   }
