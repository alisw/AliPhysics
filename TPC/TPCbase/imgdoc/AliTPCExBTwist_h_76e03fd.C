   {
   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
   TCanvas *c2 = new TCanvas("cAliTPCExBTwist","cAliTPCExBTwist",500,300);
   AliTPCExBTwist twist;
   twist.SetXTwist(0.001);  // x angle in [rad]
   twist.SetXTwist(0.0005); // y angle in [rad]
   twist.SetOmegaTauT1T2(0.32,1,1);
   twist.CreateHistoDRPhiinXY(1.)->Draw("surf2"); // A side
   return c2;
   }
