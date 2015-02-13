   {
   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
   TCanvas *c2 = new TCanvas("cAliTPCExBBShape","cAliTPCExBBShape",500,300);
   AliTPCExBBShape exb;
   AliMagF mag("mag","mag");        // 0.5 Tesla (solenoid)
   exb.SetBField(&mag);             // use Bfield from AliMagF
   exb.SetOmegaTauT1T2(-0.32,1.,1.); // values ideally from OCDB
   exb.CreateHistoDRPhiinZR(0,100,100)->Draw("surf2");
   return c2;
   }
