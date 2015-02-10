   {
   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
   TCanvas *c2 = new TCanvas("cAliTPCGGVoltError","cAliTPCGGVoltError",500,300);
   AliTPCGGVoltError gg;
   gg.SetDeltaVGGA(-40); gg.SetDeltaVGGC(-40); // 40 Volt offset
   gg.InitGGVoltErrorDistortion();
   gg.SetOmegaTauT1T2(0,1,1); // B=0
   gg.CreateHistoDRinZR(0)->Draw("surf2");
   return c2;
   }
