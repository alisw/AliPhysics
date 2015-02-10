   {
   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
   TCanvas *c2 = new TCanvas("cAliTPCVoltError3D","cAliTPCVoltError3D",500,450);
   AliTPCFCVoltError3D fc;
   fc.SetOmegaTauT1T2(0,1,1);
   fc.SetRotatedClipVoltA(0,40);
   fc.SetRodVoltShiftA(3,40);
   fc.SetCopperRodShiftA(7+18,40);
   fc.SetRodVoltShiftA(15+18,40);
   fc.CreateHistoDRPhiinXY(10)->Draw("cont4z");
   return c2;
   }
