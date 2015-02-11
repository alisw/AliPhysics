   {
   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
   TCanvas *c2 = new TCanvas("cAliTPCCorrection","cAliTPCCorrection",700,1050);  c2->Divide(2,3);
   AliTPCROCVoltError3D roc; // EXAMPLE PLOTS - SEE BELOW
   roc.SetROCDataFileName("$ALICE_ROOT/TPC/TPCcalib/maps/TPCROCdzSurvey.root");
   roc.SetOmegaTauT1T2(0,1,1); // B=0
   Float_t z0 = 1; // at +1 cm -> A side
   c2->cd(1); roc.CreateHistoDRinXY(1.,300,300)->Draw("cont4z");
   c2->cd(3);roc.CreateHistoDRPhiinXY(1.,300,300)->Draw("cont4z");
   c2->cd(5);roc.CreateHistoDZinXY(1.,300,300)->Draw("cont4z");
   Float_t phi0=0.5;
   c2->cd(2);roc.CreateHistoDRinZR(phi0)->Draw("surf2");
   c2->cd(4);roc.CreateHistoDRPhiinZR(phi0)->Draw("surf2");
   c2->cd(6);roc.CreateHistoDZinZR(phi0)->Draw("surf2");
   return c2;
   }
