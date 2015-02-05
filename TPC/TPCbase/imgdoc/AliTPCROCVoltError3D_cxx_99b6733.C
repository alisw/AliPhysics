   {
   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
   TCanvas *c2 = new TCanvas("cAliTPCROCVoltError3D","cAliTPCROCVoltError3D",500,400);
   AliTPCROCVoltError3D roc;
   roc.SetROCDataFileName("$ALICE_ROOT/TPC/TPCcalib/maps/TPCROCdzSurvey.root");
   roc.SetElectronArrivalCorrection(kFALSE);  // Correction for electron arrival offset, IROC vs OROC
   roc.SetROCDisplacement(kTRUE);   // include the chamber offset in z when calculating the dz
   roc.SetOmegaTauT1T2(0,1,1); // B=0
   roc.CreateHistoDZinXY(1.,300,300)->Draw("colz");
   return c2;
   }
