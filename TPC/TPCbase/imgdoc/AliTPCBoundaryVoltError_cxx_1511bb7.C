/// \file AliTPCBoundaryVoltError_cxx_1511bb7.C

   {
   gROOT->SetStyle("Plain"); gStyle->SetPalette(1);
   TCanvas *c2 = new TCanvas("cAliTPCBoundaryVoltError","cAliTPCBoundaryVoltError",500,300);
   AliTPCBoundaryVoltError bve;
   Float_t val = 40;// [Volt]; 40V corresponds to 1mm
   /* IFC shift, CE follows, ROC follows by factor half */
   Float_t boundA[8] = { val, val, val,0,0,0,0,val}; // voltages A-side
   Float_t boundC[6] = {-val,-val,-val,0,0,0};       // voltages C-side
   bve.SetBoundariesA(boundA);
   bve.SetBoundariesC(boundC);
   bve.SetOmegaTauT1T2(-0.32,1,1);
   bve.SetROCDisplacement(kTRUE); // include the chamber offset in z when calculating the dz distortions
   bve.CreateHistoDRinZR(0)->Draw("surf2");
   return c2;
   }
