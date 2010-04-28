void AliTPCCorrectionDemo() {

  //
  // This is a Demo function of the general class AliTPCCorrection, which is used for 
  // general space point correction due to different effects.
  // The effects used in this Demo are:
  //   1. ExB twist - general offset of the TPC axis in comparison to the B field axis
  //   2. GG error (Gating Grid volt. error) - not perfectly aligned GG voltage (in terms of voltage)
  //   3. ExBBShape - B field shape correction of the secound order
  //
  // See class descriptions for further details 
  //
  // Authors: Magnus Mager, Stefan Rossegger, Jim Thomas
  //
  //
  // omegaTau (wt) of the langevin equation
  // This is a function of the drift vel., the magnetic and electric field
  // e.g. vd=2.6 cm/usc; Ez=400 V/cm; Bz=0.5 T
  // wt =  -10.0*(Bz*10)*vd/Ez = -0.325 

  Double_t vdrift = 2.6; // [cm/us]   // to be updated: per second (ideally)
  Double_t bzField = -0.5; // [Tesla] // to be updated: per run
  Double_t ezField = 400; // [V/cm]   // to be updated: never (hopefully)
  Double_t wt = -10.0 * (bzField*10) * vdrift / ezField ; 

  // Correction Terms for effective omegaTau; obtained by a laser calibration run
  Double_t T1 = 0.9;
  Double_t T2 = 1.5;

  AliMagF mag("mag","mag");

  AliTPCExBTwist twist;
  twist.SetXTwist(0.001);
  
  AliTPCGGVoltError GGerror;
  GGerror.SetDeltaVGGA(50.);
  GGerror.SetDeltaVGGC(50.);
  GGerror.InitGGVoltErrorDistortion();

  AliTPCExBBShape exb;
  exb.SetBField(&mag);

  TObjArray cs;
  cs.Add(&twist);
  cs.Add(&GGerror);
  cs.Add(&exb);

  AliTPCComposedCorrection cc;
  cc.SetCorrections(&cs);
  cc.SetOmegaTauT1T2(wt,T1,T2);
  cc.SetMode(1);

  cc.Print("DA"); // Print used correction classes

  TCanvas *c=new TCanvas;  // Plots
  c->Divide(2,2);
  c->cd(1);twist.CreateHistoDRPhiinZR(1.,100,100)->Draw("surf2");
  c->cd(2);GGerror.CreateHistoDRPhiinZR(1.,100,100)->Draw("surf2");
  c->cd(3);exb.CreateHistoDRPhiinZR(1.)->Draw("surf2");
  c->cd(4);cc.CreateHistoDRPhiinZR(1.)->Draw("surf2");
}
