
// You have to load the class before ... ;-)
// .L Detector.cxx++

void standardPlots() {

  Detector its("ALICE","ITS");

  its.MakeAliceCurrent(1); 

  its.SetRadius("bpipe",2.1);
  its.AddLayer("spd0",2.2,0.001,0.0012,0.0012);

  its.PrintLayout();
  its.SolveViaBilloir();
 
 
  TCanvas *c1 = new TCanvas("c1","c1");
 
  c1->Divide(2,2);

  Int_t color =4;
  c1->cd(1);  gPad->SetGridx();   gPad->SetGridy(); 

  gPad->SetLogx(); 
  its.GetGraphRecoEfficiency(3,color)->Draw("AC");
  its.GetGraphRecoEfficiency(2,color)->Draw("C");
  its.GetGraphRecoEfficiency(1,color)->Draw("C");
 
  c1->cd(2); gPad->SetGridx();   gPad->SetGridy(); 
  its.GetGraphMomentumResolution(color)->Draw("AC");
 
  c1->cd(3); gPad->SetGridx();   gPad->SetGridy(); 
  gPad->SetLogx(); 
  its.GetGraphPointingResolution(0,color)->Draw("AC");
  its.GetGraphPointingResolutionTeleEqu(0,1)->Draw("C");

  c1->cd(4); gPad->SetGridx();   gPad->SetGridy(); 
  gPad->SetLogx(); 
  its.GetGraphPointingResolution(1,color)->Draw("AC");
  its.GetGraphPointingResolutionTeleEqu(1,1)->Draw("C");


}


void particleDependendResolution() { 
// particle dependency on resolution

  Detector its("ALICE","ITS");

  its.MakeAliceCurrent(); 

  its.PrintLayout();
  its.SolveViaBilloir(0);
 
  its.SetRadius("bpipe",2.1);
  its.AddLayer("spd0",2.2,0.001,0.0012,0.0012);

  TCanvas *c1 = new TCanvas("c1","c1");
 
  c1->Divide(2,1);
  c1->cd(1); gPad->SetGridx();   gPad->SetGridy(); 
  gPad->SetLogx(); //gPad->SetLogy();
  c1->cd(2); gPad->SetGridx();   gPad->SetGridy(); 
  gPad->SetLogx(); //gPad->SetLogy();

  
  // compare to telescope equation ?
  //  c1->cd(1); its.GetGraphPointingResolutionTeleEqu(0,1)->Draw("AC");
  //  c1->cd(2); its.GetGraphPointingResolutionTeleEqu(1,1)->Draw("AC");

  its.SetParticleMass(0.140); // pion  
  its.SolveViaBilloir(0,0);
  c1->cd(1); its.GetGraphPointingResolution(0,1)->Draw("AC");
  c1->cd(2); its.GetGraphPointingResolution(1,1)->Draw("AC");
 
  its.SetParticleMass(0.498); // kaon  
  its.SolveViaBilloir(0,0);
  c1->cd(1); its.GetGraphPointingResolution(0,2)->Draw("C");
  c1->cd(2); its.GetGraphPointingResolution(1,2)->Draw("C");

  its.SetParticleMass(0.00051); // electron  
  its.SolveViaBilloir(0,0);
  c1->cd(1); its.GetGraphPointingResolution(0,3)->Draw("C");
  c1->cd(2); its.GetGraphPointingResolution(1,3)->Draw("C");

  its.SetParticleMass(0.938); // proton  
  its.SolveViaBilloir(0,0);
  c1->cd(1); its.GetGraphPointingResolution(0,4)->Draw("C");
  c1->cd(2); its.GetGraphPointingResolution(1,4)->Draw("C");



}
