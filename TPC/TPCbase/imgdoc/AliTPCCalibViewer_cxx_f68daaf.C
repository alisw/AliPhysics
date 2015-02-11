  {
     Float_t mean = 0;
     Float_t sigma = 1.5;
     Float_t sigmaMax = 4;
     gROOT->SetStyle("Plain");
     TH1F *distribution = new TH1F("Distrib2", "Distribution f(x, #mu, #sigma)", 1000,-5,5);
     TRandom rand(23);
     for (Int_t i = 0; i <50000;i++) distribution->Fill(rand.Gaus(mean, sigma));
     Float_t *ar = distribution->GetArray();

     TCanvas* macro_example_canvas = new TCanvas("cAliTPCCalibViewer2", "", 350, 350);
     macro_example_canvas->Divide(0,2);
     TVirtualPad *pad1 = macro_example_canvas->cd(1);
     pad1->SetGridy();
     pad1->SetGridx();
     distribution->Draw();
     TVirtualPad *pad2 = macro_example_canvas->cd(2);
     pad2->SetGridy();
     pad2->SetGridx();
     TH1F *shist = AliTPCCalibViewer::Integrate(distribution, mean, sigma, sigmaMax);
     shist->SetNameTitle("Cumulative","Cumulative S(t, #mu, #sigma)");
     shist->Draw();

  }
