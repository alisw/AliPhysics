void GetTrackingEfficiency(const char* fileName = "PWG4_JetTasksOutput.root")
{
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libJETAN");
  gSystem->Load("libPWG4JetTasks");

  TFile::Open(fileName);
  list = (TList*) gFile->Get("PWG4_LeadingTrackUE/histosLeadingTrackUE");
  
  AliUEHistograms* corr = (AliUEHistograms*) list->FindObject("AliUEHistograms");
  corr->SetEtaRange(-0.79, 0.79);

  obj = (TH1*) corr->GetNumberDensitypT()->GetTrackingEfficiency(1)->Clone("trackingefficiency");
  
  obj->Draw();
  obj->Fit("pol2", "", "", 0.5, 1.8);
  obj->Fit("pol0", "+", "SAME", 5, 15);
  
  Printf("pol2:");
  for (Int_t i=0; i<3; i++)
    Printf("par %d: %f", i, obj->GetFunction("pol2")->GetParameter(i));
  
  Printf("pol0:");
  for (Int_t i=0; i<2; i++)
    Printf("par %d: %f", i, obj->GetFunction("pol0")->GetParameter(i));
  
  // extend up to pT 100
  for (Int_t bin=obj->GetXaxis()->FindBin(10); bin <= obj->GetNbinsX(); bin++)
    obj->SetBinContent(bin, obj->GetFunction("pol0")->Eval(obj->GetXaxis()->GetBinCenter(bin)));
    
  file = TFile::Open("trackingefficiency.root", "RECREATE");
  obj->Write();
  file->Close();
}
