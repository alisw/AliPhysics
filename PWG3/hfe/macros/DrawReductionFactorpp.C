TList *GetResults(const Char_t *testfile,const Char_t *plus="");
TObject* GetSpectrum(AliCFContainer *c, Int_t step);
TObject* GetEfficiency(AliCFContainer *c, Int_t step = 9);
TObject* GetEfficiency(AliCFContainer *c, Int_t step, Int_t step0);


void DrawReductionFactorpp(const char *testfile,const char *name) {
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);
  
  ///////////////////////////////////
  // Take the stuff
  //////////////////////////////////

  TList *results = GetResults(testfile,name);
  if(!results){
    printf("No output objects: Calculation will terminate here\n");
    return;
  }

  AliHFEcontainer *containerhfe = (AliHFEcontainer *) results->FindObject("trackContainer");
  if(!containerhfe) {
    printf("No hfe container \n");
    return;
  }

  AliCFContainer *sumcontainer = containerhfe->GetCFContainer("recTrackContReco");
  if(!sumcontainer) {
    printf("No container sum esd\n");
    return;
  }
  
  Int_t numberOfEvents = (Int_t) containerhfe->GetNumberOfEvents();
  
  printf("Number of events for  %d after Event cut\n",numberOfEvents);
  
  /////////////////////////////////////
  // Take efficiencies
  /////////////////////////////////////
  
  // ESD values

  AliCFEffGrid  *efficiencystepkineITSTPC  = (AliCFEffGrid*)  GetEfficiency(sumcontainer,AliHFEcuts::kStepRecKineITSTPC,AliHFEcuts::kStepRecNoCut);
  AliCFEffGrid  *efficiencystepPrim        = (AliCFEffGrid*)  GetEfficiency(sumcontainer,AliHFEcuts::kStepRecPrim,AliHFEcuts::kStepRecKineITSTPC);
  AliCFEffGrid  *efficiencystepHFEcutsITS  = (AliCFEffGrid*)  GetEfficiency(sumcontainer,AliHFEcuts::kStepHFEcutsITS,AliHFEcuts::kStepRecPrim);
  AliCFEffGrid  *efficiencystepHFEcutsTRD  = (AliCFEffGrid*)  GetEfficiency(sumcontainer,AliHFEcuts::kStepHFEcutsTRD,AliHFEcuts::kStepHFEcutsITS);
  

  AliCFEffGrid  *efficiencystepPIDTOF      = (AliCFEffGrid*)  GetEfficiency(sumcontainer,AliHFEcuts::kStepHFEcutsTRD + 1,AliHFEcuts::kStepHFEcutsTRD);
  AliCFEffGrid  *efficiencystepPIDTPC      = (AliCFEffGrid*)  GetEfficiency(sumcontainer,AliHFEcuts::kStepHFEcutsTRD + 2,AliHFEcuts::kStepHFEcutsTRD + 1);

  AliCFEffGrid  *efficiencystepPIDall      = (AliCFEffGrid*)  GetEfficiency(sumcontainer,AliHFEcuts::kStepHFEcutsTRD + 2,AliHFEcuts::kStepHFEcutsTRD);

  AliCFEffGrid  *efficiencystepall      = (AliCFEffGrid*)  GetEfficiency(sumcontainer,AliHFEcuts::kStepHFEcutsTRD + 2,AliHFEcuts::kStepRecNoCut);
  

  ///////////////////////////////////////
  // Plot 1D
  ///////////////////////////////////////

  TCanvas * canvascomparisonbis = new TCanvas("ITSTPCrefitStep","ITSTPCrefitStep",1000,700);
  canvascomparisonbis->Divide(2,2);
 
  canvascomparisonbis->cd(1);
  TH1D* h_effstepkineITSTPC1Donly = (TH1D *) efficiencystepkineITSTPC->Project(0);
  h_effstepkineITSTPC1Donly->SetTitle("(ITS & TPC refit) / (No cut)");
  h_effstepkineITSTPC1Donly->SetStats(0);
  h_effstepkineITSTPC1Donly->SetLineColor(kBlue);
  h_effstepkineITSTPC1Donly->SetMarkerColor(kBlue);
  h_effstepkineITSTPC1Donly->SetMarkerStyle(25);
  h_effstepkineITSTPC1Donly->GetXaxis()->SetRangeUser(0.3,6.0);
  h_effstepkineITSTPC1Donly->SetXTitle("p_{T} [GeV/c]");
  h_effstepkineITSTPC1Donly->SetYTitle("Efficiency");
  h_effstepkineITSTPC1Donly->Draw();

  canvascomparisonbis->cd(2);
  TH1D* h_effstepPrim1Donly = (TH1D *) efficiencystepPrim->Project(0);
  h_effstepPrim1Donly->SetTitle("No Kink");
  h_effstepPrim1Donly->SetStats(0);
  h_effstepPrim1Donly->SetLineColor(kBlue);
  h_effstepPrim1Donly->SetMarkerColor(kBlue);
  h_effstepPrim1Donly->SetMarkerStyle(25);
  h_effstepPrim1Donly->GetXaxis()->SetRangeUser(0.3,6.0);
  h_effstepPrim1Donly->SetXTitle("p_{T} [GeV/c]");
  h_effstepPrim1Donly->SetYTitle("Efficiency");
  h_effstepPrim1Donly->Draw();

  canvascomparisonbis->cd(3);
  TH1D* h_effstepHFEcutsITS1Donlya = (TH1D *) efficiencystepHFEcutsITS->Project(0);
  h_effstepHFEcutsITS1Donlya->SetTitle("First pixel");
  h_effstepHFEcutsITS1Donlya->SetStats(0);
  h_effstepHFEcutsITS1Donlya->SetLineColor(kBlue);
  h_effstepHFEcutsITS1Donlya->SetMarkerColor(kBlue);
  h_effstepHFEcutsITS1Donlya->SetMarkerStyle(25);
  h_effstepHFEcutsITS1Donlya->GetXaxis()->SetRangeUser(0.3,6.0);
  h_effstepHFEcutsITS1Donlya->GetYaxis()->SetRangeUser(0.0,1.0);
  h_effstepHFEcutsITS1Donlya->SetXTitle("p_{T} [GeV/c]");
  h_effstepHFEcutsITS1Donlya->SetYTitle("Efficiency");
  h_effstepHFEcutsITS1Donlya->Draw();
 
  canvascomparisonbis->cd(4);
  TH1D* h_effstepHFEcutsTRDonly = (TH1D *) efficiencystepHFEcutsTRD->Project(0);
  h_effstepHFEcutsTRDonly->SetTitle("TOF matching");
  h_effstepHFEcutsTRDonly->SetStats(0);
  h_effstepHFEcutsTRDonly->SetLineColor(kBlue);
  h_effstepHFEcutsTRDonly->SetMarkerColor(kBlue);
  h_effstepHFEcutsTRDonly->SetMarkerStyle(25);
  h_effstepHFEcutsTRDonly->GetXaxis()->SetRangeUser(0.3,6.0);
  h_effstepHFEcutsTRDonly->SetXTitle("p_{T} [GeV/c]");
  h_effstepHFEcutsTRDonly->SetYTitle("Efficiency");
  h_effstepHFEcutsTRDonly->Draw();
  
  /////////
  // PID
  /////////

  TCanvas * canvasc1Don = new TCanvas("PIDstep","PIDstep",1000,700);
  canvasc1Don->Divide(2,1);
  
  canvasc1Don->cd(1);
  TH1D* h_effstepPID1Donly = (TH1D *) efficiencystepPIDTOF->Project(0);
  h_effstepPID1Donly->SetTitle("TOF 3 sigma cut");
  h_effstepPID1Donly->SetStats(0);
  h_effstepPID1Donly->SetLineColor(kBlue);
  h_effstepPID1Donly->SetMarkerColor(kBlue);
  h_effstepPID1Donly->SetMarkerStyle(25);
  h_effstepPID1Donly->GetXaxis()->SetRangeUser(0.3,6.0);
  h_effstepPID1Donly->SetXTitle("p_{T} [GeV/c]");
  h_effstepPID1Donly->SetYTitle("Efficiency");
  h_effstepPID1Donly->Draw();
 
  canvasc1Don->cd(2);
  gPad->SetLogy();
  TH1D* h_effstepPID1Donlyy = (TH1D *) efficiencystepPIDTPC->Project(0);
  h_effstepPID1Donlyy->SetTitle("TPC PID cut");
  h_effstepPID1Donlyy->SetStats(0);
  h_effstepPID1Donlyy->SetLineColor(kBlue);
  h_effstepPID1Donlyy->SetMarkerColor(kBlue);
  h_effstepPID1Donlyy->SetMarkerStyle(25);
  h_effstepPID1Donlyy->GetXaxis()->SetRangeUser(0.3,6.0);
  h_effstepPID1Donlyy->SetXTitle("p_{T} [GeV/c]");
  h_effstepPID1Donlyy->SetYTitle("Efficiency");
  h_effstepPID1Donlyy->Draw();
  

  
  
}

// ====================================================================
TList *GetResults(const Char_t *testfile,const Char_t *plus){
  //
  // read output
  //
  TFile *f = TFile::Open(testfile);
  if(!f || f->IsZombie()){
    printf("File not readable\n");
    return 0x0;
  }
  TString name("HFE_Results");
  name += plus; 
  printf("Name of TList %s\n",(const char*)name); 
  TList *l = dynamic_cast<TList *>(f->Get((const char*)name));
  if(!l){
    printf("Output container not found\n");
    f->Close(); delete f;
    return 0x0;
  } 
  TList *returnlist = dynamic_cast<TList *>(l->Clone());
  f->Close(); delete f;
  return returnlist;
}
//_________________________________________________________________________
TObject* GetSpectrum(AliCFContainer *c, Int_t step) {
  AliCFDataGrid* data = new AliCFDataGrid("data","",*c,step);
  //data->SetMeasured(step);
  return data;
}
//__________________________________________________________________________
TObject* GetEfficiency(AliCFContainer *c, Int_t step) {
  TString name("eff");
  name += step;
  AliCFEffGrid* eff = new AliCFEffGrid((const char*)name,"",*c);
  eff->CalculateEfficiency(step,0);
  return eff;
}
//_________________________________________________________________________
TObject* GetEfficiency(AliCFContainer *c, Int_t step, Int_t step0) {
  TString name("eff");
  name += step;
  name+= step0;
  AliCFEffGrid* eff = new AliCFEffGrid((const char*)name,"",*c);
  eff->CalculateEfficiency(step,step0);
  return eff;
}
