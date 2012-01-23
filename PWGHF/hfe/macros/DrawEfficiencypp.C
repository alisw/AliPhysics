TList *GetResults(const Char_t *testfile,const Char_t *plus="");
TObject* GetSpectrum(AliCFContainer *c, Int_t step);
TObject* GetEfficiency(AliCFContainer *c, Int_t step = 9);
TObject* GetEfficiency(AliCFContainer *c, Int_t step, Int_t step0);


void DrawEfficiencypp(const char *testfile,const char *name) {
  
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

  AliCFContainer *sumcontaineresd = containerhfe->MakeMergedCFContainer("sumesd","sumesd","MCTrackCont:recTrackContReco");
  AliCFContainer *sumcontainermc = containerhfe->MakeMergedCFContainer("summc","summc","MCTrackCont:recTrackContMC");
  
  if(!sumcontaineresd) {
    printf("No container sum esd\n");
    return;
  }
  
  if(!sumcontainermc) {
    printf("No container sum mc\n");
    return;
  }

  Int_t numberOfEvents = (Int_t) containerhfe->GetNumberOfEvents();
  
  printf("Number of events for  %d after Event cut\n",numberOfEvents);
  
  /////////////////////////////////////
  // Take efficiencies
  /////////////////////////////////////
  
  // MC values

  AliCFEffGrid  *efficiencystepMCkineITSTPC  = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepRecKineITSTPC,AliHFEcuts::kStepMCGeneratedZOutNoPileUpCentralityFine);
  AliCFEffGrid  *efficiencystepMCPrim        = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepRecPrim,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepRecKineITSTPC);
  AliCFEffGrid  *efficiencystepMCHFEcutsITS  = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsITS,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepRecPrim);
  AliCFEffGrid  *efficiencystepMCHFEcutsTRD  = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsITS);
  

  AliCFEffGrid  *efficiencystepMCPIDTOF      = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 1,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD);
  AliCFEffGrid  *efficiencystepMCPIDTPC      = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 2,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 1);

  AliCFEffGrid  *efficiencystepMCPIDall      = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 2,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD);

  AliCFEffGrid  *efficiencystepMCall      = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 1,AliHFEcuts::kStepMCGeneratedZOutNoPileUpCentralityFine);
  AliCFEffGrid  *efficiencystepMCallwithTPC  = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 2,AliHFEcuts::kStepMCGeneratedZOutNoPileUpCentralityFine);
  

  ///////////////////////////////////////
  // Plot 1D
  ///////////////////////////////////////
  
  TCanvas * canvascomparisonbis = new TCanvas("ITSTPCrefitStepMC","ITSTPCrefitStepMC",1000,700);
  canvascomparisonbis->Divide(2,2);
 
  canvascomparisonbis->cd(1);
  TH1D* h_effstepMCkineITSTPC1Donly = (TH1D *) efficiencystepMCkineITSTPC->Project(0);
  h_effstepMCkineITSTPC1Donly->SetTitle("(ITS & TPC refit) / (No cut)");
  h_effstepMCkineITSTPC1Donly->SetStats(0);
  h_effstepMCkineITSTPC1Donly->SetLineColor(kBlue);
  h_effstepMCkineITSTPC1Donly->SetMarkerColor(kBlue);
  h_effstepMCkineITSTPC1Donly->SetMarkerStyle(25);
  h_effstepMCkineITSTPC1Donly->SetXTitle("p_{T} [GeV/c]");
  h_effstepMCkineITSTPC1Donly->SetYTitle("Efficiency");
  h_effstepMCkineITSTPC1Donly->SetTitleOffset(1.5,"Y");
  h_effstepMCkineITSTPC1Donly->GetXaxis()->SetRangeUser(0.3,6.0);
  h_effstepMCkineITSTPC1Donly->Draw();

  canvascomparisonbis->cd(2);
  TH1D* h_effstepMCPrim1Donly = (TH1D *) efficiencystepMCPrim->Project(0);
  h_effstepMCPrim1Donly->SetTitle("No Kink");
  h_effstepMCPrim1Donly->SetStats(0);
  h_effstepMCPrim1Donly->SetLineColor(kBlue);
  h_effstepMCPrim1Donly->SetMarkerColor(kBlue);
  h_effstepMCPrim1Donly->SetMarkerStyle(25);
  h_effstepMCPrim1Donly->SetXTitle("p_{T} [GeV/c]");
  h_effstepMCPrim1Donly->SetYTitle("Efficiency");
  h_effstepMCPrim1Donly->SetTitleOffset(1.5,"Y");
  h_effstepMCPrim1Donly->GetXaxis()->SetRangeUser(0.3,6.0);
  h_effstepMCPrim1Donly->Draw();

  canvascomparisonbis->cd(3);
  TH1D* h_effstepMCHFEcutsITS1Donlya = (TH1D *) efficiencystepMCHFEcutsITS->Project(0);
  h_effstepMCHFEcutsITS1Donlya->SetTitle("First pixel");
  h_effstepMCHFEcutsITS1Donlya->SetStats(0);
  h_effstepMCHFEcutsITS1Donlya->SetLineColor(kBlue);
  h_effstepMCHFEcutsITS1Donlya->SetMarkerColor(kBlue);
  h_effstepMCHFEcutsITS1Donlya->SetMarkerStyle(25);
  h_effstepMCHFEcutsITS1Donlya->SetXTitle("p_{T} [GeV/c]");
  h_effstepMCHFEcutsITS1Donlya->SetYTitle("Efficiency");
  h_effstepMCHFEcutsITS1Donlya->SetTitleOffset(1.5,"Y");
  h_effstepMCHFEcutsITS1Donlya->GetXaxis()->SetRangeUser(0.3,6.0);
  h_effstepMCHFEcutsITS1Donlya->Draw();
 
  canvascomparisonbis->cd(4);
  TH1D* h_effstepMCPID1Donly = (TH1D *) efficiencystepMCHFEcutsTRD->Project(0);
  h_effstepMCPID1Donly->SetTitle("TOF matching");
  h_effstepMCPID1Donly->SetStats(0);
  h_effstepMCPID1Donly->SetLineColor(kBlue);
  h_effstepMCPID1Donly->SetMarkerColor(kBlue);
  h_effstepMCPID1Donly->SetMarkerStyle(25);
  h_effstepMCPID1Donly->SetXTitle("p_{T} [GeV/c]");
  h_effstepMCPID1Donly->SetYTitle("Efficiency");
  h_effstepMCPID1Donly->SetTitleOffset(1.5,"Y");
  h_effstepMCPID1Donly->GetXaxis()->SetRangeUser(0.3,6.0);
  h_effstepMCPID1Donly->Draw();
  
  /////////
  // PID
  /////////

  TCanvas * canvasc1DonMC = new TCanvas("PIDstepMC","PIDstepMC",1000,700);
  canvasc1DonMC->Divide(3,1);
  
  canvasc1DonMC->cd(1);
  TH1D* h_effstepPID1DMConly = (TH1D *) efficiencystepMCPIDTOF->Project(0);
  h_effstepPID1DMConly->SetTitle("");
  h_effstepPID1DMConly->SetStats(0);
  h_effstepPID1DMConly->SetLineColor(kBlue);
  h_effstepPID1DMConly->SetMarkerColor(kBlue);
  h_effstepPID1DMConly->SetMarkerStyle(25);
  h_effstepPID1DMConly->SetXTitle("p_{T} [GeV/c]");
  h_effstepPID1DMConly->SetYTitle("Efficiency");
  h_effstepPID1DMConly->SetTitleOffset(1.5,"Y");
  h_effstepPID1DMConly->GetXaxis()->SetRangeUser(0.3,6.0); 
  h_effstepPID1DMConly->Draw();
 
  canvasc1DonMC->cd(2);
  TH1D* h_effstepPID1DMConlyy = (TH1D *) efficiencystepMCPIDTPC->Project(0);
  h_effstepPID1DMConlyy->SetTitle("");
  h_effstepPID1DMConlyy->SetStats(0);
  h_effstepPID1DMConlyy->SetLineColor(kBlue);
  h_effstepPID1DMConlyy->SetMarkerColor(kBlue);
  h_effstepPID1DMConlyy->SetMarkerStyle(25);
  h_effstepPID1DMConlyy->SetXTitle("p_{T} [GeV/c]");
  h_effstepPID1DMConlyy->SetYTitle("Efficiency");
  h_effstepPID1DMConlyy->SetTitleOffset(1.5,"Y");
  h_effstepPID1DMConlyy->GetXaxis()->SetRangeUser(0.3,6.0);
  h_effstepPID1DMConlyy->Draw();
  
  canvasc1DonMC->cd(3);
  TH1D* h_effstepPID1DMConlyyy = (TH1D *) efficiencystepMCPIDall->Project(0);
  h_effstepPID1DMConlyyy->SetTitle("");
  h_effstepPID1DMConlyyy->SetStats(0);
  h_effstepPID1DMConlyyy->SetLineColor(kBlue);
  h_effstepPID1DMConlyyy->SetMarkerColor(kBlue);
  h_effstepPID1DMConlyyy->SetMarkerStyle(25);
  h_effstepPID1DMConlyyy->SetXTitle("p_{T} [GeV/c]");
  h_effstepPID1DMConlyyy->SetYTitle("Efficiency");
  h_effstepPID1DMConlyyy->SetTitleOffset(1.5,"Y");
  h_effstepPID1DMConlyyy->GetXaxis()->SetRangeUser(0.3,6.0);
  h_effstepPID1DMConlyyy->Draw();

  //////////
  // all
  //////////
  
  TCanvas * canvasall = new TCanvas("AllMC","AllMC",1000,700);
  canvasall->cd(1);
  
  TH1D* h_effstepPID1Dall = (TH1D *) efficiencystepMCall->Project(0);
  h_effstepPID1Dall->SetTitle("");
  h_effstepPID1Dall->SetStats(0);
  h_effstepPID1Dall->SetLineColor(kBlue);
  h_effstepPID1Dall->SetMarkerColor(kBlue);
  h_effstepPID1Dall->SetMarkerStyle(25);
  h_effstepPID1Dall->SetXTitle("p_{T} [GeV/c]");
  h_effstepPID1Dall->SetYTitle("Efficiency");
  h_effstepPID1Dall->SetTitleOffset(1.5,"Y");
  h_effstepPID1Dall->GetXaxis()->SetRangeUser(0.3,6.0);
  h_effstepPID1Dall->Scale(0.5);


  TH1D* h_effstepPID1DallwithTPC = (TH1D *) efficiencystepMCallwithTPC->Project(0);
  h_effstepPID1DallwithTPC->SetTitle("");
  h_effstepPID1DallwithTPC->SetStats(0);
  h_effstepPID1DallwithTPC->SetLineColor(kRed);
  h_effstepPID1DallwithTPC->SetMarkerColor(kRed);
  h_effstepPID1DallwithTPC->SetMarkerStyle(26);
  h_effstepPID1DallwithTPC->SetXTitle("p_{T} [GeV/c]");
  h_effstepPID1DallwithTPC->SetYTitle("Efficiency");
  h_effstepPID1DallwithTPC->SetTitleOffset(1.5,"Y");
  h_effstepPID1DallwithTPC->GetXaxis()->SetRangeUser(0.3,6.0);
    
  
  gPad->SetLeftMargin(0.13);
  //gPad->SetLogy();
  gPad->SetTicks();
  //gPad->SetGridx();
  //gPad->SetGridy();
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  h_effstepPID1Dall->Draw();
  h_effstepPID1DallwithTPC->Draw("same");
  TLegend *leg = new TLegend(0.4,0.6,0.89,0.89);
  leg->AddEntry(h_effstepPID1Dall,"TPC PID fixed to 1/2","p");
  leg->AddEntry(h_effstepPID1DallwithTPC,"TPC PID MC efficiency","p");
  leg->Draw("same");
  
  
  
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
