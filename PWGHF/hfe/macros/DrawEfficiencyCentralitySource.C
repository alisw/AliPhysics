TH1I *GetNEvents(const Char_t *testfile,const Char_t *plus="");
void CorrectFromTheWidth(TH1D *h1);
TList *GetResults(const Char_t *testfile,const Char_t *plus="");
TObject* GetSpectrum(AliCFContainer *c, Int_t step);
TObject* GetEfficiency(AliCFContainer *c, Int_t step = 9);
TObject* GetEfficiency(AliCFContainer *c, Int_t step, Int_t step0);
AliCFContainer *GetContainerCentralitySource(AliCFContainer *container,Int_t bincentrality,Int_t source);
AliCFContainer *GetContainerSourceAsFunctionOfCentrality(AliCFContainer *container,Int_t source);

void DrawEfficiencyCentralitySource(Int_t bincentrality,Int_t source,const char *testfile = "/lustre/alice/train/V005.PbPb/2010-11-29_2209.4266/mergedPeriods/data/PbPb/LHC10h.pass1/HFEtask_PbPb.root");


void DrawEfficiencyCentralitySource(Int_t bincentrality,Int_t source,const char *testfile) {
  
  //
  // source: 0 (charm), 1 (beauty), 2 (gamma), 3 (others)
  // centrality:  from 0 to 100 with 0-5, 5-10, 10-15...
  //
  //

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);

  ///////////////////////////////////
  // Take the stuff
  //////////////////////////////////


  TList *results = GetResults(testfile,"");
  if(!results){
    printf("No output objects: Calculation will terminate here\n");
    return;
  }

  AliHFEcontainer *containerhfe = (AliHFEcontainer *) results->FindObject("trackContainer");
  if(!containerhfe) {
    printf("No hfe container \n");
    return;
  }

  // 0: pt, 1: eta, 2: phi, 3: charge, 4: source, 5: centrality
  AliCFContainer *sumcontaineresdd = containerhfe->MakeMergedCFContainer("sumesd","sumesd","MCTrackCont:recTrackContReco");
  AliCFContainer *sumcontainermcc = containerhfe->MakeMergedCFContainer("summc","summc","MCTrackCont:recTrackContMC");

  if(!sumcontaineresdd) {
    printf("No container sum esd\n");
    return;
  }

  if(!sumcontainermcc) {
    printf("No container sum mc\n");
    return;
  }

  // 0: pt, 1: eta, 2: phi for a given centrality and source
  //AliCFContainer *sumcontaineresd = GetContainerSourceCentrality(sumcontaineresdd,bincentrality,source);
  //AliCFContainer *sumcontainermc = GetContainerSourceCentrality(sumcontainermcc,bincentrality,source);

  // 0: pt, 1: eta, 2: phi, 3: centrality for a given centrality and source
  AliCFContainer *sumcontaineresd = GetContainerSourceAsFunctionOfCentrality(sumcontaineresdd,source);
  AliCFContainer *sumcontainermc = GetContainerSourceAsFunctionOfCentrality(sumcontainermcc,source);

 

  Int_t nSteps = sumcontaineresd->GetNStep();
  printf("In total %d steps\n",nSteps);

  AliCFDataGrid *dataGrida = (AliCFDataGrid *) GetSpectrum(sumcontaineresdd,nSteps-1);
 
  TH1D *spectrumcentrality = (TH1D *) dataGrida->Project(5);
  TH1D *spectrumpt = (TH1D *) dataGrida->Project(0);
  TH2D *spectrumptc = (TH2D *) dataGrida->Project(5,0);

  TAxis *xaxis = spectrumptc->GetXaxis();
  Int_t bin0 = xaxis->FindBin(0.0);
  Int_t bin5_s = xaxis->FindBin(4.99);
  Int_t bin30 = xaxis->FindBin(30.0);
  Int_t bin40_s = xaxis->FindBin(39.9);
  Int_t bin70 = xaxis->FindBin(70.0);
  Int_t bin80_s = xaxis->FindBin(79.9);

  printf("Bin 0 %d\n",bin0);
  printf("Bin 5 %d\n",bin5_s);
  printf("Bin 30 %d\n",bin30);
  printf("Bin 40 %d\n",bin40_s);
  printf("Bin 70 %d\n",bin70);
  printf("Bin 80 %d\n",bin80_s);
  
  TH1D *spectrumcentrality_0_5 = spectrumptc->ProjectionY("centrality_0_5",bin0,bin5_s);
  TH1D *spectrumcentrality_30_40 = spectrumptc->ProjectionY("centrality_30_40",bin30,bin40_s);
  TH1D *spectrumcentrality_70_80 = spectrumptc->ProjectionY("centrality_70_80",bin70,bin80_s);

  Int_t numberOfEvents = (Int_t) containerhfe->GetNumberOfEvents();
  
  printf("Number of events for a %d after Event cut\n",numberOfEvents);

  ////////////////////////////////
  // Input after ITS&TPC refit
  ///////////////////////////////
  TCanvas * canvascpt = new TCanvas("RawSpectrumCentrality","RawSpectrumCentrality",1000,700);
  canvascpt->Divide(2,1);
  canvascpt->cd(1);
  spectrumpt->SetTitle("");
  spectrumpt->SetStats(0);
  spectrumpt->SetLineColor(kBlue);
  spectrumpt->SetMarkerColor(kBlue);
  spectrumpt->SetMarkerStyle(25);
  //
  spectrumcentrality_0_5->SetTitle("");
  spectrumcentrality_0_5->SetStats(0);
  spectrumcentrality_0_5->SetLineColor(kRed);
  spectrumcentrality_0_5->SetMarkerColor(kRed);
  spectrumcentrality_0_5->SetMarkerStyle(26);
  //
  spectrumcentrality_30_40->SetTitle("");
  spectrumcentrality_30_40->SetStats(0);
  spectrumcentrality_30_40->SetLineColor(kMagenta);
  spectrumcentrality_30_40->SetMarkerColor(kBlack);
  spectrumcentrality_30_40->SetMarkerStyle(27);
  //
  spectrumcentrality_70_80->SetTitle("");
  spectrumcentrality_70_80->SetStats(0);
  spectrumcentrality_70_80->SetLineColor(kBlue);
  spectrumcentrality_70_80->SetMarkerColor(kBlue);
  spectrumcentrality_70_80->SetMarkerStyle(28);
  //
  spectrumpt->Draw();
  spectrumcentrality_0_5->Draw("same");
  spectrumcentrality_30_40->Draw("same");
  spectrumcentrality_70_80->Draw("same");
  TLegend *leg_different_centralities = new TLegend(0.4,0.6,0.89,0.89);
  leg_different_centralities->AddEntry(spectrumpt,"Minimum-bias","p");
  leg_different_centralities->AddEntry(spectrumcentrality_0_5,"0_5","p");
  leg_different_centralities->AddEntry(spectrumcentrality_30_40,"30_40","p");
  leg_different_centralities->AddEntry(spectrumcentrality_70_80,"70_80","p");
  leg_different_centralities->Draw("same");

  canvascpt->cd(2);
  spectrumptc->Draw("colz");
  
  /////////////////////////////////////
  // Take efficiencies
  /////////////////////////////////////
  
  AliCFEffGrid  *efficiencystepkineITSTPC  = (AliCFEffGrid*)  GetEfficiency(sumcontaineresd,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepRecKineITSTPC,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepRecNoCut);
  AliCFEffGrid  *efficiencystepPrim        = (AliCFEffGrid*)  GetEfficiency(sumcontaineresd,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepRecPrim,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepRecKineITSTPC);
  AliCFEffGrid  *efficiencystepHFEcutsITS  = (AliCFEffGrid*)  GetEfficiency(sumcontaineresd,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsITS,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepRecPrim);
  AliCFEffGrid  *efficiencystepHFEcutsTRD  = (AliCFEffGrid*)  GetEfficiency(sumcontaineresd,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsITS);
  AliCFEffGrid  *efficiencystepPIDTOF      = (AliCFEffGrid*)  GetEfficiency(sumcontaineresd,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD+1,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD);
  AliCFEffGrid  *efficiencystepPIDTPC      = (AliCFEffGrid*)  GetEfficiency(sumcontaineresd,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD+2,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD+1);

  AliCFEffGrid  *efficiencystepPIDall      = (AliCFEffGrid*)  GetEfficiency(sumcontaineresd,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD+2,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD);

  AliCFEffGrid  *efficiencystepall      = (AliCFEffGrid*)  GetEfficiency(sumcontaineresd,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD+2,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepRecNoCut);

  //

  AliCFEffGrid  *efficiencystepMCkineITSTPC  = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepRecKineITSTPC,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepRecNoCut);
  AliCFEffGrid  *efficiencystepMCPrim        = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepRecPrim,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepRecKineITSTPC);
  AliCFEffGrid  *efficiencystepMCHFEcutsITS  = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsITS,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepRecPrim);
  AliCFEffGrid  *efficiencystepMCHFEcutsTRD  = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsITS);
  AliCFEffGrid  *efficiencystepMCPIDTOF      = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD+1,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD);
  AliCFEffGrid  *efficiencystepMCPIDTPC      = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD+2,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD+1);

  AliCFEffGrid  *efficiencystepMCPIDall      = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD+2,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD);

  AliCFEffGrid  *efficiencystepMCall      = (AliCFEffGrid*)  GetEfficiency(sumcontainermc,AliHFEcuts::kNcutStepsMCTrack+AliHFEcuts::kStepHFEcutsTRD+2,0);

  
  ///////////////////////////////////////
  // Plot 1D
  ///////////////////////////////////////

  /////////////////////
  // Different cuts
  /////////////////////
  TCanvas * canvascomparison = new TCanvas("ITSTPCrefitStepESD","ITSTPCrefitStepESD",1000,700);
  canvascomparison->Divide(2,2);
 
  canvascomparison->cd(1);
  TH2D* h_effsteponlykineITSTPC1Donly = (TH2D *) efficiencystepkineITSTPC->Project(0,3);
  h_effsteponlykineITSTPC1Donly->SetTitle("");
  h_effsteponlykineITSTPC1Donly->SetStats(0);
  h_effsteponlykineITSTPC1Donly->SetLineColor(kBlue);
  h_effsteponlykineITSTPC1Donly->SetMarkerColor(kBlue);
  h_effsteponlykineITSTPC1Donly->SetMarkerStyle(25);
  h_effsteponlykineITSTPC1Donly->Draw("colz");
   
  canvascomparison->cd(2);
  TH2D* h_effsteponlyPrim1Donly =( TH2D *) efficiencystepPrim->Project(0,3);
  h_effsteponlyPrim1Donly->SetTitle("");
  h_effsteponlyPrim1Donly->SetStats(0);
  h_effsteponlyPrim1Donly->SetLineColor(kBlue);
  h_effsteponlyPrim1Donly->SetMarkerColor(kBlue);
  h_effsteponlyPrim1Donly->SetMarkerStyle(25);
  h_effsteponlyPrim1Donly->Draw("colz");

  
  canvascomparison->cd(3);
  TH2D* h_effsteponlyHFEcutsITS1Donly = (TH2D *) efficiencystepHFEcutsITS->Project(0,3);
  h_effsteponlyHFEcutsITS1Donly->SetTitle("");
  h_effsteponlyHFEcutsITS1Donly->SetStats(0);
  h_effsteponlyHFEcutsITS1Donly->SetLineColor(kBlue);
  h_effsteponlyHFEcutsITS1Donly->SetMarkerColor(kBlue);
  h_effsteponlyHFEcutsITS1Donly->SetMarkerStyle(25);
  h_effsteponlyHFEcutsITS1Donly->Draw("colz");
 
  
  canvascomparison->cd(4);
  TH2D* h_effsteponlyPID1Donly = (TH2D *) efficiencystepHFEcutsTRD->Project(0,3);
  h_effsteponlyPID1Donly->SetTitle("");
  h_effsteponlyPID1Donly->SetStats(0);
  h_effsteponlyPID1Donly->SetLineColor(kBlue);
  h_effsteponlyPID1Donly->SetMarkerColor(kBlue);
  h_effsteponlyPID1Donly->SetMarkerStyle(25);
  h_effsteponlyPID1Donly->Draw("colz");
 

  /////////////////////
  // Different cuts
  /////////////////////
  TCanvas * canvascomparisonbis = new TCanvas("ITSTPCrefitStepMC","ITSTPCrefitStepMC",1000,700);
  canvascomparisonbis->Divide(2,2);
 
  canvascomparisonbis->cd(1);
  TH2D* h_effstepMCkineITSTPC1Donly =  (TH2D *)efficiencystepMCkineITSTPC->Project(0,3);
  h_effstepMCkineITSTPC1Donly->SetTitle("");
  h_effstepMCkineITSTPC1Donly->SetStats(0);
  h_effstepMCkineITSTPC1Donly->SetLineColor(kBlue);
  h_effstepMCkineITSTPC1Donly->SetMarkerColor(kBlue);
  h_effstepMCkineITSTPC1Donly->SetMarkerStyle(25);
  h_effstepMCkineITSTPC1Donly->Draw("colz");
  

  canvascomparisonbis->cd(2);
  TH2D* h_effstepMCPrim1Donly = (TH2D *) efficiencystepMCPrim->Project(0,3);
  h_effstepMCPrim1Donly->SetTitle("");
  h_effstepMCPrim1Donly->SetStats(0);
  h_effstepMCPrim1Donly->SetLineColor(kBlue);
  h_effstepMCPrim1Donly->SetMarkerColor(kBlue);
  h_effstepMCPrim1Donly->SetMarkerStyle(25);
  h_effstepMCPrim1Donly->Draw("colz");
  
  canvascomparisonbis->cd(3);
  TH2D* h_effstepMCHFEcutsITS1Donly = (TH2D *) efficiencystepMCHFEcutsITS->Project(0,3);
  h_effstepMCHFEcutsITS1Donly->SetTitle("");
  h_effstepMCHFEcutsITS1Donly->SetStats(0);
  h_effstepMCHFEcutsITS1Donly->SetLineColor(kBlue);
  h_effstepMCHFEcutsITS1Donly->SetMarkerColor(kBlue);
  h_effstepMCHFEcutsITS1Donly->SetMarkerStyle(25);
  h_effstepMCHFEcutsITS1Donly->Draw("colz");
  
  canvascomparisonbis->cd(4);
  TH2D* h_effstepMCPID1Donly = (TH2D *) efficiencystepMCHFEcutsTRD->Project(0,3);
  h_effstepMCPID1Donly->SetTitle("");
  h_effstepMCPID1Donly->SetStats(0);
  h_effstepMCPID1Donly->SetLineColor(kBlue);
  h_effstepMCPID1Donly->SetMarkerColor(kBlue);
  h_effstepMCPID1Donly->SetMarkerStyle(25);
  h_effstepMCPID1Donly->Draw("colz");
  

  ///////////
  // PID
  ///////////
  TCanvas * canvasc1Don = new TCanvas("PIDstepESD","PIDstepESD",1000,700);
  canvasc1Don->Divide(3,1);
  
  canvasc1Don->cd(1);
  TH2D* h_effstepPID1Donly = (TH2D *) efficiencystepPIDTOF->Project(0,3);
  h_effstepPID1Donly->SetTitle("");
  h_effstepPID1Donly->SetStats(0);
  h_effstepPID1Donly->SetLineColor(kBlue);
  h_effstepPID1Donly->SetMarkerColor(kBlue);
  h_effstepPID1Donly->SetMarkerStyle(25);
  gPad->SetLogz();
  h_effstepPID1Donly->Draw("colz");
  
  canvasc1Don->cd(2);
  TH2D* h_effstepPID1Donlyy = (TH2D *) efficiencystepPIDTPC->Project(0,3);
  h_effstepPID1Donlyy->SetTitle("");
  h_effstepPID1Donlyy->SetStats(0);
  h_effstepPID1Donlyy->SetLineColor(kBlue);
  h_effstepPID1Donlyy->SetMarkerColor(kBlue);
  h_effstepPID1Donlyy->SetMarkerStyle(25);
  gPad->SetLogz();
  h_effstepPID1Donlyy->Draw("colz");
  
  canvasc1Don->cd(3);
  TH2D* h_effstepPID1Donlyyy = (TH2D *) efficiencystepPIDall->Project(0,3);
  h_effstepPID1Donlyyy->SetTitle("");
  h_effstepPID1Donlyyy->SetStats(0);
  h_effstepPID1Donlyyy->SetLineColor(kBlue);
  h_effstepPID1Donlyyy->SetMarkerColor(kBlue);
  h_effstepPID1Donlyyy->SetMarkerStyle(25);
  gPad->SetLogz();
  h_effstepPID1Donlyyy->Draw("colz");
  
  //

  TCanvas * canvasc1DonMC = new TCanvas("PIDstepMC","PIDstepMC",1000,700);
  canvasc1DonMC->Divide(3,1);
  
  canvasc1DonMC->cd(1);
  TH2D* h_effstepPID1DMConly = (TH2D *) efficiencystepMCPIDTOF->Project(0,3);
  h_effstepPID1DMConly->SetTitle("");
  h_effstepPID1DMConly->SetStats(0);
 
  h_effstepPID1DMConly->SetLineColor(kBlue);
 
  h_effstepPID1DMConly->SetMarkerColor(kBlue);
 
  h_effstepPID1DMConly->SetMarkerStyle(25);
  gPad->SetLogz(); 
  h_effstepPID1DMConly->Draw("colz");
 

  canvasc1DonMC->cd(2);
  TH2D* h_effstepPID1DMConlyy = (TH2D *) efficiencystepMCPIDTPC->Project(0,3);
  h_effstepPID1DMConlyy->SetTitle("");
  h_effstepPID1DMConlyy->SetStats(0);
  h_effstepPID1DMConlyy->SetLineColor(kBlue);
  h_effstepPID1DMConlyy->SetMarkerColor(kBlue);
  h_effstepPID1DMConlyy->SetMarkerStyle(25);
  gPad->SetLogz();
  h_effstepPID1DMConlyy->Draw("colz");
 

  canvasc1DonMC->cd(3);
  TH2D* h_effstepPID1DMConlyyy = (TH2D *) efficiencystepMCPIDall->Project(0,3);
  h_effstepPID1DMConlyyy->SetTitle("");
  h_effstepPID1DMConlyyy->SetStats(0);
  h_effstepPID1DMConlyyy->SetLineColor(kBlue);
  h_effstepPID1DMConlyyy->SetMarkerColor(kBlue);
  h_effstepPID1DMConlyyy->SetMarkerStyle(25);
  gPad->SetLogz();
  h_effstepPID1DMConlyyy->Draw("colz");


  //////////
  // all
  //////////

  TCanvas * canvasall = new TCanvas("AllMC","AllMC",1000,700);
  //canvasall->Divide(2,1);

  canvasall->cd(1);
  TH2D* h_effstepPID1Dall = (TH2D *) efficiencystepMCall->Project(0,3);
  h_effstepPID1Dall->SetTitle("");
  h_effstepPID1Dall->SetStats(0);
  h_effstepPID1Dall->SetLineColor(kBlue);
  h_effstepPID1Dall->SetMarkerColor(kBlue);
  h_effstepPID1Dall->SetMarkerStyle(25);
  gPad->SetLogz();
  h_effstepPID1Dall->Draw("colz");

  
 


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
  TString name(plus);
  name += "HFE_Results"; 
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
//_______________________________________________________________________
TH1I *GetNEvents(const Char_t *testfile,const Char_t *plus){
  //
  // read output
  //
  TFile *f = TFile::Open(testfile);
  if(!f || f->IsZombie()){
    printf("File not readable\n");
    return 0x0;
  }
  TString name("nEvents");
  name += plus; 
  printf("Name of nEvents %s\n",(const char*)name); 
  TH1I *l = dynamic_cast<TH1I *>(f->Get((const char*)name));
  if(!l){
    printf("nEvents not found\n");
    f->Close(); delete f;
    return 0x0;
  } 
  TH1I *returnlist = dynamic_cast<TH1I *>(l->Clone());
  if(!returnlist) return 0x0;
  returnlist->SetDirectory(0);
  f->Close(); delete f;
  return returnlist;
}
//________________________________________________________________________
void CorrectFromTheWidth(TH1D *h1) {
  //
  // Correct from the width of the bins --> dN/dp_{T} (GeV/c)^{-1}
  //

  TAxis *axis = h1->GetXaxis();
  Int_t nbinX = h1->GetNbinsX();

  for(Int_t i = 1; i <= nbinX; i++) {

    Double_t width = axis->GetBinWidth(i);
    Double_t content = h1->GetBinContent(i);
    Double_t error = h1->GetBinError(i); 
    h1->SetBinContent(i,content/width); 
    h1->SetBinError(i,error/width);
  }

}
//____________________________________________________________________________
AliCFContainer *GetContainerCentralitySource(AliCFContainer *container,Int_t bincentrality,Int_t source) {

  Int_t *vars = new Int_t[3];
  vars[0] = 0;
  vars[1] = 1;
  vars[2] = 2;
     
  Int_t nbBinPt = container->GetNBins(0);
  Int_t nbBinEta = container->GetNBins(1);
  Int_t nbBinPhi = container->GetNBins(2);
  Int_t nbBinCharge = container->GetNBins(3);
  Int_t nbBinSource = container->GetNBins(4);
  Int_t nbBinCentrality = container->GetNBins(5);
  Double_t *arrayPt = new Double_t[nbBinPt + 1];
  Double_t *arrayEta = new Double_t[nbBinEta + 1];
  Double_t *arrayPhi = new Double_t[nbBinPhi + 1];
  Double_t *arrayCharge = new Double_t[nbBinCharge + 1];
  Double_t *arraySource = new Double_t[nbBinSource + 1];
  Double_t *arrayCentrality = new Double_t[nbBinCentrality + 1];
  container->GetBinLimits(0,arrayPt);
  container->GetBinLimits(1,arrayEta);
  container->GetBinLimits(2,arrayPhi);
  container->GetBinLimits(3,arrayCharge);
  container->GetBinLimits(4,arraySource);
  container->GetBinLimits(4,arrayCentrality);

  Int_t sourcebin = source;
  Int_t centralitybin = bincentrality;
 

  Double_t *varMin = new Double_t[6];
  Double_t *varMax = new Double_t[6];
  varMin[0] = arrayPt[0];
  varMin[1] = arrayEta[0];
  varMin[2] = arrayPhi[0];
  varMin[3] = arrayCharge[0];
  if(source == -1) sourcebin = 0;
  if(bincentrality == -1) centralitybin = 0;
  varMin[4] = arraySource[sourcebin];
  varMin[5] = arraySource[centralitybin];
  varMax[0] = arrayPt[nbBinPt];
  varMax[1] = arrayEta[nbBinEta];
  varMax[2] = arrayPhi[nbBinPhi];
  varMax[3] = arrayCharge[nbBinCharge];
  if(source == -1) sourcebin = nbBinSource;
  if(bincentrality == -1) centralitybin = nbBinCentrality;
  varMax[4] = arraySource[sourcebin];
  varMax[5] = arraySource[centralitybin];

  //printf("Nb of bin charge %d\n",nbBinCharge);
  //for(Int_t nb = 0; nb <= nbBinCharge; nb++) {
  //  printf("charge %f for %d\n",arrayCharge[nb],nb);
  //}

  AliCFContainer *k = container->MakeSlice(3,vars,varMin,varMax);
  k->SetName("Other");

  return k;

}
//____________________________________________________________________________
AliCFContainer *GetContainerSourceAsFunctionOfCentrality(AliCFContainer *container,Int_t source) {

  Int_t *vars = new Int_t[4];
  vars[0] = 0;
  vars[1] = 1;
  vars[2] = 2;
  vars[3] = 5;
     
  Int_t nbBinPt = container->GetNBins(0);
  Int_t nbBinEta = container->GetNBins(1);
  Int_t nbBinPhi = container->GetNBins(2);
  Int_t nbBinCharge = container->GetNBins(3);
  Int_t nbBinSource = container->GetNBins(4);
  Int_t nbBinCentrality = container->GetNBins(5);
  Double_t *arrayPt = new Double_t[nbBinPt + 1];
  Double_t *arrayEta = new Double_t[nbBinEta + 1];
  Double_t *arrayPhi = new Double_t[nbBinPhi + 1];
  Double_t *arrayCharge = new Double_t[nbBinCharge + 1];
  Double_t *arraySource = new Double_t[nbBinSource + 1];
  Double_t *arrayCentrality = new Double_t[nbBinCentrality + 1];
  container->GetBinLimits(0,arrayPt);
  container->GetBinLimits(1,arrayEta);
  container->GetBinLimits(2,arrayPhi);
  container->GetBinLimits(3,arrayCharge);
  container->GetBinLimits(4,arraySource);
  container->GetBinLimits(5,arrayCentrality);

  Int_t sourcebin = source;
 

  Double_t *varMin = new Double_t[6];
  Double_t *varMax = new Double_t[6];
  varMin[0] = arrayPt[0];
  varMin[1] = arrayEta[0];
  varMin[2] = arrayPhi[0];
  varMin[3] = arrayCharge[0];
  if(source == -1) sourcebin = 0;
  varMin[4] = arraySource[sourcebin];
  varMin[5] = arrayCentrality[0];
  varMax[0] = arrayPt[nbBinPt];
  varMax[1] = arrayEta[nbBinEta];
  varMax[2] = arrayPhi[nbBinPhi];
  varMax[3] = arrayCharge[nbBinCharge];
  if(source == -1) sourcebin = nbBinSource;
  varMax[4] = arraySource[sourcebin];
  varMax[5] = arrayCentrality[nbBinCentrality];

  //printf("Nb of bin charge %d\n",nbBinCharge);
  //for(Int_t nb = 0; nb <= nbBinCharge; nb++) {
  //  printf("charge %f for %d\n",arrayCharge[nb],nb);
  //}

  AliCFContainer *k = container->MakeSlice(4,vars,varMin,varMax);
  k->SetName("OtherAsFunctionOfCentrality");

  return k;

}
