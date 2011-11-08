void CorrectForEfficiencypptrd(const char *filedata,const char *fileMC,const char *NameList,Bool_t withtrd=kTRUE, Float_t CorrType=0, Int_t conf=2, Bool_t smoothing=kFALSE, Bool_t hadroncontaminationsubtracted=kTRUE, Bool_t UnsetCorrelatedErrors=kFALSE);
void CorrectForEfficiencyBeauty(const char *filedata,const char *fileMC, const char *fileMCbg);
TList *GetResults(const Char_t *testfile,const Char_t *plus="");
TList *GetQA(const Char_t *testfile,const Char_t *plus="");
TObject* GetSpectrum(AliCFContainer *c, Int_t step);
THnSparseF* GetHadronEffbyIPcut(TList *hfeqa);
void CorrectFromTheWidth(TH1D *h1);

void CorrectForEfficiencypptrd(const char *filedata,const char *fileMC,const char *NameList,Bool_t withtrd, Float_t CorrType, Int_t conf, Bool_t smoothing, Bool_t hadroncontaminationsubtracted, Bool_t UnsetCorrelatedErrors) {

   ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Will produce finalSpectrum.root with results inside
  // TGraphErrors UnfoldingCorrectedSpectrum -> unfolding procedure (doesn't work in 2D with trunk)
  // TGraphErrors AlltogetherSpectrum -> corrected spectrum with other method
  // TH1D RatioUnfoldingAlltogetherSpectrum 
  // THnSparseF UnfoldingCorrectedNotNormalizedSpectrum -> unfolding procedure not yet normalized
  // AliCFDataGrid AlltogetherCorrectedNotNormalizedSpectrum -> other method not yet normalized
  //
  // NameList - possible extension for the full name of the data and MC TLists (e.g. "_PID2" for HFE_Results_PID2)
  // CorrType - 0 use fixed TPC PID 0.5
  //            1 use TPC PID efficiency from V0s 
  //
  // Conf - TRD efficiency
  // double TRDeffs[5] = {0.7,0.8,0.9,0.75,0.85};
  // conf = 2, TRDeff = 0.8...
  // (Only relevant for withtrd = kTRUE)
  //
  // f/ smoothing kFALSE
  //
  // g/ Flag for hadron contamination subtraction
  // kTRUE means we subtract the hadron contamination
  // kFALSE not
  //
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  double TRDeffs[5] = {0.7,0.8,0.9,0.75,0.85};

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);

  if (CorrType>1){
    printf("Wrong argument for the type of correction, it is maybe missing, or it has to be <=2!!! \n");
    return;
  }

  /////////////////////
  // Take the stuff
  /////////////////////
  
  printf("Get results for data \n");
  TList *resultsdata = GetResults(filedata,NameList);
  if(!resultsdata){
    printf("No output objects for data: Calculation will terminate here\n");
    return;
  }

  ///////////////////////////////////////////////////////////
  // Event container for the normalization to CINB1
  // Take the number of events as function of z and number of contributors
  // For those without primary vertex from tracks, take the fraction in the z range from MC
  // For 10cm: LHC10f6a (88.0756%), LHC10f6 (86.7461%)
  ////////////////////////////////////////////////////////////
  resultsdata->Print();

  AliCFContainer *containerdata = (AliCFContainer *) resultsdata->FindObject("eventContainer");
  if(!containerdata) {
    printf("No container \n");
    return;
  }
  AliCFDataGrid *dataGrid = (AliCFDataGrid *) GetSpectrum(containerdata,AliHFEcuts::kEventStepZRange);
  //dataGrid->PrintBinLimits();
  dataGrid->PrintNBins();
  
  // 
  TH2D *spectrum_zrange = (TH2D *) dataGrid->Project(0,3);
  TAxis *axis_0 = spectrum_zrange->GetXaxis();
  TAxis *axis_1 = spectrum_zrange->GetYaxis();
  printf("Number of bins in x %d and in y %d\n",axis_0->GetNbins(),axis_1->GetNbins());
  TH1D *spectrum1Da = (TH1D *) spectrum_zrange->ProjectionX("bin0",1,1);
  TH1D *spectrum1Db = (TH1D *) spectrum_zrange->ProjectionX("bin>0",2,12);
  TH1D *spectrum1Dc = (TH1D *) spectrum_zrange->ProjectionX("total");
  Double_t nbinbin0 = spectrum1Da->Integral();
  Double_t nbinnobin0 = spectrum1Db->Integral();
  Double_t nbintotal = spectrum1Dc->Integral();

  //TCanvas *c1teste = new TCanvas("eventcontainertest","eventcontainertest",800,800);
  //c1teste->cd(1);
  //spectrum_zrange->Draw("colz");
  

  Float_t numberofentries = nbinnobin0+nbinbin0*0.880756;
  printf("Number in bin0 %f, number out %f, total %f, normalized %f\n",nbinbin0,nbinnobin0,nbintotal,numberofentries);
  
  //////////////////////////////
  // Take more stuff
  ///////////////////////////////
 
  printf("Get results for MC \n");
  TList *resultsmc = GetResults(fileMC,NameList);
  if(!resultsmc){
    printf("No output objects for mc: Calculation will terminate here\n");
    return;
  }

  AliHFEcontainer *datahfecontainer = (AliHFEcontainer *) resultsdata->FindObject("trackContainer");
  if(!datahfecontainer) {
    printf("No container for data \n");
    return;
  }
  AliCFContainer *sumcontainer = datahfecontainer->GetCFContainer("recTrackContReco");

  ////////////////////////////////////////////
  //Plot raw spectrum for TPC TOF scenario
  ////////////////////////////////////////////
  AliCFDataGrid *dataGrida = 0x0;
  if(!withtrd) dataGrida = (AliCFDataGrid *) GetSpectrum(sumcontainer,AliHFEcuts::kStepHFEcutsTRD + 2);
  else dataGrida = (AliCFDataGrid *) GetSpectrum(sumcontainer,AliHFEcuts::kStepHFEcutsTRD + 3);
 
  TH1D *spectrumpt =  (TH1D *) dataGrida->Project(0);
  CorrectFromTheWidth(spectrumpt);
  
  TH1D *total = new TH1D("total","",1,0.3,10.3);
  total->SetMaximum(1.0e+07);
  total->SetMinimum(1000);
  total->SetXTitle("p_{T} [GeV/c]");
  total->SetYTitle("dN/dp_{T}[GeV/c]^{-1}");
  total->SetTitleOffset(1.5,"Y");
  //total->Scale(0.9);
  total->GetXaxis()->SetRangeUser(0.38,10.3);
   

  TCanvas *c1test = new TCanvas("rawspectrum","rawspectrum",800,800);
  c1test->cd(1);
  gPad->SetLeftMargin(0.13);
  gPad->SetLogy();
  gPad->SetTicks();
  //gPad->SetGridx();
  //gPad->SetGridy();
  gPad->SetFillColor(10);
  gPad->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1111);
  gStyle->SetOptTitle(0);
  total->SetStats(0);
  spectrumpt->SetStats(0);
  total->Draw();
  spectrumpt->Draw("same");
  TPaveText* t1=new TPaveText(0.49,0.45,0.79,0.52,"NDC");
  t1->SetFillStyle(0);
  t1->SetBorderSize(0);
  t1->AddText(0.,0.,"pp, #sqrt{s} =  7 TeV");
  t1->SetTextColor(kRed);
  //t1->SetTextSize(20);
  t1->SetTextFont(42);
  t1->Draw();
  TPaveText* t2=new TPaveText(0.49,0.35,0.79,0.42,"NDC");
  t2->SetFillStyle(0);
  t2->SetBorderSize(0);
  t2->AddText(0.,0.,"L = 2.6 nb^{-1}");
  t2->SetTextColor(kRed);
  //t1->SetTextSize(20);
  t2->SetTextFont(42);
  t2->Draw();
  TPaveText* t3=new TPaveText(0.49,0.35,0.79,0.42,"NDC");
  t3->SetFillStyle(0);
  t3->SetBorderSize(0);
  t3->AddText(0.,0.,"|#eta| < 0.8");
  t3->SetTextColor(kRed);
  //t1->SetTextSize(20);
  t3->SetTextFont(42);
  t3->Draw();
  
  /////////////////////////////
  // Check number of events
  /////////////////////////////


  Int_t numberOfEventsafterallcuts = (Int_t) datahfecontainer->GetNumberOfEvents();
   
  printf("Number of events not corrected %d\n",numberOfEventsafterallcuts);
  printf("Number of events corrected  %f\n",numberofentries);

  AliHFEcontainer *mchfecontainer = (AliHFEcontainer *) resultsmc->FindObject("trackContainer");
  if(!mchfecontainer) {
    printf("No mc container \n");
    return;
  }
 
  printf("Find the container V0\n");   
  AliHFEcontainer *containerhfeV0 = (AliHFEcontainer *) resultsdata->FindObject("containerV0");
  if(!containerhfeV0) {
    printf("No hfe container \n");
    return;
  }
    
  //////////////
  // Correct
  /////////////

  AliHFEspectrum *spectrum = new AliHFEspectrum("HFEspectrum");
  spectrum->SetNbDimensions(1);
  // If you want to correct positive (0) or negative (1) separately
  //spectrum->SetChargeChoosen(0);
  spectrum->SetSmoothing(smoothing);
  spectrum->SetUnSetCorrelatedErrors(UnsetCorrelatedErrors);
  spectrum->SetDebugLevel(1);
  spectrum->SetNumberOfEvents((Int_t)numberofentries);
  spectrum->SetDumpToFile(kTRUE);
  // True step in MC (events in range +- 10 cm and no pile up)
  spectrum->SetMCTruthStep(AliHFEcuts::kStepMCGeneratedZOutNoPileUpCentralityFine);

  // Step where we correct from MC (tracking + TOF)
  spectrum->SetMCEffStep(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 1);
  printf("!!! Step we correct with MC %d\n!!!",AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 1);
 

  // Step from data we correct for
  //////////****************************//////////////////////////////
  if(withtrd) spectrum->SetStepToCorrect(AliHFEcuts::kStepHFEcutsTRD + 3);
  else spectrum->SetStepToCorrect(AliHFEcuts::kStepHFEcutsTRD + 2);
  //////////****************************//////////////////////////////

  // Steps to be corrected with V0 
  if (CorrType == 1){
    printf("\n This time we correct with PID efficiency from V0s !!! \n");
    if(withtrd) {
      spectrum->SetStepBeforeCutsV0(AliHFEcuts::kStepHFEcutsTRD + 2);
      spectrum->SetStepAfterCutsV0(AliHFEcuts::kStepHFEcutsTRD + 3);
    }
    else {
      spectrum->SetStepBeforeCutsV0(AliHFEcuts::kStepHFEcutsTRD + 1);
      spectrum->SetStepAfterCutsV0(AliHFEcuts::kStepHFEcutsTRD + 2);
    }
  }

  if (CorrType < 1){
 
    Double_t effTPC = 0.5,
      effTRD = TRDeffs[conf-1]; 
    printf("effTRD = %f\n", effTRD);
    // Place for a correction as an efficiency parametrized as function of pt of pt,eta,phi or p=pt
    TF1 *hEfficiency = new TF1("efficiency", "[0]", 0., 20.);
    if(withtrd) hEfficiency->SetParameter(0, effTPC*effTRD);
    else hEfficiency->SetParameter(0, effTPC);
    spectrum->SetEfficiencyFunction(hEfficiency);
  }  

  // Give everything (data container, mc container and V0 data container)
  if (CorrType == 1){
    printf("\n This time we correct with PID efficiency from V0s !!! \n");
    spectrum->Init(datahfecontainer,mchfecontainer,containerhfeV0);
  }
  else {
    printf("\n This time we correct with PID efficiency from MC or fixed value !!! \n");
    spectrum->Init(datahfecontainer,mchfecontainer);
  }
 
  // kTRUE means subtract hadron background, kFALSE means do not subtract hadron background
  spectrum->Correct(hadroncontaminationsubtracted);
  

}

//_____________________________________________________________________
void CorrectForEfficiencyBeauty(const char *filedata,const char *fileMC, const char *fileMCbg) {
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Will produce finalSpectrum.root with results inside
  // TGraphErrors UnfoldingCorrectedSpectrum -> unfolding procedure (doesn't work in 2D with trunk)
  // TGraphErrors AlltogetherSpectrum -> corrected spectrum with other method
  // TH1D RatioUnfoldingAlltogetherSpectrum 
  // THnSparseF UnfoldingCorrectedNotNormalizedSpectrum -> unfolding procedure not yet normalized
  // AliCFDataGrid AlltogetherCorrectedNotNormalizedSpectrum -> other method not yet normalized
  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  gStyle->SetPalette(1);
  gStyle->SetOptStat(1111);
  gStyle->SetPadBorderMode(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);

  /////////////////////
  // Take the stuff
  /////////////////////

  TList *resultsdata = GetResults(filedata,"_PID2");
  //TList *resultsdata = GetResults(filedata);
  if(!resultsdata){
    printf("No output objects for data: Calculation will terminate here\n");
    return;
  }
  ///////////////////////////////////////////////////////////
  // Event container for the normalization to CINB1
  // Take the number of events as function of z and number of contributors
  // For those without primary vertex from tracks, take the fraction in the z range from MC
  // For 10cm: LHC10f6a (88.0756%), LHC10f6 (86.7461%)
  ////////////////////////////////////////////////////////////

  AliCFContainer *containerdata = (AliCFContainer *) resultsdata->FindObject("eventContainer");
  if(!containerdata) {
    printf("No container \n");
    return;
  }
  AliCFDataGrid *dataGrid = (AliCFDataGrid *) GetSpectrum(containerdata,AliHFEcuts::kEventStepZRange);
  TH2D *spectrum_zrange = (TH2D *) dataGrid->Project(0,2);
  TH1D *spectrum1Da = (TH1D *) spectrum_zrange->ProjectionX("bin0",1,1);
  TH1D *spectrum1Db = (TH1D *) spectrum_zrange->ProjectionX("bin>0",2,12);
  TH1D *spectrum1Dc = (TH1D *) spectrum_zrange->ProjectionX("total");
  Double_t nbinbin0 = spectrum1Da->Integral();
  Double_t nbinnobin0 = spectrum1Db->Integral();
  Double_t nbintotal = spectrum1Dc->Integral();

  Float_t numberofentries = nbinnobin0+nbinbin0*0.880756;
  printf("Number in bin0 %f, number out %f, total %f, normalized %f\n",nbinbin0,nbinnobin0,nbintotal,numberofentries);


  AliHFEcontainer *datahfecontainer = (AliHFEcontainer *) resultsdata->FindObject("trackContainer");
  if(!datahfecontainer) {
    printf("No container for data \n");
    return;
  }


  //////////////////////////////
  // Check MC # of events
  ///////////////////////////////

  TList *resultsmc = GetResults(fileMC,"_PID2");
  if(!resultsmc){
    printf("No output objects for mc: Calculation will terminate here\n");
    return;
  }

  AliCFContainer *containermc = (AliCFContainer *) resultsmc->FindObject("eventContainer");
  if(!containermc) {
    printf("No container \n");
    return;
  }
  AliCFDataGrid *mcGrid = (AliCFDataGrid *) GetSpectrum(containermc,AliHFEcuts::kEventStepZRange);
  TH2D *spectrum_zrangemc = (TH2D *) mcGrid->Project(0,2);
  TH1D *spectrum1Damc = (TH1D *) spectrum_zrangemc->ProjectionX("bin0",1,1);
  TH1D *spectrum1Dbmc = (TH1D *) spectrum_zrangemc->ProjectionX("bin>0",2,12);
  TH1D *spectrum1Dcmc = (TH1D *) spectrum_zrangemc->ProjectionX("total");
  Double_t nbinbin0mc = spectrum1Damc->Integral();
  Double_t nbinnobin0mc = spectrum1Dbmc->Integral();
  Double_t nbintotalmc = spectrum1Dcmc->Integral();

  Float_t numberofentriesmc = nbinnobin0mc+nbinbin0mc*0.880756;
  printf("MC: Number in bin0 %f, number out %f, total %f, normalized %f\n",nbinbin0mc,nbinnobin0mc,nbintotalmc,numberofentriesmc);



  /////////////////////////////
  // Check number of events
  /////////////////////////////

  Int_t numberOfEventsafterallcuts = (Int_t) datahfecontainer->GetNumberOfEvents();

  printf("Number of events not corrected %d\n",numberOfEventsafterallcuts);
  printf("Number of events corrected  %f\n",numberofentries);

  AliHFEcontainer *mchfecontainer = (AliHFEcontainer *) resultsmc->FindObject("trackContainer");
  if(!mchfecontainer) {
    printf("No mc container \n");
    return;
  }

  printf("Find the container V0\n");
  AliHFEcontainer *containerhfeV0 = (AliHFEcontainer *) resultsdata->FindObject("containerV0");
  if(!containerhfeV0) {
    printf("No hfe container \n");
    return;
  }

  // nonHFE backgrounds
  TList *resultsmcbg = GetResults(fileMCbg,"_PID2");
  if(!resultsmcbg){
    printf("No output objects for mc: Calculation will terminate here\n");    return;
  }  
  AliHFEcontainer *mcbghfecontainer = (AliHFEcontainer *) resultsmcbg->FindObject("trackContainer");
  if(!mcbghfecontainer) {
    printf("No mc container \n");
    return;
  }

  AliCFContainer *containermcbg = (AliCFContainer *) resultsmcbg->FindObject("eventContainer");
  if(!containermcbg) {
    printf("No container \n");
    return;
  }

  AliCFDataGrid *mcbgGrid = (AliCFDataGrid *) GetSpectrum(containermcbg,AliHFEcuts::kEventStepZRange);
  TH2D *spectrum_zrangemcbg = (TH2D *) mcbgGrid->Project(0,2);
  TH1D *spectrum1Damcbg = (TH1D *) spectrum_zrangemcbg->ProjectionX("bin0",1,1);
  TH1D *spectrum1Dbmcbg = (TH1D *) spectrum_zrangemcbg->ProjectionX("bin>0",2,12);
  TH1D *spectrum1Dcmcbg = (TH1D *) spectrum_zrangemcbg->ProjectionX("total");
  Double_t nbinbin0mcbg = spectrum1Damcbg->Integral();
  Double_t nbinnobin0mcbg = spectrum1Dbmcbg->Integral();
  Double_t nbintotalmcbg = spectrum1Dcmcbg->Integral();

  Float_t numberofentriesmcbg = nbinnobin0mcbg+nbinbin0mcbg*0.880756;
  printf("MCbg: Number in bin0 %f, number out %f, total %f, normalized %f\n",nbinbin0mcbg,nbinnobin0mcbg,nbintotalmcbg,numberofentriesmcbg);


  // hadron contamination after IP cuts
  TList *hfeqa = GetQA(fileMC,"_PID2");
  THnSparseF* hsHadronEffbyIPcut = (THnSparseF* )GetHadronEffbyIPcut(hfeqa);

  //////////////
  // Correct
  /////////////

  AliHFEspectrum *spectrum = new AliHFEspectrum("HFEspectrum");

  spectrum->SetBeautyAnalysis();
  // Enable background subtraction
  spectrum->EnableIPanaHadronBgSubtract();
  spectrum->EnableIPanaCharmBgSubtract();
  spectrum->EnableIPanaConversionBgSubtract();
  spectrum->EnableIPanaNonHFEBgSubtract();
  //spectrum->SetNonHFEBackground2ndMethod();

  spectrum->SetNbDimensions(1);
  // If you want to correct positive (0) or negative (1) separately
  //spectrum->SetChargeChoosen(0);
  spectrum->SetDebugLevel(1);
  spectrum->SetNumberOfEvents((Int_t)numberofentries);
  spectrum->SetNumberOfMCEvents((Int_t)numberofentriesmc);
  spectrum->SetNumberOfMC2Events((Int_t)numberofentriesmcbg);
  spectrum->SetDumpToFile(smoothing);
  // True step in MC (events in range +- 10 cm and no pile up)
  spectrum->SetMCTruthStep(AliHFEcuts::kStepMCGeneratedZOutNoPileUpCentralityFine);
  // Step where we correct from MC (tracking + TOF)
  spectrum->SetMCEffStep(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD + 1);
  // Step from data we correct for
  spectrum->SetStepToCorrect(AliHFEcuts::kStepHFEcutsTRD + 3); // +3 = tracking + TOF + TPC + IPcut
  // Steps to be corrected with V0 (TPC PID only)
  spectrum->SetStepBeforeCutsV0(AliHFEcuts::kStepHFEcutsTRD + 1);
  spectrum->SetStepAfterCutsV0(AliHFEcuts::kStepHFEcutsTRD + 2);

  spectrum->SetHadronEffbyIPcut(hsHadronEffbyIPcut);
  // Give everything (data container, mc container and V0 data container)
  spectrum->Init(datahfecontainer,mchfecontainer,containerhfeV0,mcbghfecontainer);

  // kTRUE means subtract hadron background, kFALSE means do not subtract hadron background
  spectrum->CorrectBeauty(kTRUE);

}

//_____________________________________________________________________
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
  printf("Name of TList %s\n",(const char*)name); 
  TList *l = dynamic_cast<TList *>(f->Get((const char*)name));
  if(!l){
    printf("Results list was not found\n");
    f->Close(); delete f;
    return 0x0;
  } 
  TList *returnlist = dynamic_cast<TList *>(l->Clone());
  f->Close(); delete f;
  return returnlist;
}

//_____________________________________________________________________
TList *GetQA(const Char_t *testfile,const Char_t *plus){
  //
  // read output
  //
  TFile *f = TFile::Open(testfile);
  if(!f || f->IsZombie()){
    printf("File not readable\n");
    return 0x0;
  }
  TString name("HFE_QA");
  name += plus;
  printf("Name of TList %s\n",(const char*)name);
  TList *l = dynamic_cast<TList *>(f->Get((const char*)name));
  if(!l){
    printf("QA list was not found\n");
    f->Close(); delete f;
    return 0x0;
  }
  TList *returnlist = dynamic_cast<TList *>(l->Clone());
  f->Close(); delete f;
  return returnlist;
}

//_____________________________________________________________________
THnSparseF* GetHadronEffbyIPcut(TList *hfeqa){


  // get hadron reduction factors due to the IP cuts
  TList* tl=(TList*)hfeqa->FindObject("list_TaskQA");
  TH1F* hbefore=(TH1F*)tl->FindObject("hadronsBeforeIPcut");
  TH1F* hafter=(TH1F*)tl->FindObject("hadronsAfterIPcut");
  TH1F* hreduction= (TH1F*)hafter->Clone("hreduction");
  hbefore->Sumw2(); 
  hafter->Sumw2();
  hreduction->Sumw2();
  hreduction->Divide(hbefore);

  Double_t* binEdges[0];
  Int_t hnBin = hreduction->GetNbinsX();
  Int_t nBin[1] = {hnBin};

  for(int i=0; i<nBin[0]; i++){
    binEdges[0][i] = hreduction->GetBinLowEdge(i+1);
  }
  binEdges[0][nBin[0]] = hreduction->GetBinLowEdge(nBin[0]) + hreduction->GetBinWidth(nBin[0]);

  THnSparseF* hsreduction = new THnSparseF("hadroncontamin", "hadroncontamin; pt[GeV/c]", 1, nBin);
  hsreduction->SetBinEdges(0, binEdges[0]);

  Double_t dataE[1];
  Double_t yval;
  for(int i=0; i<nBin[0]; i++){
    dataE[0]=hreduction->GetBinCenter(i+1);
    yval=hreduction->GetBinContent(i+1);
    hsreduction->Fill(dataE,yval);
  }

  Int_t* ibins;
  ibins = new Int_t[nBin[0] + 1];
  hsreduction->SetBinError(ibins,0);


  return hsreduction;
}

//_________________________________________________________________________
TObject* GetSpectrum(AliCFContainer *c, Int_t step) {
  AliCFDataGrid* data = new AliCFDataGrid("data","",*c,step);
  //data->SetMeasured(step);
  return data;
}
//___________________________________________________________________________
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
