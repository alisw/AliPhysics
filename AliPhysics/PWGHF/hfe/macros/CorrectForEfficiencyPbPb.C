void CorrectForEfficiencyPbPb(const char *filedata,const char *fileMC,const char *NameList,Bool_t smoothingon=kFALSE,Bool_t hadroncontaminationsubtracted=kFALSE,Bool_t unsetCorrelatedErrors=kTRUE);
TList *GetResults(const Char_t *testfile,const Char_t *plus="");
TList *GetQA(const Char_t *testfile,const Char_t *plus="");
TObject* GetSpectrum(AliCFContainer *c, Int_t step);
THnSparseF* GetHadronEffbyIPcut(TList *hfeqa);
void CorrectFromTheWidth(TH1D *h1);

void CorrectForEfficiencyPbPb(const char *filedata,const char *fileMC,const char *NameList,Bool_t smoothingon,Bool_t hadroncontaminationsubtracted,Bool_t unsetCorrelatedErrors) {
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Inputs of the macro:
  // filedata and fileMC: path of the output over data and MC
  // NameList: name of the list in this output ("HFE_Results" for example)
  // smoothingon: kTRUE means smoothing is used, should be kFALSE in general
  // hadroncontaminationsubtracted: kTRUE means hadron contamination is subtracted
  // unsetCorrelatedErrors: kTRUE means that the errors may be not properly calculated, nevertheless in 2D (centrality,pt) the proper calculation of the errors gave crazy results in the past, therefore kTRUE per default
  //
  // In the macro the centrality bins for the final corrected spectra are choosen
  // Look at: 
  // spectrum->SetNCentralityBinAtTheEnd(2);
  // spectrum->SetLowHighBoundaryCentralityBinAtTheEnd(0,2,0); // 0-10%
  // spectrum->SetLowHighBoundaryCentralityBinAtTheEnd(7,9,1); // 60-80%
  // Per definition:
  // spectrum->SetLowHighBoundaryCentralityBinAtTheEnd(0,1,0); // 0-5%
  // spectrum->SetLowHighBoundaryCentralityBinAtTheEnd(1,2,0); // 5-10%
  // spectrum->SetLowHighBoundaryCentralityBinAtTheEnd(2,3,0); // 10-20%
  // spectrum->SetLowHighBoundaryCentralityBinAtTheEnd(3,4,0); // 20-30%
  // spectrum->SetLowHighBoundaryCentralityBinAtTheEnd(4,5,0); // 30-40%
  // spectrum->SetLowHighBoundaryCentralityBinAtTheEnd(5,6,0); // 40-50%
  // spectrum->SetLowHighBoundaryCentralityBinAtTheEnd(6,7,0); // 50-60%
  // spectrum->SetLowHighBoundaryCentralityBinAtTheEnd(7,8,0); // 60-70%
  // spectrum->SetLowHighBoundaryCentralityBinAtTheEnd(8,9,0); // 70-80%
  //
  // Outputs:
  //  TFile*         finalSpectrum.root
  //
  //  Minimum-bias results not relevant for the analysis
  //
  //  TGraphErrors     UnfoldingCorrectedSpectrum;1
  //  TGraphErrors     AlltogetherSpectrum;1
  //  TH1D     RatioUnfoldingAlltogetherSpectrum;1
  //  THnSparseT<TArrayF>      UnfoldingCorrectedNotNormalizedSpectrum;1       step1 projection centralitypt
  //  AliCFDataGrid    AlltogetherCorrectedNotNormalizedSpectrum;1     dataGrid
  //
  // Results:
  //
  // Unfolded results: can be used for minimum-bias and enhanced sample
  //
  // TH1D     Unfolded_Notnormalized_centrality_bin_[0_2[;1    // Not normalized to the number of events
  // TGraphErrors     Unfolded_normalized_centrality_bin_[0_2[;1  // Normalized to the number of events ---> usual output
  //
  // Direct correction (ptesd/ptMC): can be used ONLY for minimum-bias sample
  //
  // TH1D     Dirrectcorrected_Notnormalized_centrality_bin_[0_2[;1   // Not normalized to the number of events
  // TGraphErrors     Dirrectedcorrected_normalized_centrality_bin_[0_2[;1 // Normalized to the number of events --> usual output
  //
  // Important:
  //
  // TPC pid efficiency is not included --> results have to be multiplied by 2
  //
  // For minimum-bias Dirrectcorrected_normalized_centrality_bin_[0_2[ and Unfolded_normalized_centrality_bin_[0_2[ should be similar
  // Not for enhance sample: Dirrectcorrected is wrong.
  //
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
  
  TList *resultsdata = GetResults(filedata,NameList);
  if(!resultsdata){
    printf("No output objects for data: Calculation will terminate here\n");
    return;
  }

  ///////////////////////////////////////////////////////////
  // Event container for the normalization to CINB1
  // Normalize to the number events after event selection
  ////////////////////////////////////////////////////////////
  
  AliCFContainer *containerdata = (AliCFContainer *) resultsdata->FindObject("eventContainer");
  if(!containerdata) {
    printf("No container \n");
    return;
  }
  AliCFDataGrid *dataGrid = (AliCFDataGrid *) GetSpectrum(containerdata,AliHFEcuts::kEventStepReconstructed);
  THnSparseF *eventcontainer = (THnSparseF*) dataGrid->GetGrid();
  TAxis *cenaxisb = eventcontainer->GetAxis(2);
  Int_t nbbin = cenaxisb->GetNbins();
  printf("There are %d centrality bins!!!!!!!!!!!\n",nbbin);
  Double_t *nbeventstotal = new Double_t[nbbin];
  
  for(Int_t binc = 0; binc < nbbin; binc++){
    cenaxisb->SetRange(binc+1,binc+1);
    TH1D *spectrum1Dc = (TH1D *) dataGrid->Project(0);
    nbeventstotal[binc] = spectrum1Dc->Integral();
  }

  Float_t *numberofentries = new Float_t[nbbin]; 
  Float_t numberofentriessum = 0.0;
  for(Int_t binc = 0; binc < nbbin; binc++) {
    numberofentries[binc] = nbeventstotal[binc];
    numberofentriessum += numberofentries[binc];
    printf("Number %f for centrality bin %d\n",numberofentries[binc],binc);
  }

  //////////////////////////////
  // Take more stuff
  ///////////////////////////////
 
  
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


  /////////////////////////////
  // Check number of events
  /////////////////////////////


  Int_t numberOfEventsafterallcuts = (Int_t) datahfecontainer->GetNumberOfEvents();
   
  AliHFEcontainer *mchfecontainer = (AliHFEcontainer *) resultsmc->FindObject("trackContainer");
  if(!mchfecontainer) {
    printf("No mc container \n");
    return;
  }
 
  
  //////////////
  // Correct
  /////////////

  AliHFEspectrum *spectrum = new AliHFEspectrum("HFEspectrum");
  spectrum->SetNbDimensions(1);
  // If you want to correct positive (0) or negative (1) separately
  //spectrum->SetChargeChoosen(0);
  spectrum->SetDebugLevel(1);
  // Give the number of events per centrality bins
  for(Int_t binc = 0; binc < nbbin; binc++) {
    spectrum->SetNumberOfEvents((Int_t)numberofentries[binc],binc);
  }
  // Calculation of the errors
  spectrum->SetUnSetCorrelatedErrors(unsetCorrelatedErrors);
  /////////////////////////////////////////////////////////////////////////
  // Number of centrality bin we want for the final corrected spectra
  ///////////////////////////////////////////////////////////////////////
  spectrum->SetNCentralityBinAtTheEnd(2);
  spectrum->SetLowHighBoundaryCentralityBinAtTheEnd(0,2,0); // 0-10%
  spectrum->SetLowHighBoundaryCentralityBinAtTheEnd(7,9,1); // 60-80%

  spectrum->SetDumpToFile(kTRUE);
  spectrum->SetSmoothing(smoothingon);
  // True step in MC (events in range +- 10 cm and no pile up)
  spectrum->SetMCTruthStep(AliHFEcuts::kStepMCGeneratedEventCut);
  // Step where we correct from MC (tracking + TOF)
  spectrum->SetMCEffStep(AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kStepHFEcutsTRD +1);
  // Step from data we correct for
  spectrum->SetStepToCorrect(AliHFEcuts::kStepHFEcutsTRD + 2);
  // PbPb flag
  spectrum->SetPbPbAnalysis(kTRUE);  

  // Give everything (data container, mc container and V0 data container)
  spectrum->Init(datahfecontainer,mchfecontainer);
    
  // kTRUE means subtract hadron background, kFALSE means do not subtract hadron background
  spectrum->Correct(hadroncontaminationsubtracted);



  

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
