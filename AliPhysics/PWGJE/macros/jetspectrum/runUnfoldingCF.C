void Load();
TList *GetResults(Char_t *testfile,char* listName);
TObject* GetMeasuredSpectrum(TList *fContainer);
TObject* GetGeneratedSpectrum(TList *fContainer);
TObject* GetEfficiency(TList *fContainer,Int_t iCase = 44);
THnSparse* GetResponseMatrix(TList *fContainer);
THnSparse* CreateEmptyResponseMatrix();
Float_t GetNTrials(TList *fContainer);
THnSparse* CreateGuessed(const THnSparse* h);
void PrintErrors(THnSparse *h);

Int_t fDebug = 10;

int iRebin = 8; // 

// *************************************************************************************************************
// Macro for unfolding jet spectra:
// * First Unfolding is always the MC itself 
// * if bMCtest we take a second MC which is unfolded with the first MC otherwise real data is unfolded with first MC
// * bMC2 controls wether we correct to charged particle jets or not
// ############### CHANGELOG ################
// 23.11. added method to create an empty response matrix
 
void runUnfoldingCF( int gIterations = 10, bool bMC2 = true,bool bNoPrior = true,bool bMCtest = true, int iJetF = 0) {
  Load();

  const bool bPreSmooth = false;

  // can be different depeding on the prior
  Int_t fIterations = gIterations;
  Int_t fIterations2 = gIterations;
  Int_t fIterations3 = gIterations;

  // the first one is usually MC only
  TString cJetF;
  TString listName;
  TString testfile;
  if(bMCtest) testfile = "allpt_lhc10e14_100903_Chunk1.root"; // the higher stats chunck  
  else testfile = "allpt_lhc10e14_100903.root";  
  TList *results = 0;
  if(iJetF == 0){
    cJetF = "anti-k_{T}";
    listName = "spec2_jetsAOD_FASTJET04_jetsAODMC%s_FASTJET04_0000000000";  
  }
  else if(iJetF == 1){
    cJetF = "k_{T}";
    listName = "spec2_jetsAOD_FASTKT04_jetsAODMC%s_FASTKT04_0000000000"; 
  }
  else if(iJetF == 2){
    cJetF = "UA1";
    listName = "spec2_jets_jetsAODMC%s_UA104_0000000000";
  }
  TList *results =  GetResults(testfile.Data(),Form(listName.Data(),(bMC2?"2":"")));
  
  // here comes the real data
  TString cJetF2;
  TString listName2;
  TString testfile2;
  
  if(bMCtest)testfile2 = "allpt_lhc10e14_100903_Chunk0.root";  // the smalle stats
  else testfile2 = "/Users/kleinb/alice/jets/train/100910/PWG4_JetTasksOutput_Merge_bcd.root";
  if(iJetF == 0){
    cJetF = "anti-k_{T}";
    if(bMCtest)listName2 = "spec2_jetsAOD_FASTJET04_jetsAODMC%s_FASTJET04_0000000000";  
    else listName2 = "spec2_jetsAOD_FASTJET04__0000000000";  

  }
  else if(iJetF == 1){
    cJetF = "k_{T}";
    if(bMCtest)listName2 = "spec2_jetsAOD_FASTKT04_jetsAODMC%s_FASTKT04_0000000000"; 
    else listName2 = "spec2_jetsAOD_FASTKT04__0000000000"; 

  }
  else if(iJetF == 2){
    cJetF = "k_{T}";
    if(bMCtest)listName2 = "spec2_jets_jetsAODMC%s_UA104_0000000000";
    else listName2 = "spec2_jets__0000000000";

  }
  TList *results2 = GetResults(testfile2.Data(),Form(listName2.Data(),(bMC2?"2":"")));bool bScale = false;
  


  if(!results){
    Error("No output objects: Calculation will terminate here");
    return;
  }

  bool b1D = true;
  const Float_t kMaxPt = 140; //  only for drawin... 
  const Float_t kMaxPt2 = 140;


  char cString[] = Form("unfolded/101010_unfolding%s_MC%s%s%%s_%sIter%04d_Rebin%02d_jet%d",(bNoPrior?"_NoPrior":""),(bMC2?"2":""),(bPreSmooth?"_PreSmooth":""),(bMCtest?"MCtest_":""),gIterations,iRebin,iJetF);
  char cPrintMask[] = Form("%s.png",cString);
  
  TF3 *fSmooth = new TF3("fSmooth","[0]*pow(x,[1])+[2]*y+[3]*z",5.,150,-10,10,-10,10);
  fSmooth = 0;
  TF3 *fSmooth2 = 0;
  TF3 *fSmooth3 = 0;

  //  fSmooth = 0;
  AliLog::SetGlobalLogLevel(AliLog::kDebug);
  // get the essential
  if(fSmooth){
    fSmooth->SetParameters(1,-8,0,0);
    fSmooth2 = (TF3*)fSmooth->Clone("fSmooth2");
    fSmooth3 = (TF3*)fSmooth->Clone("fSmooth3");
  }

  


  THnSparseF *efficiency  = (THnSparseF*)  GetEfficiency(results,41);
  THnSparseF *efficiencyRecMatch  = (THnSparseF*)  GetEfficiency(results,131);
  THnSparseF *response    = (THnSparseF*)  GetResponseMatrix(results);



  THnSparseF *measuredIn2    = (THnSparseF*) GetMeasuredSpectrum(results2);
  measuredIn2->Multiply(efficiencyRecMatch);
  THnSparseF *measuredIn    = (THnSparseF*) GetMeasuredSpectrum(results);
  measuredIn->Multiply(efficiencyRecMatch);
  THnSparseF *generatedOut   = (THnSparseF*) GetGeneratedSpectrum(results);
  THnSparseF *generatedOut2   = (THnSparseF*) GetGeneratedSpectrum(results2);
 
  Float_t fTrials =  GetNTrials(results);
  Float_t fTrials2 = GetNTrials(results2);
  if(bScale){
    Printf("Trials %.3E %.3E",fTrials,fTrials2);
    measuredIn2->Scale(1./fTrials2);
    measuredIn->Scale(1./fTrials);
  }

  // set the content to zero above threshold
  /*
  const float lastPt = 250;
  Int_t ibxyz[3];
  for(int ibx = 1;ibx<=measuredIn2->GetAxis(0)->GetNbins();ibx++){
    float pt = measuredIn2->GetAxis(0)->GetBinCenter(ibx);

    if(pt<lastPt)continue;
    ibxyz[0] = ibx;
    for(int iby = 1;iby<=measuredIn2->GetAxis(1)->GetNbins();iby++){
      ibxyz[1] = iby;
      for(int ibz = 1;ibz<=measuredIn2->GetAxis(2)->GetNbins();ibz++){
	ibxyz[2] = ibz;
	measuredIn2->SetBinContent(ibxyz,0);
	measuredIn2->SetBinError(ibxyz,0);
      }
    }
  }
  */


  // for testing
  const int nDim = 3;
  const int dimrec[nDim] = {0,1,2}; 
  const int dimgen[nDim] = {3,4,5}; 
 

  /* 
  THnSparseF *generated = response->Projection(nDim,dimgen);
  THnSparseF *measured = response->Projection(nDim,dimrec);
  */

  // create a guessed "a priori" distribution using binning of MC
  THnSparse* guessed2 = CreateGuessed(generatedOut) ; // can at best take the measured?
  THnSparse* guessed = CreateGuessed(generatedOut) ; // can at best take the measured?
  Printf("%s:%d %d",(char*)__FILE__,__LINE__,guessed2->GetNbins()); 
  Printf("%s:%d %d",(char*)__FILE__,__LINE__,guessed->GetNbins()); 

  //---- Dbug show the errrors
  //  PrintErrors(guessed);
  //  PrintErrors(response);
  //  PrintErrors(efficiency);
  //   PrintErrors(measuredIn);
  //  PrintErrors(generatedOut);
 
  Bool_t bCorrelatedErrors = true;
  AliCFUnfolding unfolding("unfolding","title",3,response,efficiency,measuredIn,(bNoPrior?0:guessed));
  unfolding.SetMaxNumberOfIterations(fIterations); // regulate flutuations...
  unfolding.SetMaxChi2PerDOF(0);
  if(fSmooth)unfolding.UseSmoothing(fSmooth);
  if(bCorrelatedErrors){
    unfolding.SetUseCorrelatedErrors(kTRUE);
    unfolding.SetMaxConvergencePerDOF(0.01*0.01);
  }

  unfolding.Unfold();

  THnSparse* result = unfolding.GetUnfolded();
  THnSparse* estMeasured = unfolding.GetEstMeasured();


  //----
  AliCFUnfolding unfolding2("unfolding2","title",3,response,efficiency,measuredIn2,(bNoPrior?0:guessed));
  //  AliCFUnfolding unfolding2("unfolding2","title",3,response,efficiency,measuredIn2,0);
  unfolding2.SetMaxNumberOfIterations(fIterations2);
  unfolding2.SetMaxChi2PerDOF(0);
  // carefull with 0x0 pointer neighbouring bin smoothing is switched on...
  if(fSmooth2)unfolding2.UseSmoothing(fSmooth2);
  if(bCorrelatedErrors){
    unfolding2.SetUseCorrelatedErrors(kTRUE);
    unfolding2.SetMaxConvergencePerDOF(0.01*0.01);
  }
  unfolding2.Unfold();



  THnSparse* estMeasuredTmp = unfolding2.GetEstMeasured();
  THnSparse* result2 = 0;
  THnSparse* estMeasured2 = 0;
  if(bPreSmooth){
    // use the result of the previous as smoothed input...
    AliCFUnfolding unfolding3("unfolding3","title",3,response,efficiency,estMeasuredTmp,(bNoPrior?0:guessed));
    unfolding3.SetMaxNumberOfIterations(fIterations3);
    unfolding3.SetMaxChi2PerDOF(0);
    unfolding3.UseSmoothing(fSmooth3);
    unfolding3.Unfold();
    result2 = unfolding3.GetUnfolded();
    estMeasured2 = unfolding3.GetEstMeasured();
  }
  else{
    result2 = unfolding2.GetUnfolded();
    estMeasured2 = estMeasuredTmp;
  }



  TCanvas * canvas = new TCanvas("canvas","",1200,900);
  canvas->Divide(2,3);
  TCanvas * cPrint = new TCanvas("cPrint","");



  TFile *f1 = new TFile(Form("%s.root",Form(cString,"")),"RECREATE");
  TParameter<Long_t>* fIterationsPara =  new TParameter<Long_t> ("fIterations", 0);
  fIterationsPara->SetVal((Long_t)gIterations);
  fIterationsPara->Write();
 
  TDirectory *dMC = f1->mkdir("unfoldMC");
  TDirectory *dReal = f1->mkdir("unfoldReal");

  gROOT->cd();
  // color code black is unfolded
  // blue is measured 
  // red  is generated
  // kGreen is guessed



  if(b1D){

    TH1* h_gen = generatedOut->Projection(0);

    h_gen->SetMarkerColor(kRed);
    h_gen->SetMarkerStyle(kFullSquare);
    h_gen->SetName("generated");
    h_gen->SetTitle("generated");
    TH1* h_meas = measuredIn->Projection(0);
    h_meas->SetMarkerColor(kBlue);
    h_meas->SetMarkerStyle(kFullSquare);
    h_meas->SetName("measured");
    h_meas->SetTitle("measured");
    TH1* h_guess = guessed->Projection(0);
    h_guess->SetMarkerColor(kGreen);
    h_guess->SetMarkerStyle(kFullSquare);
    h_guess->SetName("guessed");
    h_guess->SetTitle("guesse");
    TH1* h_unf = result->Projection(0);
    h_unf->SetMarkerColor(kBlack);
    h_unf->SetMarkerStyle(kFullSquare);
    h_unf->SetName("unfolded");
    h_unf->SetTitle("unfolded");
    TH1* h_estmeas = estMeasured->Projection(0);
    h_estmeas->SetMarkerColor(kGray);
    h_estmeas->SetMarkerStyle(kFullSquare);
    h_estmeas->SetName("estmeas");
    h_estmeas->SetTitle("estmeas");

    TLegend *leg = new TLegend(0.6,0.6,0.85,0.8);
    leg->SetFillColor(0);
    leg->SetTextFont(gStyle->GetTextFont());
    leg->SetBorderSize(0);
    leg->AddEntry(h_meas,"measured","P");
    leg->AddEntry(h_gen,"generated","P");
    leg->AddEntry(h_unf,"unfolded","P");
    leg->AddEntry(h_estmeas,"R * u","P");



    canvas->cd(1)->SetLogy();
    h_gen->SetAxisRange(0,kMaxPt);
    h_meas->SetAxisRange(0,kMaxPt);
    h_meas->SetXTitle("p_{T} (GeV)");
    h_meas->SetYTitle("yield (arb. units)");
    h_meas->DrawCopy("P");
    h_gen->DrawCopy("Psame");
    h_unf->DrawCopy("Psame");
    h_estmeas->DrawCopy("Psame");
    leg->Draw("");
    cPrint->cd()->SetLogy();
    h_meas->DrawCopy("P");
    h_gen->DrawCopy("Psame");
    h_unf->DrawCopy("Psame");
    h_estmeas->DrawCopy("Psame");
    leg->Draw("");
    canvas->Update();
    cPrint->Update();
    cPrint->SaveAs(Form(cPrintMask,"spectrum_mc"));
    if(!gROOT->IsBatch()){
      if(getchar()=='q')return;
    }
    // the same for the real data
    TH1* h_gen2 = generatedOut2->Projection(0);

    h_gen2->SetMarkerColor(kRed);
    h_gen2->SetMarkerStyle(kFullCircle);
    h_gen2->SetName("generated2");
    h_gen2->SetTitle("generated2");
    TH1* h_meas2 = measuredIn2->Projection(0);
    h_meas2->SetMarkerColor(kBlue);
    h_meas2->SetMarkerStyle(kFullCircle);
    h_meas2->SetName("measured2");
    h_meas2->SetTitle("measured2");
    TH1* h_guess2 = guessed2->Projection(0);
    h_guess2->SetMarkerColor(kGreen);
    h_guess2->SetMarkerStyle(kFullCircle);
    h_guess2->SetName("guessed2");
    h_guess2->SetTitle("guesse2");
    TH1* h_unf2 = result2->Projection(0);
    PrintErrors(result2);
    h_unf2->SetMarkerColor(kBlack);
    h_unf2->SetMarkerStyle(kFullCircle);
    h_unf2->SetName("unfolded2");
    h_unf2->SetTitle("unfolded2");
    TH1* h_estmeas2 = estMeasured2->Projection(0);
    h_estmeas2->SetMarkerColor(kGray);
    h_estmeas2->SetMarkerStyle(kFullCircle);
    h_estmeas2->SetName("estmeas2");
    h_estmeas2->SetTitle("estmeas2");

    TLegend *leg2 = new TLegend(0.6,0.6,0.85,0.8);
    leg2->SetFillColor(0);
    leg2->SetTextFont(gStyle->GetTextFont());
    leg2->SetBorderSize(0);
    leg2->AddEntry(h_meas2,"measured","P");
    leg2->AddEntry(h_gen2,"generated","P");
    leg2->AddEntry(h_unf2,"unfolded","P");
    leg2->AddEntry(h_estmeas2,"R * u","P");
    
    Printf("Integral Measured: %E Unfolded %E R*U %E",h_meas2->Integral(),h_unf2->Integral(),h_estmeas2->Integral());



    canvas->cd(2)->SetLogy();
    //

    h_meas2->SetXTitle("p_{T} (GeV)");
    h_meas2->SetYTitle("yield (arb. units)");
    //    h_gen2->DrawCopy("P");
    h_unf2->SetAxisRange(0,kMaxPt2);
    h_meas2->SetAxisRange(0,kMaxPt2);
    h_meas2->DrawCopy("P");
    h_unf2->DrawCopy("Psame");
    h_estmeas2->DrawCopy("Psame");
    leg2->Draw();
    cPrint->cd()->SetLogy();
    h_unf2->SetAxisRange(0,kMaxPt2);
    h_meas2->SetAxisRange(0,kMaxPt2);
    h_meas2->DrawCopy("P");
    h_unf2->DrawCopy("Psame");
    h_estmeas2->DrawCopy("Psame");
    canvas->Update();
    cPrint->Update();
    leg2->Draw();
    cPrint->SaveAs(Form(cPrintMask,"spectrum_real"));
    if(!gROOT->IsBatch()){
      if(getchar()=='q')return;
    }



    // residuals

    TH1D *hResiduals = (TH1D*) h_meas->Clone("hResiduals");
    hResiduals->Reset();
    for(int ib = 1;ib <  hResiduals->GetNbinsX();ib++){
      float val1 =  h_meas->GetBinContent(ib);
      float err1 =  h_meas->GetBinError(ib);
      float val2 =  h_estmeas->GetBinContent(ib);
      if(err1>0){
	Float_t res = (val1-val2)/err1;
	Float_t res_err = (val1-val2)/err1*0.01; // error bars of 1%
	hResiduals->SetBinContent(ib,res);
	hResiduals->SetBinError(ib,0.01);
      }
      
    }
    hResiduals->SetXTitle("p_{T} (GeV)");
    hResiduals->SetYTitle("#frac{m - R * u}{e_{m}}");
    hResiduals->SetMaximum(4.5);
    hResiduals->SetMinimum(-4.5);


    canvas->cd(3);
    //    hRatioGenUnf->DrawCopy("p");
    hResiduals->SetAxisRange(0,kMaxPt);
    hResiduals->DrawCopy("P");
    cPrint->cd();
    gPad->SetLogy(0);
    hResiduals->DrawCopy("P");
    canvas->Update();
    cPrint->Update();
    cPrint->SaveAs(Form(cPrintMask,"residuals_mc"));
    if(!gROOT->IsBatch()){
      if(getchar()=='q')return;
    }

    TH1D *hResiduals2 =  (TH1D*)h_meas2->Clone("hResiduals2");
    hResiduals2->Reset();
    for(int ib = 1;ib <  hResiduals2->GetNbinsX();ib++){
      float val1 =  h_meas2->GetBinContent(ib);
      float err1 =  h_meas2->GetBinError(ib);
      float val2 =  h_estmeas2->GetBinContent(ib);
      if(err1>0){
	Float_t res = (val1-val2)/err1;
	Float_t res_err = (val1-val2)/err1*0.01; // error bars of 1%
	hResiduals2->SetBinContent(ib,res);
	hResiduals2->SetBinError(ib,0.01);
      }
      
    }
    hResiduals2->SetMaximum(4.5);
    hResiduals2->SetMinimum(-4.5);
    hResiduals2->SetYTitle("#frac{m - R * u}{e_{m}}");
    hResiduals2->SetAxisRange(0,kMaxPt2);
    cPrint->cd();
    gPad->SetLogy(0);
    hResiduals2->DrawCopy("P");
    canvas->Update();
    cPrint->Update();
    cPrint->SaveAs(Form(cPrintMask,"residuals_real"));
    if(!gROOT->IsBatch()){
      if(getchar()=='q')return;
    }

    canvas->Update();
    cPrint->Update();

    if(!gROOT->IsBatch()){
      if(getchar()=='q')return;
    }
    hResiduals2->DrawCopy("P");
    canvas->Update();
    cPrint->Update();
    if(!gROOT->IsBatch()){
      if(getchar()=='q')return;
    }

    TH1D *hRatioGenUnf2 = (TH1D*) h_gen2->Clone("hRatioGenUnf2");
    hRatioGenUnf2->Divide(h_unf2);
    canvas->cd(4);
    //    hRatioGenUnf2->DrawCopy("p");

    hResiduals2->DrawCopy("P");

    TH1D *hRatioUnfMeas = (TH1D*) h_unf->Clone("hRatioUnfMeas");
    hRatioUnfMeas->SetXTitle("p_{T} (GeV)");
    hRatioUnfMeas->SetYTitle("unfolded/meas");
    hRatioUnfMeas->Divide(h_meas);
    hRatioUnfMeas->SetMaximum(15);
    if(bMC2)hRatioUnfMeas->SetMaximum(5);
    hRatioUnfMeas->SetMinimum(0);
    TH1D *hRatioGenMeas =  (TH1D*)h_gen->Clone("hRatioGenMeas");
    hRatioGenMeas->Divide(h_meas);

    canvas->cd(5);
 
    hRatioUnfMeas->SetAxisRange(0.,kMaxPt);
    hRatioUnfMeas->DrawCopy("p");
    hRatioGenMeas->DrawCopy("psame");




    TH1D *hRatioUnfMeas2 =  (TH1D*) h_unf2->Clone("hRatioUnfMeas2");
    hRatioUnfMeas2->Divide(h_meas2);
    canvas->cd(6);
    hRatioUnfMeas2->SetXTitle("p_{T} (GeV)");
    hRatioUnfMeas2->SetYTitle("unfolded/meas");
    hRatioUnfMeas2->SetAxisRange(0.,kMaxPt2);
    hRatioUnfMeas2->SetMaximum(15);
    if(bMC2)hRatioUnfMeas2->SetMaximum(5);
    hRatioUnfMeas2->SetMinimum(0);
    hRatioUnfMeas2->DrawCopy("p");
    hRatioUnfMeas->DrawCopy("psame");

    // plot the effieciencies 
    cPrint->cd();
    TH1* hEffGen = efficiency->Projection(0);
    hEffGen->SetName("hEffGen");
    TH1* hEffRec = efficiencyRecMatch->Projection(0);
    hEffRec->SetName("hEffRec");

    hEffGen->SetAxisRange(0,kMaxPt2);
    hEffRec->SetAxisRange(0,kMaxPt2);

    hEffRec->SetXTitle("p_{T,gen/rec} (GeV/c)");
    hEffRec->SetYTitle("Efficiency");
    hEffGen->SetXTitle("p_{T,gen/rec} (GeV/c)");
    hEffGen->SetYTitle("Efficiency");

    hEffGen->SetMarkerColor(kBlue);
    hEffGen->SetMarkerStyle(kFullCircle);
    hEffRec->SetMarkerColor(kBlue);
    hEffRec->SetMarkerStyle(kOpenCircle);
    hEffGen->SetMaximum(1.2);
    hEffGen->SetMinimum(0.3);
    hEffGen->DrawCopy("p");
    hEffRec->DrawCopy("psame");
    cPrint->Update();
    cPrint->SaveAs(Form(cPrintMask,"effs"));

    h_unf->SetDirectory(dMC);
    h_meas->SetDirectory(dMC);
    hResiduals->SetDirectory(dMC);
    h_estmeas->SetDirectory(dMC);
    h_gen->SetDirectory(dMC);

 
    // put the used effs and response to the unfolded of the MC 

    hEffGen->SetDirectory(dMC);
    hEffRec->SetDirectory(dMC);
    //    response->SetDirectory(dMC);
    dMC->cd();
    guessed->Write();
    response->Write();
    measuredIn->Write();
    gROOT->cd();
   

    h_unf2->SetDirectory(dReal);
    h_meas2->SetDirectory(dReal);
    hResiduals2->SetDirectory(dReal);
    h_estmeas2->SetDirectory(dReal);
    h_gen2->SetDirectory(dReal); // will be empty if we use real data, but can check with different data set as well...

    dReal->cd();
    guessed2->Write();
    measuredIn2->Write();
    gROOT->cd();
    // store the parameters of the unfloding with tparamter...
    

  }
  canvas->cd();
  canvas->Modified();
  canvas->SaveAs(Form(cPrintMask,"all"));

  f1->Write();
  f1->Close();


  return;
  if(fDebug)Printf("Line:%d %s Entries %f",__LINE__,h_gen->GetName(),h_gen->GetEntries());
  //printf("c1\n");
  if(b1D)h_gen->DrawCopy("P");
  else h_gen->DrawCopy("lego2");

  


  if(b1D) canvas->cd(2)->SetLogy();
  else canvas->cd(2)->SetLogz();
  TH1* h_meas = 0;
  if(b1D) h_meas = measuredIn->Projection(0);
  else h_guessed = measuredIn->Projection(0,1);
  h_meas->SetName("measured");
  h_meas->SetTitle("measured");
  if(fDebug)Printf("Line:%d %s Entries %f",__LINE__,h_meas->GetName(),h_meas->GetEntries());

  //printf("c2\n");
  if(b1D) h_meas->Draw("hist");
  else h_meas->Draw("lego2");
  

  if(b1D) canvas->cd(3)->SetLogy();
  else canvas->cd(3)->SetLogz();

  TH1* h_guessed = 0;
  if(b1D) h_guessed = guessed2->Projection(0);
  else h_guessed = guessed2->Projection(0,1);

  h_guessed->SetTitle("a priori");
  h_guessed->SetName("apriori");
  if(fDebug)Printf("Line:%d %s Entries %f",__LINE__,h_guessed->GetName(),h_guessed->GetEntries());
  //printf("c3\n");
  if(b1D) h_guessed->Draw("hist");
  else h_guessed->Draw("lego2");

  canvas->cd(4);
  TH1* h_eff = 0;
  if(b1D) h_eff = efficiency->Projection(0);
  else h_eff = efficiency->Projection(0,1);

  h_eff->SetTitle("efficiency");
  h_eff->SetName("efficiency");
  if(fDebug)Printf("Line:%d %s Entries %f",__LINE__,h_eff->GetName(),h_eff->GetEntries());
  //printf("c4\n");
  if(b1D) h_eff->Draw("hist");
  else h_eff->Draw("lego2");


  if(b1D) canvas->cd(5)->SetLogy();
  else canvas->cd(5)->SetLogz();
  TH1* h_unf = 0;
  if(b1D) h_unf = result->Projection(0);
  else h_unf = result->Projection(0,1);
    
  h_unf->SetName("unfolded");
  h_unf->SetTitle("unfolded");
  //printf("c5\n");
  if(b1D) h_unf->Draw("hist");
  else h_unf->Draw("lego2");

  canvas->cd(6);
  TH1* ratio = 0;
  if(b1D){
    ratio = (TH1D*)h_unf->Clone("ratio");
    ratio->SetTitle("unfolded / generated");
    ratio->Divide(h_unf,h_gen,1,1);
  }
  else{
    ratio = (TH2D*)h_unf->Clone("ratio");
    ratio->SetTitle("unfolded / generated");
    ratio->Divide(h_unf,h_gen,1,1);
  }
//   ratio->Draw("cont4z");
//   ratio->Draw("surf2");
  //printf("c6\n");
  if(b1D)ratio->Draw("hist");
  else ratio->Draw("colz");

  return;

  canvas->cd(7);
  TH2* orig = unfolding.GetOriginalPrior()->Projection(1,0);
  orig->SetName("originalprior");
  orig->SetTitle("original prior");
  //printf("c7\n");
  orig->Draw("lego2");

  canvas->cd(8);
  THnSparseF* corrected = (THnSparseF*)measured->Clone("corrected");
  //corrected->ApplyEffCorrection(*efficiency);
  TH2D* corr = corrected->Projection(0,1);
  corr->SetTitle("simple correction");
  corr->SetName("simplecorrection");
  //printf("c8\n");
  corr->Draw("lego2");

  canvas->cd(9);
  TH2D* ratio2 = (TH2D*) corr->Clone();
  ratio2->Divide(corr,h_gen,1,1);
  ratio2->SetTitle("simple correction / generated");
  //printf("c9\n");
  ratio2->Draw("lego2");

  return;
}

// ====================================================================


void Load(){
  gSystem->Load("libTree");
  gSystem->Load("libPhysics");
  gSystem->Load("libHist");
  gSystem->Load("libVMC");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libAOD");
  gSystem->Load("libESD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libJETAN");
  gSystem->Load("libCORRFW");
  gSystem->Load("libPWG4JetTasks");
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
}

TList *GetResults( Char_t *testfile, Char_t *cName){
  //
  // read output
  //
  TFile *f = TFile::Open(testfile);
  if(!f || f->IsZombie()){
    Error("File not readable");
    return 0x0;
  }
  TDirectory *dir = dynamic_cast<TDirectory*>(f->Get(Form("PWG4_%s",cName)));
  if(!dir){
    Printf("%s:%d: Output directory %s not found",(char*)__FILE__,__LINE__,Form("PWG4_%s",cName));
    f->ls();
    f->Close(); delete f;
    return 0x0;
  } 
  TList *l = dynamic_cast<TList *>(dir->Get(Form("pwg4%s",cName)));
  if(!l){
    Printf("%s:%d: Output list %s not found",(char*)__FILE__,__LINE__,Form("pwg4%s",cName));
    dir->ls();
    f->Close(); delete f;
    return 0x0;
  } 
  //  TList *returnlist = dynamic_cast<TList *>(l->Clone());
  // f->Close(); delete f;
  return l;
}


TObject* GetMeasuredSpectrum(TList *fContainer) {
  THnSparseF* htmp = (THnSparseF*)fContainer->FindObject(Form("fhnJetContainer%d",AliAnalysisTaskJetSpectrum2::kMaxStep+AliAnalysisTaskJetSpectrum2::kStep1));
  THnSparseF* hm = (THnSparseF*)htmp->Clone("measured");
  hm->SetName("measured");
  hm->SetTitle("measured");
  if(iRebin==1)return hm;
  Int_t fRebin[3];
  fRebin[0] = iRebin;
  fRebin[1] = 1;
  fRebin[2] = 1;
  THnSparse *hm2 = hm->Rebin(fRebin); // this creates a new histo...
  //  hm2->Scale(iRebin);
  if(fDebug)Printf("Line: %d %s Entries %f",__LINE__,hm->GetName(),hm->GetEntries()); 
  return hm2;
}

Float_t GetNTrials(TList *fContainer){
  TH1F* htmp = (TH1F*)fContainer->FindObject("fh1Trials");
  return htmp->Integral();
}

TObject* GetGeneratedSpectrum(TList *fContainer) {
  // the generated jet spectrum within eta < 0.5
  THnSparseF* htmp = (THnSparseF*)fContainer->FindObject(Form("fhnJetContainer%d",AliAnalysisTaskJetSpectrum2::kStep1));
  THnSparseF* hg = (THnSparseF*)htmp->Clone("generated");
  hg->SetName("generated");
  hg->SetTitle("generated");

  if(iRebin==1)return hg;
  Int_t fRebin[3];
  fRebin[0] = iRebin;
  fRebin[1] = 1;
  fRebin[2] = 1; 
  THnSparse* hg2 = hg->Rebin(fRebin);
  // hg2->Scale(iRebin);
  if(fDebug)Printf("Line: %d %s Entries %f",__LINE__,hg->GetName(),hg->GetEntries()); 
  return hg2;
}

TObject* GetEfficiency(TList *fContainer,Int_t iCase) {

  static int count = 0;
  THnSparseF* hn1 = 0;
  THnSparseF* hn2 = 0;
  THnSparseF* he = 0;
  if(iCase==44){
    // Unitiy efficiency...
    hn1 = (THnSparseF*)fContainer->FindObject(Form("fhnJetContainer%d",AliAnalysisTaskJetSpectrum2::kStep4));
    hn2 = (THnSparseF*)fContainer->FindObject(Form("fhnJetContainer%d",AliAnalysisTaskJetSpectrum2::kStep4));
   }
  else if(iCase==14){
    // eficiency Generated with gen < 0.5 --> gen with Rec partner in 0.5
    hn1 = (THnSparseF*)fContainer->FindObject(Form("fhnJetContainer%d",AliAnalysisTaskJetSpectrum2::kStep1));
     hn2 = (THnSparseF*)fContainer->FindObject(Form("fhnJetContainer%d",AliAnalysisTaskJetSpectrum2::kStep4));
   }
  else if(iCase==41){
    // eficiency Generated with gen < 0.5 --> gen with Rec partner in 0.5
    hn1 = (THnSparseF*)fContainer->FindObject(Form("fhnJetContainer%d",AliAnalysisTaskJetSpectrum2::kStep4));
    hn2 = (THnSparseF*)fContainer->FindObject(Form("fhnJetContainer%d",AliAnalysisTaskJetSpectrum2::kStep1));
   }
  else if(iCase==131){
    // eficiency reconstruted with < 0.5 --> with partner
    hn1 = (THnSparseF*)fContainer->FindObject(Form("fhnJetContainer%d",AliAnalysisTaskJetSpectrum2::kMaxStep+AliAnalysisTaskJetSpectrum2::kStep3));
   hn2 = (THnSparseF*)fContainer->FindObject(Form("fhnJetContainer%d",AliAnalysisTaskJetSpectrum2::kMaxStep+AliAnalysisTaskJetSpectrum2::kStep1));
   }

  if(iRebin==1){
    he = (THnSparseF*)hn1->Clone("efficiency");
    he->Divide(hn1,hn2,1.,1.,"B");
    he->SetName(Form("efficiency_%d_%d",iCase,count));
    he->SetTitle(Form("efficiency_%d_%d",iCase,count));
    count++;
    return he;
  }
  Int_t fRebin[3];
  fRebin[0] = iRebin;
  fRebin[1] = 1;
  fRebin[2] = 1;  

  THnSparseF* hn1_2 = hn1->Rebin(fRebin);  
  THnSparseF* hn2_2 = hn2->Rebin(fRebin);  
  he = (THnSparseF*)hn1_2->Clone("efficiency");
  he->Divide(hn1_2,hn2_2,1.,1.,"B");
  he->SetName(Form("efficiency_%d_%d",iCase,count));
  he->SetTitle(Form("efficiency_%d_%d",iCase,count));
  return he;
}

THnSparse* CreateEmptyResponseMatrix(){
  // Creates an empty response matrix as in the Spectrum2 Task

  const Int_t kNvar   = 3 ; //number of variables on the grid:pt,eta, phi
  const Double_t kPtmin = 0.0, kPtmax = 320.; 
  const Double_t kEtamin = -3.0, kEtamax = 3.0;
  const Double_t kPhimin = 0., kPhimax = 2. * TMath::Pi();
  const Double_t kZmin = 0., kZmax = 1;

  //arrays for the number of bins in each dimension
  Int_t iBin[kNvar];
  iBin[0] = 320; //bins in pt
  iBin[1] =  1; //bins in eta 
  iBin[2] = 1; // bins in phi

  //arrays for lower bounds :
  Double_t* binEdges[kNvar];
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    binEdges[ivar] = new Double_t[iBin[ivar] + 1];

  //values for bin lower bounds
  //  for(Int_t i=0; i<=iBin[0]; i++) binEdges[0][i]=(Double_t)TMath::Power(10,TMath::Log10(kPtmin) + (TMath::Log10(kPtmax)-TMath::Log10(kPtmin))/iBin[0]*(Double_t)i);  
  for(Int_t i=0; i<=iBin[0]; i++) binEdges[0][i]=(Double_t)kPtmin  + (kPtmax-kPtmin)/(Double_t)iBin[0]*(Double_t)i;
  for(Int_t i=0; i<=iBin[1]; i++) binEdges[1][i]=(Double_t)kEtamin  + (kEtamax-kEtamin)/iBin[1]*(Double_t)i;
  for(Int_t i=0; i<=iBin[2]; i++) binEdges[2][i]=(Double_t)kPhimin  + (kPhimax-kPhimin)/iBin[2]*(Double_t)i;

  Int_t thnDim[2*kNvar];
  for (int k=0; k<kNvar; k++) {
    //first half  : reconstructed 
    //second half : MC
    thnDim[k]      = iBin[k];
    thnDim[k+kNvar] = iBin[k];
  }

  THnSparseF* fhnCorrelation = new THnSparseF("fhnCorrelation","THnSparse with correlations",2*kNvar,thnDim);
  for (int k=0; k<kNvar; k++) {
    fhnCorrelation->SetBinEdges(k,binEdges[k]);
    fhnCorrelation->SetBinEdges(k+kNvar,binEdges[k]);
  }
  fhnCorrelation->Sumw2();
  return  fhnCorrelation;

}

THnSparse* GetResponseMatrix(TList *fContainer) {
  THnSparse* h = (THnSparse*)fContainer->FindObject("fhnCorrelation");
  if(fDebug)Printf("Line: %d %s Entries %f",__LINE__,h->GetName(),h->GetEntries());
  if(iRebin==1)return h;
  Int_t fRebin[6];
  fRebin[0] = iRebin;
  fRebin[1] = 1;
  fRebin[2] = 1; 
  fRebin[3] = iRebin;
  fRebin[4] = 1;
  fRebin[5] = 1; 

  THnSparse* h2 = h->Rebin(fRebin);
  return h2;
}

THnSparse* CreateGuessed(const THnSparse* h) {

  // is already rebinned
  THnSparse* guessed = (THnSparse*) h->Clone();
  static int icount = 0;
  guessed->SetName(Form("%s_guess_%d",h->GetName(),icount++));
  

  return guessed;
 
}

void PrintErrors(THnSparse *h){

  for(Long_t i = 0; i< h->GetNbins();i++){
    Double_t val = h->GetBinContent(i);
    Double_t err = h->GetBinError  (i);
    if(val>0){
      Printf("%s %ld %1.3E +- %1.3E rel Error: %1.3E",h->GetName(),i,val,err,err/val);
    }
  }

}
