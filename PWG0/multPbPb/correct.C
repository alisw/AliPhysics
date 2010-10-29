class AliAnalysisMultPbTrackHistoManager;

AliAnalysisMultPbTrackHistoManager * hManData = 0;
AliAnalysisMultPbTrackHistoManager * hManCorr = 0;
TH2D * hEvStatData = 0;
TH2D * hEvStatCorr = 0;

void LoadLibs();
void LoadData(TString dataFolder, TString correctionFolder);
void SetStyle();
Double_t CheckSecondaries();
void ShowAcceptanceInVzSlices() ;

void correct(TString dataFolder = "output/LHC09d_000104892_p4/", TString correctionFolder = "output/LHC10a8_104867/") {

  // Load stuff and set some styles
  LoadLibs();
  LoadData(dataFolder,correctionFolder);
  SetStyle();
  // ShowAcceptanceInVzSlices();
  // return;

  // TODO add some cool printout for cuts and centrality selection
  
  Double_t fractionWeak = CheckSecondaries();
  cout << "Rescaling weak correction: " << fractionWeak << endl;
  

  // Some shorthands
  // FIXME: Gen should be projected including overflow in z?
  TH1D * hDataPt   = (TH1D*) hManData->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRec,        -0.5,0.5,-10,10)->Clone("hDataPt");
  TH1D * hMCPtGen  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoGen,        -0.5,0.5,-10,10);
  TH1D * hMCPtRec  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRec,        -0.5,0.5,-10,10);
  TH1D * hMCPtPri  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim,    -0.5,0.5,-10,10);
  TH1D * hMCPtSeM  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat,  -0.5,0.5,-10,10);
  TH1D * hMCPtSeW  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak, -0.5,0.5,-10,10);
  TH1D * hMCPtFak  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecFake,    -0.5,0.5,-10,10);
 
  TCanvas * cdata = new TCanvas ("cData", "Data");  

  hDataPt->Draw();
  //  hMCPtRec->Draw("same");

  TCanvas * cMC = new TCanvas ("cMC", "Monte Carlo");
  hMCPtGen ->Draw();
  hMCPtRec ->Draw("same");
  hMCPtPri ->Draw("same");
  hMCPtSeM ->Draw("same");
  hMCPtSeW ->Draw("same");
  hMCPtFak ->Draw("same");

  cout << "Fake/All Rec  = " << hMCPtFak->Integral()/hMCPtRec->Integral()  << endl;
  cout << "SM/All   Rec  = " << hMCPtSeM->Integral()/hMCPtRec->Integral()  << endl;
  cout << "SW/All   Rec  = " << hMCPtSeW->Integral()/hMCPtRec->Integral()  << endl;


  // Compute efficiency and subtract secondaries and fakes after rescaling to data
  // PRIM_DATA = ALL_DATA - SEC_MC/ALL_MC*ALL_DATA - FAK_MC/ALL_MC*ALL_DATA
  // TRUE_DATA = PRIM_DATA * GEN_MC/PRIM_MC

  TH1D * hEffPt = (TH1D*) hMCPtPri->Clone("hEffPt");
  hEffPt->Divide(hMCPtPri,hMCPtGen,1,1,"B");

  TH1D * hCorSeM = (TH1D*) hMCPtSeM->Clone("hEffPt");
  hCorSeM->Divide(hMCPtSeM,hMCPtRec,1,1,"B");
  hCorSeM->Multiply(hDataPt);

  TH1D * hCorSeW = (TH1D*) hMCPtSeW->Clone("hEffPt");
  hCorSeW->Divide(hMCPtSeW,hMCPtRec,1,1,"B");
  hCorSeW->Scale(fractionWeak);// rescale weak correction
  hCorSeW->Multiply(hDataPt); 

  TH1D * hCorFak = (TH1D*) hMCPtFak->Clone("hEffPt");
  hCorFak->Divide(hMCPtFak,hMCPtRec,1,1,"B");
  hCorFak->Multiply(hDataPt);

  TH1D * hCorrected = (TH1D*) hDataPt->Clone("hCorrected");
  hCorrected->Add(hCorSeM,-1);
  hCorrected->Add(hCorSeW,-1);
  hCorrected->Add(hCorFak,-1);
  hCorrected->Divide(hEffPt);
  hCorrected->SetMarkerStyle(kOpenStar);

  cdata->cd();
  hDataPt->Draw();
  hCorrected->SetLineColor(kBlack);
  hCorrected->SetMarkerColor(kBlack);
  hCorrected->Draw("same");
  hMCPtGen->DrawClone("same");
  TF1 * f = GetLevy();
  hCorrected->Fit(f,"", "same");
  hCorrected->Fit(f,"IME", "same");
  cout << "dN/deta = " << f->Integral(0,100) << " +- " << f->IntegralError(0,100) << endl;
  cout << "Generated dN/deta (correction) = " << hMCPtGen->Integral() << endl;
  // FIXME: comment this out
  TH1D * hDataGen  = hManData->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoGen,        -0.5,0.5,-10,10);
  cout << "Generated dN/deta (data) =       " << hDataGen->Integral() << endl;
  hDataGen->Draw("same");
}

Double_t CheckSecondaries() {
  // Returns the fraction you need to rescale the secondaries from weak decays for.

  // Some shorthands
  TH1D * hDataDCA   = hManData->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec       );
  TH1D * hMCDCAGen  = hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoGen       );
  TH1D * hMCDCARec  = hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec       );
  TH1D * hMCDCAPri  = hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim   );
  TH1D * hMCDCASW  = hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak    );
  TH1D * hMCDCASM  = hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat    );
  TH1D * hMCDCAFak  = hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecFake   );
 

  TCanvas * cCumulative = new TCanvas("cDCAculumative","DCA cumulative distributions");
  GetCumulativeHisto(hMCDCAPri )->Draw();
  GetRatioIntegratedFractions(hMCDCASW, hMCDCARec  )->Draw("same");
  GetRatioIntegratedFractions(hMCDCASM, hMCDCARec  )->Draw("same");
  GetRatioIntegratedFractions(hMCDCAPri,hMCDCARec  )->Draw("same");


  TCanvas * c1 = new TCanvas("cDCAFit", "Fit to the DCA distributions");  
  c1->SetLogy();
  // Draw all
  //  hDataDCA->Draw();
  // //  hMCDCAGen ->Draw("same");
  //  hMCDCARec ->Draw("same");
  // hMCDCAPri ->Draw("same");
  // hMCDCASW ->Draw("same");
  // hMCDCASM ->Draw("same");
  // hMCDCAFak ->Draw("same");
  //  return;
  
  TH1D * hMCPrimSMFak = (TH1D*) hMCDCAPri->Clone("hMCPrimSMFak");
  hMCPrimSMFak->Add(hMCDCASM);
  hMCPrimSMFak->Add(hMCDCAFak);
  Int_t ncomp = 2;
  MyHistoFitter * fitter = new MyHistoFitter (ncomp);
  //  fitter->SetFCN(fcnhistoNoMCErr);
  fitter->SetHistoToFit(hDataDCA);
  fitter->SetMaxXValue(50);
  fitter->SetHistoComponent(0,hMCPrimSMFak);
  //fitter->SetHistoComponent(0,hMCPrimSM);
  fitter->SetHistoComponent(1,hMCDCASW);
  //  fitter->SetHistoComponent(2,hMCDCAFak);
  fitter->SetParName(0,"Primaries + Secondaries Material + Fakes");
  fitter->SetParName(1,"Secondaries (Weak)");
  //  fitter->SetParName(2,"Fakes");
  fitter->SetFactor(0,1);
  fitter->SetFactor(1,1);
  //  fitter->SetFactor(2,1);
  fitter->SetFactorLimits(0,0.5,3);
  fitter->SetFactorLimits(1,0.5,3);
  //  fitter->SetFactorLimits(2,0.5,3);
  //  fitter->SetFactorLimits(2,0.5,30);
  //  fitter->FixFactor(2);// FIXME: THIS IS BUGGY!
  //  fitter->SetFactorLimits(3,0.9,1.1);
  fitter->SetScaleHistos();
  fitter->Fit();
  fitter->PrintFitResults();
  fitter->GetHistoToFit()->Draw("");
  fitter->GetHistoFitted()->Draw("chist,same");
  for(Int_t icomp = 0; icomp < ncomp; icomp++){
    fitter->GetHistoComponent(icomp)->Draw("same");
  }
  
  
  return fitter->GetFactor(1)/fitter->GetFactor(0);

}

void LoadLibs() {

  gSystem->Load("libVMC");
  gSystem->Load("libTree");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libMinuit");
  gSystem->Load("libPWG2Spectra");
  gSystem->Load("libPWG0base"); 
   
  gROOT->ProcessLine(gSystem->ExpandPathName(".include $ALICE_ROOT/PWG0"));
  // Load helper classes
  // TODO: replace this by a list of TOBJStrings
  TString taskName("AliAnalysisTaskMultPbTracks.cxx+");
  TString histoManName("AliAnalysisMultPbTrackHistoManager.cxx+");
  TString centrName("AliAnalysisMultPbCentralitySelector.cxx+");
  TString listName("$ALICE_ROOT/PWG1/background/AliHistoListWrapper.cxx+");

  Bool_t debug=0;
  gROOT->LoadMacro(listName    +(debug?"+g":""));   
  gROOT->LoadMacro(histoManName+(debug?"+g":""));
  gROOT->LoadMacro(centrName   +(debug?"+g":""));   
  gROOT->LoadMacro(taskName    +(debug?"+g":""));   

  // Histo fitter
  gROOT->LoadMacro("/Users/mfloris/Work/ALICE/ANALYSIS/HistoFitter/fcn.cxx+g");
  gROOT->LoadMacro("/Users/mfloris/Work/ALICE/ANALYSIS/HistoFitter/MyHistoFitter.cxx+g");


}


void SetStyle() {

  hManData->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetLineColor(kBlack);    
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoGen)    ->SetLineColor(kRed  );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetLineColor(kRed  );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim)->SetLineColor(kGreen);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat) ->SetLineColor(kBlue );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak) ->SetLineColor(kBlue );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecFake)->SetLineColor(kCyan );

  hManData->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerColor(kBlack);    
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoGen)    ->SetMarkerColor(kRed  );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerColor(kRed  );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim)->SetMarkerColor(kGreen);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat) ->SetMarkerColor(kBlue );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak) ->SetMarkerColor(kBlue );
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecFake)->SetMarkerColor(kCyan );

  hManData->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerStyle(kFullCircle);    
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoGen)    ->SetMarkerStyle(kFullSquare);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim)->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat) ->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak) ->SetMarkerStyle(kOpenCircle);
  hManCorr->GetHistoPtEtaVz(AliAnalysisMultPbTrackHistoManager::kHistoRecFake)->SetMarkerStyle(kOpenSquare);

 hManData->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetLineColor(kBlack);    
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoGen)    ->SetLineColor(kRed  );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetLineColor(kRed  );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim)->SetLineColor(kGreen);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat) ->SetLineColor(kBlue );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak) ->SetLineColor(kBlue );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecFake)->SetLineColor(kCyan );

  hManData->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerColor(kBlack);    
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoGen)    ->SetMarkerColor(kRed  );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerColor(kRed  );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim)->SetMarkerColor(kGreen);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat) ->SetMarkerColor(kBlue );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak) ->SetMarkerColor(kBlue );
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecFake)->SetMarkerColor(kCyan );

  hManData->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerStyle(kFullCircle);    
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoGen)    ->SetMarkerStyle(kFullSquare);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRec)    ->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim)->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecMat) ->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecSecWeak) ->SetMarkerStyle(kOpenSquare);
  hManCorr->GetHistoDCA(AliAnalysisMultPbTrackHistoManager::kHistoRecFake)->SetMarkerStyle(kOpenSquare);


}

void LoadData(TString dataFolder, TString correctionFolder){

  // Get histo manager for data and MC + stat histos
  TFile * fData = new TFile(dataFolder+"multPbPbtracks.root");
  TFile * fCorr = new TFile(correctionFolder+"multPbPbtracks.root");
  TFile * fStatData = new TFile(dataFolder+"event_stat.root");
  TFile * fStatCorr = new TFile(correctionFolder+"event_stat.root");

  hManData = (AliAnalysisMultPbTrackHistoManager*) fData->Get("histoManager");
  hManCorr = (AliAnalysisMultPbTrackHistoManager*) fCorr->Get("histoManager");
  AliESDtrackCuts * cutsData = (AliESDtrackCuts*) fData->Get("AliESDtrackCuts");
  AliESDtrackCuts * cutsCorr = (AliESDtrackCuts*) fCorr->Get("AliESDtrackCuts");
  if (cutsData != cutsCorr) {
    cout << "ERROR: MC and data do not have the same cuts" << endl;
    // FIXME: exit here
  }
  cutsData->Print();
  hEvStatData = (TH2D*) fStatData->Get("fHistStatistics");
  hEvStatCorr = (TH2D*) fStatCorr->Get("fHistStatistics");

  // Normalize
  Int_t irowGoodTrigger = 1;
  if (hEvStatCorr && hEvStatData) {
    //  hManData->ScaleHistos(75351.36/1.015);// Nint for run 104892 estimated correcting for the trigger efficiency, multiplied for the physics selection efficiency which I'm not correcting for the time being
    hManData->ScaleHistos(hEvStatData->GetBinContent(AliPhysicsSelection::kStatAccepted,irowGoodTrigger));
    hManCorr->ScaleHistos(hEvStatCorr->GetBinContent(AliPhysicsSelection::kStatAccepted,irowGoodTrigger));
  } else {
    cout << "WARNING!!! ARBITRARY SCALING" << endl;
    hManData->ScaleHistos(1000);
    hManCorr->ScaleHistos(1000);    
  }
}

TF1 * GetHagedorn(Float_t norm=68, Float_t pt0=25, Float_t n=13) {

  TF1 * f =0;
  Double_t mass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  
  f=new TF1("fHagedorn",Form("(x/sqrt(x*x+%f*%f))*x*[0]*( 1 + x/[1] )^(-[2])",mass,mass),0,10);
  f->SetParameters(norm, pt0, n);
  f->SetParLimits(1, 0.01, 10);
  f->SetParNames("norm", "pt0", "n");
  f->SetLineWidth(1);
  return f;


}

TF1 * GetMTExp(Float_t norm=68, Float_t t=25) {

  TF1 * f =0;
  Double_t mass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  
  f=new TF1("fMTExp",Form("(x/sqrt(x*x+%f*%f))*x*[0]*exp(-sqrt(x*x+%f*%f)/[1])",mass,mass,mass,mass),0,10);
  f->SetParameters(norm, t);
  //  f->SetParLimits(1, 0.01);
  f->SetParNames("norm", "T");
  f->SetLineWidth(1);
  return f;


}

TF1 * GetLevy(Double_t temp=0.1, Double_t n=7, Double_t norm=10, const char * name="fLevy") {

  char formula[500];
  Double_t mass = TDatabasePDG::Instance()->GetParticle("pi+")->Mass();
  

  sprintf(formula,"(x/sqrt(x*x+[3]*[3]))*x*( [0]*([1]-1)*([1]-2)  )/( [1]*[2]*( [1]*[2]+[3]*([1]-2) )  ) * ( 1 + (sqrt([3]*[3]+x*x) -[3])/([1]*[2])  )^(-[1])");
  TF1* f=new TF1(name,formula,0,10);
  f->SetParameters(norm, n, temp,mass);
  f->SetParLimits(2, 0.01, 10);
  f->SetParNames("norm (dN/dy)", "n", "T", "mass");
  f->FixParameter(3,mass);
  f->SetLineWidth(1);
  return f;


}


TH1D * GetCumulativeHisto (TH1 * h) { 
  // Returns a cumulative histogram
  // FIXME: put this method in a tools class
  TH1D * hInt = h->Clone(TString(h->GetName())+"_Int");
  hInt->Reset();
  Float_t integral = h->Integral(-1,-1); // Consider under/over flow!
  Int_t nbin = h->GetNbinsX();
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    Double_t content = h->Integral(-1,ibin) / integral;
    hInt->SetBinContent(ibin, content);
  }
  return hInt;
}

TH1D * GetRatioIntegratedFractions (TH1 * hNum, TH1 * hDenum) { 
  // Returns the the ratio of integrated histograms 
  // FIXME: put this method in a tools class
  TH1D * hRatio = hNum->Clone(TString(hNum->GetName())+hDenum->GetName()+"_RatioIntegrated");
  hRatio->Reset();
  Int_t nbin = hNum->GetNbinsX();
  for(Int_t ibin = 1; ibin <= nbin; ibin++){
    Double_t content = hNum->Integral(-1,ibin) / hDenum->Integral(-1,ibin);// consider underflow
    hRatio->SetBinContent(ibin, content);
  }
  return hRatio;
}

void ShowAcceptanceInVzSlices() {
  TCanvas * cvz = new TCanvas("cvz","acc #times eff vs vz");
  for(Int_t ivz = -10; ivz < 10; ivz+=4){
    TH1D * hMCPtPri  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim ,   -0.5,0.5,ivz,ivz+4);
    TH1D * hMCPtGen  = hManCorr->GetHistoPt(AliAnalysisMultPbTrackHistoManager::kHistoGen,        -0.5,0.5,ivz,ivz+4);
    //    hEff= hMCPtGen;
    TH1D * hEff = hMCPtPri->Clone(Form("hEff_vz_%d_%d",ivz,ivz+4));
    hEff->Divide(hMCPtPri,hMCPtGen,1,1,"B");
    cout << "ivz " << ivz << endl;
    if(ivz < -9) {
      cout << "First" << endl;
      hEff->Draw();
      // hMCPtGen->Draw();
      // hMCPtPri->Draw("same");
    }
    else {
      cout << "Same" << endl;
      hEff->Draw("same");
      // hMCPtGen->Draw("");
      // hMCPtPri->Draw("same");
    }
    cvz->Update();
    //    cvz->WaitPrimitive();
  }
  
  //hManCorr->GetHistoVz(AliAnalysisMultPbTrackHistoManager::kHistoRecPrim )->Draw();
   hManCorr->GetHistoVz(AliAnalysisMultPbTrackHistoManager::kHistoRec )->Draw("");
  hManCorr->GetHistoVz(AliAnalysisMultPbTrackHistoManager::kHistoGen     )->Draw("same");      


}
