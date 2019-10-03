/************************************************************************

  V0 QA Post-processing Macro
  ---------------------------

 To be used in conjunction with AliAnalysisTaskQAV0.cxx, .h 
 Designed for visualization of the QA output from this task. 

 Creates a "summary pdf" file with several pages, each referring to 
 some aspect of the QA, including: 

  --- Event Counters
  --- Topological Variables
  --- Invariant Mass Peak Widths (resolution) 
  --- TPC dE/dx for daughters 

************************************************************************/


PostProcessQAV0(Bool_t lAttemptInvMassFit = kTRUE,
                Char_t *output            = "pdf"                          // "eps", "png" or "pdf"
               ){


  CustomGStyleSettings();  

  //==============================================================
  //Open Output File
  TFile* file = TFile::Open("AnalysisResults.root"), "READ");
  if (!file){
    cout<<"Output file not found!"<<endl;
    return;
  }
	file->cd("PWGLFQAV0_QA");
	TList* clist  = (TList*)file->FindObjectAny("clist");
  if (!clist){
    cout<<"File does not seem to hold QA list output!"<<endl;
    return;
  }  
  //==============================================================

  //==============================================================
  //Open Event Histogram: first canvas
  TH1D* fHistEvent = (TH1D*)clist->FindObject("fHistEvent");
  TCanvas *cHistEvent = new TCanvas("cHistEvent","",800,670); 
  cHistEvent->SetTopMargin(0.15);
  cHistEvent->SetGridx(); 
  cHistEvent->SetGridy();
  fHistEvent->Draw();

  fHistEvent->SetMarkerSize(1.35);
  fHistEvent->GetXaxis()->SetTitleOffset(1.2);
  fHistEvent->GetYaxis()->SetTitleOffset(1.2);
  fHistEvent->Draw("same text00");

  TLatex Tl;
  Tl.SetNDC();
  Tl.SetTextSize(0.05);
  Tl.DrawLatex(.35, .9277, "Event Counters")  ;  
  if      (output == "png") cHistEvent->SaveAs("LF_QAanalysis_V0_page1.png");  
  else if (output == "eps") cHistEvent->SaveAs("LF_QAanalysis_V0_page1.eps");
  else if (output == "pdf") cHistEvent->SaveAs("LF_QAanalysis_V0.pdf(");    

  //==============================================================

  //==============================================================
  //Invariant Mass Plots: Base, no dE/dx
  /*
  TH2D *f2dHistInvMassK0Short    = (TH2D*)clist->FindObject("f2dHistInvMassK0Short");
  TH2D *f2dHistInvMassLambda     = (TH2D*)clist->FindObject("f2dHistInvMassLambda");
  TH2D *f2dHistInvMassAntiLambda = (TH2D*)clist->FindObject("f2dHistInvMassAntiLambda");
    
  f2dHistInvMassK0Short     -> Scale ( 1. / ((Double_t)fHistEvent->GetBinContent(4) ) );
  f2dHistInvMassLambda      -> Scale ( 1. / ((Double_t)fHistEvent->GetBinContent(4) ) );
  f2dHistInvMassAntiLambda  -> Scale ( 1. / ((Double_t)fHistEvent->GetBinContent(4) ) );

  TCanvas *cInvMassK0Short = new TCanvas ( "cInvMassK0Short", "", 1200,800); 
  cInvMassK0Short->Divide(1,2); 
  cInvMassK0Short->cd(2)->Divide(3,1); 
  cInvMassK0Short->cd(1); 
  cInvMassK0Short->cd(1)->SetLogz();
  cInvMassK0Short->cd(1)->SetLeftMargin(0.065);
  cInvMassK0Short->cd(1)->SetTopMargin(0.13);
  cInvMassK0Short->cd(1)->SetBottomMargin(0.11);

  f2dHistInvMassK0Short->GetYaxis()->SetTitleSize(0.05);
  f2dHistInvMassK0Short->GetYaxis()->SetTitleOffset(0.6);
  f2dHistInvMassK0Short->GetXaxis()->SetTitleSize(0.05);
  f2dHistInvMassK0Short->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  f2dHistInvMassK0Short->GetZaxis()->SetTitle("Counts / Event");
  f2dHistInvMassK0Short->GetZaxis()->SetTitleOffset(0.7);
  f2dHistInvMassK0Short->GetZaxis()->SetRangeUser( f2dHistInvMassK0Short->GetMinimum(1e-10)*0.9, f2dHistInvMassK0Short->GetMaximum()*1.2 );
  f2dHistInvMassK0Short->Draw("colz");
  Tl.DrawLatex(.35, .9277, "Mass under K^{0}_{S} hypothesis, no TPC dE/dx")  ;   

  TH1D *fLowPtK0ShortSample = (TH1D*) f2dHistInvMassK0Short->ProjectionY( "fLowPtK0ShortSample", 
    f2dHistInvMassK0Short->GetXaxis()->FindBin(0.5001),
    f2dHistInvMassK0Short->GetXaxis()->FindBin(0.9999) );

  TH1D *fMidPtK0ShortSample = (TH1D*) f2dHistInvMassK0Short->ProjectionY( "fMidPtK0ShortSample", 
    f2dHistInvMassK0Short->GetXaxis()->FindBin(1.5001),
    f2dHistInvMassK0Short->GetXaxis()->FindBin(2.4999) );

  TH1D *fHiPtK0ShortSample = (TH1D*) f2dHistInvMassK0Short->ProjectionY( "fHiPtK0ShortSample", 
    f2dHistInvMassK0Short->GetXaxis()->FindBin(5.0001),
    f2dHistInvMassK0Short->GetXaxis()->FindBin(9.9999) );

  for(Int_t ic = 1; ic<4; ic++){ 
    cInvMassK0Short->cd(2)->cd(ic)->SetTopMargin(0.14);
    cInvMassK0Short->cd(2)->cd(ic)->SetLeftMargin(0.15);
    cInvMassK0Short->cd(2)->cd(ic)->SetBottomMargin(0.11);
    cInvMassK0Short->cd(2)->cd(ic)->SetGridx(); 
    cInvMassK0Short->cd(2)->cd(ic)->SetGridy();
  }

  cInvMassK0Short->cd(2)->cd(1); 
  fLowPtK0ShortSample->GetXaxis()->SetTitleOffset(1.1);
  fLowPtK0ShortSample->GetYaxis()->SetTitleOffset(2.33);
  fLowPtK0ShortSample->GetYaxis()->SetTitle("Counts/Event");
  fLowPtK0ShortSample->Draw(); 
  Tl.DrawLatex(.25, .9277, "0.5 < p_{T} (GeV/c^{2}) < 1.0");

  cInvMassK0Short->cd(2)->cd(2); 
  fMidPtK0ShortSample->GetXaxis()->SetTitleOffset(1.1);
  fMidPtK0ShortSample->GetYaxis()->SetTitleOffset(2.33);
  fMidPtK0ShortSample->GetYaxis()->SetTitle("Counts/Event");
  fMidPtK0ShortSample->Draw(); 
  Tl.DrawLatex(.25, .9277, "1.5 < p_{T} (GeV/c^{2}) < 2.5");

  cInvMassK0Short->cd(2)->cd(3); 
  fHiPtK0ShortSample->GetXaxis()->SetTitleOffset(1.1);
  fHiPtK0ShortSample->GetYaxis()->SetTitleOffset(2.33);
  fHiPtK0ShortSample->GetYaxis()->SetTitle("Counts/Event");
  fHiPtK0ShortSample->Draw(); 
  Tl.DrawLatex(.25, .9277, "5 < p_{T} (GeV/c^{2}) < 10");  
  if      (output == "png") cInvMassK0Short->SaveAs("LF_QAanalysis_V0_pageX.png");
  else if (output == "eps") cInvMassK0Short->SaveAs("LF_QAanalysis_V0_pageX.eps");
  else if (output == "pdf") cInvMassK0Short->SaveAs("LF_QAanalysis_V0.pdf");

  */
  //==============================================================

  //==============================================================
  //Invariant Mass Plots: WITH dE/dx
  TH2D *f2dHistInvMassWithdEdxK0Short    = (TH2D*)clist->FindObject("f2dHistInvMassWithdEdxK0Short");
  TH2D *f2dHistInvMassWithdEdxLambda     = (TH2D*)clist->FindObject("f2dHistInvMassWithdEdxLambda");
  TH2D *f2dHistInvMassWithdEdxAntiLambda = (TH2D*)clist->FindObject("f2dHistInvMassWithdEdxAntiLambda");
    
  f2dHistInvMassWithdEdxK0Short     -> Scale ( 1. / ((Double_t)fHistEvent->GetBinContent(4) ) );
  f2dHistInvMassWithdEdxLambda      -> Scale ( 1. / ((Double_t)fHistEvent->GetBinContent(4) ) );
  f2dHistInvMassWithdEdxAntiLambda  -> Scale ( 1. / ((Double_t)fHistEvent->GetBinContent(4) ) );

  TCanvas *cInvMassK0ShortWithdEdx = new TCanvas ( "cInvMassK0ShortWithdEdx", "", 1200,800); 
  cInvMassK0ShortWithdEdx->Divide(1,2); 
  cInvMassK0ShortWithdEdx->cd(2)->Divide(3,1); 
  cInvMassK0ShortWithdEdx->cd(1); 
  cInvMassK0ShortWithdEdx->cd(1)->SetLogz();
  cInvMassK0ShortWithdEdx->cd(1)->SetLeftMargin(0.065);
  cInvMassK0ShortWithdEdx->cd(1)->SetTopMargin(0.13);
  cInvMassK0ShortWithdEdx->cd(1)->SetBottomMargin(0.11);

  f2dHistInvMassWithdEdxK0Short->GetYaxis()->SetTitleSize(0.05);
  f2dHistInvMassWithdEdxK0Short->GetYaxis()->SetTitleOffset(0.6);
  f2dHistInvMassWithdEdxK0Short->GetXaxis()->SetTitleSize(0.05);
  f2dHistInvMassWithdEdxK0Short->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  f2dHistInvMassWithdEdxK0Short->GetZaxis()->SetTitle("Counts / Event");
  f2dHistInvMassWithdEdxK0Short->GetZaxis()->SetTitleOffset(0.7);
  f2dHistInvMassWithdEdxK0Short->GetZaxis()->SetRangeUser( f2dHistInvMassWithdEdxK0Short->GetMinimum(1e-10)*0.9, f2dHistInvMassWithdEdxK0Short->GetMaximum()*1.2 );
  f2dHistInvMassWithdEdxK0Short->Draw("colz");
  Tl.DrawLatex(.35, .9277, "Mass under K^{0}_{S} hypothesis, With TPC dE/dx")  ;   

  TH1D *fLowPtK0ShortSampleWithdEdx = (TH1D*) f2dHistInvMassWithdEdxK0Short->ProjectionY( "fLowPtK0ShortSampleWithdEdx", 
    f2dHistInvMassWithdEdxK0Short->GetXaxis()->FindBin(0.5001),
    f2dHistInvMassWithdEdxK0Short->GetXaxis()->FindBin(0.9999) );

  TH1D *fMidPtK0ShortSampleWithdEdx = (TH1D*) f2dHistInvMassWithdEdxK0Short->ProjectionY( "fMidPtK0ShortSampleWithdEdx", 
    f2dHistInvMassWithdEdxK0Short->GetXaxis()->FindBin(1.5001),
    f2dHistInvMassWithdEdxK0Short->GetXaxis()->FindBin(2.4999) );

  TH1D *fHiPtK0ShortSampleWithdEdx = (TH1D*) f2dHistInvMassWithdEdxK0Short->ProjectionY( "fHiPtK0ShortSampleWithdEdx", 
    f2dHistInvMassWithdEdxK0Short->GetXaxis()->FindBin(5.0001),
    f2dHistInvMassWithdEdxK0Short->GetXaxis()->FindBin(9.9999) );

  for(Int_t ic = 1; ic<4; ic++){ 
    cInvMassK0ShortWithdEdx->cd(2)->cd(ic)->SetTopMargin(0.14);
    cInvMassK0ShortWithdEdx->cd(2)->cd(ic)->SetLeftMargin(0.15);
    cInvMassK0ShortWithdEdx->cd(2)->cd(ic)->SetBottomMargin(0.11);
    cInvMassK0ShortWithdEdx->cd(2)->cd(ic)->SetGridx(); 
    cInvMassK0ShortWithdEdx->cd(2)->cd(ic)->SetGridy();
  }

  cInvMassK0ShortWithdEdx->cd(2)->cd(1); 
  fLowPtK0ShortSampleWithdEdx->GetXaxis()->SetTitleOffset(1.1);
  fLowPtK0ShortSampleWithdEdx->GetYaxis()->SetTitleOffset(2.33);
  fLowPtK0ShortSampleWithdEdx->GetYaxis()->SetTitle("Counts/Event");
  fLowPtK0ShortSampleWithdEdx->Draw(); 
  Tl.DrawLatex(.25, .9277, "0.5 < p_{T} (GeV/c^{2}) < 1.0");

  cInvMassK0ShortWithdEdx->cd(2)->cd(2); 
  fMidPtK0ShortSampleWithdEdx->GetXaxis()->SetTitleOffset(1.1);
  fMidPtK0ShortSampleWithdEdx->GetYaxis()->SetTitleOffset(2.33);
  fMidPtK0ShortSampleWithdEdx->GetYaxis()->SetTitle("Counts/Event");
  fMidPtK0ShortSampleWithdEdx->Draw(); 
  Tl.DrawLatex(.25, .9277, "1.5 < p_{T} (GeV/c^{2}) < 2.5");

  cInvMassK0ShortWithdEdx->cd(2)->cd(3); 
  fHiPtK0ShortSampleWithdEdx->GetXaxis()->SetTitleOffset(1.1);
  fHiPtK0ShortSampleWithdEdx->GetYaxis()->SetTitleOffset(2.33);
  fHiPtK0ShortSampleWithdEdx->GetYaxis()->SetTitle("Counts/Event");
  fHiPtK0ShortSampleWithdEdx->Draw(); 
  Tl.DrawLatex(.25, .9277, "5 < p_{T} (GeV/c^{2}) < 10");  

  TLatex Tli;
  Tli.SetNDC();
  Tli.SetTextSize(0.05);  

  if( lAttemptInvMassFit ){ 
    //Attempt rough signal extraction   
    TF1 *f1 = new TF1("f1","[0]+[1]*x+[2]*x*x+[3]*TMath::Gaus(x,[4],[5])",0.435, 0.565);
    //Reasonable first guess
    f1->SetParameter(0, fLowPtK0ShortSampleWithdEdx -> GetBinContent( fLowPtK0ShortSampleWithdEdx->FindBin(0.497-0.01) ) );
    f1->SetParameter(1, 0 );
    f1->SetParameter(2, 0 );
    f1->SetParameter(3, 0.5*fLowPtK0ShortSampleWithdEdx->GetMaximum() ); 
    f1->SetParameter(4, 0.497);
    f1->SetParameter(5, 0.0045);
    cout<<"Will now fit, please wait!"<<endl;
    fLowPtK0ShortSampleWithdEdx->Fit("f1","REM0");
    cInvMassK0ShortWithdEdx->cd(2)->cd(1); 
    f1->Draw("same"); 
    //Peak Position, Width: printout 
    fLowPtK0ShortSampleWithdEdx->GetYaxis()->SetRangeUser(0,fLowPtK0ShortSampleWithdEdx->GetMaximum() * 1.2);
    Tli.DrawLatex(0.2,0.8,Form("peak = %.4f #pm %.4f MeV/c^{2}",f1->GetParameter(4),f1->GetParameter(5)) ) ;
    
    //Attempt rough signal extraction   
    TF1 *f2 = new TF1("f2","[0]+[1]*x+[2]*x*x+[3]*TMath::Gaus(x,[4],[5])",0.375, 0.635);
    //Reasonable first guess
    f2->SetParameter(0, fMidPtK0ShortSampleWithdEdx -> GetBinContent( fMidPtK0ShortSampleWithdEdx->FindBin(0.497-0.01) ) );
    f2->SetParameter(1, 0 );
    f2->SetParameter(2, 0 );
    f2->SetParameter(3, 0.5*fMidPtK0ShortSampleWithdEdx->GetMaximum() ); 
    f2->SetParameter(4, 0.497);
    f2->SetParameter(5, 0.0045);
    cout<<"Will now fit, please wait!"<<endl;
    fMidPtK0ShortSampleWithdEdx->Fit("f2","REM0");
    cInvMassK0ShortWithdEdx->cd(2)->cd(2); 
    f2->Draw("same"); 
    //Peak Position, Width: printout 
    fMidPtK0ShortSampleWithdEdx->GetYaxis()->SetRangeUser(0,fMidPtK0ShortSampleWithdEdx->GetMaximum() * 1.2);
    Tli.DrawLatex(0.2,0.8,Form("peak = %.4f #pm %.4f MeV/c^{2}",f2->GetParameter(4),f2->GetParameter(5)) ) ;

    //Attempt rough signal extraction   
    TF1 *f3 = new TF1("f3","[0]+[1]*x+[2]*x*x+[3]*TMath::Gaus(x,[4],[5])",0.375, 0.665);
    //Reasonable first guess
    f3->SetParameter(0, fHiPtK0ShortSampleWithdEdx -> GetBinContent( fHiPtK0ShortSampleWithdEdx->FindBin(0.497-0.01) ) );
    f3->SetParameter(1, 0 );
    f3->SetParameter(2, 0 );
    f3->SetParameter(3, 0.5*fHiPtK0ShortSampleWithdEdx->GetMaximum() ); 
    f3->SetParameter(4, 0.497);
    f3->SetParameter(5, 0.0045);
    cout<<"Will now fit, please wait!"<<endl;
    fHiPtK0ShortSampleWithdEdx->Fit("f3","REM0");
    cInvMassK0ShortWithdEdx->cd(2)->cd(3); 
    f3->Draw("same"); 
    //Peak Position, Width: printout 
    fHiPtK0ShortSampleWithdEdx->GetYaxis()->SetRangeUser(0,fHiPtK0ShortSampleWithdEdx->GetMaximum() * 1.2);
    Tli.DrawLatex(0.2,0.8,Form("peak = %.4f #pm %.4f MeV/c^{2}",f3->GetParameter(4),f3->GetParameter(5)) ) ;
  }

  if      (output == "png") cInvMassK0ShortWithdEdx->SaveAs("LF_QAanalysis_V0_page2.png");      
  else if (output == "eps") cInvMassK0ShortWithdEdx->SaveAs("LF_QAanalysis_V0_page2.eps");
  else if (output == "pdf") cInvMassK0ShortWithdEdx->SaveAs("LF_QAanalysis_V0.pdf"); 

  //==============================================================  


  //==============================================================
  //Invariant Mass Plots: Base, no dE/dx
  /*
  TCanvas *cInvMassLambda = new TCanvas ( "cInvMassLambda", "", 1200,800); 
  cInvMassLambda->Divide(1,2); 
  cInvMassLambda->cd(2)->Divide(3,1); 
  cInvMassLambda->cd(1); 
  cInvMassLambda->cd(1)->SetLogz();
  cInvMassLambda->cd(1)->SetLeftMargin(0.065);
  cInvMassLambda->cd(1)->SetTopMargin(0.13);
  cInvMassLambda->cd(1)->SetBottomMargin(0.11);

  f2dHistInvMassLambda->GetYaxis()->SetTitleSize(0.05);
  f2dHistInvMassLambda->GetYaxis()->SetTitleOffset(0.6);
  f2dHistInvMassLambda->GetXaxis()->SetTitleSize(0.05);
  f2dHistInvMassLambda->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  f2dHistInvMassLambda->GetZaxis()->SetTitle("Counts / Event");
  f2dHistInvMassLambda->GetZaxis()->SetTitleOffset(0.7);
  f2dHistInvMassLambda->GetZaxis()->SetRangeUser( f2dHistInvMassLambda->GetMinimum(1e-10)*0.9, f2dHistInvMassLambda->GetMaximum()*1.2 );
  f2dHistInvMassLambda->Draw("colz");
  Tl.DrawLatex(.35, .9277, "Mass under #Lambda hypothesis, no TPC dE/dx")  ;   

  TH1D *fLowPtLambdaSample = (TH1D*) f2dHistInvMassLambda->ProjectionY( "fLowPtLambdaSample", 
    f2dHistInvMassLambda->GetXaxis()->FindBin(0.5001),
    f2dHistInvMassLambda->GetXaxis()->FindBin(0.9999) );

  TH1D *fMidPtLambdaSample = (TH1D*) f2dHistInvMassLambda->ProjectionY( "fMidPtLambdaSample", 
    f2dHistInvMassLambda->GetXaxis()->FindBin(1.5001),
    f2dHistInvMassLambda->GetXaxis()->FindBin(2.4999) );

  TH1D *fHiPtLambdaSample = (TH1D*) f2dHistInvMassLambda->ProjectionY( "fHiPtLambdaSample", 
    f2dHistInvMassLambda->GetXaxis()->FindBin(5.0001),
    f2dHistInvMassLambda->GetXaxis()->FindBin(9.9999) );

  for(Int_t ic = 1; ic<4; ic++){ 
    cInvMassLambda->cd(2)->cd(ic)->SetTopMargin(0.14);
    cInvMassLambda->cd(2)->cd(ic)->SetLeftMargin(0.15);
    cInvMassLambda->cd(2)->cd(ic)->SetBottomMargin(0.11);
    cInvMassLambda->cd(2)->cd(ic)->SetGridx(); 
    cInvMassLambda->cd(2)->cd(ic)->SetGridy();
  }

  cInvMassLambda->cd(2)->cd(1); 
  fLowPtLambdaSample->GetXaxis()->SetTitleOffset(1.1);
  fLowPtLambdaSample->GetYaxis()->SetTitleOffset(2.33);
  fLowPtLambdaSample->GetYaxis()->SetTitle("Counts/Event");
  fLowPtLambdaSample->Draw(); 
  Tl.DrawLatex(.25, .9277, "0.5 < p_{T} (GeV/c^{2}) < 1.0");

  cInvMassLambda->cd(2)->cd(2); 
  fMidPtLambdaSample->GetXaxis()->SetTitleOffset(1.1);
  fMidPtLambdaSample->GetYaxis()->SetTitleOffset(2.33);
  fMidPtLambdaSample->GetYaxis()->SetTitle("Counts/Event");
  fMidPtLambdaSample->Draw(); 
  Tl.DrawLatex(.25, .9277, "1.5 < p_{T} (GeV/c^{2}) < 2.5");

  cInvMassLambda->cd(2)->cd(3); 
  fHiPtLambdaSample->GetXaxis()->SetTitleOffset(1.1);
  fHiPtLambdaSample->GetYaxis()->SetTitleOffset(2.33);
  fHiPtLambdaSample->GetYaxis()->SetTitle("Counts/Event");
  fHiPtLambdaSample->Draw(); 
  Tl.DrawLatex(.25, .9277, "5 < p_{T} (GeV/c^{2}) < 10");  
  if      (output == "png") cInvMassLambda->SaveAs("LF_QAanalysis_V0_pageY.png");
  else if (output == "eps") cInvMassLambda->SaveAs("LF_QAanalysis_V0_pageY.eps");
  else if (output == "pdf") cInvMassLambda->SaveAs("LF_QAanalysis_V0.pdf");

  */
  //==============================================================

  //==============================================================
  //Invariant Mass Plots: WITH dE/dx
  TCanvas *cInvMassLambdaWithdEdx = new TCanvas ( "cInvMassLambdaWithdEdx", "", 1200,800); 
  cInvMassLambdaWithdEdx->Divide(1,2); 
  cInvMassLambdaWithdEdx->cd(2)->Divide(3,1); 
  cInvMassLambdaWithdEdx->cd(1); 
  cInvMassLambdaWithdEdx->cd(1)->SetLogz();
  cInvMassLambdaWithdEdx->cd(1)->SetLeftMargin(0.065);
  cInvMassLambdaWithdEdx->cd(1)->SetTopMargin(0.13);
  cInvMassLambdaWithdEdx->cd(1)->SetBottomMargin(0.11);

  f2dHistInvMassWithdEdxLambda->GetYaxis()->SetTitleSize(0.05);
  f2dHistInvMassWithdEdxLambda->GetYaxis()->SetTitleOffset(0.6);
  f2dHistInvMassWithdEdxLambda->GetXaxis()->SetTitleSize(0.05);
  f2dHistInvMassWithdEdxLambda->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  f2dHistInvMassWithdEdxLambda->GetZaxis()->SetTitle("Counts / Event");
  f2dHistInvMassWithdEdxLambda->GetZaxis()->SetTitleOffset(0.7);
  f2dHistInvMassWithdEdxLambda->GetZaxis()->SetRangeUser( f2dHistInvMassWithdEdxLambda->GetMinimum(1e-10)*0.9, f2dHistInvMassWithdEdxLambda->GetMaximum()*1.2 );
  f2dHistInvMassWithdEdxLambda->Draw("colz");
  Tl.DrawLatex(.35, .9277, "Mass under #Lambda hypothesis, With TPC dE/dx")  ;   

  TH1D *fLowPtLambdaSampleWithdEdx = (TH1D*) f2dHistInvMassWithdEdxLambda->ProjectionY( "fLowPtLambdaSampleWithdEdx", 
    f2dHistInvMassWithdEdxLambda->GetXaxis()->FindBin(0.5001),
    f2dHistInvMassWithdEdxLambda->GetXaxis()->FindBin(0.9999) );

  TH1D *fMidPtLambdaSampleWithdEdx = (TH1D*) f2dHistInvMassWithdEdxLambda->ProjectionY( "fMidPtLambdaSampleWithdEdx", 
    f2dHistInvMassWithdEdxLambda->GetXaxis()->FindBin(1.5001),
    f2dHistInvMassWithdEdxLambda->GetXaxis()->FindBin(2.4999) );

  TH1D *fHiPtLambdaSampleWithdEdx = (TH1D*) f2dHistInvMassWithdEdxLambda->ProjectionY( "fHiPtLambdaSampleWithdEdx", 
    f2dHistInvMassWithdEdxLambda->GetXaxis()->FindBin(5.0001),
    f2dHistInvMassWithdEdxLambda->GetXaxis()->FindBin(9.9999) );

  for(Int_t ic = 1; ic<4; ic++){ 
    cInvMassLambdaWithdEdx->cd(2)->cd(ic)->SetTopMargin(0.14);
    cInvMassLambdaWithdEdx->cd(2)->cd(ic)->SetLeftMargin(0.15);
    cInvMassLambdaWithdEdx->cd(2)->cd(ic)->SetBottomMargin(0.11);
    cInvMassLambdaWithdEdx->cd(2)->cd(ic)->SetGridx(); 
    cInvMassLambdaWithdEdx->cd(2)->cd(ic)->SetGridy();
  }

  cInvMassLambdaWithdEdx->cd(2)->cd(1); 
  fLowPtLambdaSampleWithdEdx->GetXaxis()->SetTitleOffset(1.1);
  fLowPtLambdaSampleWithdEdx->GetYaxis()->SetTitleOffset(2.33);
  fLowPtLambdaSampleWithdEdx->GetYaxis()->SetTitle("Counts/Event");
  fLowPtLambdaSampleWithdEdx->Draw(); 
  Tl.DrawLatex(.25, .9277, "0.5 < p_{T} (GeV/c^{2}) < 1.0");

  cInvMassLambdaWithdEdx->cd(2)->cd(2); 
  fMidPtLambdaSampleWithdEdx->GetXaxis()->SetTitleOffset(1.1);
  fMidPtLambdaSampleWithdEdx->GetYaxis()->SetTitleOffset(2.33);
  fMidPtLambdaSampleWithdEdx->GetYaxis()->SetTitle("Counts/Event");
  fMidPtLambdaSampleWithdEdx->Draw(); 
  Tl.DrawLatex(.25, .9277, "1.5 < p_{T} (GeV/c^{2}) < 2.5");

  cInvMassLambdaWithdEdx->cd(2)->cd(3); 
  fHiPtLambdaSampleWithdEdx->GetXaxis()->SetTitleOffset(1.1);
  fHiPtLambdaSampleWithdEdx->GetYaxis()->SetTitleOffset(2.33);
  fHiPtLambdaSampleWithdEdx->GetYaxis()->SetTitle("Counts/Event");
  fHiPtLambdaSampleWithdEdx->Draw(); 
  Tl.DrawLatex(.25, .9277, "5 < p_{T} (GeV/c^{2}) < 10"); 

  if( lAttemptInvMassFit ){ 
    //Attempt rough signal extraction   
    TF1 *fl1 = new TF1("fl1","[0]+[1]*x+[2]*x*x+[3]*TMath::Gaus(x,[4],[5])",1.095, 1.14);
    //Reasonable first guess
    fl1->SetParameter(0, fLowPtLambdaSampleWithdEdx -> GetBinContent( fLowPtLambdaSampleWithdEdx->FindBin(1.116+0.01) ) );
    fl1->SetParameter(1, 0 );
    fl1->SetParameter(2, 0 );
    fl1->SetParameter(3, 0.35*fLowPtLambdaSampleWithdEdx->GetMaximum() ); 
    fl1->SetParameter(4, 1.115683);
    fl1->SetParLimits(4,1.116-0.01,1.116+0.01);
    fl1->SetParameter(5, 0.002);
    fl1->SetParLimits(5,0.0005,0.01);
    cout<<"Will now fit, please wait!"<<endl;
    fLowPtLambdaSampleWithdEdx->Fit("fl1","REM0");
    cInvMassLambdaWithdEdx->cd(2)->cd(1); 
    fl1->Draw("same"); 
    //Peak Position, Width: printout 
    fLowPtLambdaSampleWithdEdx->GetYaxis()->SetRangeUser(0,fLowPtLambdaSampleWithdEdx->GetMaximum() * 1.2);
    Tli.DrawLatex(0.2,0.8,Form("peak = %.4f #pm %.4f MeV/c^{2}",fl1->GetParameter(4),fl1->GetParameter(5)) ) ;
    
    //Attempt rough signal extraction   
    TF1 *fl2 = new TF1("fl2","[0]+[1]*x+[2]*x*x+[3]*TMath::Gaus(x,[4],[5])",1.085, 1.15);
    //Reasonable first guess
    fl2->SetParameter(0, fMidPtLambdaSampleWithdEdx -> GetBinContent( fMidPtLambdaSampleWithdEdx->FindBin(1.116+0.01) ) );
    fl2->SetParameter(1, 0 );
    fl2->SetParameter(2, 0 );
    fl2->SetParameter(3, 0.6*fMidPtLambdaSampleWithdEdx->GetMaximum() ); 
    fl2->SetParLimits(4,1.116-0.01,1.116+0.01);
    fl2->SetParameter(4, 1.116);
    fl2->SetParameter(5, 0.0025);
    fl2->SetParLimits(5,0.0005,0.01);
    cout<<"Will now fit, please wait!"<<endl;
    fMidPtLambdaSampleWithdEdx->Fit("fl2","REM0");
    cInvMassLambdaWithdEdx->cd(2)->cd(2); 
    fl2->Draw("same"); 
    //Peak Position, Width: printout 
    fMidPtLambdaSampleWithdEdx->GetYaxis()->SetRangeUser(0,fMidPtLambdaSampleWithdEdx->GetMaximum() * 1.2);
    Tli.DrawLatex(0.2,0.8,Form("peak = %.4f #pm %.4f MeV/c^{2}",fl2->GetParameter(4),fl2->GetParameter(5)) ) ;

    //Attempt rough signal extraction   
    TF1 *fl3 = new TF1("fl3","[0]+[1]*x+[2]*x*x+[3]*TMath::Gaus(x,[4],[5])",1.095, 1.15);
    //Reasonable first guess
    fl3->SetParameter(0, fHiPtLambdaSampleWithdEdx -> GetBinContent( fHiPtLambdaSampleWithdEdx->FindBin(1.116+0.01) ) );
    fl3->SetParameter(1, 0 );
    fl3->SetParameter(2, 0 );
    fl3->SetParameter(3, 0.6*fHiPtLambdaSampleWithdEdx->GetMaximum() ); 
    fl3->SetParameter(4, 1.116);
    fl3->SetParLimits(4,1.116-0.005,1.116+0.005);
    fl3->SetParameter(5, 0.0035);
    fl3->SetParLimits(5,0.0005,0.01);
    cout<<"Will now fit, please wait!"<<endl;
    fHiPtLambdaSampleWithdEdx->Fit("fl3","REM0");
    cInvMassLambdaWithdEdx->cd(2)->cd(3); 
    fl3->Draw("same"); 
    //Peak Position, Width: printout 
    fHiPtLambdaSampleWithdEdx->GetYaxis()->SetRangeUser(0,fHiPtLambdaSampleWithdEdx->GetMaximum() * 1.2);
    Tli.DrawLatex(0.2,0.8,Form("peak = %.4f #pm %.4f MeV/c^{2}",fl3->GetParameter(4),fl3->GetParameter(5)) ) ;
  }

 
  if      (output == "png") cInvMassLambdaWithdEdx->SaveAs("LF_QAanalysis_V0_page3.png");
  else if (output == "eps") cInvMassLambdaWithdEdx->SaveAs("LF_QAanalysis_V0_page3.eps");
  else if (output == "pdf") cInvMassLambdaWithdEdx->SaveAs("LF_QAanalysis_V0.pdf");
  //==============================================================  

  //==============================================================
  //Invariant Mass Plots: Base, no dE/dx
  /*
  TCanvas *cInvMassAntiLambda = new TCanvas ( "cInvMassAntiLambda", "", 1200,800); 
  cInvMassAntiLambda->Divide(1,2); 
  cInvMassAntiLambda->cd(2)->Divide(3,1); 
  cInvMassAntiLambda->cd(1); 
  cInvMassAntiLambda->cd(1)->SetLogz();
  cInvMassAntiLambda->cd(1)->SetLeftMargin(0.065);
  cInvMassAntiLambda->cd(1)->SetTopMargin(0.13);
  cInvMassAntiLambda->cd(1)->SetBottomMargin(0.11);

  f2dHistInvMassAntiLambda->GetYaxis()->SetTitleSize(0.05);
  f2dHistInvMassAntiLambda->GetYaxis()->SetTitleOffset(0.6);
  f2dHistInvMassAntiLambda->GetXaxis()->SetTitleSize(0.05);
  f2dHistInvMassAntiLambda->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  f2dHistInvMassAntiLambda->GetZaxis()->SetTitle("Counts / Event");
  f2dHistInvMassAntiLambda->GetZaxis()->SetTitleOffset(0.7);
  f2dHistInvMassAntiLambda->GetZaxis()->SetRangeUser( f2dHistInvMassAntiLambda->GetMinimum(1e-10)*0.9, f2dHistInvMassAntiLambda->GetMaximum()*1.2 );
  f2dHistInvMassAntiLambda->Draw("colz");
  Tl.DrawLatex(.35, .9277, "Mass under #bar{#Lambda} hypothesis, no TPC dE/dx")  ;   

  TH1D *fLowPtAntiLambdaSample = (TH1D*) f2dHistInvMassAntiLambda->ProjectionY( "fLowPtAntiLambdaSample", 
    f2dHistInvMassAntiLambda->GetXaxis()->FindBin(0.5001),
    f2dHistInvMassAntiLambda->GetXaxis()->FindBin(0.9999) );

  TH1D *fMidPtAntiLambdaSample = (TH1D*) f2dHistInvMassAntiLambda->ProjectionY( "fMidPtAntiLambdaSample", 
    f2dHistInvMassAntiLambda->GetXaxis()->FindBin(1.5001),
    f2dHistInvMassAntiLambda->GetXaxis()->FindBin(2.4999) );

  TH1D *fHiPtAntiLambdaSample = (TH1D*) f2dHistInvMassAntiLambda->ProjectionY( "fHiPtAntiLambdaSample", 
    f2dHistInvMassAntiLambda->GetXaxis()->FindBin(5.0001),
    f2dHistInvMassAntiLambda->GetXaxis()->FindBin(9.9999) );

  for(Int_t ic = 1; ic<4; ic++){ 
    cInvMassAntiLambda->cd(2)->cd(ic)->SetTopMargin(0.14);
    cInvMassAntiLambda->cd(2)->cd(ic)->SetLeftMargin(0.15);
    cInvMassAntiLambda->cd(2)->cd(ic)->SetBottomMargin(0.11);
    cInvMassAntiLambda->cd(2)->cd(ic)->SetGridx(); 
    cInvMassAntiLambda->cd(2)->cd(ic)->SetGridy();
  }

  cInvMassAntiLambda->cd(2)->cd(1); 
  fLowPtAntiLambdaSample->GetXaxis()->SetTitleOffset(1.1);
  fLowPtAntiLambdaSample->GetYaxis()->SetTitleOffset(2.33);
  fLowPtAntiLambdaSample->GetYaxis()->SetTitle("Counts/Event");
  fLowPtAntiLambdaSample->Draw(); 
  Tl.DrawLatex(.25, .9277, "0.5 < p_{T} (GeV/c^{2}) < 1.0");

  cInvMassAntiLambda->cd(2)->cd(2); 
  fMidPtAntiLambdaSample->GetXaxis()->SetTitleOffset(1.1);
  fMidPtAntiLambdaSample->GetYaxis()->SetTitleOffset(2.33);
  fMidPtAntiLambdaSample->GetYaxis()->SetTitle("Counts/Event");
  fMidPtAntiLambdaSample->Draw(); 
  Tl.DrawLatex(.25, .9277, "1.5 < p_{T} (GeV/c^{2}) < 2.5");

  cInvMassAntiLambda->cd(2)->cd(3); 
  fHiPtAntiLambdaSample->GetXaxis()->SetTitleOffset(1.1);
  fHiPtAntiLambdaSample->GetYaxis()->SetTitleOffset(2.33);
  fHiPtAntiLambdaSample->GetYaxis()->SetTitle("Counts/Event");
  fHiPtAntiLambdaSample->Draw(); 
  Tl.DrawLatex(.25, .9277, "5 < p_{T} (GeV/c^{2}) < 10");  
  if      (output == "png") cInvMassAntiLambda->SaveAs("LF_QAanalysis_V0_pageZ.png");
  else if (output == "eps") cInvMassAntiLambda->SaveAs("LF_QAanalysis_V0_pageZ.eps");
  else if (output == "pdf") cInvMassAntiLambda->SaveAs("LF_QAanalysis_V0.pdf");

  */
  //==============================================================

  //==============================================================
  //Invariant Mass Plots: WITH dE/dx
  TCanvas *cInvMassAntiLambdaWithdEdx = new TCanvas ( "cInvMassAntiLambdaWithdEdx", "", 1200,800); 
  cInvMassAntiLambdaWithdEdx->Divide(1,2); 
  cInvMassAntiLambdaWithdEdx->cd(2)->Divide(3,1); 
  cInvMassAntiLambdaWithdEdx->cd(1); 
  cInvMassAntiLambdaWithdEdx->cd(1)->SetLogz();
  cInvMassAntiLambdaWithdEdx->cd(1)->SetLeftMargin(0.065);
  cInvMassAntiLambdaWithdEdx->cd(1)->SetTopMargin(0.13);
  cInvMassAntiLambdaWithdEdx->cd(1)->SetBottomMargin(0.11);

  f2dHistInvMassWithdEdxAntiLambda->GetYaxis()->SetTitleSize(0.05);
  f2dHistInvMassWithdEdxAntiLambda->GetYaxis()->SetTitleOffset(0.6);
  f2dHistInvMassWithdEdxAntiLambda->GetXaxis()->SetTitleSize(0.05);
  f2dHistInvMassWithdEdxAntiLambda->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  f2dHistInvMassWithdEdxAntiLambda->GetZaxis()->SetTitle("Counts / Event");
  f2dHistInvMassWithdEdxAntiLambda->GetZaxis()->SetTitleOffset(0.7);
  f2dHistInvMassWithdEdxAntiLambda->GetZaxis()->SetRangeUser( f2dHistInvMassWithdEdxAntiLambda->GetMinimum(1e-10)*0.9, f2dHistInvMassWithdEdxAntiLambda->GetMaximum()*1.2 );
  f2dHistInvMassWithdEdxAntiLambda->Draw("colz");
  Tl.DrawLatex(.35, .9277, "Mass under #bar{#Lambda} hypothesis, With TPC dE/dx")  ;   

  TH1D *fLowPtAntiLambdaSampleWithdEdx = (TH1D*) f2dHistInvMassWithdEdxAntiLambda->ProjectionY( "fLowPtAntiLambdaSampleWithdEdx", 
    f2dHistInvMassWithdEdxAntiLambda->GetXaxis()->FindBin(0.5001),
    f2dHistInvMassWithdEdxAntiLambda->GetXaxis()->FindBin(0.9999) );

  TH1D *fMidPtAntiLambdaSampleWithdEdx = (TH1D*) f2dHistInvMassWithdEdxAntiLambda->ProjectionY( "fMidPtAntiLambdaSampleWithdEdx", 
    f2dHistInvMassWithdEdxAntiLambda->GetXaxis()->FindBin(1.5001),
    f2dHistInvMassWithdEdxAntiLambda->GetXaxis()->FindBin(2.4999) );

  TH1D *fHiPtAntiLambdaSampleWithdEdx = (TH1D*) f2dHistInvMassWithdEdxAntiLambda->ProjectionY( "fHiPtAntiLambdaSampleWithdEdx", 
    f2dHistInvMassWithdEdxAntiLambda->GetXaxis()->FindBin(5.0001),
    f2dHistInvMassWithdEdxAntiLambda->GetXaxis()->FindBin(9.9999) );

  for(Int_t ic = 1; ic<4; ic++){ 
    cInvMassAntiLambdaWithdEdx->cd(2)->cd(ic)->SetTopMargin(0.14);
    cInvMassAntiLambdaWithdEdx->cd(2)->cd(ic)->SetLeftMargin(0.15);
    cInvMassAntiLambdaWithdEdx->cd(2)->cd(ic)->SetBottomMargin(0.11);
    cInvMassAntiLambdaWithdEdx->cd(2)->cd(ic)->SetGridx(); 
    cInvMassAntiLambdaWithdEdx->cd(2)->cd(ic)->SetGridy();
  }

  cInvMassAntiLambdaWithdEdx->cd(2)->cd(1); 
  fLowPtAntiLambdaSampleWithdEdx->GetXaxis()->SetTitleOffset(1.1);
  fLowPtAntiLambdaSampleWithdEdx->GetYaxis()->SetTitleOffset(2.33);
  fLowPtAntiLambdaSampleWithdEdx->GetYaxis()->SetTitle("Counts/Event");
  fLowPtAntiLambdaSampleWithdEdx->Draw(); 
  Tl.DrawLatex(.25, .9277, "0.5 < p_{T} (GeV/c^{2}) < 1.0");

  cInvMassAntiLambdaWithdEdx->cd(2)->cd(2); 
  fMidPtAntiLambdaSampleWithdEdx->GetXaxis()->SetTitleOffset(1.1);
  fMidPtAntiLambdaSampleWithdEdx->GetYaxis()->SetTitleOffset(2.33);
  fMidPtAntiLambdaSampleWithdEdx->GetYaxis()->SetTitle("Counts/Event");
  fMidPtAntiLambdaSampleWithdEdx->Draw(); 
  Tl.DrawLatex(.25, .9277, "1.5 < p_{T} (GeV/c^{2}) < 2.5");

  cInvMassAntiLambdaWithdEdx->cd(2)->cd(3); 
  fHiPtAntiLambdaSampleWithdEdx->GetXaxis()->SetTitleOffset(1.1);
  fHiPtAntiLambdaSampleWithdEdx->GetYaxis()->SetTitleOffset(2.33);
  fHiPtAntiLambdaSampleWithdEdx->GetYaxis()->SetTitle("Counts/Event");
  fHiPtAntiLambdaSampleWithdEdx->Draw(); 
  Tl.DrawLatex(.25, .9277, "5 < p_{T} (GeV/c^{2}) < 10"); 

  if( lAttemptInvMassFit ){ 
    //Attempt rough signal extraction   
    TF1 *fal1 = new TF1("fal1","[0]+[1]*x+[2]*x*x+[3]*TMath::Gaus(x,[4],[5])",1.095, 1.14);
    //Reasonable first guess
    fal1->SetParameter(0, fLowPtAntiLambdaSampleWithdEdx -> GetBinContent( fLowPtAntiLambdaSampleWithdEdx->FindBin(1.116+0.01) ) );
    fal1->SetParameter(1, 0 );
    fal1->SetParameter(2, 0 );
    fal1->SetParameter(3, 0.35*fLowPtAntiLambdaSampleWithdEdx->GetMaximum() ); 
    fal1->SetParameter(4, 1.115683);
    fal1->SetParLimits(4,1.116-0.01,1.116+0.01);
    fal1->SetParameter(5, 0.002);
    fal1->SetParLimits(5,0.0005,0.01);
    cout<<"Will now fit, please wait!"<<endl;
    fLowPtAntiLambdaSampleWithdEdx->Fit("fal1","REM0");
    cInvMassAntiLambdaWithdEdx->cd(2)->cd(1); 
    fal1->Draw("same"); 
    //Peak Position, Width: printout 
    fLowPtAntiLambdaSampleWithdEdx->GetYaxis()->SetRangeUser(0,fLowPtAntiLambdaSampleWithdEdx->GetMaximum() * 1.2);
    Tli.DrawLatex(0.2,0.8,Form("peak = %.4f #pm %.4f MeV/c^{2}",fal1->GetParameter(4),fal1->GetParameter(5)) ) ;
    
    //Attempt rough signal extraction   
    TF1 *fal2 = new TF1("fal2","[0]+[1]*x+[2]*x*x+[3]*TMath::Gaus(x,[4],[5])",1.085, 1.15);
    //Reasonable first guess
    fal2->SetParameter(0, fMidPtAntiLambdaSampleWithdEdx -> GetBinContent( fMidPtAntiLambdaSampleWithdEdx->FindBin(1.116+0.01) ) );
    fal2->SetParameter(1, 0 );
    fal2->SetParameter(2, 0 );
    fal2->SetParameter(3, 0.6*fMidPtAntiLambdaSampleWithdEdx->GetMaximum() ); 
    fal2->SetParLimits(4,1.116-0.01,1.116+0.01);
    fal2->SetParameter(4, 1.116);
    fal2->SetParameter(5, 0.0025);
    fal2->SetParLimits(5,0.0005,0.01);
    cout<<"Will now fit, please wait!"<<endl;
    fMidPtAntiLambdaSampleWithdEdx->Fit("fal2","REM0");
    cInvMassAntiLambdaWithdEdx->cd(2)->cd(2); 
    fal2->Draw("same"); 
    //Peak Position, Width: printout 
    fMidPtAntiLambdaSampleWithdEdx->GetYaxis()->SetRangeUser(0,fMidPtAntiLambdaSampleWithdEdx->GetMaximum() * 1.2);
    Tli.DrawLatex(0.2,0.8,Form("peak = %.4f #pm %.4f MeV/c^{2}",fal2->GetParameter(4),fal2->GetParameter(5)) ) ;

    //Attempt rough signal extraction   
    TF1 *fal3 = new TF1("fal3","[0]+[1]*x+[2]*x*x+[3]*TMath::Gaus(x,[4],[5])",1.095, 1.15);
    //Reasonable first guess
    fal3->SetParameter(0, fHiPtAntiLambdaSampleWithdEdx -> GetBinContent( fHiPtAntiLambdaSampleWithdEdx->FindBin(1.116+0.01) ) );
    fal3->SetParameter(1, 0 );
    fal3->SetParameter(2, 0 );
    fal3->SetParameter(3, 0.6*fHiPtAntiLambdaSampleWithdEdx->GetMaximum() ); 
    fal3->SetParameter(4, 1.116);
    fal3->SetParLimits(4,1.116-0.005,1.116+0.005);
    fal3->SetParameter(5, 0.0035);
    fal3->SetParLimits(5,0.0005,0.01);
    cout<<"Will now fit, please wait!"<<endl;
    fHiPtAntiLambdaSampleWithdEdx->Fit("fal3","REM0");
    cInvMassAntiLambdaWithdEdx->cd(2)->cd(3); 
    fal3->Draw("same"); 
    //Peak Position, Width: printout 
    fHiPtAntiLambdaSampleWithdEdx->GetYaxis()->SetRangeUser(0,fHiPtAntiLambdaSampleWithdEdx->GetMaximum() * 1.2);
    Tli.DrawLatex(0.2,0.8,Form("peak = %.4f #pm %.4f MeV/c^{2}",fal3->GetParameter(4),fal3->GetParameter(5)) ) ;
  }

 
  if      (output == "png") cInvMassAntiLambdaWithdEdx->SaveAs("LF_QAanalysis_V0_page4.png");
  else if (output == "eps") cInvMassAntiLambdaWithdEdx->SaveAs("LF_QAanalysis_V0_page4.eps");
  else if (output == "pdf") cInvMassAntiLambdaWithdEdx->SaveAs("LF_QAanalysis_V0.pdf");

  //==============================================================  

  //==============================================================  
  // Strict Lambda Analysis for dE/dx Calibration Check 
  TH2D *f2dHistdEdxSignalPionFromLambda    = (TH2D*)clist->FindObject("f2dHistdEdxSignalPionFromLambda");
  TH2D *f2dHistdEdxSignalProtonFromLambda  = (TH2D*)clist->FindObject("f2dHistdEdxSignalProtonFromLambda");
  TH2D *f2dHistResponsePionFromLambda      = (TH2D*)clist->FindObject("f2dHistResponsePionFromLambda");
  TH2D *f2dHistResponseProtonFromLambda    = (TH2D*)clist->FindObject("f2dHistResponseProtonFromLambda");

  f2dHistdEdxSignalPionFromLambda->Rebin2D(2,10);
  f2dHistdEdxSignalProtonFromLambda->Rebin2D(2,10);
    
  
  TH1D *fLowPtPionResponse = (TH1D*) f2dHistResponsePionFromLambda->ProjectionY( "fLowPtPionResponse", 
    f2dHistResponsePionFromLambda->GetXaxis()->FindBin(0.5001),
    f2dHistResponsePionFromLambda->GetXaxis()->FindBin(0.9999) );

  TH1D *fMidPtPionResponse = (TH1D*) f2dHistResponsePionFromLambda->ProjectionY( "fMidPtPionResponse", 
    f2dHistResponsePionFromLambda->GetXaxis()->FindBin(1.5001),
    f2dHistResponsePionFromLambda->GetXaxis()->FindBin(2.4999) );

  TH1D *fLowPtProtonResponse = (TH1D*) f2dHistResponseProtonFromLambda->ProjectionY( "fLowPtProtonResponse", 
    f2dHistResponseProtonFromLambda->GetXaxis()->FindBin(0.5001),
    f2dHistResponseProtonFromLambda->GetXaxis()->FindBin(0.9999) );

  TH1D *fMidPtProtonResponse = (TH1D*) f2dHistResponseProtonFromLambda->ProjectionY( "fMidPtProtonResponse", 
    f2dHistResponseProtonFromLambda->GetXaxis()->FindBin(1.5001),
    f2dHistResponseProtonFromLambda->GetXaxis()->FindBin(2.4999) );



  TCanvas *cdEdxPure = new TCanvas("cdEdxPure","",1500,800);
  cdEdxPure->Divide(4,2); 

  for(Int_t ic = 1; ic<9; ic++){ 
    cdEdxPure->SetLeftMargin(0.15);
    cdEdxPure->cd(ic)->SetLogz();
    //cdEdxPure->cd(ic)->SetTopMargin(0.133);
    if( ic%4 == 1 || ic%4 == 2){
      cdEdxPure->cd(ic)->SetRightMargin(0.133);
    }
    if( ic%4 != 1 && ic%4 != 2){
      cdEdxPure->cd(ic)->SetGridx();
      cdEdxPure->cd(ic)->SetGridy();
    }
  }



  cdEdxPure->cd(1); 
  f2dHistdEdxSignalPionFromLambda->Draw("colz"); 
  cdEdxPure->cd(2); 
  f2dHistResponsePionFromLambda->Draw("colz"); 
  cdEdxPure->cd(3); 
  fLowPtPionResponse->Draw();
  cdEdxPure->cd(4); 
  fMidPtPionResponse->Draw();

  cdEdxPure->cd(5); 
  f2dHistdEdxSignalProtonFromLambda->Draw("colz"); 
  cdEdxPure->cd(6); 
  f2dHistResponseProtonFromLambda->Draw("colz"); 
  cdEdxPure->cd(7); 
  fLowPtProtonResponse->Draw();
  cdEdxPure->cd(8); 
  fMidPtProtonResponse->Draw();



  //Write explanations on canvases
  cdEdxPure->cd(1); Tl.DrawLatex(.25, .9277, "#pi^{-} from #Lambda: TPC Signal");
  cdEdxPure->cd(2); Tl.DrawLatex(.15, .9277, "#pi^{-} from #Lambda: AliPIDResponse Value");
  cdEdxPure->cd(3); Tl.DrawLatex(.21, .9277, "#pi^{-} N#sigma, 0.5 < p_{T} (GeV/c) < 1.0");
  cdEdxPure->cd(4); Tl.DrawLatex(.21, .9277, "#pi^{-} N#sigma, 1.5 < p_{T} (GeV/c) < 2.5");

  cdEdxPure->cd(5); Tl.DrawLatex(.25, .9277, "p from #Lambda: TPC Signal");
  cdEdxPure->cd(6); Tl.DrawLatex(.15, .9277, "p from #Lambda: AliPIDResponse Value");
  cdEdxPure->cd(7); Tl.DrawLatex(.21, .9277, "p N#sigma, 0.5 < p_{T} (GeV/c) < 1.0");
  cdEdxPure->cd(8); Tl.DrawLatex(.21, .9277, "p N#sigma, 1.5 < p_{T} (GeV/c) < 2.5");

  Double_t lLowPtPeakPion         = fLowPtPionResponse->GetBinCenter( fLowPtPionResponse->GetMaximumBin() );
  Double_t lMidPtPeakPion         = fMidPtPionResponse->GetBinCenter( fMidPtPionResponse->GetMaximumBin() );
  Double_t lLowPtPeakProton       = fLowPtProtonResponse->GetBinCenter( fLowPtProtonResponse->GetMaximumBin() );
  Double_t lMidPtPeakProton       = fMidPtProtonResponse->GetBinCenter( fMidPtProtonResponse->GetMaximumBin() );

  //List Maximal Values
  cout<<"Maximal Value for pion from Lambda at low pt...............: " <<lLowPtPeakPion<<endl;
  cout<<"Maximal Value for pion from Lambda at mid pt...............: " <<lMidPtPeakPion<<endl;
  cout<<"Maximal Value for proton from Lambda at low pt.............: " <<lLowPtPeakProton<<endl;
  cout<<"Maximal Value for proton from Lambda at mid pt.............: " <<lMidPtPeakProton<<endl;

  if( TMath::Abs( lLowPtPeakPion ) > 0.3)   cout<<"*** WARNING: Check Low Pt pion PID Response! ***"<<endl;
  if( TMath::Abs( lMidPtPeakPion ) > 0.3)   cout<<"*** WARNING: Check Mid Pt pion PID Response! ***"<<endl;
  if( TMath::Abs( lLowPtPeakProton ) > 0.3) cout<<"*** WARNING: Check Low Pt proton PID Response! ***"<<endl;
  if( TMath::Abs( lMidPtPeakProton ) > 0.3) cout<<"*** WARNING: Check Mid Pt proton PID Response! ***"<<endl;

  TLatex Tlq;
  Tlq.SetNDC();
  Tlq.SetTextSize(0.06);

  //Draw Arrows to be sure!
  Double_t lFractionHeight = 0.33;
  cdEdxPure->cd(3);  
  TArrow *ar1 = new TArrow(lLowPtPeakPion,lFractionHeight * fLowPtPionResponse->GetMaximum(),lLowPtPeakPion,0.0, 0.008 ,"|>");
  ar1->SetLineWidth(2);
  ar1->Draw();
  if( TMath::Abs( lLowPtPeakPion ) < 0.3) {
    Tlq.SetTextColor(8); 
    Tlq.DrawLatex ( 0.15,0.8, "OK!");
  } else { 
    Tlq.SetTextColor(kRed); 
    Tlq.DrawLatex ( 0.15,0.8, "Prob.! Check dEd/x!");
  }

  cdEdxPure->cd(4);  
  TArrow *ar2 = new TArrow(lMidPtPeakPion,lFractionHeight * fMidPtPionResponse->GetMaximum(),lLowPtPeakPion,0.0, 0.008 ,"|>");
  ar2->SetLineWidth(2);
  ar2->Draw();
  if( TMath::Abs( lMidPtPeakPion ) < 0.3) {
    Tlq.SetTextColor(8); 
    Tlq.DrawLatex ( 0.15,0.8, "OK!");
  } else { 
    Tlq.SetTextColor(kRed); 
    Tlq.DrawLatex ( 0.15,0.8, "Prob.! Check dEd/x!");
  }
  
  cdEdxPure->cd(7);  
  TArrow *ar3 = new TArrow(lLowPtPeakProton,lFractionHeight * fLowPtProtonResponse->GetMaximum(),lLowPtPeakPion,0.0, 0.008 ,"|>");
  ar3->SetLineWidth(2);
  ar3->Draw();
  if( TMath::Abs( lLowPtPeakProton ) < 0.3) {
    Tlq.SetTextColor(8); 
    Tlq.DrawLatex ( 0.15,0.8, "OK!");
  } else { 
    Tlq.SetTextColor(kRed); 
    Tlq.DrawLatex ( 0.15,0.8, "Prob.! Check dEd/x!");
  }

  cdEdxPure->cd(8);  
  TArrow *ar4 = new TArrow(lMidPtPeakProton,lFractionHeight * fMidPtProtonResponse->GetMaximum(),lLowPtPeakPion,0.0, 0.008 ,"|>");
  ar4->SetLineWidth(2);
  ar4->Draw();
  if( TMath::Abs( lMidPtPeakProton ) < 0.3) {
    Tlq.SetTextColor(8); 
    Tlq.DrawLatex ( 0.15,0.8, "OK!");
  } else { 
    Tlq.SetTextColor(kRed); 
    Tlq.DrawLatex ( 0.15,0.8, "Prob.! Check dEd/x!");
  }

  if      (output == "png") cdEdxPure->SaveAs("LF_QAanalysis_V0_page5.png");
  else if (output == "eps") cdEdxPure->SaveAs("LF_QAanalysis_V0_page5.eps");
  else if (output == "pdf") cdEdxPure->SaveAs("LF_QAanalysis_V0.pdf");


  //==============================================================

  //==============================================================    
  //Topological variable QA
  // FIXME: This is still only a rough first version. 
  // Adjustments will be required for easy / long-term operation.
  TH1D *fHistTopDCAPosToPV      = (TH1D*)clist->FindObject("fHistSelectedTopDCAPosToPV");
  TH1D *fHistTopDCANegToPV      = (TH1D*)clist->FindObject("fHistSelectedTopDCANegToPV");
  TH1D *fHistTopDCAV0Daughters  = (TH1D*)clist->FindObject("fHistSelectedTopDCAV0Daughters");
  TH1D *fHistTopCosinePA        = (TH1D*)clist->FindObject("fHistSelectedTopCosinePA");
  TH1D *fHistTopV0Radius        = (TH1D*)clist->FindObject("fHistSelectedTopV0Radius");

  //Zoom in on selection in Cosine of pointing angle...
  Int_t iLowBin=-1; 
  Int_t iLowBin2 = -1; 
  Int_t iLowBin3 = -1; 

  //Normalize to per-event 
  fHistTopDCAPosToPV      -> Scale ( 1. / ((Double_t)fHistEvent->GetBinContent(4) ) );
  fHistTopDCANegToPV      -> Scale ( 1. / ((Double_t)fHistEvent->GetBinContent(4) ) );
  fHistTopDCAV0Daughters  -> Scale ( 1. / ((Double_t)fHistEvent->GetBinContent(4) ) );
  fHistTopCosinePA        -> Scale ( 1. / ((Double_t)fHistEvent->GetBinContent(4) ) );
  fHistTopV0Radius        -> Scale ( 1. / ((Double_t)fHistEvent->GetBinContent(4) ) );

  fHistTopDCAPosToPV      -> GetYaxis() -> SetTitle("Counts / Event");
  fHistTopDCANegToPV      -> GetYaxis() -> SetTitle("Counts / Event");
  fHistTopDCAV0Daughters  -> GetYaxis() -> SetTitle("Counts / Event");
  fHistTopCosinePA        -> GetYaxis() -> SetTitle("Counts / Event");
  fHistTopV0Radius        -> GetYaxis() -> SetTitle("Counts / Event");

  fHistTopDCAPosToPV      -> GetYaxis() -> SetTitleSize(0.05);
  fHistTopDCANegToPV      -> GetYaxis() -> SetTitleSize(0.05);
  fHistTopDCAV0Daughters  -> GetYaxis() -> SetTitleSize(0.05);
  fHistTopCosinePA        -> GetYaxis() -> SetTitleSize(0.05);
  fHistTopV0Radius        -> GetYaxis() -> SetTitleSize(0.05);

  fHistTopDCAPosToPV      -> GetXaxis() -> SetTitleSize(0.05);
  fHistTopDCANegToPV      -> GetXaxis() -> SetTitleSize(0.05);
  fHistTopDCAV0Daughters  -> GetXaxis() -> SetTitleSize(0.05);
  fHistTopCosinePA        -> GetXaxis() -> SetTitleSize(0.05);
  fHistTopV0Radius        -> GetXaxis() -> SetTitleSize(0.05);

  Double_t lMinimumCosPADraw = 1.-1.25*(1.-GetXForMinValue ( fHistTopCosinePA ) ); //assumes monotonically increasing dist
  cout<<"Function test: "<< lMinimumCosPADraw <<endl;
  fHistTopCosinePA->GetXaxis()->SetRangeUser(lMinimumCosPADraw,1.001);  

  Double_t lmin[5]; 
  Double_t lminPrecision[5]; 
  lmin[3] = GetXForMinValue ( fHistTopCosinePA );
  lmin[4] = GetXForMinValue ( fHistTopV0Radius );
  lmin[0] = GetXForMaxValue ( fHistTopDCAPosToPV );
  lmin[1] = GetXForMaxValue ( fHistTopDCANegToPV );
  lmin[2] = GetXForMinValue ( fHistTopDCAV0Daughters );

  lminPrecision[3] = fHistTopCosinePA        -> GetBinWidth(1);
  lminPrecision[4] = fHistTopV0Radius        -> GetBinWidth(1);
  lminPrecision[0] = fHistTopDCAPosToPV      -> GetBinWidth(1);
  lminPrecision[1] = fHistTopDCANegToPV      -> GetBinWidth(1);
  lminPrecision[2] = fHistTopDCAV0Daughters  -> GetBinWidth(1);

  cout<<"Minimum Values Found: "<<endl;
  cout<<"Cosine of Pointing Angle...........: "<<GetXForMinValue (fHistTopCosinePA )<<" precision = "<<lminPrecision[3]<<endl;
  cout<<"DCA Neg Daughter to PV.............: "<<GetXForMaxValue (fHistTopDCAPosToPV )<<" precision = "<<lminPrecision[1]<<endl;
  cout<<"DCA Pos Daughter to PV.............: "<<GetXForMaxValue (fHistTopDCANegToPV )<<" precision = "<<lminPrecision[0]<<endl;
  cout<<"DCA V0 Daughers....................: "<<GetXForMinValue (fHistTopDCAV0Daughters )<<" precision = "<<lminPrecision[2]<<endl;
  cout<<"V0 Decay Radius....................: "<<GetXForMinValue (fHistTopV0Radius )<<" precision = "<<lminPrecision[4]<<endl;


  TCanvas *cTopo = new TCanvas ("cTopo","",1200,800);
  cTopo->Divide(3,2); 
  for(Int_t ic = 1; ic<7; ic++){ 
    cTopo->cd(ic)->SetLeftMargin(0.15);
    cTopo->cd(ic)->SetGridx(); 
    cTopo->cd(ic)->SetGridy(); 
  }

  cTopo->cd(1);   
  fHistTopDCAPosToPV->Draw(); 
  cTopo->cd(2);   
  fHistTopDCANegToPV->Draw(); 
  cTopo->cd(3);   
  fHistTopDCAV0Daughters->Draw(); 
  cTopo->cd(4);   
  fHistTopCosinePA->Draw(); 
  cTopo->cd(5);   
  fHistTopV0Radius->Draw(); 

  TLatex Tlt;
  Tlt.SetNDC();
  Tlt.SetTextSize(0.05);
  cTopo->cd(6);
    
  Tlt.DrawLatex(.22, .9, "Boundary Checks")  ;


        
  TString lCut[5];
  lCut [ 0 ] = "Min DCA Pos D. To PV (cm)";
  lCut [ 1 ] = "Min DCA Neg D. To PV (cm)";
  lCut [ 2 ] = "Max DCA V0 Daughters (#sigma)";
  lCut [ 3 ] = "Min Cosine PA";
  lCut [ 4 ] = "Min 2D Decay Radius (cm)";

  TString lCutVal[5]; 
  TString lCutValPrec[5]; 
    
  Tlt.SetTextSize(0.04);
  Tlt.SetTextFont(42);

  Tlt.DrawLatex(.01, .80, "Topological Var.")  ;
  Tlt.DrawLatex(.6, .80, "Value")  ;
  Tlt.DrawLatex(.75, .80, "Precision")  ;

  for (Int_t il=0;il<5;il++){ 
    Tlt.DrawLatex(.01,0.72-((double)il)*0.075, lCut[il].Data() );
    lCutVal[il] = Form( "%.4f", lmin[il] );
    Tlt.DrawLatex(.5925,0.72-((double)il)*0.075, lCutVal[il].Data() );
    lCutValPrec[il] = Form( "%.4f", lminPrecision[il] );
    Tlt.DrawLatex(.7675,0.72-((double)il)*0.075, lCutValPrec[il].Data() );
  }

  Tlt.SetTextSize(0.05);
  Tlt.SetTextFont(42);

  //Try to make a wild guess... 
  if ( TMath::Abs( lmin[0] - 0.100 ) < 2*lminPrecision[0] &&
       TMath::Abs( lmin[1] - 0.100 ) < 2*lminPrecision[1] &&
       TMath::Abs( lmin[2] - 1.000 ) < 2*lminPrecision[2] &&
       TMath::Abs( lmin[3] - 0.998 ) < 2*lminPrecision[3] &&
       TMath::Abs( lmin[4] - 0.500 ) < 2*lminPrecision[4] ){
    Tlt.DrawLatex( 0.17, 0.275, "#bf{Autodetect Interpretation}: "); 
    Tlt.DrawLatex( 0.1, 0.2, "Typical #bf{Pb-Pb} Reconstruction Cuts");
  }
  if ( TMath::Abs( lmin[0] - 0.020 ) < 2*lminPrecision[0] &&
       TMath::Abs( lmin[1] - 0.020 ) < 2*lminPrecision[1] &&
       TMath::Abs( lmin[2] - 1.500 ) < 2*lminPrecision[2] &&
       TMath::Abs( lmin[3] - 0.98 )  < 2*lminPrecision[3] &&
       TMath::Abs( lmin[4] - 0.500 ) < 2*lminPrecision[4] ){
    Tlt.DrawLatex( 0.15, 0.29, "#bf{Autodetect Interpretation}: ");  
    Tlt.DrawLatex( 0.1, 0.2, "Typical pp Reconstruction Cuts"); //To be checked
  }

  if      (output == "png") cTopo->SaveAs("LF_QAanalysis_V0_page6.png");
  else if (output == "eps") cTopo->SaveAs("LF_QAanalysis_V0_page6.eps");
  else if (output == "pdf") cTopo->SaveAs("LF_QAanalysis_V0.pdf)");


}

CustomGStyleSettings(){ 
  gStyle->SetOptStat(0); 
  gStyle->SetGridColor(kGray); 

}

Double_t GetXForMaxValue( TH1D *lHist){ 
  Double_t lMaximum = 1e-10; 
  Double_t lMaximumX = 1000; 
  Double_t lThisBin = -1; 
  for(Long_t ib=1; ib<lHist->GetNbinsX()+1; ib++){ 
    lThisBin = lHist->GetBinContent(ib);
    if ( lThisBin > lMaximum ){ 
      lMaximum = lHist->GetBinContent(ib); 
      lMaximumX = lHist->GetBinLowEdge(ib); 
    }
  }
  return lMaximumX; 
}

Double_t GetXForMinValue( TH1D *lHist){ 
  Double_t lThreshold = 1e-7;
  Double_t lMinimum = 1e+10; 
  Double_t lMinimumX = 1000; 
  Double_t lThisBin = -1; 
  for(Long_t ib=1; ib<lHist->GetNbinsX()+1; ib++){ 
    lThisBin = lHist->GetBinContent(ib);
    if ( lThisBin < lMinimum && lThisBin>lThreshold ){ 
      lMinimum = lHist->GetBinContent(ib); 
      lMinimumX = lHist->GetBinLowEdge(ib); 
    }
  }
  return lMinimumX; 
}

//Trick for background fitting (pol2) 
Double_t MyBgPol2(const Double_t *x, const Double_t *par)
{
        //par[0] and par[1] -> linear and angular coefficients
        //par[2] and par[3] -> lower and upper peak boundaries
        if ( x[0] > par[2] && x[0] < par[3]) {
                TF1::RejectPoint();
                return 0;
        }
        return par[0] + par[1]*x[0];
}



