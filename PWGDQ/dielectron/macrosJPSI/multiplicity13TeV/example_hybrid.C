void example_hybrid(TString filename = "dstAnalysisHistograms.root",  TString cutName = "relaxedPID" ){
  
  TFile file( filename.Data() );
  TList* qa = (TList *)file.Get( "jpsi2eeHistos" );
  if(!qa){
    printf("QA list not found!");
    return;
  }
  file.Close();
  
  
  THashList* listJpsiPM = (THashList*) qa->FindObject( Form("PairSEPM_%s", cutName.Data() ) );
  THnF* hPM = (THnF*) listJpsiPM->FindObject("PairInvMass");
    
  THashList* listJpsiTR = (THashList*) qa->FindObject( Form("PairTRPM_%s", cutName.Data() ) );
  THnF* hTR = (THnF*) listJpsiTR->FindObject("PairInvMass");
   
  
  TH1F* hMassPM = (TH1F*) hPM->Projection(4);
  TH1F* hMassTR = (TH1F*) hTR->Projection(4);
  
  
  cout << hMassPM << endl;
  cout << hMassTR << endl;
  
  AliDielectronSignalFunc* sig = new AliDielectronSignalFunc;
  sig->SetIntegralRange(2.92, 3.16);
  
  
  sig->SetMethod( AliDielectronSignalBase::kCombinatorialPlusFit );
  sig->SetScaleRawToBackground( 4.5, 5.9,0.,0.);

  double fitFrom = 1.5;
  double fitUntil = 4.9;

    
  TF1* fB = new TF1( "fitBg"   , "1+expo", fitFrom, fitUntil);
  fB->SetParameter(0, 2.);
  fB->SetParLimits(0, 0., 10.);
  fB->SetParameter(1, -1.);
  fB->SetParLimits(1, -10., 0.);


  sig->SetFunctions( fB, fB, fB );
  sig->SetFitRange( fitFrom, fitUntil );
               

  TObjArray *arr = new TObjArray;
  arr->SetOwner();
  arr->AddAt( hMassPM, AliDielectron::kEv1PM );
  arr->AddAt( hMassTR, AliDielectron::kEv1PMRot );

  sig->Process(arr);
  plot(sig);
    
}
    
 void plot ( AliDielectronSignalFunc* sig) {
    
    
  TCanvas *c = new TCanvas( "portrait",  "portrait",  600, 800 );
  c->SetTopMargin(0.01);
  c->SetRightMargin(0.01);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.01);
  c->Clear();
  c->Divide(1,2,0,0);
  
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  
  TH1F* hUS = sig->GetUnlikeSignHistogram();
  TH1F* hBackground = sig->GetBackgroundHistogram();
  TH1F* hSignal = sig->GetSignalHistogram();
        
  
  Double_t sigN = sig->GetSignal();    
  Double_t sigEr = sig->GetSignalError();
  Double_t sigS2B = sig->GetSB();
  Double_t sigS2Ber = sig->GetSBError();
  Double_t sigSignif = sig->GetSignificance();
  Double_t sigSignifEr = sig->GetSignificanceError();


  TH1* combinatorialBg  = sig->GetCombinatorialBackgroundHistogram();
  combinatorialBg->SetLineColor(kAzure+1);
        
  
  hUS->SetMarkerStyle(20);
  hUS->SetMarkerSize(0.7);
  hUS->SetMarkerColor(kRed);
  hUS->SetLineColor(kRed);
  hUS->SetTitle(Form("Unlike sign spectrum M inv. ee;counts per %.4g MeV", 1000*hUS->GetXaxis()->GetBinWidth(1)));

  
  hBackground->SetMarkerStyle(21);
  hBackground->SetMarkerSize(0.7);
  hBackground->SetMarkerColor(kBlack);
  hBackground->SetLineColor(kBlack);

    
  hSignal->GetXaxis()->SetTitle("#it{m}_{e^{+}e^{-}} (GeV/#it{c}^{2})");
  hSignal->GetYaxis()->SetTitle(Form("Counts per %.4g MeV/#it{c}^{2}", 1000*hSignal->GetXaxis()->GetBinWidth(1)));
    
  hSignal->SetMarkerSize(.7);
  hSignal->SetMarkerStyle(20); 
  hSignal->SetMarkerColor(2); 
  hSignal->SetLineColor(2);
    
        
  c->cd(1);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.0);
  gPad->SetTickx();
  gPad->SetTicky();

  hUS->GetXaxis()->SetRangeUser( .5, 5. );
  
  hUS->GetYaxis()->SetRangeUser(-.1, 1.5 * hUS->GetMaximum() );
  
  hUS->GetYaxis()->SetTitle(Form("Counts per %.4g MeV/#it{c}^{2}", 1000*hUS->GetXaxis()->GetBinWidth(1)));
  hBackground->GetXaxis()->SetRangeUser( .5, 5. );
  
  hUS->SetBinErrorOption(TH1::kPoisson);
  hUS->Draw("e");
  hBackground->Draw("same");
  combinatorialBg->Draw("same"); 

  
  
  TLatex *lat = new TLatex;
  lat->SetNDC(kTRUE);
  lat->SetTextFont(42);
  lat->SetTextSize(.03);
   
  
  lat->DrawLatex( 0.68, 0.70, Form("#bf{Signal:\t%3.1f #pm %2.1f}", sigN, sigEr) );
  lat->DrawLatex( 0.68, 0.63, Form("S/B:\t\t%1.2f #pm %1.2f",sigS2B,sigS2Ber ) );
  lat->DrawLatex( 0.68, 0.56, Form("S/#sqrt{S+B}:\t%1.2f #pm %1.2f ",  sigSignif, sigSignifEr) );
  

    
  c->cd(2);                 
  gPad->SetRightMargin(0.01);
  gPad->SetTopMargin(0);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(.2);
  gPad->SetTickx();
  gPad->SetTicky();


  hSignal->GetXaxis()->SetRangeUser( .5, 5. );
  hSignal->GetYaxis()->SetRangeUser( 2.*hSignal->GetMinimum() , 1.4*hSignal->GetMaximum() );
  hSignal->Draw();

  c->SaveAs("invMass.pdf");
  
}