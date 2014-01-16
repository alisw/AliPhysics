void MakeFinalSpectrum()
{
  //-----------------------------------------------------------------------------
  // This macro takes the raw pi0 spectrum from LHC11a_FitResult_20130314.root,
  // correct it for feed down, then for efficiency,
  // adds all systematic errors and produces the production invariant spectrum
  //--
  // Last modification: 06.07.2012 by Yuri Kharlov
  //-----------------------------------------------------------------------------

  TFile *f  = new TFile("LHC11a_FitResult_20130913.root");
  TH1D * nr1CB    = (TH1D*)f->Get("Mix_CB_yr1")  ;
  TH1D * nr1intCB = (TH1D*)f->Get("Mix_CB_yr1int")  ;
  TH1D * nr2CB    = (TH1D*)f->Get("Mix_CB_yr2")  ;
  TH1D * nr2intCB = (TH1D*)f->Get("Mix_CB_yr2int")  ;

  TH1D * nr1GS    = (TH1D*)f->Get("Mix_yr1")  ;
  TH1D * nr1intGS = (TH1D*)f->Get("Mix_yr1int")  ;
  TH1D * nr2GS    = (TH1D*)f->Get("Mix_yr2")  ;
  TH1D * nr2intGS = (TH1D*)f->Get("Mix_yr2int")  ;

  //Divide by bin width
  for(Int_t i=1;i<= nr1CB->GetNbinsX();i++){
    nr1CB   ->SetBinContent(i,nr1CB->GetBinContent(i)/nr1CB->GetXaxis()->GetBinWidth(i)) ;
    nr1CB   ->SetBinError  (i,nr1CB->GetBinError(i)/nr1CB->GetXaxis()->GetBinWidth(i)) ;
    nr1intCB->SetBinContent(i,nr1intCB->GetBinContent(i)/nr1intCB->GetXaxis()->GetBinWidth(i)) ;
    nr1intCB->SetBinError  (i,nr1intCB->GetBinError(i)/nr1intCB->GetXaxis()->GetBinWidth(i)) ;
    nr2CB   ->SetBinContent(i,nr2CB->GetBinContent(i)/nr2CB->GetXaxis()->GetBinWidth(i)) ;
    nr2CB   ->SetBinError  (i,nr2CB->GetBinError(i)/nr2CB->GetXaxis()->GetBinWidth(i)) ;
    nr2intCB->SetBinContent(i,nr2intCB->GetBinContent(i)/nr2intCB->GetXaxis()->GetBinWidth(i)) ;
    nr2intCB->SetBinError  (i,nr2intCB->GetBinError(i)/nr2intCB->GetXaxis()->GetBinWidth(i)) ;

    nr1GS   ->SetBinContent(i,nr1GS->GetBinContent(i)/nr1GS->GetXaxis()->GetBinWidth(i)) ;
    nr1GS   ->SetBinError  (i,nr1GS->GetBinError(i)/nr1GS->GetXaxis()->GetBinWidth(i)) ;
    nr1intGS->SetBinContent(i,nr1intGS->GetBinContent(i)/nr1intGS->GetXaxis()->GetBinWidth(i)) ;
    nr1intGS->SetBinError  (i,nr1intGS->GetBinError(i)/nr1intGS->GetXaxis()->GetBinWidth(i)) ;
    nr2GS   ->SetBinContent(i,nr2GS->GetBinContent(i)/nr2GS->GetXaxis()->GetBinWidth(i)) ;
    nr2GS   ->SetBinError  (i,nr2GS->GetBinError(i)/nr2GS->GetXaxis()->GetBinWidth(i)) ;
    nr2intGS->SetBinContent(i,nr2intGS->GetBinContent(i)/nr2intGS->GetXaxis()->GetBinWidth(i)) ;
    nr2intGS->SetBinError  (i,nr2intGS->GetBinError(i)/nr2intGS->GetXaxis()->GetBinWidth(i)) ;
  }

  // feed down correction

  TF1 *fKaonContaminationToPi0 = new TF1("kaonCont","1./(1.-1.33*1.2*exp(-2.95-0.16*x))",0.,30.);
  nr1GS   ->Divide(fKaonContaminationToPi0);
  nr1intGS->Divide(fKaonContaminationToPi0);
  nr2GS   ->Divide(fKaonContaminationToPi0);
  nr2intGS->Divide(fKaonContaminationToPi0);
  nr1CB   ->Divide(fKaonContaminationToPi0);
  nr2CB   ->Divide(fKaonContaminationToPi0);
  nr1intCB->Divide(fKaonContaminationToPi0);
  nr2intCB->Divide(fKaonContaminationToPi0);

  // SPD pileup correction

  TF1 *fSPDpileup = new TF1("SPDpileup","0.988",0.,30.);
  nr1GS   ->Multiply(fSPDpileup);
  nr1intGS->Multiply(fSPDpileup);
  nr2GS   ->Multiply(fSPDpileup);
  nr2intGS->Multiply(fSPDpileup);
  nr1CB   ->Multiply(fSPDpileup);
  nr2CB   ->Multiply(fSPDpileup);
  nr1intCB->Multiply(fSPDpileup);
  nr2intCB->Multiply(fSPDpileup);

  //correct for efficiency
  TFile *fEff = new TFile("Pi0_efficiency_LHC11a__20131029_Mall.root") ;

  TF1 * effGS=fEff->Get("eff_Pi0_Gaus_2760GeV") ;
  TF1 * effCB=fEff->Get("eff_Pi0_CB_2760GeV") ;
  effGS   ->SetRange(0.,25.) ;
  effCB   ->SetRange(0.,25.) ;

  nr1GS   ->Divide(effGS) ;
  nr1intGS->Divide(effGS) ;
  nr2GS   ->Divide(effGS) ;
  nr2intGS->Divide(effGS) ;
  nr1CB   ->Divide(effCB) ;
  nr2CB   ->Divide(effCB) ;
  nr1intCB->Divide(effCB) ;
  nr2intCB->Divide(effCB) ;

  //make 1/pt
  for(Int_t i=1;i<=nr1CB->GetNbinsX();i++){
    Double_t pt = TMath::TwoPi()*nr1CB->GetXaxis()->GetBinCenter(i);
    nr1CB   ->SetBinContent(i,nr1CB   ->GetBinContent(i)/pt) ;
    nr1CB   ->SetBinError  (i,nr1CB   ->GetBinError(i)  /pt) ;
    nr1intCB->SetBinContent(i,nr1intCB->GetBinContent(i)/pt) ;
    nr1intCB->SetBinError  (i,nr1intCB->GetBinError(i)  /pt) ;
    nr2CB   ->SetBinContent(i,nr2CB   ->GetBinContent(i)/pt) ;
    nr2CB   ->SetBinError  (i,nr2CB   ->GetBinError(i)  /pt) ;
    nr2intCB->SetBinContent(i,nr2intCB->GetBinContent(i)/pt) ;
    nr2intCB->SetBinError  (i,nr2intCB->GetBinError(i)  /pt) ;

    nr1GS   ->SetBinContent(i,nr1GS   ->GetBinContent(i)/pt) ;
    nr1GS   ->SetBinError  (i,nr1GS   ->GetBinError(i)  /pt) ;
    nr1intGS->SetBinContent(i,nr1intGS->GetBinContent(i)/pt) ;
    nr1intGS->SetBinError  (i,nr1intGS->GetBinError(i)  /pt) ;
    nr2GS   ->SetBinContent(i,nr2GS   ->GetBinContent(i)/pt) ;
    nr2GS   ->SetBinError  (i,nr2GS   ->GetBinError(i)  /pt) ;
    nr2intGS->SetBinContent(i,nr2intGS->GetBinContent(i)/pt) ;
    nr2intGS->SetBinError  (i,nr2intGS->GetBinError(i)  /pt) ;
  }

  //For the final spectrum we take average of fits
  //with numerical integration of entries in signal
  
  TH1D * hStat = (TH1D*)nr2intCB->Clone("hPi02760GeVStat") ;
  TH1D * hSys  = (TH1D*)hStat   ->Clone("hPi02760GeVSys") ;
  TH1D * hSys2 = (TH1D*)hStat   ->Clone("hPi02760GeVSysTypeB") ;
  TH1D * hSys3 = (TH1D*)hStat   ->Clone("hPi02760GeVSysTypeC") ;
  hStat->SetAxisRange(0.,14.9,"X");
  hSys ->SetAxisRange(0.,14.9,"X");

  //For systematic error estimate take largest deviation
  //of integrated yeilds (note, they are efficiency corrected)
  for(Int_t i=1;i<=nr1CB->GetNbinsX();i++){
    Double_t mean= hStat->GetBinContent(i) ;
    Double_t dev = TMath::Max(
                   TMath::Max(TMath::Abs(nr1intCB->GetBinContent(i)-mean),
                              TMath::Abs(nr2intCB->GetBinContent(i)-mean)),
                   TMath::Max(TMath::Abs(nr1intGS->GetBinContent(i)-mean),
                              TMath::Abs(nr2intGS->GetBinContent(i)-mean))
			      );
    hSys ->SetBinError(i,dev) ;
    hSys2->SetBinError(i,dev) ;
  }

  //Add other sys errors
  TF1 * globalE = new TF1("globalE","1.-((x+1.354)/(x*1.002+1.354))^6.18 ",1.,30.) ; 
  TF1 * conv    = new TF1("conversion","0.035",0.,30.) ;
  TF1 * accept  = new TF1("accept"    ,"0.01" ,0.,30.) ;
  TF1 * pileup  = new TF1("pileup"    ,"0.004",0.,30.) ;
  TF1 * calib   = new TF1("calib"     ,"0.005",0.,30.) ;
  TF1 * modDiff = new TF1("modDiff"   ,"0.04",0.,30.) ;
  // TF1 * modDiff = new TF1("modDiff"   ,"16.9*exp(-4.5*x)+0.033",0.,30.) ;
  TF1 * tofCut  = new TF1("tofCut"    ,"0.0105" ,0.,30.) ;

  //Borya's estimate of non-linearity (found for pp @ 7 TeV)
  TF1 * nonlin= new TF1("nl","0.015+7.38*exp(-x/0.24)",0.,30.) ;

  //Draw sys errors
  TH1D * hRelSysRaw = (TH1D*)hSys->Clone("RelSysRaw") ;
  hRelSysRaw->SetTitle("Summary of systematic errors");
  for(Int_t i=1;i<=hSys->GetNbinsX();i++){
    Double_t mean= hSys->GetBinContent(i) ;
    Double_t a=hSys->GetBinError(i) ;
    if(mean>0)
      hRelSysRaw->SetBinContent(i,a/mean) ;
    else
      hRelSysRaw->SetBinContent(i,0.) ;
      hRelSysRaw->SetBinError(i,0.) ;
  }

  //Add errors in sys errors

  TH1D * hRelSysTot = (TH1D*)hSys->Clone("RelSys") ;
  hRelSysTot->SetTitle("Summary of systematic uncertainties");

  for(Int_t i=1;i<=hSys->GetNbinsX();i++){
    Double_t pt   = hSys->GetXaxis()->GetBinCenter(i) ;
    Double_t mean = hSys->GetBinContent(i) ;
    Double_t a    = hSys->GetBinError(i) ;
    // Double_t b    = mean * hSysErrModules->GetBinContent(i) ;
    Double_t tot= mean*mean*nonlin ->Eval(pt)*nonlin ->Eval(pt) 
                 +mean*mean*conv   ->Eval(pt)*conv   ->Eval(pt)
                 +mean*mean*accept ->Eval(pt)*accept ->Eval(pt)
                 +mean*mean*pileup ->Eval(pt)*pileup ->Eval(pt)
                 +mean*mean*calib  ->Eval(pt)*calib  ->Eval(pt)
                 +mean*mean*modDiff->Eval(pt)*modDiff->Eval(pt)
                 +mean*mean*tofCut ->Eval(pt)*tofCut ->Eval(pt)
                 +mean*mean*globalE->Eval(pt)*globalE->Eval(pt); 
    Double_t raa= mean*mean*nonlin->Eval(pt)*nonlin->Eval(pt) 
                 +mean*mean*pileup->Eval(pt)*pileup->Eval(pt)
                 +mean*mean*calib ->Eval(pt)*calib ->Eval(pt); 
    hSys3->SetBinError(i,TMath::Sqrt(tot)) ;
    // hSys->SetBinError(i,TMath::Sqrt(tot + a*a + b*b)) ;
    hSys2->SetBinError(i,TMath::Sqrt(raa + a*a)) ;
    hSys ->SetBinError(i,TMath::Sqrt(tot + a*a)) ;
    
    a = hSys->GetBinError(i) ;
    printf("i=%d, %g+-%g\n",i,mean,a);
    if(mean>0)
      hRelSysTot->SetBinContent(i,a/mean) ;
    else {
      hRelSysTot->SetBinContent(i,0.) ;
      hRelSysTot->SetBinError(i,0.) ;
    }
  }

  // TFile *fSysErrModules = TFile::Open("PHOS_sysErr_modules.root");
  // TH1D* hSysErrModules = (TH1D*)fSysErrModules->Get("hSyserr");

  gStyle->SetOptStat(0);
  TCanvas * c = new TCanvas("c","SysErrors") ;
  c->SetLogy();
  c->cd() ;
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.07);
  hRelSysTot->SetAxisRange(0.8    ,11.9,"X");
  hRelSysTot->SetAxisRange(0.0031,0.45,"Y");
  hRelSysTot->GetYaxis()->SetMoreLogLabels();
  hRelSysTot->GetYaxis()->SetNoExponent();
  hRelSysTot->SetNdivisions(520,"X");
  hRelSysTot->SetLineColor(1) ;
  hRelSysTot->SetLineWidth(2) ;
  hRelSysRaw->SetLineColor(2) ;
  hRelSysRaw->SetLineWidth(2) ;
  globalE->SetLineColor(kGreen+2) ;
  nonlin ->SetLineColor(4) ;
  conv   ->SetLineColor(6) ;
  accept ->SetLineColor(kOrange) ;
  pileup ->SetLineColor(42);
  calib  ->SetLineColor(44);
  calib  ->SetLineStyle(2);
  modDiff->SetLineColor(52);
  modDiff->SetLineStyle(2);
  tofCut ->SetLineColor(53);
  tofCut ->SetLineStyle(3);
  // hSysErrModules->SetLineColor(kOrange+2);
  // hSysErrModules->SetLineStyle(2);
  hRelSysTot->SetXTitle("p_{T} (GeV/c)") ;
  hRelSysTot->SetYTitle("Rel.syst.error") ;
  hRelSysTot->GetYaxis()->SetTitleOffset(1.2) ;
  hRelSysTot->Draw("hist") ;
  hRelSysRaw->Draw("h same") ;
  globalE->Draw("same") ; 
  nonlin ->Draw("same") ;
  conv   ->Draw("same") ;
  accept ->Draw("same") ;
  pileup ->Draw("same");
  calib  ->Draw("same");
  modDiff->Draw("same");
  tofCut ->Draw("same");
  // hSysErrModules->Draw("same ][");
  TLegend * l = new TLegend(0.57,0.62,0.85,0.925) ;
  l->SetFillColor(kWhite);
  l->SetBorderSize(0);
  l->AddEntry(hRelSysTot,"Total uncertainty","l") ;
  l->AddEntry(hRelSysRaw,"Raw extraction","l") ;
  l->AddEntry(conv   ,"Conversion","l") ;
  l->AddEntry(nonlin ,"Non-linearity","l") ;
  l->AddEntry(accept ,"Acceptance","l");
  l->AddEntry(pileup ,"Pileup","l");
  l->AddEntry(calib  ,"Rel.calib.","l");
  l->AddEntry(modDiff,"Per-module yield","l");
  l->AddEntry(tofCut ,"Timing cut","l");
  // l->AddEntry(hSysErrModules,"Intermodule spectra","l");
  l->AddEntry(globalE,"Global E scale","l") ;
  l->Draw() ;

  c->Print("LHC11a_SysErrors.eps");

  hStat->SetTitle("Normalized production #pi^{0} yield, pp @ 2.76 TeV, stat.err.");
  hStat->SetXTitle("p_{T}, GeV/c");
  hStat->SetYTitle("1/N_{ev} 1/(2#pi p_{T}) d^{2}N/dydp_{T} (GeV^{-2}c^{2})");
  hSys ->SetTitle("Normalized production #pi^{0} yield, pp @ 2.76 TeV, total syst.err.");
  hSys ->SetXTitle("p_{T}, GeV/c");
  hSys ->SetYTitle("1/N_{ev} 1/(2#pi p_{T}) d^{2}N/dydp_{T} (GeV^{-2}c^{2})");
  hSys2->SetTitle("Normalized production #pi^{0} yield, pp @ 2.76 TeV, syst.err. for R_{AA}");
  hSys3->SetTitle("Normalized production #pi^{0} yield, pp @ 2.76 TeV, apparatus syst.err.");
  hSys ->SetFillColor(kBlue-10) ;
  hSys ->SetNdivisions(520,"X");

  hSysFinal = (TH1F*)hSys->Clone("hSysFinal");
  hSysFinal ->SetTitle("Production #pi^{0} yield, pp @ 2.76 TeV, 07.11.2013");
  hSysFinal ->GetYaxis()->SetTitleOffset(1.2);
  hStat->GetYaxis()->SetTitleOffset(1.2);
  hSysFinal ->SetAxisRange(0.8,9.9,"X");

  TCanvas *c2 = new TCanvas("c2","Production spectrum");
  gPad->SetRightMargin(0.02);
  gPad->SetTopMargin(0.07);
  gPad->SetGridx();
  gPad->SetGridy();
  c2->SetLogy();
  hSysFinal->Draw("E2") ;
  hStat->SetMarkerStyle(20) ;
  hStat->SetMarkerColor(4) ;
  hStat->SetLineColor(4) ;
  hStat->Draw("same") ;
  
  // //Apply bin width correction
  // TH1D * hBWcorr = BinWidthCorrection(hStat) ;
  // hSys->Divide(hBWcorr) ;
  // hSys2->Divide(hBWcorr) ;
  // hSys3->Divide(hBWcorr) ;
  // hStat->Divide(hBWcorr) ;

  c2->Print("LHC11a_pi0Spectrum.eps");

  TF1 *tsallis = FitTsallis(hStat,hSys);
  
  TFile fout("PHOS_pi0_2760eV_noBinWidthCorr_20131112.root","recreate") ;
  hSys   ->Write() ;
  hSys2  ->Write() ;
  hSys3  ->Write() ;
  hStat  ->Write() ;
  tsallis->Write() ;
  fout.Close() ;


}

//-----------------------------------------------------------------------------
TF1 *FitTsallis(TH1D* hStat, TH1D* hSys)
{
  TF1 * fit = new TF1("Tsalis","[0]/2./3.1415*([2]-1.)*([2]-2.)/([2]*[1]*([2]*[1]+0.135*([2]-2.)))*(1.+(sqrt(x*x+0.135*0.135)-0.135)/([2]*[1]))^-[2]",0.5,25.) ;
  fit->SetParameters(10.,0.2,8.) ;
  fit->SetLineColor(kBlack) ;
  fit->SetLineWidth(2) ;
  fit->SetParName(0,"dN/dy") ;
  fit->SetParName(1,"T") ;
  fit->SetParName(2,"n") ;
  
  TH1D * hSysRatio  = (TH1D*)hSys ->Clone("RatioSys") ;
  TH1D * hStatRatio = (TH1D*)hStat->Clone("RatioStat") ;

  TH1D * hsum = (TH1D*)hStat->Clone("sum") ;
  for(Int_t i=1; i<=hsum->GetNbinsX();i++){
    Double_t a=hSys->GetBinError(i) ;
    Double_t b=hStat->GetBinError(i) ;
    hsum->SetBinError(i,TMath::Sqrt(a*a+b*b)) ;
  }
  hsum->SetStats(1) ;
  hsum->Fit(fit,"Q") ;
  Double_t meanPt=fit->Moment(1,1.,10.) ;
  printf("<pt>=%f \n",meanPt) ;

  hSysRatio ->Divide(fit) ;
  hStatRatio->Divide(fit) ;

  return fit;
}
//-----------------------------------------------------------------------------
TH1D * BinWidthCorrection(TH1D * h){
  //We apply bin width a-la PHENIX 
  //Use Tsalis fit to calculate shift in y direction
 
  TF1 * fit = new TF1("hag","[0]*(([2]*[1]+sqrt(x*x+0.135*0.135))/([2]*[1]+0.135))^-[2]",0.5,25.) ;
  fit->SetParameters(10.,0.2,8.) ;
  TCanvas * corr = new TCanvas("BWcorr","Bin width correction") ;
  Int_t col[6]={kRed,kOrange,kMagenta,kGreen,kCyan,kBlue} ;
  TH1D * hcorr[20] ;
  char key[55] ;
  Double_t rMax=10 ;
  Int_t iteration=0 ;
  TH1D * htmp = (TH1D*)h->Clone("tmp") ;
  while(iteration<6){
    printf(" Iteration %d: rMax=%f \n",iteration, rMax) ;
    htmp->Fit(fit,"N") ;
    sprintf(key,"Ineration%d",iteration) ;
    hcorr[iteration]=(TH1D*)h->Clone(key);
    rMax= 0; 
    for(Int_t i=1;i<=h->GetNbinsX();i++){
      Double_t a=h->GetXaxis()->GetBinLowEdge(i) ;
      Double_t b=h->GetXaxis()->GetBinUpEdge(i) ;
      Double_t r=fit->Integral(a,b)/(b-a)/fit->Eval(0.5*(a+b)) ;
      hcorr[iteration]->SetBinContent(i,r) ;
      hcorr[iteration]->SetBinError(i,0.) ;
      if(rMax<r)rMax=r ;
    }
    delete htmp ;
    htmp = (TH1D*)h->Clone("tmp") ;
    htmp->Divide(hcorr[iteration]) ;
    corr->cd() ;
    hcorr[iteration]->SetLineColor(col[iteration]);
    if(iteration==0)
      hcorr[iteration]->Draw() ;
    else
      hcorr[iteration]->Draw("same") ;
    corr->Update() ;
    iteration++ ;
  } 

  hcorr[5]->SetTitle("Bin-width correction for #pi^{0} spectrum");
  hcorr[5]->SetYTitle("Bin-width corrected / uncorrected");
  TFile fout("PHOS_pi0_7TeV_BinWidthCorrection.root","recreate") ;
  hcorr[5]->Write();
  fout.Close();

  return hcorr[5] ;
}
