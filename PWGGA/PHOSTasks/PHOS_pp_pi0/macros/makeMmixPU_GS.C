//-----------------------------------------------------------------------------
void makeMmixPU_GS(const TString histoFile="LHC11a_pass4_20130913.root",
		   const Int_t nSigma=2,
		   const char* module="A10")
{
  //Fit Real/Mixed ratio, normalize Mixed and subtract it from Real.
  // The pi0 peak if fitted by the Gaussian function, 
  // the background is fitted by pol1 or pol2

  TString hMassName;
  TFile * file = new TFile(histoFile) ;
  THashList *hList = (THashList*)file->Get("histESD");
  char key[125] ;
  hMassName = "hMassPt";
  hMassName += module;
  TH2F * h   = (TH2F*)hList->FindObject(hMassName) ;

  hMassName = "hMiMassPt";
  hMassName += module;
  TH2F * hm  = (TH2F*)hList->FindObject(hMassName) ;

  TH1F * hev = (TH1F*)hList->FindObject("hSelEvents") ;

  // Array of pt bins
  Int_t nbin=18 ;
  Double_t xa[]={0.6,0.8,1.,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,3.0,3.5,4.0,5.,6.,8.,10.,12.} ;
  PPRstyle();
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadBottomMargin(0.11);
  gStyle->SetTitleX(0.80);
  gStyle->SetTitleY(0.98);

  //Fit real only 
  //Linear Bg
  char key2[155];
  sprintf(key,"Mix") ;
  sprintf(key2,"%s_mr1",key) ;
  TH1D * mr1 = new TH1D(key2,"Mass",nbin,xa) ;
  sprintf(key2,"%s_sr1",key) ;
  TH1D * sr1 = new TH1D(key2,"Width",nbin,xa) ;
  sprintf(key2,"%s_ar1",key) ;
  TH1D * ar1 = new TH1D(key2,"a",nbin,xa) ;
  sprintf(key2,"%s_br1",key) ;
  TH1D * br1 = new TH1D(key2,"a",nbin,xa) ;
  sprintf(key2,"%s_yr1",key) ;
  TH1D * nr1 = new TH1D(key2,"Raw yield",nbin,xa) ;
  sprintf(key2,"%s_yr1int",key) ;
  TH1D * nr1int = new TH1D(key2,"Raw yield, integrated",nbin,xa) ;

  //Quadratic Bg
  sprintf(key2,"%s_mr2",key) ;
  TH1D * mr2 = new TH1D(key2,"Mass",nbin,xa) ;
  sprintf(key2,"%s_sr2",key) ;
  TH1D * sr2 = new TH1D(key2,"Width",nbin,xa) ;
  sprintf(key2,"%s_ar2",key) ;
  TH1D * ar2 = new TH1D(key2,"a",nbin,xa) ;
  sprintf(key2,"%s_br2",key) ;
  TH1D * br2 = new TH1D(key2,"a",nbin,xa) ;
  sprintf(key2,"%s_cr2",key) ;
  TH1D * cr2 = new TH1D(key2,"a",nbin,xa) ;
  sprintf(key2,"%s_yr2",key) ;
  TH1D * nr2 = new TH1D(key2,"Raw yield",nbin,xa) ;
  sprintf(key2,"%s_yr2int",key) ;
  TH1D * nr2int = new TH1D(key2,"Raw yield, integrated",nbin,xa) ;

  TF1 * fit1 = new TF1("fit",CB,0.,1.,5) ;
  fit->SetParName(0,"A") ;
  fit->SetParName(1,"m_{0}") ;
  fit->SetParName(2,"#sigma") ;
  fit->SetParName(3,"a_{0}") ;
  fit->SetParName(4,"a_{1}") ;
  fit->SetLineWidth(2) ;
  fit->SetLineColor(2) ;
  TF1 * fgs = new TF1("gs",CBs,0.,1.,3) ;
  fgs->SetLineColor(2) ;
  fgs->SetLineWidth(1) ;

  TF1 * fit2 = new TF1("fit2",CB2,0.,1.,6) ;
  fit2->SetParName(0,"A") ;
  fit2->SetParName(1,"m_{0}") ;
  fit2->SetParName(2,"#sigma") ;
  fit2->SetParName(3,"a_{0}") ;
  fit2->SetParName(4,"a_{1}") ;
  fit2->SetLineWidth(2) ;
  fit2->SetLineColor(4) ;
  fit2->SetLineStyle(2) ;
  TF1 * fbg1 = new TF1("bg1",BG1,0.,1.,2) ;
  TF1 * fbg2 = new TF1("bg2",BG2,0.,1.,3) ;

  TCanvas * c3 = new TCanvas("mggFit1_Signal","mggFitCB",10,10,1400,800) ;
  c3->Divide(6,3) ;

  TCanvas * cReal = new TCanvas("cReal","Mgg real events",10,10,1400,800) ;
  cReal->Divide(6,3) ;

  TCanvas * c1 = new TCanvas("mggFit1","mggFit1",10,10,1400,800) ;
  c1->Divide(6,3) ;
  c1->cd(0) ;
  TCanvas * c2=0,*c4=0,*c5=0,*c6=0 ; 

  TAxis * pta=h->GetYaxis() ;
  TAxis * ma=h->GetXaxis() ;
  for(Int_t i=1;i<=nbin;i++){
    c1->cd(i) ;
    Int_t imin=pta->FindBin(xa[i-1]+0.0001);
    Int_t imax=pta->FindBin(xa[i]-0.0001) ;
    Double_t pt=(xa[i]+xa[i-1])/2. ;
    TH1D * hp = h->ProjectionX(Form("re%d",i),imin,imax) ;
    hp->Sumw2() ;
    TH1D * hpm= hm->ProjectionX(Form("mi%d",i),imin,imax) ;
    hpm->Sumw2() ;
    if(i<15){
      hp ->Rebin(2) ;
      hpm->Rebin(2) ;
    }
    else{
      hp ->Rebin(5) ;
      hpm->Rebin(5) ;
    }
    hp ->SetNdivisions(506);
    hpm->SetNdivisions(506);
    hp ->SetLabelSize(0.05,"X");
    hpm->SetLabelSize(0.05,"X");
    hp ->SetTitleSize(0.00,"X");
    hpm->SetTitleSize(0.00,"X");
    hp ->SetLabelSize(0.05,"Y");
    hpm->SetLabelSize(0.05,"Y");
    hp ->SetTitleSize(0.05,"Y");
    hpm->SetTitleSize(0.05,"Y");
    hp ->SetTitleOffset(1.0,"X");
    hpm->SetTitleOffset(1.0,"X");
    if (i>12) {
      hp ->SetLabelSize(0.05,"X");
      hpm->SetLabelSize(0.05,"X");
      hp ->SetTitleSize(0.05,"X");
      hpm->SetTitleSize(0.05,"X");
      hp ->SetXTitle("m_{#gamma#gamma} (GeV/c^{2})");
      hpm->SetXTitle("m_{#gamma#gamma} (GeV/c^{2})");
    }

    // for(Int_t ib=1; ib<=hp->GetNbinsX();ib++){if(hp ->GetBinContent(ib)==0)hp ->SetBinError(ib,1.);}
    // for(Int_t ib=1; ib<=hp->GetNbinsX();ib++){if(hpm->GetBinContent(ib)==0)hpm->SetBinError(ib,1.);}
    TH1D * hpm2   = (TH1D*)hpm->Clone("Bg1") ;
    TH1D * hpcopy = (TH1D*)hp ->Clone("hpcopy") ;
    TH1D * hp2    = (TH1D*)hp ->Clone("hp2") ;
    TH1D * hpReal = (TH1D*)hp ->Clone("hpReal") ;
    hpcopy->Divide(hpm) ;
    sprintf(key,"%3.1f<p_{t}<%3.1f GeV/c",xa[i-1],xa[i]) ;
    hpcopy->SetTitle(key) ;
    hpReal->SetTitle(key) ;
    hpReal->SetLineWidth(2);
    hpcopy->SetMarkerStyle(20) ;
    hpcopy->SetMarkerSize(0.7) ;
    
    fit1->SetParameters(TMath::Min(0.3,0.001*i*i),0.136,0.0077,0.0013,-0.0007) ;
    fit1->SetParLimits(2,0.003,0.010) ;
    fit1->SetParLimits(1,0.132,0.139) ;
    fit1->SetParLimits(0,0.,1.) ;

    Double_t rangeMin=0.06 ;
    Double_t rangeMax=0.22 ;
    if(i>=16)rangeMax=0.3 ; //More points to fix background
    hpcopy->Fit(fit1,"NQ","",rangeMin,rangeMax) ;
    hpcopy->Fit(fit1,"MQ","",rangeMin,rangeMax) ;
    // hpcopy->Fit(fit1,"LQ","",rangeMin,rangeMax) ;
    ar1->SetBinContent(i,fit1->GetParameter(3)) ;
    ar1->SetBinError  (i,fit1->GetParError (3)) ;
    br1->SetBinContent(i,fit1->GetParameter(4)) ;
    br1->SetBinError  (i,fit1->GetParError (4)) ;

    fit2->SetParameters(fit1->GetParameters()) ;
    fit2->SetParLimits(2,0.003,0.010) ;
    fit2->SetParLimits(1,0.132,0.139) ;
    fit2->SetParLimits(0,0.,1.) ;
    fit2->SetParameter(5,0.) ;

    hpcopy->Fit(fit2,"+NQ","",rangeMin,rangeMax) ;
    hpcopy->Fit(fit2,"+MQ","",rangeMin,rangeMax) ;
    // hpcopy->Fit(fit2,"+LQ","",rangeMin,rangeMax) ;
    ar2->SetBinContent(i,fit2->GetParameter(3)) ;
    ar2->SetBinError  (i,fit2->GetParError (3)) ;
    br2->SetBinContent(i,fit2->GetParameter(4)) ;
    br2->SetBinError  (i,fit2->GetParError (4)) ;
    cr2->SetBinContent(i,fit2->GetParameter(5)) ;
    cr2->SetBinError  (i,fit2->GetParError (5)) ;

    c1->cd(i) ;
    hpcopy->SetAxisRange(0.065,0.229,"X") ;
    hpcopy->Draw() ;

    cReal->cd(i) ;
    hpReal->SetAxisRange(0.065,0.229,"X") ;
    hpReal->Draw("ehist") ;

    fbg1->SetParameters(fit1->GetParameter(3),fit1->GetParameter(4)); 
    fbg2->SetParameters(fit2->GetParameter(3),fit2->GetParameter(4),fit2->GetParameter(5)); 

    Double_t intRangeMin = PeakPosition(pt)-nSigma*PeakWidth(pt) ;
    Double_t intRangeMax = PeakPosition(pt)+nSigma*PeakWidth(pt) ;
    Int_t    intBinMin   = hp->GetXaxis()->FindBin(intRangeMin) ;
    Int_t    intBinMax   = hp->GetXaxis()->FindBin(intRangeMax) ;
    Double_t errStat     = hpm->Integral(intBinMin,intBinMax); 

    hpm ->Multiply(fbg1) ;
    hpm2->Multiply(fbg2) ;
    hp  ->Add(hpm,-1.) ;
    hp2 ->Add(hpm2,-1.) ;

    c3->cd(i) ;
    
    if(i<15)
      fgs->SetParameters(hp->Integral(32,36)/5.,fit1->GetParameter(1),fit1->GetParameter(2)) ;
    else
      fgs->SetParameters(hp->Integral(13,15)/3.,0.136,0.005) ;
    fgs->SetParLimits(2,0.003,0.010) ;
    fgs->SetParLimits(1,0.132,0.142) ;
    fgs->SetParLimits(0,0.,1.e+6) ;
    
    hp->Fit(fgs,"Q","",rangeMin,rangeMax) ;   
    hp->SetMaximum(hp2->GetMaximum()*1.1) ;
    hp->SetMinimum(hp2->GetMinimum()*1.1) ;
    hp->SetMarkerStyle(20) ;
    hp->SetMarkerSize(0.7) ;
    mr1->SetBinContent(i,fgs->GetParameter(1)) ;
    mr1->SetBinError  (i,fgs->GetParError(1)) ;
    sr1->SetBinContent(i,TMath::Abs(fgs->GetParameter(2))) ;
    sr1->SetBinError  (i,fgs->GetParError(2)) ;

    Double_t y=fgs->Integral(intRangeMin,intRangeMax)/hp->GetXaxis()->GetBinWidth(1) ;
    nr1->SetBinContent(i,y) ;
    Double_t ey=0 ;
    if(fgs->GetParameter(0)!=0. && fgs->GetParameter(2)!=0.){
      Double_t en=fgs->GetParError(0)/fgs->GetParameter(0) ;
      Double_t es=fgs->GetParError(2)/fgs->GetParameter(2) ;
      ey=y*TMath::Sqrt(en*en+es*es) ;
    }
    nr1->SetBinError(i,ey) ;

    Double_t npiInt=hp->Integral(intBinMin,intBinMax) ;
    Double_t norm=fbg1->GetParameter(0) ;
    Double_t normErr=fbg1->GetParError(0) ;
    if(npiInt>0.){
      nr1int->SetBinContent(i,npiInt) ;
      nr1int->SetBinError(i,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
    }
    hp2->GetXaxis()->SetRangeUser(0.06,0.229) ;
    hp2->SetMaximum(hp2->GetMaximum()*1.1) ;
    hp2->SetMinimum(hp2->GetMinimum()*1.1) ;
    hp2->SetMarkerStyle(20) ;
    hp2->SetMarkerSize(0.7) ;

    hp2->Fit(fgs,"Q","",rangeMin,rangeMax) ;
    mr2->SetBinContent(i,fgs->GetParameter(1)) ;
    mr2->SetBinError  (i,fgs->GetParError(1)) ;
    sr2->SetBinContent(i,TMath::Abs(fgs->GetParameter(2))) ;
    sr2->SetBinError  (i,fgs->GetParError(2)) ;
    y=fgs->Integral(intRangeMin,intRangeMax)/hp->GetXaxis()->GetBinWidth(1) ;
    nr2->SetBinContent(i,y) ;
    ey=0 ;
    if(fgs->GetParameter(0)!=0. && fgs->GetParameter(2)!=0.){
      Double_t en=fgs->GetParError(0)/fgs->GetParameter(0) ;
      Double_t es=fgs->GetParError(2)/fgs->GetParameter(2) ;
      ey=y*TMath::Sqrt(en*en+es*es) ;
    }
    nr2->SetBinError(i,ey) ;
    npiInt=hp2->Integral(intBinMin,intBinMax) ;
    norm=fbg2->GetParameter(0) ;
    normErr=fbg2->GetParError(0) ;
    if(npiInt>0.){
      nr2int->SetBinContent(i,npiInt) ;
      nr2int->SetBinError(i,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
    } 
    hp2->SetTitle(key) ;
    hp2->SetAxisRange(0.065,0.229,"X") ;
    hp2->Draw() ;
    hp ->Draw("same") ;
    delete hpm ;
    delete hpm2 ;
    c1->Update() ;
    c3->Update() ;
    cReal->Update() ;
  }

  if (c3) c3->Print("Pi0_Signal_Gaus.eps") ;
  if (c1) c1->Print("Pi0_Ratio_Gaus.eps") ;
  if (cReal) cReal->Print("Pi0_InvMass_Gaus.eps") ;

  //Normalize by the number of events
  Double_t nMBOR   = hev->GetBinContent(2); // MBOR events pass4
  Double_t nevents = nMBOR;
  printf("==============\nN events = %d\n==============\n",nevents);

  nr1   ->Scale(1./nevents) ;
  nr1int->Scale(1./nevents) ;
  nr2   ->Scale(1./nevents) ;
  nr2int->Scale(1./nevents) ;


  TFile fout("LHC11a_FitResult.root","update");
  mr1->Write(0,TObject::kOverwrite) ;
  sr1->Write(0,TObject::kOverwrite) ;
  ar1->Write(0,TObject::kOverwrite) ;
  br1->Write(0,TObject::kOverwrite) ;
  nr1->Write(0,TObject::kOverwrite) ;
  nr1int->Write(0,TObject::kOverwrite) ;
  ar2->Write(0,TObject::kOverwrite) ;
  br2->Write(0,TObject::kOverwrite) ;
  cr2->Write(0,TObject::kOverwrite) ;
  mr2->Write(0,TObject::kOverwrite) ;
  sr2->Write(0,TObject::kOverwrite) ;
  nr2->Write(0,TObject::kOverwrite) ;
  nr2int->Write(0,TObject::kOverwrite) ;
  fout.Close() ;

}

//-----------------------------------------------------------------------------
const Double_t kMean=0.135 ; //Approximate peak position to facilitate error estimate

Double_t PeakPosition(Double_t pt){
  //Fit to the measured peak position
  return 0.454*exp(-pt*7.13)+0.13569 ;
}
Double_t PeakWidth(Double_t pt){
  //fit to the measured peak width
  return 1.69e-02*exp(-pt*2.82)+4.93e-03 ;
}
 
Double_t CB(Double_t * x, Double_t * par){
  //Parameterization of Real/Mixed ratio
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+par[3]+par[4]*(x[0]-kMean) ;
}
Double_t CB2(Double_t * x, Double_t * par){
  //Another parameterization of Real/Mixed ratio
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+par[3]+par[4]*(x[0]-kMean)+par[5]*(x[0]-kMean)*(x[0]-kMean) ;
}
Double_t CBs(Double_t * x, Double_t * par){
  //Parameterizatin of signal
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.) ;
}
Double_t BG1(Double_t * x, Double_t * par){
  //Normalizatino of Mixed
  return par[0]+par[1]*(x[0]-kMean) ;
}
Double_t BG2(Double_t * x, Double_t * par){
  //Another normalization of  Mixed
  return par[0]+par[1]*(x[0]-kMean)+par[2]*(x[0]-kMean)*(x[0]-kMean) ;
}


//-----------------------------------------------------------------------------
PPRstyle()
{

  //////////////////////////////////////////////////////////////////////
  //
  // ROOT style macro for the TRD TDR
  //
  //////////////////////////////////////////////////////////////////////

  gStyle->SetPalette(1);
  gStyle->SetCanvasBorderMode(-1);
  gStyle->SetCanvasBorderSize(1);
  gStyle->SetCanvasColor(10);

  gStyle->SetFrameFillColor(10);
  gStyle->SetFrameBorderSize(1);
  gStyle->SetFrameBorderMode(-1);
  gStyle->SetFrameLineWidth(1.2);
  gStyle->SetFrameLineColor(1);

  gStyle->SetHistFillColor(0);
  gStyle->SetHistLineWidth(1);
  gStyle->SetHistLineColor(1);

  gStyle->SetPadColor(10);
  gStyle->SetPadBorderSize(1);
  gStyle->SetPadBorderMode(-1);

  gStyle->SetStatColor(10);
  gStyle->SetTitleColor(kBlack,"X");
  gStyle->SetTitleColor(kBlack,"Y");

  gStyle->SetLabelSize(0.04,"X");
  gStyle->SetLabelSize(0.04,"Y");
  gStyle->SetLabelSize(0.04,"Z");
  gStyle->SetTitleSize(0.04,"X");
  gStyle->SetTitleSize(0.04,"Y");
  gStyle->SetTitleSize(0.04,"Z");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelFont(42,"Z");
  gStyle->SetStatFont(42);

  gStyle->SetTitleOffset(1.0,"X");
  gStyle->SetTitleOffset(1.4,"Y");

  gStyle->SetFillColor(kWhite);
  gStyle->SetTitleFillColor(kWhite);

  gStyle->SetOptDate(0);
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

}

