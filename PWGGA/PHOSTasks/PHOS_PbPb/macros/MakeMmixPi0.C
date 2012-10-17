/* $Id$ */

const Double_t kMean=0.136 ; //Approximate peak position to facilitate error estimate

//-----------------------------------------------------------------------------
void MakeMmixPi0(const TString filename = "Pi0Flow_000167920.root",
		 const TString listPath = "AliPHOSPi0Flow/PHOSPi0FlowCoutput1",
		 const Int_t centrality=0, 
		 const char* pid="CPV")
{
  //Fit Real/Mixed ratio, normalize Mixed and subtract it from Real

  TFile * f = new TFile(filename) ;
  //TList *histoList = (TList*)f->Get("PHOSPi0Flow");
  TList *histoList = (TList*)f->Get("PHOSPi0Flow/PHOSPi0FlowCoutput1"); // lego train

  char key[125] ;

  TH1F * hev = (TH1F*)histoList->FindObject("hTotSelEvents") ;
  TH2F * hCentrality  = (TH2F*)histoList->FindObject("hCenPHOSCells") ;
  TH1D * hCentrality1 = hCentrality->ProjectionX();  
  
  printf("TotSelEvents: %f \n",hev->GetBinContent(4)) ;
  printf("Centrality:   %f \n",hCentrality1->Integral()) ;
  
  TString inputKey;
  TString outputKey = Form("%s_cent%d",pid,centrality);

  TH2F *h , *hAdd;
  TH2F *hm, *hmAdd;

  if (centrality >= 0 && centrality < 4) {
    printf("\tCentrality %d\n",centrality);
    inputKey = Form("hPi0%s_cen%d"  ,pid,centrality);
    TH2F *h = (TH2F*)histoList->FindObject(inputKey) ;
    inputKey = Form("hMiPi0%s_cen%d",pid,centrality);
    TH2F *hm= (TH2F*)histoList->FindObject(inputKey) ;
    if (h==0) {
      printf("Missing histogram %s\n",inputKey);
      return;
    }
  }
  else {
    printf("Wrong centrality %d. Allowed values are 0,1,2,3,10.\n",centrality);
    return;
  }

  Int_t nPadX = 3, nPadY = 2;
  Int_t nPtBins=6 ;
  Double_t ptBinEdges[21]={1., 2., 3., 4., 5., 7., 10.} ;

  PPRstyle();
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadBottomMargin(0.08);

  //Fit real only 
  //Linear Bg
  char kkey[55];
  sprintf(kkey,outputKey.Data()) ;
  char key2[155];
  sprintf(key,"Mix%s",kkey) ;
  sprintf(key2,"%s_mr1",key) ;
  TH1D * mr1 = new TH1D(key2,"Mass",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_sr1",key) ;
  TH1D * sr1 = new TH1D(key2,"Width",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_ar1",key) ;
  TH1D * ar1 = new TH1D(key2,"a",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_br1",key) ;
  TH1D * br1 = new TH1D(key2,"a",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_yr1",key) ;
  TH1D * nr1 = new TH1D(key2,"Raw yield",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_yr1int",key) ;
  TH1D * nr1int = new TH1D(key2,"Raw yield, integrated",nPtBins,ptBinEdges) ;

  //Quadratic Bg
  sprintf(key2,"%s_mr2",key) ;
  TH1D * mr2 = new TH1D(key2,"Mass",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_sr2",key) ;
  TH1D * sr2 = new TH1D(key2,"Width",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_ar2",key) ;
  TH1D * ar2 = new TH1D(key2,"a",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_br2",key) ;
  TH1D * br2 = new TH1D(key2,"a",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_cr2",key) ;
  TH1D * cr2 = new TH1D(key2,"a",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_yr2",key) ;
  TH1D * nr2 = new TH1D(key2,"Raw yield",nPtBins,ptBinEdges) ;
  sprintf(key2,"%s_yr2int",key) ;
  TH1D * nr2int = new TH1D(key2,"Raw yield, integrated",nPtBins,ptBinEdges) ;

  TF1 * fit1 = new TF1("fit",CB,0.,1.,6) ;
  fit1->SetParName(0,"A") ;
  fit1->SetParName(1,"m_{0}") ;
  fit1->SetParName(2,"#sigma") ;
  fit1->SetParName(3,"a_{0}") ;
  fit1->SetParName(4,"a_{1}") ;
  fit1->SetParName(5,"a_{2}") ;
  fit1->SetLineWidth(2) ;
  fit1->SetLineColor(2) ;
  TF1 * fgs = new TF1("gs",CBs,0.,1.,4) ;
  fgs->SetParName(0,"A") ;
  fgs->SetParName(1,"m_{0}") ;
  fgs->SetParName(2,"#sigma") ;
  fgs->SetParName(3,"B") ;
  fgs->SetLineColor(2) ;
  fgs->SetLineWidth(1) ;

  TF1 * fit2 = new TF1("fit2",CB2,0.,1.,7) ;
  fit2->SetParName(0,"A") ;
  fit2->SetParName(1,"m_{0}") ;
  fit2->SetParName(2,"#sigma") ;
  fit2->SetParName(3,"a_{0}") ;
  fit2->SetParName(4,"a_{1}") ;
  fit2->SetParName(5,"a_{2}") ;
  fit2->SetParName(6,"a_{3}") ;
  fit2->SetLineWidth(2) ;
  fit2->SetLineColor(4) ;
  fit2->SetLineStyle(2) ;

  TF1 * fbg1 = new TF1("bg1",BG1,0.,1.,3) ;
  TF1 * fbg2 = new TF1("bg2",BG2,0.,1.,4) ;

  TCanvas * c3 = new TCanvas("mggFit1_Signal","mggFitCB",10,10,1200,800) ;
  c3->Divide(nPadX,nPadY) ;

  TCanvas * c1 = new TCanvas("mggFit1","mggFit1",10,10,1200,800) ;
  c1->Divide(nPadX,nPadY) ;
  c1->cd(0) ;

  TCanvas * rawCanvas = new TCanvas("rawCanvas","rawCanvas",10,10,1200,800);
  rawCanvas->Divide(nPadX, nPadY);

  TAxis * pta=h->GetYaxis() ;
  TAxis * ma=h->GetXaxis() ;
  for(Int_t ptBin=1; ptBin<=nPtBins; ptBin++){
    c1->cd(ptBin) ;
    printf("\n\t%.1f<pt<%.1f GeV/c\n",ptBinEdges[ptBin-1],ptBinEdges[ptBin]);
    Int_t imin=pta->FindBin(ptBinEdges[ptBin-1]+0.0001);
    Int_t imax=pta->FindBin(ptBinEdges[ptBin]-0.0001) ;
    Double_t pt=(ptBinEdges[ptBin]+ptBinEdges[ptBin-1])/2. ;
    TH1D * hp = h->ProjectionX(Form("re_%d",ptBin),imin,imax) ;
    hp->Sumw2() ;
    TH1D * hpm= hm->ProjectionX("mi",imin,imax) ;
    hpm->Sumw2() ;
    if(ptBin<=17){
      hp ->Rebin(4) ;
      hpm->Rebin(4) ;
    }
    else{
      hp ->Rebin(5) ;
      hpm->Rebin(5) ;
    }
    for(Int_t ib=1; ib<=hp->GetNbinsX();ib++){if(hp ->GetBinContent(ib)==0)hp ->SetBinError(ib,1.);}
    for(Int_t ib=1; ib<=hp->GetNbinsX();ib++){if(hpm->GetBinContent(ib)==0)hpm->SetBinError(ib,1.);}
    TH1D * hpm2   = (TH1D*)hpm->Clone("Bg1") ;
    TH1D * hpmScaled = (TH1D*)hpm->Clone("hpmScaled") ;
    TH1D * hpcopy = (TH1D*)hp ->Clone("hpcopy") ;
    TH1D * hp2    = (TH1D*)hp ->Clone("hp2") ;
    hpcopy->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
    hp2   ->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
    hpcopy->Divide(hpm) ;
    sprintf(key,"PID=%s, %3.1f<p_{T}<%3.1f GeV/c",pid,ptBinEdges[ptBin-1],ptBinEdges[ptBin]) ;
    hpcopy->SetTitle(key) ;
    hpcopy->SetMarkerStyle(20) ;
    hpcopy->SetMarkerSize(0.7) ;
    
    fit1->SetParameters(0.001,0.136,0.0055,0.0002,-0.002,0.0) ;
    fit1->SetParLimits(0,0.000,1.000) ;
    fit1->SetParLimits(1,0.120,0.145) ;
    fit1->SetParLimits(2,0.005,0.012) ;

    Double_t rangeMin=0.05 ;
    Double_t rangeMax=0.3 ;
    if(centrality==0) rangeMax=0.4 ;
    if(ptBin==1){
      rangeMin=0.06 ;
      rangeMax=0.25 ;
    }
    hpcopy->Fit(fit1,"Q" ,"",rangeMin,rangeMax) ;
    hpcopy->Fit(fit1,"MQ","",rangeMin,rangeMax) ;

    ar1->SetBinContent(ptBin,fit1->GetParameter(3)) ;
    ar1->SetBinError  (ptBin,fit1->GetParError(3)) ;
    br1->SetBinContent(ptBin,fit1->GetParameter(4)) ;
    br1->SetBinError  (ptBin,fit1->GetParError(4)) ;

    fit2->SetParameters(fit1->GetParameters()) ;
    fit2->SetParLimits(0,0.000,1.000) ;
    fit2->SetParLimits(1,0.120,0.145) ;
    fit2->SetParLimits(2,0.005,0.012) ;

    hpcopy->Fit(fit2,"+NQ","",rangeMin,rangeMax) ;
    hpcopy->Fit(fit2,"+MQ","",rangeMin,rangeMax) ;

    ar2->SetBinContent(ptBin,fit2->GetParameter(3)) ;
    ar2->SetBinError  (ptBin,fit2->GetParError(3)) ;
    br2->SetBinContent(ptBin,fit2->GetParameter(4)) ;
    br2->SetBinError  (ptBin,fit2->GetParError(4)) ;
    cr2->SetBinContent(ptBin,fit2->GetParameter(5)) ;
    cr2->SetBinError  (ptBin,fit2->GetParError(5)) ;
    hpcopy->GetXaxis()->SetRangeUser(0.05,0.30) ;
    hpcopy->Draw() ;
    c1->Update() ;

    fbg1->SetParameters(fit1->GetParameter(3),fit1->GetParameter(4),fit1->GetParameter(5)); 
    fbg2->SetParameters(fit2->GetParameter(3),fit2->GetParameter(4),fit2->GetParameter(5),
			fit2->GetParameter(6)); 

    Double_t intRangeMin = PeakPosition(pt)-3.*PeakWidth(pt) ;
    Double_t intRangeMax = PeakPosition(pt)+3.*PeakWidth(pt) ;
    Int_t    intBinMin   = hp->GetXaxis()->FindBin(intRangeMin) ;
    Int_t    intBinMax   = hp->GetXaxis()->FindBin(intRangeMax) ;
    Double_t errStat     = hpm->Integral(intBinMin,intBinMax); 

    hpm ->Multiply(fbg1) ;
    hpm2->Multiply(fbg2) ;
    hp  ->Add(hpm ,-1.) ;
    hp2 ->Add(hpm2,-1.) ;

    c3->cd(ptBin) ;

    Int_t binPi0 = hp->FindBin(kMean);
    fgs->SetParameters(hp->Integral(binPi0-1,binPi0+1)/3.,fit1->GetParameter(1),fit1->GetParameter(2)) ;
    fgs->SetParLimits(0,0.000,1.e+5) ;
    fgs->SetParLimits(1,0.120,0.145) ;
    fgs->SetParLimits(2,0.004,0.02) ;
    hp->Fit(fgs,"Q","",rangeMin,rangeMax) ;   
    hp->SetMaximum(hp2->GetMaximum()*1.5) ;
    hp->SetMinimum(hp2->GetMinimum()*1.1) ;
    hp->SetMarkerStyle(20) ;
    hp->SetMarkerSize(0.7) ;
    mr1->SetBinContent(ptBin,fgs->GetParameter(1)) ;
    mr1->SetBinError  (ptBin,fgs->GetParError(1) ) ;
    sr1->SetBinContent(ptBin,TMath::Abs(fgs->GetParameter(2))) ;
    sr1->SetBinError  (ptBin,fgs->GetParError(2) ) ;

    Double_t y=fgs->GetParameter(0)/hp->GetXaxis()->GetBinWidth(1) ;
    nr1->SetBinContent(ptBin,y) ;
    Double_t ey=fgs->GetParError(0)/hp->GetXaxis()->GetBinWidth(1) ;
    nr1->SetBinError(ptBin,ey) ;

    Double_t npiInt = hp->Integral(intBinMin,intBinMax) ;
    Double_t norm   = fbg1->GetParameter(0) ;
    Double_t normErr= fbg1->GetParError(0) ;
    if(npiInt>0.){
      nr1int->SetBinContent(ptBin,npiInt) ;
      nr1int->SetBinError(ptBin,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
    }
    hp2->GetXaxis()->SetRangeUser(0.05,0.3) ;
    hp2->SetMaximum(hp2->GetMaximum()*1.5) ;
    hp2->SetMinimum(hp2->GetMinimum()*1.1) ;
    hp2->SetMarkerStyle(20) ;
    hp2->SetMarkerSize(0.7) ;
    hp2->Fit(fgs,"Q","",rangeMin,rangeMax) ;
    mr2->SetBinContent(ptBin,fgs->GetParameter(1)) ;
    mr2->SetBinError  (ptBin,fgs->GetParError(1)) ;
    sr2->SetBinContent(ptBin,TMath::Abs(fgs->GetParameter(2))) ;
    sr2->SetBinError  (ptBin,fgs->GetParError(2)) ;
    y=fgs->GetParameter(0)/hp->GetXaxis()->GetBinWidth(1) ;
    nr2->SetBinContent(ptBin,y) ;
    ey= fgs->GetParError(0)/hp->GetXaxis()->GetBinWidth(1) ;
    nr2->SetBinError(ptBin,ey) ;
    npiInt=hp2->Integral(intBinMin,intBinMax) ;
    norm=fbg2->GetParameter(0) ;
    normErr=fbg2->GetParError(0) ;
    if(npiInt>0.){
      nr2int->SetBinContent(ptBin,npiInt) ;
      nr2int->SetBinError(ptBin,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
    } 
    hp2->SetTitle(key) ;
    hp2->Draw() ;
    c3->Update() ;

    delete hp ;
    delete hpm ;
    delete hpm2 ;
  }
  char name[55] ;
  
  sprintf(name,"Pi0_ratio_%s.eps",kkey) ;
  c1->Print(name) ;
  sprintf(name,"Pi0_signal_%s.eps",kkey) ;
  c3->Print(name) ;

  //Normalize by the number of events
  Int_t cMin,cMax;
  if      (centrality == 0) {
    cMin=1;
    cMax=10;
  }
  else if (centrality == 1) {
    cMin=21;
    cMax=40;
  }
  else if (centrality == 2) {
    cMin=41;
    cMax=60;
  }
  else if (centrality == 3) {
    cMin=61;
    cMax=80;
  }
  Double_t nevents = hCentrality1->Integral(cMin,cMax);
  nr1   ->Scale(1./nevents) ;
  nr1int->Scale(1./nevents) ;
  nr2   ->Scale(1./nevents) ;
  nr2int->Scale(1./nevents) ;
  printf("Nevents=%f \n",nevents) ;
  
  
  TFile fout("LHC10h_Pi0_FitResult.root","update");
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
Double_t PeakPosition(Double_t pt){
  //Fit to the measured peak position
  return 4.99292e-003*exp(-pt*9.32300e-001)+1.34944e-001 ;
}
//-----------------------------------------------------------------------------
Double_t PeakWidth(Double_t pt){
  //fit to the measured peak width
  Double_t a=0.0068 ;
  Double_t b=0.0025 ;
  Double_t c=0.000319 ;
  return TMath::Sqrt(a*a+b*b/pt/pt+c*c*pt*pt) ;
}
 
//-----------------------------------------------------------------------------
Double_t CB(Double_t * x, Double_t * par){
  //Parameterization of Real/Mixed ratio
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+par[3]+par[4]*(x[0]-kMean) ;
}

//-----------------------------------------------------------------------------
Double_t CB2(Double_t * x, Double_t * par){
  //Another parameterization of Real/Mixed ratio7TeV/README
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+par[3]+par[4]*(x[0]-kMean)+par[5]*(x[0]-kMean)*(x[0]-kMean) ;
}
//-----------------------------------------------------------------------------
Double_t CBs(Double_t * x, Double_t * par){
  //Parameterizatin of signal
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)/TMath::Sqrt(TMath::TwoPi())/s+par[3] ;
}
//-----------------------------------------------------------------------------
Double_t BG1(Double_t * x, Double_t * par){
  //Normalizatino of Mixed
  return par[0]+par[1]*(x[0]-kMean) ;
}
//-----------------------------------------------------------------------------
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
  gStyle->SetOptFit(111);

}

