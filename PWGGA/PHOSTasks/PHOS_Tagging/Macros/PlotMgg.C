const Double_t kMean=0.136 ; //Approximate peak position to facilitate error estimate
//-----------------------------------------------------------------------------
Double_t CB1(Double_t * x, Double_t * par){
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t n=par[3] ;
  Double_t a=par[4] ;
  Double_t dx=(x[0]-m)/s ;
  if(dx>-a)
    return par[0]*exp(-dx*dx/2.)+
    par[5]+
    par[6]*(x[0]-kMean)+
    par[7]*(x[0]-kMean)*(x[0]-kMean) ;
  else{
    Double_t A=TMath::Power((n/TMath::Abs(a)),n)*TMath::Exp(-a*a/2) ;
    Double_t B=n/TMath::Abs(a)-TMath::Abs(a) ;
    return par[0]*A*TMath::Power((B-dx),-n)+
    par[5]+
    par[6]*(x[0]-kMean)+
    par[7]*(x[0]-kMean)*(x[0]-kMean) ;
  }

}
Double_t CB2(Double_t * x, Double_t * par){
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t n=par[3] ;
  Double_t a=par[4] ;
  Double_t dx=(x[0]-m)/s ;
  if(dx>-a)
    return par[0]*exp(-dx*dx/2.)+
    par[5]+
    par[6]*(x[0]-kMean)+
    par[7]*(x[0]-kMean)*(x[0]-kMean) +
    par[8]*(x[0]-kMean)*(x[0]-kMean)*(x[0]-kMean) ;
  else{
    Double_t A=TMath::Power((n/TMath::Abs(a)),n)*TMath::Exp(-a*a/2) ;
    Double_t B=n/TMath::Abs(a)-TMath::Abs(a) ;
    return par[0]*A*TMath::Power((B-dx),-n)+
    par[5]+
    par[6]*(x[0]-kMean)+
    par[7]*(x[0]-kMean)*(x[0]-kMean) +
    par[8]*(x[0]-kMean)*(x[0]-kMean)*(x[0]-kMean) ;
  }

}

Double_t CBs(Double_t * x, Double_t * par){
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t n=par[3] ;
  Double_t a=par[4] ;
  Double_t dx=(x[0]-m)/s ;
  if(dx>-a)
    return par[0]*exp(-dx*dx/2.) + par[5] ;
  else{
    Double_t A=TMath::Power((n/TMath::Abs(a)),n)*TMath::Exp(-a*a/2) ;
    Double_t B=n/TMath::Abs(a)-TMath::Abs(a) ;
    return par[0]*A*TMath::Power((B-dx),-n)  + par[5] ;
  }
}
//-----------------------------------------------------------------------------
Double_t GS1(Double_t * x, Double_t * par){
  //Parameterization of Real/Mixed ratio
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+
  par[3]+
  par[4]*(x[0]-kMean)+
  par[5]*(x[0]-kMean)*(x[0]-kMean);
}

//-----------------------------------------------------------------------------
Double_t GS2(Double_t * x, Double_t * par){
  //Another parameterization of Real/Mixed ratio7TeV/README
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+
  par[3]+
  par[4]*(x[0]-kMean)+
  par[5]*(x[0]-kMean)*(x[0]-kMean)+
  par[6]*(x[0]-kMean)*(x[0]-kMean)*(x[0]-kMean)  ;
}
//-----------------------------------------------------------------------------
Double_t GSs(Double_t * x, Double_t * par){
  //Parameterizatin of signal
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)/TMath::Sqrt(TMath::TwoPi())/s+par[3] ;
}
//-----------------------------------------------------------------------------
Double_t BG(Double_t * x, Double_t * par){
  //Normalizatino of Mixed
  return par[0]+par[1]*(x[0]-kMean) ;
}
//-----------------------------------------------------------------------------
Double_t BG1(Double_t * x, Double_t * par){
  //Normalizatino of Mixed
  return par[0]+par[1]*(x[0]-kMean)+par[2]*(x[0]-kMean)*(x[0]-kMean) ;
}
//-----------------------------------------------------------------------------
Double_t BG2(Double_t * x, Double_t * par){
  //Another normalization of  Mixed
  return par[0]+par[1]*(x[0]-kMean)+par[2]*(x[0]-kMean)*(x[0]-kMean)+par[3]*(x[0]-kMean)*(x[0]-kMean)*(x[0]-kMean) ;
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



void PlotMgg(Double_t ptMin=2.8, Double_t ptMax=3.0){

    
  TFile fout("mgg.root","update") ;
  TFile * f = new TFile("LHC13bc_MB.root") ; 
  

  TH1F * hev = (TH1F*)f->Get("hSelEvents") ;
  TH1F * hCentrality1  = (TH1F*)f->Get("hCentrality") ;
  
  printf("TotSelEvents: %f \n",hev->GetBinContent(7)) ;
  printf("Centrality:   %f \n",hCentrality1->Integral()) ;


  TH2F * hRe ;
  TH2F * hMi ;
  TH2F * tmp = 0x0 ;
	hRe = (TH2F*)f->Get(Form("hInvM_Re_Emin3_All_cent0")) ;
        tmp = (TH2F*)f->Get(Form("hInvM_Re_Emin3_All_cent1")) ;
        hRe->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Re_Emin3_All_cent2")) ;
        hRe->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Re_Emin3_All_cent3")) ;
        hRe->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Re_Emin3_All_cent4")) ;
        hRe->Add(tmp) ;delete tmp ;
        hRe->Sumw2() ;	
	
	hMi = (TH2F*)f->Get(Form("hInvM_Mi_Emin3_All_cent0")) ;
        tmp = (TH2F*)f->Get(Form("hInvM_Mi_Emin3_All_cent1")) ;
        hMi->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Mi_Emin3_All_cent2")) ;
        hMi->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Mi_Emin3_All_cent3")) ;
        hMi->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Mi_Emin3_All_cent4")) ;
        hMi->Add(tmp) ;delete tmp ;
        hMi->Sumw2() ;	
      
  PPRstyle();
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadBottomMargin(0.08);


  TF1 * fun1 = new TF1("ft1",CB1,0.,1.,8) ;
  fun1->SetParName(0,"A") ;
  fun1->SetParName(1,"m_{0}") ;
  fun1->SetParName(2,"#sigma") ;
  fun1->SetParName(3,"n") ;
  fun1->SetParName(4,"#alpha") ;
  fun1->SetParName(5,"a_{0}") ;
  fun1->SetParName(6,"a_{1}") ;
  fun1->SetParName(7,"a_{2}") ;
  fun1->FixParameter(3,3.28) ;
  fun1->FixParameter(4,1.56) ;
  fun1->SetLineWidth(2) ;
  fun1->SetLineColor(4) ;
  fun1->SetLineStyle(2) ;

  TF1 * fun2 = new TF1("ft2",CB2,0.,1.,9) ;
  fun2->SetParName(0,"A") ;
  fun2->SetParName(1,"m_{0}") ;
  fun2->SetParName(2,"#sigma") ;
  fun2->SetParName(3,"n") ;
  fun2->SetParName(4,"#alpha") ;
  fun2->SetParName(5,"a_{0}") ;
  fun2->SetParName(6,"a_{1}") ;
  fun2->SetParName(7,"a_{2}") ;
  fun2->SetParName(8,"a_{3}") ;
  fun2->FixParameter(3,3.28) ;
  fun2->FixParameter(4,1.56) ;
  fun2->SetLineWidth(2) ;
  fun2->SetLineColor(4) ;
  fun2->SetLineStyle(2) ;
  
  TF1 * funGS1 = new TF1("ft1",GS1,0.,1.,6) ;
  funGS1->SetParName(0,"A") ;
  funGS1->SetParName(1,"m_{0}") ;
  funGS1->SetParName(2,"#sigma") ;
  funGS1->SetParName(3,"a_{0}") ;
  funGS1->SetParName(4,"a_{1}") ;
  funGS1->SetParName(5,"a_{2}") ;
  funGS1->SetLineWidth(2) ;
  funGS1->SetLineColor(8) ;
//   funGS1->SetLineStyle(2) ;

  TF1 * funGS2 = new TF1("ft2",GS2,0.,1.,7) ;
  funGS2->SetParName(0,"A") ;
  funGS2->SetParName(1,"m_{0}") ;
  funGS2->SetParName(2,"#sigma") ;
  funGS2->SetParName(3,"a_{0}") ;
  funGS2->SetParName(4,"a_{1}") ;
  funGS2->SetParName(5,"a_{2}") ;
  funGS2->SetParName(6,"a_{3}") ;
  funGS2->SetLineWidth(2) ;
  funGS2->SetLineColor(8) ;
//   funGS2->SetLineStyle(2) ;
  
  
  
  TF1 * fit1=0x0 ;
  TF1 * fit2=0x0 ;
  
  TF1 * fbgP0 = new TF1("bg",BG,0.,1.,2) ;
  TF1 * fbgP1 = new TF1("bg1",BG1,0.,1.,3) ;
  TF1 * fbgP2 = new TF1("bg2",BG2,0.,1.,4) ;
  TF1 * fbg1=0x0 ;
  TF1 * fbg2=0x0 ;
  TF1 * fgs = new TF1("gs",GSs,0.,1.,4) ;
  TF1 * fcb = new TF1("cb",CBs,0.,1.,6) ;
  fgs->SetLineColor(8) ;
  fgs->SetLineWidth(2) ;
  fcb->SetLineColor(4) ;
  fcb->SetLineWidth(2) ;
  fcb->SetLineStyle(2) ;
  fcb->FixParameter(3,3.28) ;
  fcb->FixParameter(4,1.56) ;

  TCanvas * c1 = new TCanvas(Form("Ratio"),Form("Ratio"),10,10,1200,800) ;

  TAxis * pta=hRe->GetYaxis() ;
  TAxis * ma=hRe->GetXaxis() ;
    
  Int_t imin=pta->FindBin(ptMin+0.0001);
  Int_t imax=pta->FindBin(ptMax-0.0001) ;
  TH1D *hReal = hRe->ProjectionX(Form("real"),imin,imax) ;
  TH1D *hBg= hMi->ProjectionX(Form("mixed"),imin,imax) ;
//      hpm->Sumw2() ;
//     hp ->Rebin(2) ;
//     hpm->Rebin(2) ;

    if(ptMax<2.){ //2pol,3pol
      fit1=funGS2 ;
      fit2=fun2 ;
      fbg1=fbgP2 ;
      fbg2=fbgP2 ;
      
    }
    else{ //pol1,pol2
      fit1=funGS1 ;
      fit2=fun1 ;
      fbg1=fbgP1 ;
      fbg2=fbgP1 ;
      
    }
    
    TH1D * hBgNoPol   = (TH1D*)hBg->Clone(Form("BgPolCorrection")) ;
    TH1D * hpcopy = (TH1D*)hReal ->Clone(Form("hpcopy")) ;
    TH1D * hSignal    = (TH1D*)hReal ->Clone(Form("Signal")) ;
    hpcopy->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
    hSignal   ->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
    hpcopy->Divide(hBg) ;
    hpcopy->SetTitle(Form("%3.1f<p_{T}<%3.1f GeV/c",ptMin,ptMax)) ;
    hpcopy->SetMarkerStyle(20) ;
    hpcopy->SetMarkerSize(0.7) ;
    hpcopy->GetXaxis()->SetRangeUser(0.05,0.25) ;
      
      fit1->SetParameters(0.0002,0.136,0.005,0.0002,-0.002,0.0) ;
      fit1->SetParLimits(0,0.000,2.*hpcopy->GetMaximum()) ;
      fit1->SetParLimits(1,0.125,0.145) ;
      fit1->SetParLimits(2,0.0045,0.010) ;

      Double_t rangeMin=0.05 ;//TMath::Max(0.06,0.11-0.01*i) ;
      Double_t rangeMax=0.25; //TMath::Min(0.25,0.18+0.01*i) ;
//      Double_t rangeMin=TMath::Max(0.06,0.12-0.01*i) ;
//      Double_t rangeMax=TMath::Min(0.25,0.16+0.01*i) ;

      hpcopy->Fit(fit1,"Q","",rangeMin,rangeMax) ;
      hpcopy->Fit(fit1,"MQ","",rangeMin,rangeMax) ;
      
      fbg1->SetParameters(fit1->GetParameter(3),fit1->GetParameter(4),fit1->GetParameter(5),fit1->GetParameter(6)); 
      
      hBg ->Multiply(fbg1) ;
      hBgNoPol->Scale(fbg1->Eval(0.135)) ;
      hSignal  ->Add(hBg ,-1.) ;
                
      
//   Double_t nevents = hCentrality1->Integral(1,100);
//   printf("Nevents=%f \n",nevents) ;
//   hSignal->Scale(1./nevents) ;
//   hReal->Scale(1./nevents) ;
//   hBg->Scale(1./nevents) ;
//   hBgNoPol->Scale(1./nevents) ;
  
  fgs->SetParameters(60.,0.135,0.005,1.) ;
  fgs->FixParameter(2,0.0056) ;
  fgs->SetLineColor(2) ;
  hSignal->Fit(fgs,"L","",0.08,0.22);
  
  hReal->SetTitle(Form("%3.1f<p_{T}<%3.1f GeV/c",ptMin,ptMax)) ;
  hBg->SetTitle(Form("%3.1f<p_{T}<%3.1f GeV/c",ptMin,ptMax)) ;
  hBgNoPol->SetTitle(Form("%3.1f<p_{T}<%3.1f GeV/c",ptMin,ptMax)) ;
  hSignal->SetTitle(Form("%3.1f<p_{T}<%3.1f GeV/c",ptMin,ptMax)) ;
  hReal->GetXaxis()->SetRangeUser(0.05,0.25) ;
  hReal->SetLineColor(1) ;
  hReal->SetLineWidth(2) ;
  hReal->Draw("h") ;
  hBg->SetLineColor(4) ;
  hBg->SetLineWidth(2) ;
  hBg->Draw("same") ;
//   hBgNoPol->Draw("same") ;
  hSignal->SetLineWidth(2) ;
  hSignal->SetLineColor(2) ;
  hSignal->Draw("same") ;
  
  
  TLegend * ll = new TLegend(0.5,0.6,0.9,0.9) ;
  ll->AddEntry(hReal,"Real","l") ;
  ll->AddEntry(hSignal,"Signal","l") ;
  ll->AddEntry(hBg,"Normalized Mixed","l") ;
  ll->Draw() ;
  
  fout.cd() ;
  hReal->Write(0,TObject::kOverwrite) ;
  hSignal->Write(0,TObject::kOverwrite) ;
  hBg->Write(0,TObject::kOverwrite) ;
  hBgNoPol->Write(0,TObject::kOverwrite) ;
  fout.Close() ;



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

