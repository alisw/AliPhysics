const Double_t kMean=0.135 ; //Approximate peak position to facilitate error estimate

//-----------------------------------------------------------------------------
void MakeMmixPi0(const Int_t centrality=0,
		 const char* cModule="")
{
  //---------------------------------------------------------------------------
  // This macro processes PWG1 QA output of the analysis task PHOSPbPbQA
  // (see analysis code in the class AliAnalysisTaskPHOSPbPbQA).
  // It fits Real/Mixed ratio, normalize Mixed and subtract it from Real
  // Arguments:
  // centrality: 0 or 1 for centralities 0-20% or 20-100%
  // cModule: string with the PHOS module to analyze: "", "SM1", "SM2" or "SM3"
  //-
  // Yuri Kharlov. 19.11.2011
  //---------------------------------------------------------------------------
/* $Id$ */

  TFile * f = new TFile("PHOSPbPb_all.root") ;
  TList *histESD = (TList*) f->Get("PHOSPbPbQAResults");
  char key[125] ;

  const char* pid="All";
  const char* txtCent[] = {"0-20%","20-100%"};
  TH2F * hCentrality  = (TH2F*)f->Get("hCenPHOS") ;
  TH1D * hCentrality1 = hCentrality->ProjectionX();
  TString inputKey;
  TString outputKey = Form("%s10_cent%d",pid,centrality);

  TH2F *h , *hAdd;
  TH2F *hm, *hmAdd;
  if (centrality == 0 || centrality == 1) {// centrality 0-20%, 20-100%
    inputKey = Form("hPi0%s%s_cen%d"  ,pid,cModule,centrality);
    TH2F *h = (TH2F*)f->Get(inputKey) ;
    inputKey = Form("hMiPi0%s%s_cen%d",pid,cModule,centrality);
    TH2F *hm= (TH2F*)f->Get(inputKey) ;
    if (h==0) {
      printf("Missing histogram %s\n",inputKey);
      return;
    }
  }
  else {
    printf("Wrong centrality %d. Allowed values are 0,1,2,3,10.\n",centrality);
    return;
  }

  // Int_t nPadX=2,nPadY=2;
  // Int_t nbin=4 ;
  // Double_t xa[]={0.8, 2.0, 3.5, 6.0, 20.} ;

  Int_t nPadX=1,nPadY=1;
  Int_t nbin=1 ;
  Double_t xa[]={2.0, 20.} ;

  PPRstyle();
  gStyle->SetOptFit(111);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.08);

  TString txtModule;
  if      (strcmp(cModule,"")   ==0) txtModule="All PHOS modules";
  else if (strcmp(cModule,"SM1")==0) txtModule="PHOS module 4";
  else if (strcmp(cModule,"SM2")==0) txtModule="PHOS module 3";
  else if (strcmp(cModule,"SM3")==0) txtModule="PHOS module 2";

  TPaveText *label = new TPaveText(0.15,0.80,0.40,0.90,"NDC");
  label->SetFillColor(kWhite);
  label->SetBorderSize(1);
  label->AddText(Form("Centrality %s" ,txtCent[centrality]));
  label->AddText(Form("%s",txtModule.Data()));

  //Fit real only 
  //Linear Bg
  char kkey[55];
  sprintf(kkey,outputKey.Data()) ;
  char key2[155];
  sprintf(key,"Mix%s",kkey) ;
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

  TF1 * fit1 = new TF1("fit",CB,0.,1.,6) ;
  fit1->SetParName(0,"A") ;
  fit1->SetParName(1,"m_{0}") ;
  fit1->SetParName(2,"#sigma") ;
  fit1->SetParName(3,"a_{0}") ;
  fit1->SetParName(4,"a_{1}") ;
  fit1->SetParName(5,"a_{2}") ;
  fit1->SetLineWidth(2) ;
  fit1->SetLineColor(2) ;
  TF1 * fgs = new TF1("gs",CBs,0.,1.,3) ;
  fgs->SetParName(0,"A") ;
  fgs->SetParName(1,"m_{0}") ;
  fgs->SetParName(2,"#sigma") ;
  fgs->SetLineColor(2) ;
  fgs->SetLineWidth(2) ;

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

  TCanvas *c4, *c2;
  TAxis * pta=h->GetYaxis() ;
  TAxis * ma=h->GetXaxis() ;
  for(Int_t i=1;i<=nbin;i++){
    if(i<=15)
       c1->cd(i) ;
    else{
      if(c2==0){
	c2 = new TCanvas("mggFit2","mggFit2",10,10,1200,800) ;   
	c2->Divide(nPadX,nPadY) ;
      }
      c2->cd(i-16) ;
    }
    printf("\t%.1f<pt<%.1f GeV/c\n",xa[i-1],xa[i]);
    Int_t imin=pta->FindBin(xa[i-1]+0.0001);
    Int_t imax=pta->FindBin(xa[i]-0.0001) ;
    Double_t pt=(xa[i]+xa[i-1])/2. ;
    TH1D * hp = h->ProjectionX("re",imin,imax) ;
    hp->Sumw2() ;
    TH1D * hpm= hm->ProjectionX("mi",imin,imax) ;
    hpm->Sumw2() ;
    // if(i<=11){
    if(i<=10){
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
    TH1D * hpcopy = (TH1D*)hp ->Clone("hpcopy") ;
    TH1D * hp2    = (TH1D*)hp ->Clone("hp2") ;
    hpcopy->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
    hp2   ->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
    hpcopy->Divide(hpm) ;
    sprintf(key,"%3.1f<p_{T}<%3.1f GeV/c",xa[i-1],xa[i]) ;
    hpcopy->SetTitle(key) ;
    hpcopy->SetMarkerStyle(20) ;
    hpcopy->SetMarkerSize(0.7) ;
    
    fit1->SetParameters(0.001,0.136,0.005,0.0002,-0.002,0.0) ;
    fit1->SetParLimits(0,0.000,1.000) ;
    fit1->SetParLimits(1,0.125,0.145) ;
    fit1->SetParLimits(2,0.002,0.015) ;

    Double_t rangeMin=0.06;
    Double_t rangeMax=0.26;

    hpcopy->Fit(fit1,"Q" ,"",rangeMin,rangeMax) ;
    hpcopy->Fit(fit1,"MQ","",rangeMin,rangeMax) ;

    ar1->SetBinContent(i,fit1->GetParameter(3)) ;
    ar1->SetBinError  (i,fit1->GetParError(3)) ;
    br1->SetBinContent(i,fit1->GetParameter(4)) ;
    br1->SetBinError  (i,fit1->GetParError(4)) ;

    fit2->SetParameters(fit1->GetParameters()) ;
    fit2->SetParLimits(0,0.000,1.000) ;
    fit2->SetParLimits(1,0.125,0.145) ;
    fit2->SetParLimits(2,0.002,0.015) ;

    hpcopy->Fit(fit2,"+NQ","",rangeMin,rangeMax) ;
    hpcopy->Fit(fit2,"+MQ","",rangeMin,rangeMax) ;

    ar2->SetBinContent(i,fit2->GetParameter(3)) ;
    ar2->SetBinError  (i,fit2->GetParError(3)) ;
    br2->SetBinContent(i,fit2->GetParameter(4)) ;
    br2->SetBinError  (i,fit2->GetParError(4)) ;
    cr2->SetBinContent(i,fit2->GetParameter(5)) ;
    cr2->SetBinError  (i,fit2->GetParError(5)) ;
    hpcopy->SetAxisRange(0.04,0.3,"X") ;
    hpcopy->Draw() ;
    label ->Draw();
    if(c2)
      c2->Update() ;
    else
      c1->Update() ;
//    if(getchar()=='q')return ;

    fbg1->SetParameters(fit1->GetParameter(3),fit1->GetParameter(4),fit1->GetParameter(5)); 
    fbg2->SetParameters(fit2->GetParameter(3),fit2->GetParameter(4),fit2->GetParameter(5),
			fit2->GetParameter(6)); 

    Double_t intRangeMin = PeakPosition(pt,centrality)-3.*PeakWidth(pt,centrality) ;
    Double_t intRangeMax = PeakPosition(pt,centrality)+3.*PeakWidth(pt,centrality) ;
    // printf("pt=%f, %f < M < %f\n",pt,intRangeMin,intRangeMax);
    Int_t    intBinMin   = hp->GetXaxis()->FindBin(intRangeMin) ;
    Int_t    intBinMax   = hp->GetXaxis()->FindBin(intRangeMax) ;
    Double_t errStat     = hpm->Integral(intBinMin,intBinMax); 

    hpm ->Multiply(fbg1) ;
    hpm2->Multiply(fbg2) ;
    hp  ->Add(hpm ,-1.) ;
    hp2 ->Add(hpm2,-1.) ;

    if(i<=15)
       c3->cd(i) ;
    else{
      if(c4==0){
	c4 = new TCanvas("mggFit2_Signal","mggFit2",10,10,1200,800) ;   
	c4->Divide(nPadX,nPadY) ;
      }
      c4->cd(i-15) ;
    }

    Int_t binPi0 = hp->FindBin(kMean);
    fgs->SetParameters(hp->Integral(binPi0-1,binPi0+1)/3.,fit1->GetParameter(1),fit1->GetParameter(2)) ;
    fgs->SetParLimits(0,0.000,1.e+5) ;
    fgs->SetParLimits(1,0.120,0.145) ;
    fgs->SetParLimits(2,0.002,0.015) ;
    hp->Fit(fgs,"Q","",rangeMin,rangeMax) ;   
    hp->SetMaximum(hp2->GetMaximum()*1.4) ;
    hp->SetMinimum(hp2->GetMinimum()*1.1) ;
    hp->SetMarkerStyle(20) ;
    hp->SetMarkerSize(0.7) ;
    mr1->SetBinContent(i,fgs->GetParameter(1)) ;
    mr1->SetBinError  (i,fgs->GetParError(1) ) ;
    sr1->SetBinContent(i,TMath::Abs(fgs->GetParameter(2))) ;
    sr1->SetBinError  (i,fgs->GetParError(2) ) ;

    Double_t y=fgs->Integral(intRangeMin,intRangeMax)/hp->GetXaxis()->GetBinWidth(1) ;
    nr1->SetBinContent(i,y) ;
    Double_t ey=0 ;
    if(fgs->GetParameter(0)!=0. && fgs->GetParameter(2)!=0.){
      Double_t en=fgs->GetParError(0)/fgs->GetParameter(0) ;
      Double_t es=fgs->GetParError(2)/fgs->GetParameter(2) ;
      ey=y*TMath::Sqrt(en*en+es*es) ;
    }
    nr1->SetBinError(i,ey) ;

    Double_t npiInt = hp->Integral(intBinMin,intBinMax) ;
    Double_t norm   = fbg1->GetParameter(0) ;
    Double_t normErr= fbg1->GetParError(0) ;
    if(npiInt>0.){
      nr1int->SetBinContent(i,npiInt) ;
      nr1int->SetBinError(i,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
    }
    hp2->SetAxisRange(0.04,0.3,"X") ;
    hp2->SetMaximum(hp2->GetMaximum()*1.4) ;
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
    hp2->Draw() ;
    label ->Draw();
    if(c4)
      c4->Update() ;
    else
      c3->Update() ;

    delete hp ;
//    delete hp2 ;
//    delete hpcopy ;
    delete hpm ;
    delete hpm2 ;
  }
  char name[55] ;

  sprintf(name,"Pi0_ratio_cent%d%s.eps" ,centrality,cModule) ;
  c1->Print(name) ;
  sprintf(name,"Pi0_signal_cent%d%s.eps",centrality,cModule) ;
  c3->Print(name) ;

  //Normalize by the number of events
  Int_t cMin,cMax;
  if      (centrality == 0) {
    cMin=1;
    cMax=20;
  }
  else if (centrality == 1) {
    cMin=21;
    cMax=80;
  }
  Double_t nevents = hCentrality1->Integral(cMin,cMax);
  nr1   ->Scale(1./nevents) ;
  nr1int->Scale(1./nevents) ;
  nr2   ->Scale(1./nevents) ;
  nr2int->Scale(1./nevents) ;

  printf("\t==============================\n");
  printf("\t|       N events = %i8    |\n",nevents);
  printf("\t==============================\n");

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
Double_t PeakPosition(Double_t pt, Int_t centrality = 0){
  //Fit to the measured peak position
  Double_t pi0Peak = 0.135;
  if      (centrality == 0) pi0Peak=0.1413;
  else if (centrality == 1) pi0Peak=0.1380;
  else if (centrality == 2) pi0Peak=0.1367;
  else if (centrality == 3) pi0Peak=0.1359;
  return pi0Peak;
  // return 4.99292e-003*exp(-pt*9.32300e-001)+1.34944e-001 ;
}
//-----------------------------------------------------------------------------
Double_t PeakWidth(Double_t pt, Int_t centrality = 0){
  //fit to the measured peak width
  Double_t pi0Sigma = 0.07;
  if      (centrality == 0) pi0Sigma=0.0084;
  else if (centrality == 1) pi0Sigma=0.0071;
  else if (centrality == 2) pi0Sigma=0.0068;
  else if (centrality == 3) pi0Sigma=0.0064;
  return pi0Sigma;

  // Double_t a=0.0068 ;
  // Double_t b=0.0025 ;
  // Double_t c=0.000319 ;
  // return TMath::Sqrt(a*a+b*b/pt/pt+c*c*pt*pt) ;
}
 
//-----------------------------------------------------------------------------
Double_t CB(Double_t * x, Double_t * par){
  //Parameterization of Real/Mixed ratio
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  Double_t dMpi0 = x[0]-kMean;
  Double_t peak = par[0] * TMath::Exp(-dx*dx/2.);
  Double_t bg   = par[3] + par[4]*dMpi0 + par[5]*dMpi0*dMpi0;
  return peak+bg ;
}

//-----------------------------------------------------------------------------
Double_t CB2(Double_t * x, Double_t * par){
  //Another parameterization of Real/Mixed ratio7TeV/README
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  Double_t dMpi0 = x[0]-kMean;
  Double_t peak = par[0] * TMath::Exp(-dx*dx/2.);
  Double_t bg   = par[3] + par[4]*dMpi0 + par[5]*dMpi0*dMpi0 + par[6]*dMpi0*dMpi0*dMpi0;
  return peak+bg ;
}
//-----------------------------------------------------------------------------
Double_t CBs(Double_t * x, Double_t * par){
  //Parameterizatin of signal
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  Double_t peak = par[0] * TMath::Exp(-dx*dx/2.);
  return peak ;
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

