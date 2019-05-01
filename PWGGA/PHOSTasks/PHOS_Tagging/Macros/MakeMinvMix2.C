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
    par[6]*(x[0]-kMean);
//    par[7]*(x[0]-kMean)*(x[0]-kMean) ;
  else{
    Double_t A=TMath::Power((n/TMath::Abs(a)),n)*TMath::Exp(-a*a/2) ;
    Double_t B=n/TMath::Abs(a)-TMath::Abs(a) ;
    return par[0]*A*TMath::Power((B-dx),-n)+
    par[5]+
    par[6]*(x[0]-kMean) ;
//    par[7]*(x[0]-kMean)*(x[0]-kMean) ;
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
    par[7]*(x[0]-kMean)*(x[0]-kMean) ;
//    par[8]*(x[0]-kMean)*(x[0]-kMean)*(x[0]-kMean) ;
  else{
    Double_t A=TMath::Power((n/TMath::Abs(a)),n)*TMath::Exp(-a*a/2) ;
    Double_t B=n/TMath::Abs(a)-TMath::Abs(a) ;
    return par[0]*A*TMath::Power((B-dx),-n)+
    par[5]+
    par[6]*(x[0]-kMean)+
    par[7]*(x[0]-kMean)*(x[0]-kMean) ;
//    par[8]*(x[0]-kMean)*(x[0]-kMean)*(x[0]-kMean) ;
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
  par[4]*(x[0]-kMean) ;
//  par[5]*(x[0]-kMean)*(x[0]-kMean);
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
  par[5]*(x[0]-kMean)*(x[0]-kMean);
//  par[6]*(x[0]-kMean)*(x[0]-kMean)*(x[0]-kMean)  ;
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
  return par[0]+par[1]*(x[0]-kMean) ;//+par[2]*(x[0]-kMean)*(x[0]-kMean) ;
}
//-----------------------------------------------------------------------------
Double_t BG2(Double_t * x, Double_t * par){
  //Another normalization of  Mixed
  return par[0]+par[1]*(x[0]-kMean)+par[2]*(x[0]-kMean)*(x[0]-kMean); //+par[3]*(x[0]-kMean)*(x[0]-kMean)*(x[0]-kMean) ;
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



void MakeMinvMix2(Int_t cen=6){
  TFile fout("raw_MBmix2.root","update") ;
  TFile * f = new TFile("LHC13bc_MB.root") ; //Result of data scan
  

  TH1F * hev = (TH1F*)f->Get("hSelEvents") ;
  TH1F * hCentrality1  = (TH1F*)f->Get("hCentrality") ;
  
  printf("TotSelEvents: %f \n",hev->GetBinContent(7)) ;
  printf("Centrality:   %f \n",hCentrality1->Integral()) ;

   Int_t nbin=33 ;
   Double_t xa[34] ={0.8,1.0,1.2,1.4,1.6, 1.8,2.0,2.2,2.4,2.6, 2.8,3.0,3.2,3.4,3.6, 3.8,4.0,4.5,5.0,5.5, 6.,7.,8.,10.,12.,16.,20.,22.,  24.,26.,28.,30.,35.,  40.};
/*   Int_t nbin=26 ;
   Double_t xa[27] ={0.8,1.0,1.2,1.4,1.6, 1.8,2.0,2.2,2.4,2.6, 2.8,3.0,3.2,3.4,3.6, 3.8,4.0,4.5,5.0,5.5, 6.,7.,8.,10.,12.,16.,20.};
 */  
  const Int_t nPID=4 ;
  char cPID[14][15] ;
  snprintf(cPID[0],15,"Emin3_All") ;
  snprintf(cPID[1],15,"Emin3_Disp");
  snprintf(cPID[2],15,"Emin3_CPV") ;
  snprintf(cPID[3],15,"Emin3_Both"); 

  TH2F * hRe[20] ;
  TH2F * hMi[20] ;
  TH2F * tmp = 0x0 ;
  if(cen==6){ //6:0-20, 7:0-10
      for(Int_t iPID=0;iPID<nPID;iPID++){
	hRe[iPID] = (TH2F*)f->Get(Form("hInvM_Re_%s_cent0",cPID[iPID])) ;
        tmp = (TH2F*)f->Get(Form("hInvM_Re_%s_cent1",cPID[iPID])) ;
        hRe[iPID]->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Re_%s_cent2",cPID[iPID])) ;
        hRe[iPID]->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Re_%s_cent3",cPID[iPID])) ;
        hRe[iPID]->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Re_%s_cent4",cPID[iPID])) ;
        hRe[iPID]->Add(tmp) ;delete tmp ;
        hRe[iPID]->Sumw2() ;	
	
	hMi[iPID] = (TH2F*)f->Get(Form("hInvM_Mi_%s_cent0",cPID[iPID])) ;
        tmp = (TH2F*)f->Get(Form("hInvM_Mi_%s_cent1",cPID[iPID])) ;
        hMi[iPID]->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Mi_%s_cent2",cPID[iPID])) ;
        hMi[iPID]->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Mi_%s_cent3",cPID[iPID])) ;
        hMi[iPID]->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Mi_%s_cent4",cPID[iPID])) ;
        hMi[iPID]->Add(tmp) ;delete tmp ;
        hMi[iPID]->Sumw2() ;	
      }
  }
  else{
      for(Int_t iPID=0;iPID<nPID;iPID++){
	hRe[iPID] = (TH2F*)f->Get(Form("hInvM_Re_%s_cent%d",cPID[iPID],cen)) ;	
        hRe[iPID]->Sumw2() ;	
	hMi[iPID] = (TH2F*)f->Get(Form("hInvM_Mi_%s_cent%d",cPID[iPID],cen)) ;
        hMi[iPID]->Sumw2() ;	
      }    
  }
      
  PPRstyle();
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadBottomMargin(0.08);


  //Gaus
  TH1D *mr1[nPID],*sr1[nPID],*nr1[nPID],*nr1int[nPID] ;
  for(Int_t iPID=0;iPID<nPID;iPID++){
    mr1[iPID] = new TH1D(Form("mass1_GS_%s_cen%d",cPID[iPID],cen),"Mass",nbin,xa) ;
    sr1[iPID] = new TH1D(Form("width1_GS_%s_cen%d",cPID[iPID],cen),"Width",nbin,xa) ;
    nr1[iPID] = new TH1D(Form("yeild1_GS_%s_cen%d",cPID[iPID],cen),"Raw yield",nbin,xa) ;
    nr1int[iPID] = new TH1D(Form("yeild1_int_GS_%s_cen%d",cPID[iPID],cen),"Raw yield, integrated",nbin,xa) ;
  }

  //CB
  TH1D *mr2[nPID],*sr2[nPID],*nr2[nPID],*nr2int[nPID] ;
  for(Int_t iPID=0;iPID<nPID;iPID++){  
    mr2[iPID] = new TH1D(Form("mass2_CB_%s_cen%d",cPID[iPID],cen),"Mass",nbin,xa) ;
    sr2[iPID] = new TH1D(Form("width2_CB_%s_cen%d",cPID[iPID],cen),"Width",nbin,xa) ;
    nr2[iPID] = new TH1D(Form("yeild2_CB_%s_cen%d",cPID[iPID],cen),"Raw yield",nbin,xa) ;
    nr2int[iPID] = new TH1D(Form("yeild2_int_CB_%s_cen%d",cPID[iPID],cen),"Raw yield, integrated",nbin,xa) ;
  }

  TF1 * fun1 = new TF1("ft1",CB1,0.,1.,7) ;
  fun1->SetParName(0,"A") ;
  fun1->SetParName(1,"m_{0}") ;
  fun1->SetParName(2,"#sigma") ;
  fun1->SetParName(3,"n") ;
  fun1->SetParName(4,"#alpha") ;
  fun1->SetParName(5,"a_{0}") ;
  fun1->SetParName(6,"a_{1}") ;
//   fun1->SetParName(7,"a_{2}") ;
  fun1->FixParameter(3,3.28) ;
  fun1->FixParameter(4,1.56) ;
  fun1->SetLineWidth(2) ;
  fun1->SetLineColor(4) ;
  fun1->SetLineStyle(2) ;

  TF1 * fun2 = new TF1("ft2",CB2,0.,1.,8) ;
  fun2->SetParName(0,"A") ;
  fun2->SetParName(1,"m_{0}") ;
  fun2->SetParName(2,"#sigma") ;
  fun2->SetParName(3,"n") ;
  fun2->SetParName(4,"#alpha") ;
  fun2->SetParName(5,"a_{0}") ;
  fun2->SetParName(6,"a_{1}") ;
  fun2->SetParName(7,"a_{2}") ;
//   fun2->SetParName(8,"a_{3}") ;
  fun2->FixParameter(3,3.28) ;
  fun2->FixParameter(4,1.56) ;
  fun2->SetLineWidth(2) ;
  fun2->SetLineColor(4) ;
  fun2->SetLineStyle(2) ;
  
  TF1 * funGS1 = new TF1("ft1",GS1,0.,1.,5) ;
  funGS1->SetParName(0,"A") ;
  funGS1->SetParName(1,"m_{0}") ;
  funGS1->SetParName(2,"#sigma") ;
  funGS1->SetParName(3,"a_{0}") ;
  funGS1->SetParName(4,"a_{1}") ;
//   funGS1->SetParName(5,"a_{2}") ;
  funGS1->SetLineWidth(2) ;
  funGS1->SetLineColor(8) ;
//   funGS1->SetLineStyle(2) ;

  TF1 * funGS2 = new TF1("ft2",GS2,0.,1.,6) ;
  funGS2->SetParName(0,"A") ;
  funGS2->SetParName(1,"m_{0}") ;
  funGS2->SetParName(2,"#sigma") ;
  funGS2->SetParName(3,"a_{0}") ;
  funGS2->SetParName(4,"a_{1}") ;
  funGS2->SetParName(5,"a_{2}") ;
//   funGS2->SetParName(6,"a_{3}") ;
  funGS2->SetLineWidth(2) ;
  funGS2->SetLineColor(8) ;
//   funGS2->SetLineStyle(2) ;
  
  
  
  TF1 * fit1=0x0 ;
  TF1 * fit2=0x0 ;
  
  TF1 * fbgP0 = new TF1("bg",BG,0.,1.,2) ;
  TF1 * fbgP1 = new TF1("bg1",BG1,0.,1.,2) ;
  TF1 * fbgP2 = new TF1("bg2",BG2,0.,1.,3) ;
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

  TCanvas * c1[nPID] ;
  TCanvas * c2[nPID] ;
  TCanvas * c1b[nPID] ;
  TCanvas * c2b[nPID] ;
  for(Int_t iPID=0;iPID<nPID;iPID++){
    c1[iPID] = new TCanvas(Form("Ratio_%s",cPID[iPID]),Form("Ratio_%s",cPID[iPID]),10+10*iPID,10,1200+10*iPID,800) ;
    c1[iPID]->Divide(4,4) ;
    c2[iPID] = new TCanvas(Form("Signal_%s",cPID[iPID]),Form("Signal_%s",cPID[iPID]),10+10*iPID,10,1200+10*iPID,800) ;
    c2[iPID]->Divide(4,4) ;
    c1b[iPID] = new TCanvas(Form("RatioB_%s",cPID[iPID]),Form("Ratio_%s",cPID[iPID]),10+10*iPID,10,1200+10*iPID,800) ;
    c1b[iPID]->Divide(4,4) ;
    c2b[iPID] = new TCanvas(Form("SignalB_%s",cPID[iPID]),Form("Signal_%s",cPID[iPID]),10+10*iPID,10,1200+10*iPID,800) ;
    c2b[iPID]->Divide(4,4) ;
  }
  c1[0]->cd(0) ;

  TAxis * pta=hRe[0]->GetYaxis() ;
  TAxis * ma=hRe[0]->GetXaxis() ;
    
  for(Int_t i=1;i<=nbin;i++){
    Int_t imin=pta->FindBin(xa[i-1]+0.0001);
    Int_t imax=pta->FindBin(xa[i]-0.0001) ;
    TH1D *hp;
    Double_t pt=(xa[i]+xa[i-1])/2. ;
        
    for(Int_t iPID=0;iPID<nPID;iPID++){
      if(i<17)
        c1[iPID]->cd(i) ;
      else
        c1b[iPID]->cd(i-16) ;
      hp = hRe[iPID]->ProjectionX(Form("re%d_%d",i,iPID),imin,imax) ;
//      hp->Sumw2() ;
      hpm= hMi[iPID]->ProjectionX(Form("mi%d_%d",i,iPID),imin,imax) ;
//      hpm->Sumw2() ;
      if(pt>12.){
        hp ->Rebin(2) ;
        hpm->Rebin(2) ;
      }

    if(pt<2.){ //2pol,3pol
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
    
//      for(Int_t ib=1; ib<=hp->GetNbinsX();ib++){if(hp ->GetBinContent(ib)==0)hp ->SetBinError(ib,1.);}
//      for(Int_t ib=1; ib<=hp->GetNbinsX();ib++){if(hpm->GetBinContent(ib)==0)hpm->SetBinError(ib,1.);}
      TH1D * hpm2   = (TH1D*)hpm->Clone(Form("Bg1_%d",iPID)) ;
      TH1D * hpcopy = (TH1D*)hp ->Clone(Form("hpcopy_%d",iPID)) ;
      TH1D * hp2    = (TH1D*)hp ->Clone(Form("hp2_%d",iPID)) ;
      hpcopy->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
      hp2   ->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
      hpcopy->Divide(hpm) ;
      hpcopy->SetTitle(Form("%3.1f<p_{T}<%3.1f GeV/c",xa[i-1],xa[i])) ;
      hpcopy->SetMarkerStyle(20) ;
      hpcopy->SetMarkerSize(0.7) ;
      hpcopy->GetXaxis()->SetRangeUser(0.05,0.25) ;
      
//      fit1->SetParameters(0.0002+0.0001*i*i,0.136,0.011,0.0002,-0.002,0.0) ;
      fit1->SetParameters(0.0002+0.0001*i*i,0.136,0.005,0.0002,-0.002,0.0) ;
//      fit1->SetParameters(0.0002,0.136,0.011,0.0002,-0.002,0.0) ;
      if(cen==0)
        fit1->SetParameters(0.001,0.146,0.005,0.12,-0.002,0.0) ;
      fit1->SetParLimits(0,0.000,2.*hpcopy->GetMaximum()) ;
      fit1->SetParLimits(1,0.125,0.145) ;
      fit1->SetParLimits(2,0.0045,0.010) ;

      Double_t rangeMin=0.09 ;//TMath::Max(0.06,0.11-0.01*i) ;
      Double_t rangeMax=0.18; //TMath::Min(0.25,0.18+0.01*i) ;
      if(pt>7.){
         rangeMin=0.07;
         rangeMax=0.22;
      }
//      Double_t rangeMin=TMath::Max(0.06,0.12-0.01*i) ;
//      Double_t rangeMax=TMath::Min(0.25,0.16+0.01*i) ;

      hpcopy->Fit(fit1,"Q","",rangeMin,rangeMax) ;
      hpcopy->Fit(fit1,"MQ","",rangeMin,rangeMax) ;

      fit2->SetParameters(fit1->GetParameters()) ;
      fit2->SetParameter(3,3.2) ;
      fit2->SetParameter(4,1.56) ;
      fit2->SetParameter(5,fit1->GetParameter(3)) ;
      fit2->SetParameter(6,fit1->GetParameter(4)) ;
      fit2->SetParameter(7,fit1->GetParameter(5)) ;      
      fit2->SetParLimits(0,0.000,2.*hpcopy->GetMaximum()) ;
      fit2->SetParLimits(1,0.125,0.145) ;
      fit2->SetParLimits(2,0.0045,0.010) ;
      
      hpcopy->Fit(fit2,"+QN","",rangeMin,rangeMax) ;
      hpcopy->Fit(fit2,"+MQ","",rangeMin,rangeMax) ;
      c1[iPID]->cd(i) ;
      hpcopy->Draw() ;
      if(i<17)
        c1[iPID]->Update() ;
      else
        c1b[iPID]->Update() ;

      if(i<17)
        c2[iPID]->cd(i) ;
      else
        c2b[iPID]->cd(i-16) ;
      
      fbg1->SetParameters(fit1->GetParameter(3),fit1->GetParameter(4),fit1->GetParameter(5),fit1->GetParameter(6)); 
      fbg2->SetParameters(fit2->GetParameter(5),fit2->GetParameter(6),fit2->GetParameter(7),fit2->GetParameter(8)); 
      
      Double_t intRangeMin = PeakPosition(pt)-3.*PeakWidth(pt) ;
      Double_t intRangeMax = PeakPosition(pt)+3.*PeakWidth(pt) ;
      Int_t    intBinMin   = hp->GetXaxis()->FindBin(intRangeMin) ;
      Int_t    intBinMax   = hp->GetXaxis()->FindBin(intRangeMax) ;
      Double_t errStat     = hpm->Integral(intBinMin,intBinMax); 

      hpm ->Multiply(fbg1) ;
      hpm2->Multiply(fbg2) ;
      if(pt<10.){ //Mixed is poor at pT>10 GeV
        hp  ->Add(hpm ,-1.) ;
        hp2 ->Add(hpm2,-1.) ;
      }          
      
      Int_t binPi0 = hp->FindBin(kMean);
      fgs->SetParameters(hp->Integral(binPi0-1,binPi0+1)/3.,fit1->GetParameter(1),0.006) ;
      fgs->SetParLimits(0,0.000,10.*hp->GetMaximum()) ;
      fgs->SetParLimits(1,0.120,0.175) ;
      fgs->SetParLimits(2,0.0055,0.010) ;
      hp->Fit(fgs,"QL","",rangeMin,rangeMax) ;   
      hp->Fit(fgs,"QML","",rangeMin,rangeMax) ;   
      hp->SetMaximum(hp->GetMaximum()*1.4) ;
      hp->SetMinimum(hp->GetMinimum()*1.1) ;
      hp->SetMarkerStyle(20) ;
      hp->SetMarkerSize(0.7) ;
      mr1[iPID]->SetBinContent(i,fgs->GetParameter(1)) ;
      mr1[iPID]->SetBinError  (i,fgs->GetParError(1) ) ;
      sr1[iPID]->SetBinContent(i,TMath::Abs(fgs->GetParameter(2))) ;
      sr1[iPID]->SetBinError  (i,fgs->GetParError(2) ) ;

      Double_t y=fgs->GetParameter(0)/hp->GetXaxis()->GetBinWidth(1) ;
      nr1[iPID]->SetBinContent(i,y) ;
      Double_t ey=fgs->GetParError(0)/hp->GetXaxis()->GetBinWidth(1) ;
      nr1[iPID]->SetBinError(i,ey) ;

      Double_t npiInt = hp->Integral(intBinMin,intBinMax)-(intBinMax-intBinMin)*fgs->GetParameter(3) ;
      Double_t norm   = fbg1->GetParameter(0) ;
      Double_t normErr= fbg1->GetParError(0) ;
      if(npiInt>0.){
        nr1int[iPID]->SetBinContent(i,npiInt) ;
        nr1int[iPID]->SetBinError(i,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
      }
      
//printf(" Nint1 =%f+-%f \n",npiInt,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
      hp2->GetXaxis()->SetRangeUser(rangeMin,rangeMax) ;
      hp2->SetMaximum(hp2->GetMaximum()*1.4) ;
      hp2->SetMinimum(hp2->GetMinimum()*1.1) ;
      hp2->SetMarkerStyle(24) ;
      hp2->SetMarkerSize(0.8) ;

      fcb->SetParameters(hp->Integral(binPi0-1,binPi0+1)/3.,fit2->GetParameter(1),0.006,fit2->GetParameter(3),fit2->GetParameter(4)) ;
      fcb->SetParLimits(0,0.000,10.*hp->GetMaximum()) ;
      fcb->SetParLimits(1,0.120,0.175) ;
      fcb->SetParLimits(2,0.0055,0.010) ;
      hp2->Fit(fcb,"QL","",rangeMin,rangeMax) ;
      hp2->Fit(fcb,"QML","",rangeMin,rangeMax) ;
      mr2[iPID]->SetBinContent(i,fcb->GetParameter(1)) ;
      mr2[iPID]->SetBinError  (i,fcb->GetParError(1)) ;
      sr2[iPID]->SetBinContent(i,TMath::Abs(fcb->GetParameter(2))) ;
      sr2[iPID]->SetBinError  (i,fcb->GetParError(2)) ;
      
      y=(fcb->Integral(0.05,0.25)-fcb->GetParameter(5)*(0.25-0.05))/hp->GetXaxis()->GetBinWidth(1) ;
      nr2[iPID]->SetBinContent(i,y) ;
      Double_t ey=fcb->IntegralError(0.05,0.25)/hp->GetXaxis()->GetBinWidth(1) ;
      nr2[iPID]->SetBinError(i,ey) ;
   
      npiInt=hp2->Integral(intBinMin,intBinMax)-(intBinMax-intBinMin)*fcb->GetParameter(5) ;
      norm=fbg2->GetParameter(0) ;
      normErr=fbg2->GetParError(0) ;
      if(npiInt>0.){
        nr2int[iPID]->SetBinContent(i,npiInt) ;
        nr2int[iPID]->SetBinError(i,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
      } 
//printf(" Nint2 =%f+-%f \n",npiInt,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
      hp2->SetTitle(Form("%3.1f<p_{T}<%3.1f GeV/c",xa[i-1],xa[i])) ;
      hp2->Draw() ;
      hp->SetMarkerColor(6) ;
      hp->Draw("same") ;
      hp2->Draw("same") ;
      c2[iPID]->Update() ;

//      delete hp ;
//    delete hp2 ;
//    delete hpcopy ;
      delete hpm ;
//      delete hpm2 ;
      
    }
  }
/*  
  for(Int_t iPID=0; iPID<4; iPID++){
   c1[iPID]->Print(Form("Ratio_%s_cen%d.eps",cPID[iPID],cen)) ;
   c2[iPID]->Print(Form("Signal_%s_cen%d.eps",cPID[iPID],cen)) ;
  }
*/  
  //Normalize by the number of events
  Int_t cMin,cMax;
  if      (cen == 0) {
    cMin=1;
    cMax=20;
  }
  else if (cen == 1) {
    cMin=21;
    cMax=40;
  }
  else if (cen == 2) {
    cMin=41;
    cMax=60;
  }
  else if (cen == 3) {
    cMin=61;
    cMax=80;
  }
  else if (cen == 4) {
    cMin=81;
    cMax=100;
  }
  else if (cen == 5) {
    cMin=61;
    cMax=80;
  }
  else if (cen == 6) {
    cMin=1;
    cMax=100;
  }
  else if (cen == 7) {
    cMin=1;
    cMax=10;
  }
  else if (cen == 8) {
    cMin=41;
    cMax=80;
  }
  else if (cen == 9) {
    cMin=1;
    cMax=40;
  }
  Double_t nevents = hCentrality1->Integral(cMin,cMax);
  printf("Nevents=%f \n",nevents) ;
  for(Int_t iPID=0;iPID<nPID;iPID++){
    nr1[iPID]   ->Scale(1./nevents) ;
    nr1int[iPID]->Scale(1./nevents) ;
    nr2[iPID]   ->Scale(1./nevents) ;
    nr2int[iPID]->Scale(1./nevents) ;
    nr1[iPID]->SetMarkerStyle(20) ;
    nr1[iPID]->SetMarkerColor(2) ;
    nr1[iPID]->SetTitle("#pi^{0} raw yield per event") ;
    nr2[iPID]->SetTitle("#pi^{0} raw yield per event") ;
  }
  
  fout.cd() ;
  for(Int_t iPID=0;iPID<nPID;iPID++){
    nr1[iPID]->Write(0,TObject::kOverwrite) ;
    nr2[iPID]->Write(0,TObject::kOverwrite) ;
    nr1int[iPID]->Write(0,TObject::kOverwrite) ;
    nr2int[iPID]->Write(0,TObject::kOverwrite) ;
    mr1[iPID]->Write(0,TObject::kOverwrite) ;
    mr2[iPID]->Write(0,TObject::kOverwrite) ;
    sr1[iPID]->Write(0,TObject::kOverwrite) ;
    sr2[iPID]->Write(0,TObject::kOverwrite) ;
  }
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

