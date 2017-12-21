const Double_t kMean=0.136 ; //Approximate peak position to facilitate error estimate
//-----------------------------------------------------------------------------
Double_t CB(Double_t * x, Double_t * par){
  //Parameterization of Real/Mixed ratio
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+par[3]+par[4]*(x[0]-kMean);
}
//-----------------------------------------------------------------------------
Double_t CB1(Double_t * x, Double_t * par){
  //Parameterization of Real/Mixed ratio
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+par[3]+par[4]*(x[0]-kMean)+par[5]*(x[0]-kMean)*(x[0]-kMean);
}

//-----------------------------------------------------------------------------
Double_t CB2(Double_t * x, Double_t * par){
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
Double_t CBs(Double_t * x, Double_t * par){
  //Parameterizatin of signal
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)/TMath::Sqrt(TMath::TwoPi())/s+par[3] ;
}
//-----------------------------------------------------------------------------
Double_t BG0(Double_t * x, Double_t * par){
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



void MakePi0Eff(Int_t cen=6){


  TFile * f = new TFile("MB.root") ; //result of MC scan
  TFile *fout = new TFile("Eff_MB.root","update") ;
  TFile *fmass= new TFile("mass_MB.root","update") ;
  
  char key[55] ;
  TH1F * hPrim = 0; 
  if(cen==6){ //0-20%
    hPrim =(TH1F*)f->Get("hMC_unitEta_pi0_cent0") ;
    TH1F * tmp1 = (TH1F*)f->Get("hMC_unitEta_pi0_cent1") ;
    hPrim->Add(tmp1) ;delete tmp1 ;
    tmp1 = (TH1F*)f->Get("hMC_unitEta_pi0_cent2") ;
    hPrim->Add(tmp1) ;delete tmp1 ;
    tmp1 = (TH1F*)f->Get("hMC_unitEta_pi0_cent3") ;
    hPrim->Add(tmp1) ;delete tmp1 ;
    tmp1 = (TH1F*)f->Get("hMC_unitEta_pi0_cent4") ;
    hPrim->Add(tmp1) ;delete tmp1 ;
    
//     hPrim->Sumw2() ;
  }

    
  hPrim->Scale(0.5/0.12*360./60.) ; //to have 2pi and 1 in rapidity

// hPrim->Draw(); return ;  
 
//    Int_t nbin=35 ;
//    Double_t xa[36] ={0.8,1.0,1.2,1.4,1.6, 1.8,2.0,2.2,2.4,2.6, 2.8,3.0,3.2,3.4,3.6, 3.8,4.0,4.5,5.0,5.5, 6.,7.,8.,10.,12., 14.,16.,18.,20.,22.,  24.,26.,28.,30.,35.,  40.};
//   Int_t nbin=21 ;
//   Double_t xa[22] ={0.0,0.3,0.5,0.7,0.9,1.1,1.4,1.8,2.2,2.6,3.0,3.5,4.0,5.0,6.0,8.0,10.0,12.0,15.0,20.0,25.,30.};

//   Int_t nbin=67;
//   Double_t xa[68]={0.,0.1,0.2,0.3,0.4, 0.5,0.6,0.7,0.8,0.9, 1.,1.1,1.2,1.3,1.4, 1.5,1.6,1.7,1.8,1.9, 2.0,2.2,2.4,2.6,2.8, 
//                      3.,3.2,3.4,3.6,3.8, 4.0,4.2,4.4,4.6,4.8, 5.,5.5,6.0,6.5,7.0, 7.5,8.0,8.5,9.0,9.5, 10.,11.,12.,13.,14.,
//                      15.,16.,17.,18.,19.,20.,22.,24.,26.,28., 30.,35.,40.,45.,50.,55.,60.,65.};

   Int_t nbin=33 ;
   Double_t xa[34] ={0.8,1.0,1.2,1.4,1.6, 1.8,2.0,2.2,2.4,2.6, 2.8,3.0,3.2,3.4,3.6, 3.8,4.0,4.5,5.0,5.5, 6.,7.,8.,10.,12.,16.,20.,22.,  24.,26.,28.,30.,35.,  40.};
  
   
   
   
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
//         hRe[iPID]->Sumw2() ;	
	
	hMi[iPID] = (TH2F*)f->Get(Form("hInvM_Mi_%s_cent0",cPID[iPID])) ;
        tmp = (TH2F*)f->Get(Form("hInvM_Mi_%s_cent1",cPID[iPID])) ;
        hMi[iPID]->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Mi_%s_cent2",cPID[iPID])) ;
        hMi[iPID]->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Mi_%s_cent3",cPID[iPID])) ;
        hMi[iPID]->Add(tmp) ;delete tmp ;
        tmp = (TH2F*)f->Get(Form("hInvM_Mi_%s_cent4",cPID[iPID])) ;
        hMi[iPID]->Add(tmp) ;delete tmp ;
//         hMi[iPID]->Sumw2() ;	
      }
  }

  
  PPRstyle();
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetPadRightMargin(0.01);
  gStyle->SetPadTopMargin(0.01);
  gStyle->SetPadBottomMargin(0.08);


  //Linear Bg
  TH1D *mr1[nPID],*sr1[nPID],*nr1[nPID],*nr1int[nPID] ;
  for(Int_t iPID=0;iPID<nPID;iPID++){
    sprintf(key,"mass1_GS_%s_cen%d",cPID[iPID],cen) ;
    mr1[iPID] = new TH1D(key,"Mass",nbin,xa) ;
    sprintf(key,"width1_GS_%s_cen%d",cPID[iPID],cen) ;
    sr1[iPID] = new TH1D(key,"Width",nbin,xa) ;
    sprintf(key,"yeild1_GS_%s_cen%d",cPID[iPID],cen) ;
    nr1[iPID] = new TH1D(key,"Raw yield",nbin,xa) ;
    sprintf(key,"yeild1_int_GS_%s_cen%d",cPID[iPID],cen) ;
    nr1int[iPID] = new TH1D(key,"Raw yield, integrated",nbin,xa) ;
  }

  //Quadratic Bg
  TH1D *mr2[nPID],*sr2[nPID],*nr2[nPID],*nr2int[nPID] ;
  for(Int_t iPID=0;iPID<nPID;iPID++){  
    sprintf(key,"mass2_GS_%s_cen%d",cPID[iPID],cen) ;
    mr2[iPID] = new TH1D(key,"Mass",nbin,xa) ;
    sprintf(key,"width2_GS_%s_cen%d",cPID[iPID],cen) ;
    sr2[iPID] = new TH1D(key,"Width",nbin,xa) ;
    sprintf(key,"yeild2_GS_%s_cen%d",cPID[iPID],cen) ;
    nr2[iPID] = new TH1D(key,"Raw yield",nbin,xa) ;
    sprintf(key,"yeild2_int_GS_%s_cen%d",cPID[iPID],cen) ;
    nr2int[iPID] = new TH1D(key,"Raw yield, integrated",nbin,xa) ;
  }
  TF1 * fun = new TF1("ft0",CB,0.,1.,5) ;
  fun->SetParName(0,"A") ;
  fun->SetParName(1,"m_{0}") ;
  fun->SetParName(2,"#sigma") ;
  fun->SetParName(3,"a_{0}") ;
  fun->SetParName(4,"a_{1}") ;
  fun->SetLineWidth(2) ;
  fun->SetLineColor(2) ;

  TF1 * fun1 = new TF1("ft1",CB1,0.,1.,6) ;
  fun1->SetParName(0,"A") ;
  fun1->SetParName(1,"m_{0}") ;
  fun1->SetParName(2,"#sigma") ;
  fun1->SetParName(3,"a_{0}") ;
  fun1->SetParName(4,"a_{1}") ;
  fun1->SetLineWidth(2) ;
  fun1->SetLineColor(4) ;

  TF1 * fun2 = new TF1("ft2",CB2,0.,1.,7) ;
  fun2->SetParName(0,"A") ;
  fun2->SetParName(1,"m_{0}") ;
  fun2->SetParName(2,"#sigma") ;
  fun2->SetParName(3,"a_{0}") ;
  fun2->SetParName(4,"a_{1}") ;
  fun2->SetLineWidth(2) ;
  fun2->SetLineColor(8) ;
  fun2->SetLineStyle(2) ;
  
  TF1 * fbgP0 = new TF1("bg",BG0,0.,1.,2) ;
  TF1 * fbgP1 = new TF1("bg1",BG1,0.,1.,3) ;
  TF1 * fbgP2 = new TF1("bg2",BG2,0.,1.,4) ;
  TF1 * fgs = new TF1("gs",CBs,0.,1.,4) ;
  fgs->SetLineColor(2) ;
  fgs->SetLineWidth(2) ;

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
  
  TH1D * aPrim = (TH1D*)hPrim->Rebin(nbin,"Primary",xa) ;
  
  for(Int_t i=1;i<=nbin;i++){
    Int_t imin=pta->FindBin(xa[i-1]+0.0001);
    Int_t imax=pta->FindBin(xa[i]-0.0001) ;
    TH1D *hp;
    Double_t pt=(xa[i]+xa[i-1])/2. ;
    
    Double_t nPrim=hPrim->Integral(hPrim->GetXaxis()->FindBin(xa[i-1]+0.0001),hPrim->GetXaxis()->FindBin(xa[i]-0.0001)) ;
//    Double_t nPrim=aPrim->GetBinContent(i) ;
   if(nPrim==0.) {
     printf("i=%d, no Primary \n",i) ;
      continue ;
   }
   
    if(pt<3.){ //2pol,3pol
      fit1=fun1 ;
      fit2=fun2 ;
      fbg1=fbgP1 ;
      fbg2=fbgP2 ;
      
    }
    else{ //pol1,pol2
      fit1=fun ;
      fit2=fun1 ;
      fbg1=fbgP0 ;
      fbg2=fbgP1 ;
      
    }
   

//    printf("pt(%d)=%f-%f: nPrim=%f, Old=%f \n",i,xa[i-1],xa[i],nPrim,noldPrim) ;
    
    for(Int_t iPID=0;iPID<nPID;iPID++){
      if(i<17)
        c1[iPID]->cd(i) ;
      else
        c1b[iPID]->cd(i-16) ;
      hp = hRe[iPID]->ProjectionX(Form("re%d_%d",i,iPID),imin,imax) ;
//       hp->Sumw2() ;
      hpm= hMi[iPID]->ProjectionX(Form("mi%d_%d",i,iPID),imin,imax) ;
//       hpm->Sumw2() ;

     for(Int_t ib=1; ib<=hp->GetNbinsX();ib++){if(hp ->GetBinError(ib)==0)hp ->SetBinError(ib,hp ->GetBinError(25));}
     for(Int_t ib=1; ib<=hpm->GetNbinsX();ib++){if(hpm->GetBinError(ib)==0)hpm->SetBinError(ib,hpm->GetBinError(25));}
      
    hp ->Rebin(2) ;
    hpm->Rebin(2) ;
      hp->SetBinError(20,hp->GetBinError(13));
      hp->SetBinError(21,hp->GetBinError(13));
      hpm->SetBinError(20,hpm->GetBinError(13));
      hpm->SetBinError(21,hpm->GetBinError(13));
     
      
      TH1D * hpm2   = (TH1D*)hpm->Clone(Form("Bg1_%d",iPID)) ;
      TH1D * hpcopy = (TH1D*)hp ->Clone(Form("hpcopy_%d",iPID)) ;
      TH1D * hp2    = (TH1D*)hp ->Clone(Form("hp2_%d",iPID)) ;
      hpcopy->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
      hp2   ->SetXTitle("M_{#gamma#gamma} (GeV/c^{2})");
      hpcopy->Divide(hpm) ;
      hpcopy->SetTitle(Form("%3.1f<p_{T}<%3.1f GeV/c",xa[i-1],xa[i])) ;
      hpcopy->SetMarkerStyle(20) ;
      hpcopy->SetMarkerSize(0.7) ;
      hpcopy->GetXaxis()->SetRangeUser(0.05,0.24) ;
      
//      fit1->SetParameters(0.0002+0.0001*i*i,0.136,0.011,0.0002,-0.002,0.0) ;
      fit1->SetParameters(10.,0.135,0.006,0.1,0.,0.0) ;
//      fit1->SetParLimits(0,0.000,5.*hpcopy->GetMaximum()) ;
      fit1->SetParLimits(1,0.125,0.145) ;
      fit1->SetParLimits(2,0.0035,0.010) ;

      Double_t rangeMin=0.06 ;//TMath::Max(0.06,0.11-0.01*i) ;
      Double_t rangeMax=0.20; //TMath::Min(0.25,0.18+0.01*i) ;
//      Double_t rangeMin=TMath::Max(0.06,0.12-0.01*i) ;
//      Double_t rangeMax=TMath::Min(0.25,0.16+0.01*i) ;
      hpcopy->Fit(fit1,"Q","",rangeMin,rangeMax) ;
      hpcopy->Fit(fit1,"MQ","",rangeMin,rangeMax) ;

      fit2->SetParameters(fit1->GetParameters()) ;
      fit2->SetParameter(7,0.) ;
//       fit2->SetParLimits(0,0.000,5.*hpcopy->GetMaximum()) ;
      fit2->SetParLimits(1,0.125,0.145) ;
      fit2->SetParLimits(2,0.0035,0.010) ;
      
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
      
      fbg1->SetParameters(fit1->GetParameter(3),fit1->GetParameter(4)); 
      fbg2->SetParameters(fit2->GetParameter(3),fit2->GetParameter(4),fit2->GetParameter(5)); 
      
      Double_t intRangeMin = PeakPosition(pt)-3.*PeakWidth(pt) ;
      Double_t intRangeMax = PeakPosition(pt)+3.*PeakWidth(pt) ;
      Int_t    intBinMin   = hp->GetXaxis()->FindBin(intRangeMin) ;
      Int_t    intBinMax   = hp->GetXaxis()->FindBin(intRangeMax) ;
      Double_t errStat =0;
      hpm->IntegralAndError(intBinMin,intBinMax,errStat); 

      hpm ->Multiply(fbg1) ;
      hpm2->Multiply(fbg2) ;
      hp  ->Add(hpm ,-1.) ;
      hp2 ->Add(hpm2,-1.) ;
                
      
      Int_t binPi0 = hp->FindBin(kMean);
//      fgs->SetParameters(hp->Integral(binPi0-1,binPi0+1)/3.,fun->GetParameter(1),fun->GetParameter(2)) ;
      fgs->SetParameters(1,0.135,0.005) ;
      fgs->SetParLimits(0,0.000,50.*hp->GetMaximum()) ;
      fgs->SetParLimits(1,0.120,0.175) ;
      fgs->SetParLimits(2,0.0035,0.010) ;
      hp->Fit(fgs,"Q","",rangeMin,rangeMax) ;   
      hp->Fit(fgs,"QM","",rangeMin,rangeMax) ;   
      hp->SetMaximum(hp->GetMaximum()*1.4) ;
      hp->SetMinimum(hp->GetMinimum()*1.1) ;
      hp->SetMarkerStyle(20) ;
      hp->SetMarkerSize(0.7) ;
      mr1[iPID]->SetBinContent(i,fgs->GetParameter(1)) ;
      mr1[iPID]->SetBinError  (i,fgs->GetParError(1) ) ;
      sr1[iPID]->SetBinContent(i,TMath::Abs(fgs->GetParameter(2))) ;
      sr1[iPID]->SetBinError  (i,fgs->GetParError(2) ) ;

      Double_t y=fgs->Integral(intRangeMin,intRangeMax)/hp->GetXaxis()->GetBinWidth(1) ;
      nr1[iPID]->SetBinContent(i,y/nPrim) ;
      Double_t ey=0 ;
      if(fgs->GetParameter(0)!=0. && fgs->GetParameter(2)!=0.){
        Double_t en=fgs->GetParError(0)/fgs->GetParameter(0) ;
        Double_t es=fgs->GetParError(2)/fgs->GetParameter(2) ;
        ey=y*TMath::Sqrt(en*en+es*es) ;
      }
      nr1[iPID]->SetBinError(i,ey/nPrim) ;
      
      
      
      //printf("   Na1 =%f+-%f \n",y,ey) ;
      Double_t enpiInt=0 ;
      Double_t npiInt = hp->IntegralAndError(intBinMin,intBinMax,enpiInt) ;
      Double_t norm   = fbgP0->GetParameter(0) ;
      Double_t normErr= fbgP0->GetParError(0) ;
      if(npiInt>0.){
        nr1int[iPID]->SetBinContent(i,npiInt/nPrim) ;
        nr1int[iPID]->SetBinError(i,TMath::Sqrt(enpiInt*enpiInt + 0*norm*errStat + 0*normErr*normErr*errStat*errStat + 0*norm*norm*errStat)/nPrim) ;
      }
      
//printf(" Nint1 =%f+-%f \n",npiInt,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
      hp2->GetXaxis()->SetRangeUser(0.05,0.24) ;
      fgs->SetParLimits(0,0.000,50.*hp2->GetMaximum()) ;
      hp2->SetMaximum(hp2->GetMaximum()*1.4) ;
      hp2->SetMinimum(hp2->GetMinimum()*1.1) ;
      hp2->SetMarkerStyle(20) ;
      hp2->SetMarkerSize(0.7) ;
      hp2->Fit(fgs,"Q","",rangeMin,rangeMax) ;
      hp2->Fit(fgs,"QM","",rangeMin,rangeMax) ;
      mr2[iPID]->SetBinContent(i,fgs->GetParameter(1)) ;
      mr2[iPID]->SetBinError  (i,fgs->GetParError(1)) ;
      sr2[iPID]->SetBinContent(i,TMath::Abs(fgs->GetParameter(2))) ;
      sr2[iPID]->SetBinError  (i,fgs->GetParError(2)) ;

      y=fgs->Integral(intRangeMin,intRangeMax)/hp->GetXaxis()->GetBinWidth(1) ;
      nr2[iPID]->SetBinContent(i,y/nPrim) ;
      ey=0 ;
      if(fgs->GetParameter(0)!=0. && fgs->GetParameter(2)!=0.){
        Double_t en=fgs->GetParError(0)/fgs->GetParameter(0) ;
        Double_t es=fgs->GetParError(2)/fgs->GetParameter(2) ;
        ey=y*TMath::Sqrt(en*en+es*es) ;
      }
      nr2[iPID]->SetBinError(i,ey/nPrim) ;
//printf("   Na2 =%f+-%f \n",y,ey) ;
      npiInt = hp2->IntegralAndError(intBinMin,intBinMax,enpiInt) ;
      norm=fbgP1->GetParameter(0) ;
      normErr=fbgP1->GetParError(0) ;
      if(npiInt>0.){
        nr2int[iPID]->SetBinContent(i,npiInt/nPrim) ;
        nr2int[iPID]->SetBinError(i,TMath::Sqrt(enpiInt*enpiInt + 0*norm*errStat + 0*normErr*normErr*errStat*errStat + 0*norm*norm*errStat)/nPrim) ;
      } 
//printf(" Nint2 =%f+-%f \n",npiInt,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
      hp2->SetTitle(Form("%3.1f<p_{T}<%3.1f GeV/c",xa[i-1],xa[i])) ;
      hp2->Draw() ;
      hp->SetMarkerColor(6) ;
      hp->Draw("same") ;
      c2[iPID]->Update() ;

//      delete hp ;
//    delete hp2 ;
//    delete hpcopy ;
//      delete hpm ;
//      delete hpm2 ;
      
    }
  }
  
  
  fout->cd() ;
  for(Int_t iPID=0;iPID<nPID;iPID++){
    nr1[iPID]->Write(0,TObject::kOverwrite) ;
    nr2[iPID]->Write(0,TObject::kOverwrite) ;
    nr1int[iPID]->Write(0,TObject::kOverwrite) ;
    nr2int[iPID]->Write(0,TObject::kOverwrite) ;
  }
//   fout.Close() ;

  fmass->cd() ;
  for(Int_t iPID=0;iPID<nPID;iPID++){
    mr1[iPID]->Write(0,TObject::kOverwrite) ;
    mr2[iPID]->Write(0,TObject::kOverwrite) ;
    sr1[iPID]->Write(0,TObject::kOverwrite) ;
    sr2[iPID]->Write(0,TObject::kOverwrite) ;
  }
//   fmass.Close() ;


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

