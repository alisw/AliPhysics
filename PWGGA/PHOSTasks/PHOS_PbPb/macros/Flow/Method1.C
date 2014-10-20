#include "Methods.h"

const Double_t kMean=0.135 ; //Approximate peak position to facilitate error estimate
Int_t effcor=1; //correct on efficiency

Double_t PeakPosition(Double_t pt){
  //Fit to the measured peak position
  return 4.99292e-003*exp(-pt*9.32300e-001)+1.34944e-001 ;
}
Double_t PeakWidth(Double_t pt){
  //fit to the measured peak width
  Double_t a=0.0068 ;
  Double_t b=0.0025 ;
  Double_t c=0.000319 ;
  return TMath::Sqrt(a*a+b*b/pt/pt+c*c*pt*pt) ;
}
 
Double_t CB(Double_t * x, Double_t * par){
  //Parameterization of Real/Mixed ratio
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+par[3] + par[4]*(x[0]-kMean) /*+ par[5]*(x[0]-kMean)*(x[0]-kMean) */;
}
Double_t CB2(Double_t * x, Double_t * par){
  //Another parameterization of Real/Mixed ratio
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+par[3]+par[4]*(x[0]-kMean) + par[5]*(x[0]-kMean)*(x[0]-kMean) /*+ par[6]*(x[0]-kMean)*(x[0]-kMean)*(x[0]-kMean)*/;
}
Double_t CBs(Double_t * x, Double_t * par){
  //Parameterizatin of signal
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)/TMath::Sqrt(TMath::TwoPi())/s+par[3] ;
}
Double_t BG1(Double_t * x, Double_t * par){
  //Normalizatino of Mixed
  return par[0]  + par[1]*(x[0]-kMean) /*+ par[2]*(x[0]-kMean)*(x[0]-kMean)*/;
}
Double_t BG2(Double_t * x, Double_t * par){
  //Another normalization of  Mixed
  return par[0]+par[1]*(x[0]-kMean) + par[2]*(x[0]-kMean)*(x[0]-kMean) /*+ par[3]*(x[0]-kMean)*(x[0]-kMean)*(x[0]-kMean)*/;
}

void Method1(Int_t cen=1,Int_t side=1){
//side = 1 - VOC, side = 0 - V0A
//side = 2 - TPC 3 sub method, 3 - TPC 2 sub method

  //Fit Real/Mixed ratio, normalize Mixed and subtract it from Real

//  TH1F * hCen = GetCen() ;
//  TH2F * hCenPHOS = GetCenPHOS() ;


  char kind[15], w[25];
  sprintf(w,""); //Weight: no weight - "NW", with weight - ""

  char PID[25] ;
  sprintf(PID,"Both2core") ;
  char key[55];
  char kind[15],kind2[15] ; //Reaction plane
  char d[15]; //detector
  sprintf(kind,"") ; //"" for real 
  sprintf(kind2,"") ; //"" for mixed

  if(side==1)sprintf(d,"V0C");
  else if(side==0)sprintf(d,"V0A");
  else sprintf(d,"TPC");

  TH3D *h3DR, *h3DM;
cout<<"Getting histos from data..."<<endl;
  h3DR = (TH3D*)GetRealMixed(cen,w,d,kind,kind2,PID,1);
  h3DM = (TH3D*)GetRealMixed(cen,w,d,kind,kind2,PID,0);

  TH1F * hCen = GetCen() ;
  TH2F * hCenPHOS = GetCenPHOS() ;
cout<<"Finished"<<endl;

  //  sprintf(kind,"") ; //"" for all
  const Double_t nphi=5;
  Double_t xphi[6];
  for(Int_t i=0;i<=nphi;i++){
	xphi[i]=TMath::PiOver2()*i/nphi ;
//	cout<<xphi[i]<<" ";
  }
//cout<<endl;

  if(cen==21)
  sprintf(inname,"../flow11h_Apr14.root");
  //Read resolution
  TFile * f = new TFile(inname) ;

  TH2F * hcos2AC = (TH2F*)f->Get("cos2V0AC");
  TH2F * hcos2FAC= (TH2F*)f->Get("cos2FAC");  
  Double_t resAC = GetCos(hcos2AC,cen);
  resAC=resAC ;

  hcos2AC = (TH2F*)f->Get("cos2V0ATPC");
  Double_t resAT = GetCos(hcos2AC,cen);
  resAT=resAT ;

  hcos2AC = (TH2F*)f->Get("cos2V0CTPC");
  Double_t resCT = GetCos(hcos2AC,cen);
  resCT=resCT ;

  hcos2AC = (TH2F*)f->Get("cos2AC");
  Double_t resT = GetRes(hcos2AC,cen);

  //side = 1 - VOC, side = 0 - V0A
  //side = 2 - TPC 3 sub method, 3 - TPC 2 sub method

  Double_t res=1;
  if(side==1) res = TMath::Sqrt(resAC*resCT/resAT);
  else if(side==0)res = TMath::Sqrt(resAC*resAT/resCT);
  else if(side==2)res = TMath::Sqrt(resAT*resCT/resAC);
  else res = resT;
  
//  res=1./res ;

  cout<<"RP resolution = "<< res <<endl;

  PPRstyle();
cout<<"Init... ";  
  //Fit real only 
  //Linear Bg
  char key2[155];
  sprintf(key,"PID%s_cen%d",kind,cen) ;
  sprintf(key2,"%s_mr1",key) ;
  TH1D * mr1 = new TH1D(key2,"Mass",nbin,xa) ;
  sprintf(key2,"%s_sr1",key) ;
  TH1D * sr1 = new TH1D(key2,"Width",nbin,xa) ;
  sprintf(key2,"%s_ar1",key) ;
  TH2D * ar1 = new TH2D(key2,"a",nbin,xa,nphi,xphi) ;
  sprintf(key2,"%s_br1",key) ;
  TH2D * br1 = new TH2D(key2,"b",nbin,xa,nphi,xphi) ;
  sprintf(key2,"%s_yr1",key) ;
  TH2D * nr1 = new TH2D(key2,"Raw yield",nbin,xa,nphi,xphi) ;
  sprintf(key2,"%s_yr1int",key) ;
  TH2D * nr1int = new TH2D(key2,"Raw yield, integrated",nbin,xa,nphi,xphi) ;

  //Quadratic Bg
  sprintf(key2,"%s_mr2",key) ;
  TH1D * mr2 = new TH1D(key2,"Mass",nbin,xa) ;
  sprintf(key2,"%s_sr2",key) ;
  TH1D * sr2 = new TH1D(key2,"Width",nbin,xa) ;
  sprintf(key2,"%s_ar2",key) ;
  TH2D * ar2 = new TH2D(key2,"a",nbin,xa,nphi,xphi) ;
  sprintf(key2,"%s_br2",key) ;
  TH2D * br2 = new TH2D(key2,"b",nbin,xa,nphi,xphi) ;
  sprintf(key2,"%s_cr2",key) ;
  TH2D * cr2 = new TH2D(key2,"c",nbin,xa,nphi,xphi) ;
  sprintf(key2,"%s_yr2",key) ;
  TH2D * nr2 = new TH2D(key2,"Raw yield",nbin,xa,nphi,xphi) ;
  sprintf(key2,"%s_yr2int",key) ;

  TH2D * nr2int = new TH2D(key2,"Raw yield, integrated",nbin,xa,nphi,xphi) ;

  TF1 * fit1 = new TF1("fit1",CB,0.,1.,5) ;
  fit1->SetParName(0,"A") ;
  fit1->SetParName(1,"m_{0}") ;
  fit1->SetParName(2,"#sigma") ;
  fit1->SetParName(3,"a_{0}") ;
  fit1->SetParName(4,"a_{1}") ;
  fit1->SetLineWidth(2) ;
  fit1->SetLineColor(2) ;
  TF1 * fgs = new TF1("gs",CBs,0.,1.,4) ;
  fgs->SetLineColor(2) ;
  fgs->SetLineWidth(1) ;

  TF1 * fit2 = new TF1("fit2",CB2,0.,1.,6) ;
  fit2->SetParName(0,"A") ;
  fit2->SetParName(1,"m_{0}") ;
  fit2->SetParName(2,"#sigma") ;
  fit2->SetParName(3,"a_{0}") ;
  fit2->SetParName(4,"a_{1}") ;
  fit2->SetParName(5,"a_{2}") ;

  fit2->SetLineWidth(2) ;
  fit2->SetLineColor(4) ;
  fit2->SetLineStyle(2) ;
  TF1 * fbg1 = new TF1("bg1",BG1,0.,1.,2) ;
  TF1 * fbg2 = new TF1("bg2",BG2,0.,1.,3) ;
  
  //calculate average mass and width for all pt bins
  TCanvas * cav = new TCanvas("Average","Average",10,10,800,750) ;
  cav->Divide(3,4) ;

cout<<"Finished"<<endl;

cout<<"Do avarage minv... ";  
  for(Int_t i=1;i<=nbin;i++){
    h3DR->GetYaxis()->SetRangeUser(xa[i-1],xa[i]) ;
    h3DM->GetYaxis()->SetRangeUser(xa[i-1],xa[i]) ;
    Double_t pt=(xa[i]+xa[i-1])/2. ;    
    TH2D * hpav = (TH2D*)h3DR->Project3D("x");
    TH2D * hpmav= (TH2D*)h3DM->Project3D("x") ;
    hpav->Sumw2() ;
    hpmav->Sumw2() ;
    if(i<7){
      hpav->Rebin(2) ;
      hpmav->Rebin(2) ;
    }
    else{
      hpav->Rebin(2) ;
      hpmav->Rebin(2) ;
    }
    for(Int_t ib=1; ib<=hpav->GetNbinsX();ib++){if(hpav->GetBinContent(ib)==0)hpav->SetBinError(ib,1.);}
    for(Int_t ib=1; ib<=hpav->GetNbinsX();ib++){if(hpmav->GetBinContent(ib)==0)hpmav->SetBinError(ib,1.);}
    TH1D * hpm2av = (TH1D*)hpmav->Clone("Bg1av") ;
    TH1D * hpcopyav = (TH1D*)hpav->Clone("hpcopyav") ;
    TH1D * hp2av = (TH1D*)hpav->Clone("hp2av") ;
    TH1D * hpm3av = (TH1D*)hpmav->Clone("Bg2av") ;
    TH1D * hp3av = (TH1D*)hpav->Clone("hp3av") ;

    hpcopyav->Divide(hpmav) ;
    sprintf(key,"%3.1f<p_{t}<%3.1f GeV/c",xa[i-1],xa[i]) ;
    hpcopyav->SetTitle(key) ;
    hpcopyav->SetMarkerStyle(20) ;
    hpcopyav->SetMarkerSize(0.7) ;
    
    fit1->SetParameters(TMath::Min(0.3,0.001*i),0.14,0.01,0.008,-0.0002,0) ;

    fit1->SetParLimits(2,0.003,0.03) ; //peak width
    fit1->SetParLimits(1,0.1,0.2); //peak mean
    fit1->SetParLimits(0,0.,10.) ; //peak height

    cav->cd(i) ;
    Double_t rangeMin=0.1 ;
    Double_t rangeMax=0.25 ; //TMath::Min(0.15+i*0.05,0.3) ;
//    if(cen==0)rangeMax=0.30; //TMath::Min(0.15+i*0.05,0.3) ;
    hpcopyav->Fit("fit1","Q","",rangeMin,rangeMax) ;
    hpcopyav->Fit("fit1","MQ","",rangeMin,rangeMax) ;

    mr1->SetBinContent(i,fit1->GetParameter(1)) ;
    mr1->SetBinError(i,fit1->GetParError(1)) ;
    sr1->SetBinContent(i,TMath::Abs(fit1->GetParameter(2))) ;
    sr1->SetBinError(i,fit1->GetParError(2)) ;

    fit2->SetParameters(fit1->GetParameters()) ;
    fit2->SetParameter(5,0) ;  
    hpcopyav->Fit("fit2","+Q","",rangeMin,rangeMax) ;
    hpcopyav->Fit("fit2","+MQ","",rangeMin,rangeMax) ;
    hpcopyav->GetXaxis()->SetRangeUser(0.05,0.3) ;

    cav->Update() ;

    mr2->SetBinContent(i,fit2->GetParameter(1)) ;
    mr2->SetBinError(i,fit2->GetParError(1)) ;
    sr2->SetBinContent(i,TMath::Abs(fit2->GetParameter(2))) ;
    sr2->SetBinError(i,fit2->GetParError(2)) ;
        

    delete hpav ;
//    delete hpcopyav ;
    delete hp2av ;
    delete hpmav; 
    delete hpm2av; 
  }
cout<<"Finished"<<endl;
      
  h3DR->GetYaxis()->SetRangeUser(0.,19.) ;
  h3DM->GetYaxis()->SetRangeUser(0.,19.) ;

cout<<"Do dNdPhi bins... ";  
//!!!!!!!!!!!!!!!!!DNDPHI bins:

  for(Int_t iphi=1;iphi<=nphi;iphi++){

    h3DR->GetZaxis()->SetRange(iphi,iphi) ;
    h3DM->GetZaxis()->SetRange(iphi,iphi) ;

    TH2D * h = (TH2D*)h3DR->Project3D("yx");
    TH2D * hm= (TH2D*)h3DM->Project3D("yx") ;
    h->SetName("sliceRe") ;
    hm->SetName("sliceMi") ;

    h3DR->GetZaxis()->SetRange(11-iphi,11-iphi) ;
    h3DM->GetZaxis()->SetRange(11-iphi,11-iphi) ;

    TH2D * ha = (TH2D*)h3DR->Project3D("yx");
    TH2D * hma= (TH2D*)h3DM->Project3D("yx") ;

    h->Add(ha) ; delete ha ;
    hm->Add(hma); delete hma ;

    sprintf(key2,"mggFit%d_Signal",iphi) ;
    TCanvas * c3 = new TCanvas(key2,key2,10,10,800,750) ;
    c3->Divide(3,4) ;
    
    sprintf(key2,"mggFit%d",iphi) ;
    TCanvas * c1 = new TCanvas(key2,key2,10,10,800,750) ;
    c1->Divide(3,4) ;
    c1->cd(0) ;

//    TCanvas * c2=0,*c4=0,*c5=0,*c6=0 ; 
    TAxis * pta=h->GetYaxis() ;
    TAxis * ma=h->GetXaxis() ;

//loop over pT bins
    for(Int_t i=1;i<=nbin;i++){
      c1->cd(i) ;
      Int_t imin=pta->FindBin(xa[i-1]+0.0001);
      Int_t imax=pta->FindBin(xa[i]-0.0001) ;
      Double_t pt=(xa[i]+xa[i-1])/2. ;
      TH1D * hp = h->ProjectionX("re",imin,imax) ;
      hp->Sumw2() ;
      TH1D * hpm= hm->ProjectionX("mi",imin,imax) ;
      hpm->Sumw2() ;

      if(i<7){
        hp->Rebin(2) ;
        hpm->Rebin(2) ;
      }
      else{
        hp->Rebin(2) ;
        hpm->Rebin(2) ;
      }

      for(Int_t ib=1; ib<=hp->GetNbinsX();ib++){if(hp->GetBinContent(ib)==0)hp->SetBinError(ib,1.);}
      for(Int_t ib=1; ib<=hp->GetNbinsX();ib++){if(hpm->GetBinContent(ib)==0)hpm->SetBinError(ib,1.);}
      TH1D * hpm2 = (TH1D*)hpm->Clone("Bg1") ;
      TH1D * hpcopy = (TH1D*)hp->Clone("hpcopy") ;
      TH1D * hp2 = (TH1D*)hp->Clone("hp2") ;
      hpcopy->Divide(hpm) ;
      sprintf(key,"%3.1f<p_{t}<%3.1f GeV/c",xa[i-1],xa[i]) ;
      hpcopy->SetTitle(key) ;
      hpcopy->SetMarkerStyle(20) ;
      hpcopy->SetMarkerSize(0.7) ;
      
      fit1->SetParameters(TMath::Min(0.3,0.0001*i),0.135,0.008,0.008,-0.0002,0) ;
      fit1->FixParameter(1,mr1->GetBinContent(i)) ;
      fit1->FixParameter(2,sr1->GetBinContent(i)) ;

//      fit1->SetParLimits(2,0.005,0.015) ;
//      fit1->SetParLimits(0,0.,1.) ;
      
      Double_t rangeMin=0.08 ;
      Double_t rangeMax=0.25 ; //TMath::Min(0.15+i*0.05,0.3) ;
//      if(cen==0)rangeMax=0.25 ;

      hpcopy->Fit("fit1","NQ","",rangeMin,rangeMax) ;
      hpcopy->Fit("fit1","MQ","",rangeMin,rangeMax) ;

      ar1->SetBinContent(i,iphi,fit1->GetParameter(3)) ;
      ar1->SetBinError(i,iphi,fit1->GetParError(3)) ;
      br1->SetBinContent(i,iphi,fit1->GetParameter(4)) ;
      br1->SetBinError(i,iphi,fit1->GetParError(4)) ;
      
      fit2->SetParameters(fit1->GetParameters()) ;
      fit2->FixParameter(1,mr2->GetBinContent(i)) ;
      fit2->FixParameter(2,sr2->GetBinContent(i)) ;
//      fit2->SetParLimits(2,0.005,0.015) ;
//      fit2->SetParLimits(0,0.,1.) ;
      fit2->SetParameter(5,0) ;
      
      hpcopy->Fit("fit2","+NQ","",rangeMin,rangeMax) ;
      hpcopy->Fit("fit2","+MQ","",rangeMin,rangeMax) ;

      ar2->SetBinContent(i,iphi,fit2->GetParameter(3)) ;
      ar2->SetBinError(i,iphi,fit2->GetParError(3)) ;
      br2->SetBinContent(i,iphi,fit2->GetParameter(4)) ;
      br2->SetBinError(i,iphi,fit2->GetParError(4)) ;
      cr2->SetBinContent(i,iphi,fit2->GetParameter(5)) ;
      cr2->SetBinError(i,iphi,fit2->GetParError(5)) ;

      hpcopy->GetXaxis()->SetRangeUser(0.05,0.3) ;
//      hpcopy->Draw() ;
/*      if(c5)
	c5->Update() ;
      else
	if(c2)
	  c2->Update() ;
	else
*/
	  c1->Update() ;
      //    if(getchar()=='q')return ;
      
      
      fbg1->SetParameters(fit1->GetParameter(3),fit1->GetParameter(4),fit1->GetParameter(5)); 
      fbg2->SetParameters(fit2->GetParameter(3),fit2->GetParameter(4),fit2->GetParameter(5),fit2->GetParameter(6)); 
            
      hpm->Multiply(fbg1) ;
      hpm2->Multiply(fbg2) ;
      hp->Add(hpm,-1.) ;
      hp2->Add(hpm2,-1.) ;
      
      c3->cd(i) ;
      
      fgs->SetParameters(hp->Integral(13,15)/3.,fit1->GetParameter(1),fit1->GetParameter(2),0.) ;
      fgs->SetParLimits(1,0.05,0.15);
      fgs->FixParameter(1,mr1->GetBinContent(i)) ;
      fgs->FixParameter(2,sr1->GetBinContent(i)) ;
      
      hp->Fit(fgs,"Q","",rangeMin,rangeMax) ;   
      hp->SetMaximum(hp2->GetMaximum()*1.1) ;
      hp->SetMinimum(hp2->GetMinimum()*1.1) ;
      hp->SetMarkerStyle(20) ;
      hp->SetMarkerSize(0.7) ;
      
      //0.96 - due to 2 sigma
      nr1->SetBinContent(i,iphi,fgs->GetParameter(0)) ;
      nr1->SetBinError(i,iphi,fgs->GetParError(0)) ;      
      
      Int_t intBinMin=hp->GetXaxis()->FindBin(fgs->GetParameter(1)-4.*TMath::Abs(fgs->GetParameter(2))) ;
      Int_t intBinMax=hp->GetXaxis()->FindBin(fgs->GetParameter(1)+4.*TMath::Abs(fgs->GetParameter(2))) ;
      Double_t errStat=hpm->Integral(intBinMin,intBinMax); 
      Double_t npiInt=hp->Integral(intBinMin,intBinMax) ;
      Double_t norm=fbg1->GetParameter(0) ;
      Double_t normErr=fbg1->GetParError(0) ;
      if(npiInt>0.){
	nr1int->SetBinContent(i,iphi,npiInt) ;
	nr1int->SetBinError(i,iphi,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
      }
      hp2->GetXaxis()->SetRangeUser(0.05,0.3) ;
      hp2->SetMaximum(hp2->GetMaximum()*1.1) ;
      hp2->SetMinimum(hp2->GetMinimum()*1.1) ;
      hp2->SetMarkerStyle(20) ;
      hp2->SetMarkerSize(0.7) ;
      
      fgs->FixParameter(1,mr2->GetBinContent(i)) ;
      fgs->FixParameter(2,sr2->GetBinContent(i)) ;
      hp2->Fit(fgs,"Q","",rangeMin,rangeMax) ;
      nr2->SetBinContent(i,iphi,fgs->GetParameter(0)) ;
      nr2->SetBinError(i,iphi,fgs->GetParError(0)) ;      
      npiInt=hp2->Integral(intBinMin,intBinMax) ;
      norm=fbg2->GetParameter(0) ;
      normErr=fbg2->GetParError(0) ;
      if(npiInt>0.){
	nr2int->SetBinContent(i,iphi,npiInt) ;
	nr2int->SetBinError(i,iphi,TMath::Sqrt(npiInt + norm*errStat + normErr*normErr*errStat*errStat + norm*norm*errStat)) ;
      } 
      hp2->SetTitle(key) ;
      hp2->Draw() ;

      c3->Update() ;
      
      delete hp ;
//      delete hp2 ;
//      delete hpcopy ;
      delete hpm ;
      delete hpm2 ;
    }

    delete h ;
    delete hm ;
    
  }
  
cout<<"Fit histos init... ";
  TH1D * v2A1 = new TH1D("v2all1","v_{2} two harmonics",nbin,xa) ;
  TH1D * v2A2 = new TH1D("v2all2","v_{2} two harmonics",nbin,xa) ;
  TH1D * v2A1int = new TH1D("v2all1int","v_{2} two harmonics",nbin,xa) ;
  TH1D * v2A2int = new TH1D("v2all2int","v_{2} two harmonics",nbin,xa) ;

  TH1D * v2F1 = new TH1D("v2gap1","v_{2} two harmonics",nbin,xa) ;
  TH1D * v2F2 = new TH1D("v2gap2","v_{2} two harmonics",nbin,xa) ;
  TH1D * v2F1int = new TH1D("v2gap1int","v_{2} two harmonics",nbin,xa) ;
  TH1D * v2F2int = new TH1D("v2gap2int","v_{2} two harmonics",nbin,xa) ;

  TH1D * v2A1s= new TH1D("v2all1S","v_{2} single harmonic",nbin,xa) ;
  TH1D * v2A2s= new TH1D("v2all2S","v_{2} single harmonic",nbin,xa) ;
  TH1D * v2A1ints= new TH1D("v2all1intS","v_{2} single harmonic",nbin,xa) ;
  TH1D * v2A2ints= new TH1D("v2all2intS","v_{2} single harmonic",nbin,xa) ;

  TH1D * v2F1s= new TH1D("v2gap1S","v_{2} single harmonic",nbin,xa) ;
  TH1D * v2F2s= new TH1D("v2gap2S","v_{2} single harmonic",nbin,xa) ;
  TH1D * v2F1ints= new TH1D("v2gap1intS","v_{2} single harmonic",nbin,xa) ;
  TH1D * v2F2ints= new TH1D("v2gap2intS","v_{2} single harmonic",nbin,xa) ;


  TH2D * hA1 = nr1->Clone("hA1");
  TH2D * hA2 = nr2->Clone("hA2");
  TH2D * hA1int = nr1int->Clone("hA1int");
  TH2D * hA2int = nr2int->Clone("hA2int");

cout<<"Finished"<<endl;

  Int_t col[4]={kRed,kBlue,kGreen+3,kMagenta} ;
  Int_t sym[4]={20,21,24,25} ;

  //  TF1 * fit = new TF1("fit","[0]*(1+2.*[1]*cos(2*(x-[2])))",0.,5.) ;
  TF1 * v2fit = new TF1("v2fit","[0]*(1+2.*[1]*cos(2.*x)+2.*[2]*cos(4.*x) )",0.,5.) ;
  v2fit->SetLineColor(2) ;
  TF1 * v2fit2 = new TF1("v2fit","[0]*(1+2.*[1]*cos(2.*x) )",0.,5.) ;
  TCanvas * c = new TCanvas("fit","fit") ;
  c->Divide(3,4);

cout<<"Fit dNdPhi... ";  
  for(Int_t i=1; i<=v2A1->GetNbinsX() ;i++){ //over pt
    c->cd(i);
    char hname[100];
    sprintf(hname,"A1%d",i);
    TH1D * tmp0 = hA1->ProjectionY(hname,i,i) ;
    tmp0->SetTitle(Form("%3.1f<p_{t}<%3.1f GeV",xa[i-1],xa[i])) ;
    tmp0->SetXTitle("#Delta#phi (rad)") ;
    tmp0->SetYTitle("dN/d#phi (a.u.)") ;
    tmp0->SetMarkerStyle(sym[0]) ;
    tmp0->SetMarkerColor(col[0]) ;
    
    v2fit2->SetLineStyle(1) ;
    v2fit2->SetLineColor(col[0]) ;
    v2fit2->SetParameters(tmp0->GetBinContent(5),0.1) ;
    tmp0->Fit(v2fit2,"Qi") ;
    v2A1s->SetBinContent(i,v2fit2->GetParameter(1)) ;
    v2A1s->SetBinError(i,v2fit2->GetParError(1)) ;
    
    
    sprintf(hname,"A2%d",i);
    //tmp2d = (TH2D*)hA2->Clone(hname);
    TH1D * tmp1 = hA2->ProjectionY(hname,i,i) ;
    tmp1->SetMarkerStyle(sym[1]) ;
    tmp1->SetMarkerColor(col[1]) ;
    v2fit2->SetLineColor(col[1]) ;
    v2fit2->SetParameters(tmp1->GetBinContent(5),0.1) ;
    v2fit2->SetLineColor(col[1]) ;
    tmp1->Fit(v2fit2,"Qi") ;
    v2A2s->SetBinContent(i,v2fit2->GetParameter(1)) ;
    v2A2s->SetBinError(i,v2fit2->GetParError(1)) ;

    sprintf(hname,"A1int%d",i);
    TH1D * tmp2 = hA1int->ProjectionY(hname,i,i) ;
    tmp2->SetMarkerStyle(sym[2]) ;
    tmp2->SetMarkerColor(col[2]) ;
    v2fit2->SetParameters(tmp2->GetBinContent(5),0.1) ;
    v2fit2->SetLineStyle(5) ;
    v2fit2->SetLineColor(col[2]) ;
    tmp2->Fit(v2fit2,"i") ;
    v2A1ints->SetBinContent(i,v2fit2->GetParameter(1)) ;
    v2A1ints->SetBinError(i,v2fit2->GetParError(1)) ;

    sprintf(hname,"A2int%d",i);
    //tmp2d = (TH2D*)hA2int->Clone(hname);
    TH1D * tmp3 = hA2int->ProjectionY(hname,i,i) ;
    tmp3->SetMarkerStyle(sym[3]) ;
    tmp3->SetMarkerColor(col[3]) ;
    v2fit2->SetParameters(tmp3->GetBinContent(5),0.1) ;
    v2fit2->SetLineColor(col[3]) ;
    tmp3->Fit(v2fit2,"Qi") ;
    v2A2ints->SetBinContent(i,v2fit2->GetParameter(1)) ;
    v2A2ints->SetBinError(i,v2fit2->GetParError(1)) ;
    
    tmp0->Draw() ;
    tmp1->Draw("same") ;
//    tmp2->Scale(6.28) ;
//    tmp3->Scale(6.28) ;
    tmp2->Draw("same") ;
    tmp3->Draw("same") ;
  }
cout<<"Finished"<<endl;

  TH1D *effcorrect = v2A1s->Clone("effcorrect");
  TH1D *effcorrect2 = v2A1s->Clone("effcorrect2");

  Double_t eff=0, err=0;
  if(strcmp(PID,"Both2core")==0){
        if(cen==0) {eff=.03; err=.04;}
        if(cen==1) {eff=.08; err=.04;}
        if(cen==2) {eff=.03; err=.03;}
        if(cen==3) {eff=.02; err=.02;}
        if(cen==4) {eff=.02; err=.02;}
        if(cen==5) {eff=.01; err=.02;}
        if(cen==10) {eff=.05; err=.02;}
        if(cen==11) {eff=.02; err=.02;}
        if(cen==20) {eff=.02; err=.02;}
        if(cen==21) {eff=.05; err=.02;}

        eff=eff*2.*TMath::Pi()/5./sin(2.*TMath::Pi()/5.)/4.;
        err=err*2.*TMath::Pi()/5./sin(2.*TMath::Pi()/5.)/4.;

cout<<"res: T="<<resT<<", resolution="<<res<<endl;

	eff=eff*res/resT;
	err=err*res/resT;
        cout<<"Pi0 efficiency correction before resolution correction: "<<eff<<"+-"<<err<<endl;
  }

  Double_t x[100],y[100],ex[100],ey[100] ;

  for(Int_t i=0;i<effcorrect->GetNbinsX();i++){
        effcorrect->SetBinContent(i,eff);
        effcorrect->SetBinError(i,err);
        effcorrect2->SetBinContent(i,eff);
        effcorrect2->SetBinError(i,0);
  }

/*
  if(effcor){
	v2A1s->Add(effcorrect);
        v2A2s->Add(effcorrect);
        v2A1ints->Add(effcorrect);
        v2A2ints->Add(effcorrect);
  }
  
  v2A1s->Scale(res) ;
  v2A2s->Scale(res) ;
  v2A1ints->Scale(res) ;
  v2A2ints->Scale(res) ;
*/

//return 1;
  
  v2A1s->SetMarkerStyle(sym[0]) ;
  v2A1s->SetMarkerColor(col[0]) ;
  v2A2s->SetMarkerStyle(sym[1]) ;
  v2A2s->SetMarkerColor(col[1]) ;
  v2A1ints->SetMarkerStyle(sym[2]) ;
  v2A1ints->SetMarkerColor(col[2]) ;
  v2A2ints->SetMarkerStyle(sym[3]) ;
  v2A2ints->SetMarkerColor(col[3]) ;

  TCanvas * c2 = new TCanvas("v2","v2") ;
  TH1D * box = new TH1D("box","V_{2}{RP} with different fit methods",150,0.,20.) ;
  box->SetMinimum(0.) ;
  box->SetMaximum(0.5) ;
  box->SetXTitle("p_{t}^{#pi} (GeV/c)") ;
  box->SetYTitle("v{2}") ;
  box->Draw() ;
  v2A1s->Draw("same") ;
  v2A1ints->Draw("same") ;
  v2A2s->Draw("same") ;
  v2A2ints->Draw("same") ;
    
  TLegend * l = new TLegend(0.6,0.7,0.9,0.9) ;
  l->AddEntry(v2A1s,"Fit, pol1","p") ;
  l->AddEntry(v2A2s,"Fit, pol2","p") ;
  l->AddEntry(v2A1ints,"Int, pol1","p") ;
  l->AddEntry(v2A2ints,"Int, pol2","p") ;

  l->Draw() ;
  
//  return ;

  TCanvas * cmass = new TCanvas("Mass","Pi0 mass vs pT") ;

  mr1->SetMarkerStyle(20);
  mr2->SetMarkerStyle(21);
  mr1->Draw();
  mr2->SetLineColor(2);
  mr2->SetMarkerColor(2);
  mr2->Draw("same");

  TLegend * l = new TLegend(0.6,0.7,0.9,0.9) ;
  l->AddEntry(mr1,"pol1","p") ;
  l->AddEntry(mr2,"pol2","p") ;

  l->Draw() ;

  TCanvas * cwidth = new TCanvas("Width","Pi0 width vs pT") ;

  sr1->SetMarkerStyle(20);
  sr2->SetMarkerStyle(21);
  sr1->Draw();
  sr2->SetLineColor(2);
  sr2->SetMarkerColor(2);
  sr2->Draw("same");

  TLegend * l = new TLegend(0.6,0.7,0.9,0.9) ;
  l->AddEntry(mr1,"pol1","p") ;
  l->AddEntry(mr2,"pol2","p") ;

  l->Draw() ;

//!!!!!!!!!!!!!!!!!!!!!!final v2 plot
  TCanvas * c5 = new TCanvas("v2_res","v2 res") ;


  TH1D * v2Astat = new TH1D("v2allStat","v_{2}",nbin,xa) ;
  TH1D * v2Asys = new TH1D("v2allSys","v_{2}",nbin,xa) ;
  
  for(Int_t i=1; i<=v2A1->GetNbinsX() ;i++){
    Double_t mean=
       v2A1s->GetBinContent(i)/v2A1s->GetBinError(i)/v2A1s->GetBinError(i)
      +v2A2s->GetBinContent(i)/v2A2s->GetBinError(i)/v2A2s->GetBinError(i)
      +v2A1ints->GetBinContent(i)/v2A1ints->GetBinError(i)/v2A1ints->GetBinError(i) 
      +v2A2ints->GetBinContent(i)/v2A2ints->GetBinError(i)/v2A2ints->GetBinError(i) ;
    Double_t weight=
       1./v2A1s->GetBinError(i)/v2A1s->GetBinError(i)
      +1./v2A2s->GetBinError(i)/v2A2s->GetBinError(i) 
      +1./v2A1ints->GetBinError(i)/v2A1ints->GetBinError(i) 
      +1./v2A2ints->GetBinError(i)/v2A2ints->GetBinError(i) ;

    mean/=weight ;

    Double_t rms = (v2A1s->GetBinContent(i) - mean)*(v2A1s->GetBinContent(i) - mean)/v2A1s->GetBinError(i)/v2A1s->GetBinError(i) + (v2A2s->GetBinContent(i) - mean)*(v2A2s->GetBinContent(i) - mean)/v2A2s->GetBinError(i)/v2A2s->GetBinError(i) + (v2A1ints->GetBinContent(i) - mean)*(v2A1ints->GetBinContent(i) - mean)/v2A1ints->GetBinError(i)/v2A1ints->GetBinError(i) + (v2A2ints->GetBinContent(i) - mean)*(v2A2ints->GetBinContent(i) - mean)/v2A2ints->GetBinError(i)/v2A2ints->GetBinError(i);
     Double_t statErr = TMath::Min(v2A1s->GetBinError(i),v2A2s->GetBinError(i)) ;
//				   TMath::Min(v2A1ints->GetBinError(i),v2A2ints->GetBinError(i))) ;
cout<<rms<<endl;

    rms/=weight;

    Double_t ptbin = v2A1->GetBinCenter(i);
    Double_t sys=TMath::Sqrt(rms);

    cout<<"Sys error for pT="<<ptbin<<" is "<<sys<<" in percent: "<<sys/mean<<endl;
    cout<<"Sys error due to eff corr: "<<eff/res<<" in percent: "<<eff/res/mean<<endl;

    Double_t mean_corr = mean ;// - v2ch * Ntrack_cen7080 / Ntrack_cen;
    
    v2Astat->SetBinContent(i,mean_corr) ;
    v2Astat->SetBinError(i,statErr) ;
    v2Asys->SetBinContent(i,mean_corr) ;
    v2Asys->SetBinError(i,sys) ;

  }

  if(effcor){
        v2Astat->Add(effcorrect2);
        v2Asys->Add(effcorrect);
  }

  v2Asys->Scale(1/res) ;
  v2Astat->Scale(1/res) ;

  v2Asys->SetFillStyle(1001) ;
  v2Asys->SetFillColor(kYellow) ;
  
  v2Astat->SetMarkerColor(kOrange+7) ;
  v2Astat->SetMarkerStyle(20) ;
  v2Astat->SetLineColor(kOrange+7) ;

  /*
  v2Fsys->SetFillStyle(0) ;
  v2Fsys->SetFillColor(0) ;
  v2Fsys->SetLineColor(kBlue+3) ;
  */  

  v2Astat->SetMarkerColor(kOrange+7) ;
  v2Astat->SetMarkerStyle(20) ;
  v2Astat->SetLineColor(kOrange+7) ;

  /*
  v2Fstat->SetMarkerColor(kBlue+3) ;
  v2Fstat->SetMarkerStyle(21) ;
  v2Fstat->SetLineColor(kBlue+3) ;
  */

  v2Asys->SetXTitle("p_{t} (GeV/c)") ;
  v2Asys->SetYTitle("v_{2}") ;
  
  box->Draw() ;
  v2Asys->Draw("sameE1") ;
  //  v2Fsys->Draw("E2same") ;
  v2Astat->Draw("same") ;
  //  v2Fstat->Draw("same") ;

  char nname[255];
  TFile fout("v2_method1_QA.root","update");
  sprintf(nname,"v2sys_m1_%s_%s_%d_%d",w,PID,cen,side) ;
  v2Asys->SetName(nname) ;
  sprintf(nname,"v2stat_m1_%s_%s_%d_%d",w,PID,cen,side) ;
  v2Astat->SetName(nname) ;

  v2Asys->Write(0,TObject::kOverwrite) ;
  v2Astat->Write(0,TObject::kOverwrite) ;
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

