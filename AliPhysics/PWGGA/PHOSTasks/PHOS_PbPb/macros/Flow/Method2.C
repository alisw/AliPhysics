//AA analysis - V2 calculation with Kentaro Miki method 
//27 march 2011. First attempt 
//28 march 2011. DP's fixes

#include "Methods.h"

Int_t effcor = 1;
Int_t drawopt=0;
Bool_t reject=kTRUE;
Double_t mean=0.137;
Double_t width=2*0.007;

//======================================================
Double_t fbg2(Double_t *x, Double_t *par)
{ //pol 4 - to fit v2 bg
  if (reject && x[0] > mean-width && x[0] < mean+width) {
    TF1::RejectPoint();
  }
  return par[0] + par[1]*(x[0] - mean)  + par[2]*(x[0] - mean)*(x[0] - mean) + par[3]*(x[0] - mean)*(x[0] - mean)*(x[0] - mean) + par[4]*(x[0] - mean)*(x[0] - mean)*(x[0] - mean)*(x[0] - mean);
}
//======================================================
Double_t fbg(Double_t *x, Double_t *par)
{ //pol 3 - to fit v2 bg
  if (reject && x[0] > mean-width && x[0] < mean+width) {
    TF1::RejectPoint();
  }
  return par[0] + par[1]*(x[0] - mean) + par[2]*(x[0] - mean)*(x[0] - mean)  + par[3]*(x[0] - mean)*(x[0] - mean)*(x[0] - mean);
}
//======================================================
Double_t fitMgg1(Double_t *x, Double_t *par){
// gauss + pol 1
  Double_t fitval = par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2])+par[3]+par[4]*(x[0]-mean);
  return fitval;
}
//======================================================
Double_t fitBg1(Double_t *x, Double_t *par)
{ // pol 1
  Double_t fitval = par[0]+par[1]*(x[0]-mean);
  return fitval;
}
//======================================================
Double_t fitMgg2(Double_t *x, Double_t *par)
{ //gauss + pol 2
 Double_t fitval = par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2])+par[3]+par[4]*(x[0]-mean)+par[5]*(x[0]-mean)*(x[0]-mean);
  return fitval;
}
//======================================================
Double_t fitBg2(Double_t *x, Double_t *par)
{
  Double_t fitval = par[0]+par[1]*(x[0]-mean)+par[2]*(x[0]-mean)*(x[0]-mean);
  return fitval;
}
//======================================================
Double_t fitV2p1(Double_t *x, Double_t *par){
// gauss + pol 2
//  Double_t fitval = par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2])+
//  par[3]+par[4]*(x[0]-mean)+par[5]*(x[0]-mean)*(x[0]-mean);
Double_t fitval = par[0]*TMath::Exp(-(x[0]-mean)*(x[0]-mean)/2./width/width)+
  par[3]+par[4]*(x[0]-mean)+par[5]*(x[0]-mean)*(x[0]-mean);

  return fitval;
}
//======================================================
Double_t fitV2p2(Double_t *x, Double_t *par){
// gauss + pol 3
//  Double_t fitval = par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2])+
//  par[3]+par[4]*(x[0]-mean)+par[5]*(x[0]-mean)*(x[0]-mean)+par[6]*(x[0]-mean)*(x[0]-mean)*(x[0]-mean);

Double_t fitval = par[0]*TMath::Exp(-(x[0]-mean)*(x[0]-mean)/2./width/width)+
  par[3]+par[4]*(x[0]-mean)+par[5]*(x[0]-mean)*(x[0]-mean)+par[6]*(x[0]-mean)*(x[0]-mean)*(x[0]-mean);

  return fitval;
}

//======================================================
Double_t fitV2mix(Double_t *x, Double_t *par){
// gauss + pol1*mixed
  TH1D * hm = (TH1D*)gROOT->FindObjectAny("v2m") ;
  Double_t mix = hm->GetBinContent(hm->GetXaxis()->FindBin(x[0])) ;
  Double_t fitval = par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2])+
  (par[3]+par[4]*(x[0]-par[1]))*mix;
  return fitval;
}

//======================================================
Double_t fitGauss(Double_t *x, Double_t *par)
{
  Double_t fitval = par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/2./par[2]/par[2]);
  return fitval;
}
//======================================================
void CalcV2(TH1D * h, Double_t &v2, Double_t &v2err){

  v2=0.;
  v2err=0;
  Double_t N=0 ;
  Double_t C=0 ;
  TAxis * x=h->GetXaxis() ;
  for(Int_t i=1;i<=h->GetNbinsX();i++){
    Double_t phi=x->GetBinCenter(i) ;
    Double_t cj=TMath::Cos(2.*phi) ;
    N+=h->GetBinContent(i) ;
    C+=cj*h->GetBinContent(i) ;
  }
  if(N>0.)
    v2=C/N ;
  else 
    return ;
  for(Int_t i=1;i<=h->GetNbinsX();i++){
    Double_t phi=x->GetBinCenter(i) ;
    Double_t cj=TMath::Cos(2.*phi) ;
    v2err+=(N*cj-C)*(N*cj-C)*h->GetBinContent(i)/N/N/N/N ;
  }
  if(v2err<0)v2err=0 ;
  v2err=TMath::Sqrt(v2err) ;

}
//======================================================
void ScaleMixed(TH1D * v2,TH1D * v2m,TH1D * v2mScaled){
//Scale v2m to reproduce v2 around the pi0 peak
  //left side
  Double_t leftMin=0.04 ;
  Double_t leftMax=0.1 ;
  Int_t iMinLeft=v2->GetXaxis()->FindBin(leftMin+0.00001) ;
  Int_t iMaxLeft=v2->GetXaxis()->FindBin(leftMax-0.00001) ;
  
  //right side
  Double_t rightMin=0.16;
  Double_t rightMax=0.4 ;
  Int_t iMinRight=v2->GetXaxis()->FindBin(rightMin+0.00001) ;
  Int_t iMaxRight=v2->GetXaxis()->FindBin(rightMax-0.00001) ;

  Double_t a=0.,b=0.,c=0.,d=0.,e=0. ;
  
  for(Int_t i=iMinLeft; i<=iMaxLeft;i++){
    Double_t rY=v2->GetBinContent(i) ; 
    Double_t eY=v2->GetBinError(i) ; 
    Double_t mY=v2m->GetBinContent(i) ;
    if(eY==0.)continue ;
    a+=rY/eY/eY ;
    b+=1./eY/eY ;
    c+=mY/eY/eY ;
    d+=mY*rY/eY/eY ;
    e+=mY*mY/eY/eY ;
  }  

  for(Int_t i=iMinRight; i<=iMaxRight;i++){
    Double_t rY=v2->GetBinContent(i) ; 
    Double_t eY=v2->GetBinError(i) ; 
    Double_t mY=v2m->GetBinContent(i) ;
    if(eY==0.)continue ;
    a+=rY/eY/eY ;
    b+=1./eY/eY ;
    c+=mY/eY/eY ;
    d+=mY*rY/eY/eY ;
    e+=mY*mY/eY/eY ;
  }  
  
  Double_t aSlope = 1.;//(a*c-d*b)/(c*c-e*b) ;
//  Double_t aSlope = (a*c-d*b)/(c*c-e*b) ;

  Double_t aConst = (a - aSlope*c)/b ;
 
  for(Int_t i=1; i<=v2mScaled->GetNbinsX();i++){
//    v2mScaled->SetBinContent(i+1,aConst+aSlope*v2m->GetBinContent(i)) ;
//    v2mScaled->SetBinError(i+1,aSlope*v2m->GetBinError(i)) ;
    v2mScaled->SetBinContent(i,aConst+aSlope*v2m->GetBinContent(i)) ;
    v2mScaled->SetBinError(i,aSlope*v2m->GetBinError(i)) ;
  }
}
//================================================//MAIN
void Method2(Int_t cen=1,Int_t side=1){
//side = 1 - VOC, side = 0 - V0A
//side = 2 - TPC 3 sub method, 3 - TPC 2 sub method
  
  gStyle->SetFillStyle(1)   ;  
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
//  gStyle->SetOptTitle(0);
  
//  TH2F* hCenTrack = (TH2F*)f->Get("hCenPHOS") ;

  char PID[25] ;
  sprintf(PID,"Both2core") ;
  char key[55];
  char kind[15],kind2[15] ; //Reaction plain
  char d[15] ; //detector
  sprintf(kind,"") ; //"" for real 
  sprintf(kind2,"") ; //"" for mixed
  char w[15]; //weight
  sprintf(w,""); //NW

  if(side==1)sprintf(d,"V0C");
  else if(side==0)sprintf(d,"V0A");
  else sprintf(d,"TPC");

  TH3F *h3DR, *h3DM;
  h3DR = GetRealMixed(cen,w,d,kind,kind2,PID,1);
  h3DM = GetRealMixed(cen,w,d,kind,kind2,PID,0);

  TH1F * hCen = GetCen() ;
  TH2F * hCenPHOS = GetCenPHOS() ;

  if(cen==21)
  sprintf(inname,"../flow11h_Apr14.root");

  TFile * f = new TFile(inname) ;

  //Read resolution
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
  resT=resT ;

  Double_t res=1;
  if(side==1) res = TMath::Sqrt(resAC*resCT/resAT);
  else if(side==0)res = TMath::Sqrt(resAC*resAT/resCT);
  else if(side==2)res = TMath::Sqrt(resAT*resCT/resAC);
  else res = resT;

  cout<<"res = "<<res <<endl;

//--------------------------------------------------------------------
  //invariant mass rebin
  Int_t nRebin=2 ;

  sprintf(key,"v2Pi_mixed%d",cen) ;
  TH1D * v2Pi_mix = new TH1D(key,"V2 pi0 using mixed",nbin,xa) ;
  sprintf(key,"v2Pi_pol1%d",cen) ;
  TH1D * v2Pi_pol1 = new TH1D(key,"V2 pi0 using pol1",nbin,xa) ;
  sprintf(key,"v2Pi_pol2%d",cen) ;
  TH1D * v2Pi_pol2 = new TH1D(key,"V2 pi0 using pol2",nbin,xa) ;

  sprintf(key,"nSig%d",cen) ;
  TH1D * nSig = new TH1D(key,"n Pi0",nbin,xa) ;
  sprintf(key,"nSig2%d",cen) ;
  TH1D * nSig2 = new TH1D(key,"n Pi0",nbin,xa) ;


  sprintf(key,"v2sys_m2_%s_%s_%d_%d",w,PID,cen,side) ;
  TH1D * v2Asys = new TH1D(key,"V2 sys method2",nbin,xa) ;
  sprintf(key,"v2stat_m2_%s_%s_%d_%d",w,PID,cen,side) ;
  TH1D * v2Astat = new TH1D(key,"V2 stat method2",nbin,xa) ;

  TCanvas * cf = new TCanvas("cf","yield vs phi",10,10,400,400) ;
  cf->SetFillColor(0) ;
  cf->SetFillStyle(0) ;
  cf->Range(0,0,1,1);
  cf->SetBorderSize(0);

  TCanvas * cinv = new TCanvas("cinv","inv mass",10,10,400,400) ;

  TCanvas * c4= new TCanvas("c4","v2 vs minv",10,10,800,800);
  c4->SetFillColor(0) ;
  c4->SetFillStyle(0) ;
  c4->Range(0,0,1,1);
  c4->SetBorderSize(0);
  c4->Divide(3,2) ;
  
  TCanvas * c3= new TCanvas("c3","v2 vs minv",10,10,800,800);
  c3->SetFillColor(0) ;
  c3->SetFillStyle(0) ;
  c3->Range(0,0,1,1);
  c3->SetBorderSize(0);
  c3->Divide(3,2) ;

//========================================================== minv
  TF1 *fitMgg1 = new TF1("fitMgg1",fitMgg1,0.,0.4,5);
  fitMgg1->SetLineWidth(2) ;
  fitMgg1->SetLineColor(2) ;
  TF1 *fitMgg2 = new TF1("fitMgg2",fitMgg2,0.,0.4,6);
  fitMgg2->SetLineWidth(2) ;
  fitMgg2->SetLineColor(4) ;
  fitMgg2->SetLineStyle(7) ;
  TF1 *fitGauss = new TF1("fitGauss",fitGauss,0.,0.4,3);
  fitGauss->SetLineWidth(2) ;
  fitGauss->SetLineColor(2) ;

  TF1 *fitV2p1 = new TF1("fitV2p1",fitV2p1,0.,0.4,6);
  fitV2p1->SetLineWidth(2) ;
  fitV2p1->SetLineColor(3) ;
  TF1 *fitV2p2 = new TF1("fitV2p2",fitV2p2,0.,0.4,7);
  fitV2p2->SetLineWidth(2) ;
  fitV2p2->SetLineColor(6) ;
  fitV2p2->SetLineStyle(7) ;
    
  TF1 *fitV2mix = new TF1("fitV2mix",fitV2mix,0.,0.4,5);
  fitV2mix->SetLineWidth(2) ;
  fitV2mix->SetLineColor(kRed) ;
  fitV2mix->SetLineStyle(1) ;

  TF1 *fitBg1 = new TF1("fitBg1",fitBg1,0.,0.4,2);
  TF1 *fitBg2 = new TF1("fitBg2",fitBg2,0.,0.4,3);

  TF1 * fit= new TF1("fit","[0]*(1+2*[1]*cos(2.*(x)))",0.,10.) ;
  TF1 * fitv4= new TF1("fitv4","[0]*(1+2*[1]*cos(2.*(x))+2*[2]*cos(4.*(x)))",0.,10.) ;

  TF1 *f1=new TF1("f1",fbg,0.0,0.4,4);
  f1->SetLineColor(2) ;
  TF1 *f2=new TF1("f2",fbg2,0.0,0.4,5);
  f2->SetLineColor(4) ;
  f2->SetLineStyle(4) ;

  TF1 * unit = new TF1("unit","1.",0.,10.) ;
  
  TF1 * fconst = new TF1("fconst","[0]",0.,1.) ;
  //============================

  char name[255] ;
  TAxis * axis = v2Pi_mix->GetXaxis() ;
  for(Int_t ii=1;ii<=nbin;ii++){ //over pT bins
    if(ii<7)
      c3->cd(ii) ;
    else
      c4->cd(ii-6) ;
    if(ii>4)
       nRebin=2 ;
    cinv->cd() ;    

    Double_t ptmin=axis->GetBinLowEdge(ii);
    Double_t ptmax=axis->GetBinUpEdge(ii);

    h3DR->GetYaxis()->SetRangeUser(ptmin+0.0001,ptmax-0.0001) ;
    h3DR->GetXaxis()->SetRangeUser(0.0,0.5) ;
    
    h3DM->GetYaxis()->SetRangeUser(ptmin+0.0001,ptmax-0.0001) ;
    h3DM->GetXaxis()->SetRangeUser(0.0,0.5) ;

//    h3DR->Rebin(nRebin,"");
//    h3DM->Rebin(nRebin,"");
    
    //V2 vs minv for real
    TH1D* v2 = (TH1D*)h3DR->Project3D("x") ;
    sprintf(name,"v2_%d",ii) ;
    v2->SetName(name) ;
    sprintf(name,"%3.1f<p_{t}<%3.1f GeV/c",ptmin,ptmax) ;
    v2->SetTitle(name) ;
    //V2 vs minv for mixed
    TH1D* v2m = (TH1D*)h3DM->Project3D("x");
//    sprintf(name,"v2m_%d",ii) ;
    sprintf(name,"v2m") ;
    v2m->SetName(name) ;

    TH1D* hRatio =(TH1D*) h3DR->Project3D("x");
    sprintf(name,"Ratio_%d",ii) ;
    hRatio->SetName(name) ;
    TH1D* hInvmassM1d = (TH1D*)h3DM->Project3D("x");
    sprintf(name,"Mixed_%d",ii) ;
    hInvmassM1d->SetName(name) ;
 
    v2->Rebin(nRebin) ;
    v2m->Rebin(nRebin) ;

    sprintf(name,"v2v4_%d",ii) ;
    TH1D * v2v4=(TH1D*)v2->Clone(name);
    sprintf(name,"v2mv4_%d",ii) ;
    TH1D * v2mv4=(TH1D*)v2m->Clone(name);

    hRatio->Rebin(nRebin) ;
    hInvmassM1d->Rebin(nRebin) ;
    
    hRatio->Sumw2();
    hRatio->Divide(hInvmassM1d);
    hRatio->Fit(fconst,"q0","",0.18,0.2) ;
    hRatio->Scale(1./fconst->GetParameter(0)) ;


    if(cen>0)
      fitMgg1->SetParameters(TMath::Min(0.3,0.0003*ii),0.137,0.0107,0.008,-0.0002,0) ;
    else
      fitMgg1->SetParameters(TMath::Min(0.3,0.0003*ii),0.137,0.0077,0.01,-0.0002,0) ;

    fitMgg1->SetParLimits(2,0.005,0.015) ;
    fitMgg1->SetParLimits(1,0.125,0.145);
//    fitMgg1->SetParLimits(0,0.,1.) ;

    cinv->cd() ;
    hRatio->Fit(fitMgg1,"qr0","",0.07,0.25); //fit and draw InvMass real/mixed
    hRatio->Fit(fitMgg1,"qrM0","",0.07,0.25); //fit and draw InvMass real/mixed
    fitMgg2->SetParameters(fitMgg1->GetParameters());
    fitMgg2->SetParameter(5,0.) ;
    hRatio->Fit(fitMgg2,"qr0+","",0.07,0.25); //fit and draw InvMass real/mixed
    cinv->Update(); 
 
    Double_t mean1 = fitMgg1->GetParameter(1) ;
    Double_t mean1_e = fitMgg1->GetParError(1);
    Double_t sigma1 =TMath::Abs(fitMgg1->GetParameter(2)) ;

    mean=mean1;
    width=2.*sigma1;
    cout<< "Bin " << ii <<":  sigma: "<<sigma1<<endl;
    
    fitBg1->SetParameters(fitMgg1->GetParameter(3),fitMgg1->GetParameter(4)) ;
    fitBg2->SetParameters(fitMgg2->GetParameter(3),fitMgg2->GetParameter(4),fitMgg2->GetParameter(5)) ;


    //=========================================================== get v2 real vs minv (200 bins)
    v2->Reset() ;
    v2->Sumw2();
    
    v2m->Reset() ;
    v2m->Sumw2();
    
    v2v4->Reset() ;
    v2v4->Sumw2();

    v2mv4->Reset() ;
    v2mv4->Sumw2();

  
    Int_t mMean = v2->FindBin(mean1) ;
    
    Int_t mMin2 = v2->FindBin(mean1-sigma1) ;
    Int_t mMax2 = v2->FindBin(mean1+sigma1) ;

    //========================================  get v2 Real vs minv (200 bins)
    for(Int_t i=1;i<=v2->GetNbinsX();i++){      
//      h3DR->GetXaxis()->SetRange(1+(i-1)*nRebin,(i)*nRebin) ;
//      TH1D * tmp =(TH1D*) h3DR->Project3D("z") ;
//      sprintf(name,"real_%d_%d",ii,i);
//      tmp->SetName(name);
      Int_t iptMin=h3DM->GetYaxis()->FindBin(ptmin+0.0001) ;
      Int_t iptMax=h3DM->GetYaxis()->FindBin(ptmax-0.0001) ;
      
      TH1D * tmp =(TH1D*) h3DR->ProjectionZ(Form("real_%d_%d",ii,i),1+(i-1)*nRebin,(i)*nRebin,iptMin,iptMax) ;
      Double_t v2mean=0,v2err=0 ;
      CalcV2(tmp,v2mean,v2err) ;
      v2->SetBinContent(i,v2mean);
      v2->SetBinError(i,v2err);      

fit->SetParameters(tmp->GetBinContent(1),0.1);
tmp->Fit(fit,"qr0");
v2->SetBinContent(i,fit->GetParameter(1));
v2->SetBinError(i,fit->GetParError(1));
tmp->Draw();

cout<<"v2{CalcV2}: "<<v2mean<<"+-"<<v2err<<", v2{Fit}: "<<fit->GetParameter(1)<<"+-"<<fit->GetParError(1)<<endl;

      delete tmp ;
}    


    //========================================  get v2 mixed vs minv (200 bins)
    for(Int_t i=1;i<=v2m->GetNbinsX();i++){      
      h3DM->GetXaxis()->SetRange(1+(i-1)*nRebin,(i)*nRebin) ;
      TH1D * tmp2 = (TH1D*)h3DM->Project3D("z") ;
      sprintf(name,"mixed_%d_%d",ii,i);
      tmp2->SetName(name);

/*      Double_t v2mean=0,v2err=0 ;
      CalcV2(tmp2,v2mean,v2err) ;
      v2m->SetBinContent(i,v2mean);
      v2m->SetBinError(i,v2err);      
*/
tmp2->Fit(fit,"qr0");
v2m->SetBinContent(i,fit->GetParameter(1));
v2m->SetBinError(i,fit->GetParError(1));

      delete tmp2 ;
    }    

    
    ///////////////////////////////////////////////with single harmonic v2
    //====================================== fit v2 real to find background

    //This is to find v2Bg with v2Real fit with polynomial functions with gap
/*
    reject=kTRUE;
    Double_t mMin_BgFit=0.05 ;
    Double_t mMax_BgFit=0.30 ;
    Double_t mMax_BgFit2=0.40 ;

    f1->SetParameters(0.1,0.,0.,0.);
    v2->Fit(f1,"qM","",mMin_BgFit,mMax_BgFit);

    f2->SetParameters(0.1,0.,0.,0.,0.);
    v2->Fit(f2,"qM+","",mMin_BgFit,mMax_BgFit2);
*/

    //Scale Mixed V2 to reproduce real V2
    TH1D * v2mScaled = (TH1D*)v2m->Clone(name) ;
    ScaleMixed(v2,v2m,v2mScaled) ;

    sprintf(name,"%_bg",v2->GetName()) ;
    TH1D * v2bg = (TH1D*)v2->Clone(name) ;
    
    fitV2p1->SetParameters(0.1,0.135,0.05,0.1,0.) ;
    fitV2p1->FixParameter(1,fitMgg1->GetParameter(1)) ;
    fitV2p1->FixParameter(2,fitMgg1->GetParameter(2)) ;
    v2->Fit(fitV2p1,"q","",0.07,0.25) ;
    fitV2p2->SetParameters(0.1,0.135,0.05,0.1,0.) ;
    fitV2p2->FixParameter(1,fitMgg2->GetParameter(1)) ;
    fitV2p2->FixParameter(2,fitMgg2->GetParameter(2)) ;
    v2->Fit(fitV2p2,"q+","",0.07,0.25) ;
    fitV2mix->FixParameter(1,fitMgg1->GetParameter(1)) ;
    fitV2mix->FixParameter(2,fitMgg1->GetParameter(2)) ;
    v2->Fit(fitV2mix,"q+","",0.07,0.25) ;
    
    
    
    //==========Draw===================
    if(ii<7)
      c3->cd(ii) ;
    else
      c4->cd(ii-6) ;


    spectrum_1 = new TPad("1", "1",0.001,0.32,0.99,0.99);
    spectrum_1->Draw();
    spectrum_1->cd();
    spectrum_1->Range(0,0,1,1);
    spectrum_1->SetFillColor(0);
    spectrum_1->SetFillStyle(1);
    spectrum_1->SetBorderSize(1);
    spectrum_1->SetBottomMargin(0.0);
    spectrum_1->SetTopMargin(0.03);
    spectrum_1->SetLeftMargin(0.1);
    spectrum_1->SetRightMargin(0.05);
    //  spectrum_1->SetLogx();
    //  spectrum_1->SetLogy();
    
    v2->GetXaxis()->SetRangeUser(0.05,0.3) ;
    v2->SetMarkerStyle(20) ;
    v2->SetMarkerSize(0.8) ;
    v2->SetYTitle("v_{2}^{raw}") ;
    v2->GetYaxis()->SetTitleSize(0.08) ;
    v2->GetYaxis()->SetTitleOffset(0.40) ;
//    v2->SetMinimum(0.) ;
//    v2->SetMaximum(0.3);
    v2->Draw(); 
    v2m->SetMarkerStyle(24) ;
    v2m->SetMarkerSize(0.8) ;
    v2m->SetMarkerColor(6) ;
    v2m->SetLineColor(6) ;
//    v2m->Draw("same");
    v2mScaled->SetLineColor(kOrange) ;
    v2mScaled->SetLineWidth(3) ;
    v2mScaled->SetFillColor(kOrange) ;
    v2mScaled->SetFillStyle(1001) ;
//    v2mScaled->Draw("sameL");
    
    
    if(ii==1){
     TLegend * lv2 = new TLegend(0.54,0.59,0.94,0.9) ;
     lv2->AddEntry(v2,"Data","p") ;
     lv2->AddEntry(fitV2mix,"Fit, mixed + G","l") ;
     lv2->AddEntry(fitV2p1,"Fit, pol2 + G","l") ;
     lv2->AddEntry(fitV2p2,"Fit, pol3 + G","l") ;
     lv2->Draw();            
    }
    
    TLine * linea = new TLine(mean-width,0,mean-width,1.) ;
    linea->Draw() ;
    TLine * linea2 = new TLine(mean+width,0,mean+width,1.) ;
    linea2->Draw() ;

    //Evalulate V2
    //
    //   v2^S = v2^Total - Bg/Signal (v2^Total - v2^Bg)
    //
    //Bg/signal
    Double_t bgS = 0 ;

  TH1D * hm = (TH1D*)gROOT->FindObjectAny("v2m") ;
  Double_t mix = hm->GetBinContent(hm->GetXaxis()->FindBin(fitMgg1->GetParameter(1))) ;

    if(fitMgg1->GetParameter(0)>0.) bgS = fitMgg1->Eval(fitMgg1->GetParameter(1))/fitMgg1->GetParameter(0)-1. ;
    Double_t v2M = fitV2mix->Eval(fitMgg1->GetParameter(1)) + bgS*fitV2mix->GetParameter(0) ;
    Double_t errA = fitV2mix->GetParError(0) + fitV2mix->GetParError(3)*mix ;
    Double_t sA = fitMgg1->GetParError(0)/fitMgg1->GetParameter(0) ;
    Double_t sB = fitMgg1->GetParError(3)/fitMgg1->GetParameter(3) ;
    Double_t errB = bgS*fitV2mix->GetParError(0) ; ;
    Double_t errC = bgS*TMath::Sqrt(sA*sA+sB*sB)*fitV2mix->GetParameter(0) ;
    Double_t v2Merr = TMath::Sqrt(errA*errA+errB*errB+errC*errC) ;
    
    v2Pi_mix->SetBinContent(ii,v2M) ;
    v2Pi_mix->SetBinError(ii,v2Merr) ;

    v2M = fitV2p1->Eval(fitMgg1->GetParameter(1)) + bgS*fitV2p1->GetParameter(0) ;
//    v2Merr = bgS*fitV2p1->GetParError(0) ;
    v2Pi_pol1->SetBinContent(ii,v2M) ;
    v2Pi_pol1->SetBinError(ii,v2Merr) ;
    
    if(fitMgg2->GetParameter(0)>0.) bgS = fitMgg2->Eval(fitMgg2->GetParameter(1))/fitMgg2->GetParameter(0)-1. ;
    Double_t v2M = fitV2p2->Eval(fitMgg1->GetParameter(1)) + bgS*fitV2p2->GetParameter(0) ;
//    Double_t v2Merr = bgS*fitV2p2->GetParError(0) ;
    v2Pi_pol2->SetBinContent(ii,v2M) ;
    v2Pi_pol2->SetBinError(ii,v2Merr) ;
        
    if(ii<7)
      c3->cd(ii) ;
    else
      c4->cd(ii-6) ;

    TPad *spectrum_2 = new TPad("2", "2",0.001,0.01,0.99,0.32);
    spectrum_2->SetFillColor(0) ;
    spectrum_2->SetFillStyle(0) ;
    spectrum_2->SetLogy(0) ;
    //    spectrum_2->SetGridy() ;
    spectrum_2->Draw();
    spectrum_2->Range(0,0,1,1);
    spectrum_2->SetFillColor(0);
    spectrum_2->SetBorderSize(1);
    spectrum_2->SetTopMargin(0.0);
    spectrum_2->SetBottomMargin(0.25);
    spectrum_2->SetLeftMargin(0.1);
    spectrum_2->SetRightMargin(0.05);
    //  spectrum_2->SetLogx();
    spectrum_2->cd() ;
    
        

    //================================================= S/Bg
    
    hRatio->SetMarkerStyle(20) ;
    hRatio->SetMarkerSize(0.8) ;
    hRatio->GetXaxis()->SetRangeUser(0.05,0.3) ;
//    hRatio->SetMaximum(hRatio->GetMaximum()*1.6) ;
//    hRatio->GetMinimum(0.0014) ;
//    hRatio->GetMinimum(0.002) ;
    hRatio->GetXaxis()->SetLabelSize(0.1) ;
    hRatio->SetLabelOffset(0.05) ;
    hRatio->SetLabelSize(0.1) ;
    hRatio->GetXaxis()->SetTitleOffset(0.8);
    hRatio->SetXTitle("m_{#gamma#gamma} (GeV/c^{2})");
    hRatio->GetXaxis()->SetTitleSize(0.14) ;
    hRatio->SetYTitle("Real/Mixed");
    hRatio->GetYaxis()->SetTitleSize(0.12) ;
    hRatio->GetYaxis()->SetTitleOffset(0.40) ;
    hRatio->GetYaxis()->SetLabelOffset(0.01) ;
    hRatio->GetYaxis()->SetNdivisions(503) ;
    hRatio->GetYaxis()->SetLabelSize(0.1) ;
    hRatio->Draw();
    TLine * line = new TLine(mean-width,0,mean-width,0.1) ;
    line->Draw() ;
    TLine * line2 = new TLine(mean+width,0,mean+width,0.1) ;
    line2->Draw() ;

   // Get Signal and Background if pi0 yield
    
    hRatio->Fit(fconst,"+0rq","",mean-width,mean+width) ; // fit Real2mixed ratio with const
    Double_t nRe_const=fconst->GetParameter(0) ;
    Double_t nReErr_const=fconst->GetParError(0) ;

    Double_t nBg_lin=fitBg1->Integral(mean-width,mean+width)/2./width ; // N in bg = fitbg1 integral under pi0 peak
    Double_t nBgErr_lin=fitMgg1->GetParError(3); //fitBg1->GetParError(0);

    Double_t nBg_quad=fitBg2->Integral(mean-width,mean+width)/2./width ; // N in bg = fitbg1 integral under pi0 peak
    Double_t nBgErr_quad=fitMgg2->GetParError(3); //fitBg2->GetParError(0);

    nSig->SetBinContent(ii,(nRe_const - nBg_lin)/nBg_lin) ;
    nSig->SetBinError(ii,TMath::Sqrt(nReErr_const*nReErr_const + nBgErr_lin*nBgErr_lin)) ;
    nSig2->SetBinContent(ii,(nRe_const - nBg_quad)/nBg_quad) ;
    nSig2->SetBinError(ii,TMath::Sqrt(nReErr_const*nReErr_const + nBgErr_quad*nBgErr_quad)) ;
 
    c3->Update() ;

  }


  TCanvas * csig = new TCanvas("signal","signal") ;

  TLegend *leg_n = new TLegend(0.1,0.9,0.5,0.7);
  leg_n->AddEntry(nSig,"N real with bg linear fit","p");
  leg_n->AddEntry(nSig2,"N real with bg quad fit","p");

  nSig->Draw("p");
  nSig->SetLineWidth(6);
  nSig2->SetMarkerStyle(21);
  nSig2->SetMarkerColor(2);
  nSig2->SetLineColor(2);
  nSig2->SetLineWidth(5);
  nSig2->Draw("psame");

  leg_n->Draw();

//============================================= Draw every v2

  v2Pi_mix->SetMarkerStyle(21) ;
  v2Pi_mix->SetLineColor(kOrange) ;
  v2Pi_mix->SetMarkerColor(kOrange) ;
  v2Pi_mix->Draw() ;
 
  v2Pi_pol1->SetMarkerStyle(25) ;
  v2Pi_pol1->SetLineColor(kOrange) ;
  v2Pi_pol1->SetMarkerColor(kOrange) ;
  v2Pi_pol1->Draw("same") ;

  v2Pi_pol2->SetMarkerStyle(23) ;
  v2Pi_pol2->SetLineColor(6) ;
  v2Pi_pol2->SetMarkerColor(6) ;
  v2Pi_pol2->Draw("same");

  
  TLegend *leg_v2 = new TLegend(0.1,0.9,0.5,0.7);
  leg_v2->AddEntry(v2Pi_mix,"Corrected Mixed","p");
  leg_v2->AddEntry(v2Pi_pol1,"Bg fitted with pol1","p");
  leg_v2->AddEntry(v2Pi_pol2,"Bg fitted with pol2","p");  
  leg_v2->Draw();

  //calculate sys errors
  TCanvas * csys = new TCanvas("V2sys","V2sys") ;
  for(Int_t i=1; i<=v2Pi_mix->GetNbinsX() ;i++){
//    Double_t mean=v2Pi_mix->GetBinContent(i);

  Double_t mean = v2Pi_pol1->GetBinContent(i)/v2Pi_pol1->GetBinError(i)/v2Pi_pol1->GetBinError(i) + v2Pi_pol2->GetBinContent(i)/v2Pi_pol2->GetBinError(i)/v2Pi_pol2->GetBinError(i);

    Double_t weight=
     1./v2Pi_pol1->GetBinError(i)/v2Pi_pol1->GetBinError(i)
     +1./v2Pi_pol2->GetBinError(i)/v2Pi_pol2->GetBinError(i) ;

    mean/=weight;

    Double_t rms = 
      (v2Pi_pol1->GetBinContent(i)-mean)*(v2Pi_pol1->GetBinContent(i)-mean)/v2Pi_pol1->GetBinError(i)/v2Pi_pol1->GetBinError(i)
     +(v2Pi_pol2->GetBinContent(i)-mean)*(v2Pi_pol2->GetBinContent(i)-mean)/v2Pi_pol2->GetBinError(i)/v2Pi_pol2->GetBinError(i) 
     +(v2Pi_mix->GetBinContent(i)-mean)*(v2Pi_mix->GetBinContent(i)-mean)/v2Pi_mix->GetBinError(i)/v2Pi_mix->GetBinError(i) ;

    Double_t weight=
     1./v2Pi_pol1->GetBinError(i)/v2Pi_pol1->GetBinError(i)
     +1./v2Pi_pol2->GetBinError(i)/v2Pi_pol2->GetBinError(i)
     +1./v2Pi_mix->GetBinError(i)/v2Pi_mix->GetBinError(i) ;

    
    rms/=weight ;

    v2Astat->SetBinContent(i,mean) ;
    v2Astat->SetBinError(i,v2Pi_mix->GetBinError(i)) ;
    v2Asys->SetBinContent(i,mean) ;
    v2Asys->SetBinError(i,TMath::Sqrt(rms)) ;
  }
  
  v2Asys->SetFillStyle(1001) ;
  v2Asys->SetFillColor(kYellow) ;

  v2Astat->SetMarkerColor(kOrange+7) ;
  v2Astat->SetMarkerStyle(20) ;
  v2Astat->SetLineColor(kOrange+7) ;

  switch(cen){
  case 0:
    sprintf(key,"Centrality 0-5%%") ; break ;
  case 1:
    sprintf(key,"Centrality 5-10%%") ; break ;
  case 2:
    sprintf(key,"Centrality 10-20%%") ; break ;
  case 3:
    sprintf(key,"Centrality 20-30%%") ; break ;
  case 4:
    sprintf(key,"Centrality 30-40%%") ; break ;
  case 5:
    sprintf(key,"Centrality 40-50%%") ; break ;
  case 10:
    sprintf(key,"Centrality 0-10%%") ; break ;
  case 11:
    sprintf(key,"Centrality 20-40%%") ; break ;
  }

  v2Asys->SetTitle(key) ;
  v2Asys->SetXTitle("p_{t} (GeV/c)") ;
  v2Asys->SetYTitle("v_{2}") ;
  v2Asys->SetMinimum(-0.05) ;
  v2Asys->SetMaximum(0.6) ;

  TH1D *effcorrect = v2Asys->Clone("effcorrect");
  TH1D *effcorrect2 = v2Asys->Clone("effcorrect2");

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

        eff=eff*res/resT;
        err=err*res/resT;
        cout<<"Pi0 efficiency correction: "<<eff<<"+-"<<err<<endl;
  }

  Double_t x[100],y[100],ex[100],ey[100] ;

  for(Int_t i=0;i<effcorrect->GetNbinsX();i++){
        effcorrect->SetBinContent(i,eff);
        effcorrect->SetBinError(i,err);
        effcorrect2->SetBinContent(i,eff);
        effcorrect2->SetBinError(i,0);
  }

  if(effcor){
        v2Astat->Add(effcorrect2);
        v2Asys->Add(effcorrect);
  }

  v2Asys->Scale(1/res) ;
  v2Astat->Scale(1/res) ;
  //Error due to reaction plane into sys


  v2Asys->Draw("E2") ;
  v2Astat->Draw("same") ;



//=================================================== save to file
  TFile fout("v2_method2_QA.root","update") ;
  v2Asys->Write(0,TObject::kOverwrite) ;
  v2Astat->Write(0,TObject::kOverwrite) ;
  

}
