#include "MethodsP.h"

const Double_t kMean=0.135 ; //Approximate peak position to facilitate error estimate
Int_t effcor=0; //correct on efficiency

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
  return par[0]*exp(-dx*dx/2.)+par[3] + par[4]*(x[0]-kMean) ;
}
Double_t CB2(Double_t * x, Double_t * par){
  //Another parameterization of Real/Mixed ratio
  Double_t m=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-m)/s ;
  return par[0]*exp(-dx*dx/2.)+par[3]+par[4]*(x[0]-kMean) + par[5]*(x[0]-kMean)*(x[0]-kMean);
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
  return par[0] ;
}
Double_t BG2(Double_t * x, Double_t * par){
  //Another normalization of  Mixed
  return par[0]+par[1]*(x[0]-kMean) ;
}

void Photon(Int_t cen=0,Int_t side=0){
//side = 1 - VOC, side = 0 - V0A
//side = 2 - TPC 3 sub method, 3 - TPC 2 sub method

  gStyle->SetOptFit(1);

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

  TH2F *h2DR, *h2DM;
  h2DR = GetRealMixed(cen,w,d,kind,kind2,PID,1);
  h2DM = GetRealMixed(cen,w,d,kind,kind2,PID,0);

  TH1F * hCen = GetCen() ;
  TH2F * hCenPHOS = GetCenPHOS() ;

  //  sprintf(kind,"") ; //"" for all
  const Int_t nphi=5;
  Double_t xphi[6];
  for(Int_t i=0;i<=nphi;i++)xphi[i]=TMath::PiOver2()*i/nphi ;

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

  cout<<"res = "<< res <<endl;

//now to the v2 calculation
  
  TH1D * v2A1s= new TH1D("v2 with 2nd harmonic fit","v_{2} 1",nbin,xa) ;
  TH1D * v2A2s= new TH1D("v2 with 4rd harmonic fit","v_{2} 2",nbin,xa) ;
  TH1D * v2Asys= new TH1D("v2 sys","v_{2} sys",nbin,xa) ;
  TH1D * v2Astat= new TH1D("v2 stat","v_{2} stat",nbin,xa) ;

  char key2[155];

  Int_t col[4]={kRed,kBlue,kGreen+3,kMagenta} ;
  Int_t sym[4]={20,21,24,25} ;

  //  TF1 * fit = new TF1("fit","[0]*(1+2.*[1]*cos(2*(x-[2])))",0.,5.) ;
  TF1 * v2fit4 = new TF1("v2fit","[0]*(1+2.*[1]*cos(2.*x)+2.*[2]*cos(4.*x) )",0.,5.) ;
  v2fit->SetLineColor(2) ;
  TF1 * v2fit2 = new TF1("v2fit","[0]*(1+2.*[1]*cos(2.*x) )",0.,5.) ;
  TCanvas * c = new TCanvas("fit","fit") ;
  c->Divide(4,3);

  
  for(Int_t i=1;i<=nbin;i++){
//    h3DR->GetYaxis()->SetRangeUser(xa[i-1],xa[i]) ;
//    Double_t pt=(xa[i]+xa[i-1])/2. ;    
    sprintf(key2,"dNdPhi_%d",i);    
    Int_t iMin=h2DR->GetXaxis()->FindBin(xa[i-1]+0.0001) ;
    Int_t iMax=h2DR->GetXaxis()->FindBin(xa[i]-0.0001) ;

    TH1D * hdNdPhi = (TH1D*)h2DR->ProjectionY(key2,iMin,iMax);
    hdNdPhi->Sumw2();
    c->cd(i);

    hdNdPhi->SetXTitle("#Delta#phi (rad)") ;
    hdNdPhi->SetYTitle("dN/d#phi (a.u.)") ;
    hdNdPhi->SetTitle(Form("pT: %f-%f",xa[i-1],xa[i]));
    hdNdPhi->SetMarkerStyle(sym[0]) ;
    hdNdPhi->SetMarkerColor(col[0]) ;
    
    v2fit2->SetLineStyle(1) ;
    v2fit2->SetLineColor(col[0]) ;
    v2fit2->SetParameters(hdNdPhi->GetBinContent(5),0.1) ;
    hdNdPhi->Fit(v2fit2,"Qi") ;

    v2A1s->SetBinContent(i,v2fit2->GetParameter(1)) ;
    v2A1s->SetBinError(i,v2fit2->GetParError(1)) ;

    sprintf(key2,"dNdPhi2_%d",i);
    TH1D * hdNdPhi2 = hdNdPhi->Clone(key2);

    v2fit4->SetLineStyle(2) ;
    v2fit4->SetLineColor(col[1]) ;
    v2fit4->SetParameters(hdNdPhi2->GetBinContent(5),0.1) ;
    hdNdPhi2->Fit(v2fit4,"Qi") ;

    v2A2s->SetBinContent(i,v2fit4->GetParameter(1)) ;
    v2A2s->SetBinError(i,v2fit4->GetParError(1)) ;

   
    hdNdPhi->Draw() ;
    hdNdPhi2->Draw("same") ;
  }

  TH1D *effcorrect = v2A1s->Clone("effcorrect");
  TH1D *effcorrect2 = v2A1s->Clone("effcorrect2");

  Double_t eff=0, err=0;
  if(strcmp(PID,"Both2core")==0){
        if(cen==0) {eff=.05; err=.04;}
        if(cen==1) {eff=.02; err=.04;}
        if(cen==2) {eff=.02; err=.03;}
        if(cen==3) {eff=.01; err=.02;}
        if(cen==4) {eff=.004; err=.02;}
        if(cen==5) {eff=.004; err=.02;}
        if(cen==10) {eff=.03; err=.02;}
        if(cen==11) {eff=.01; err=.02;}
        if(cen==20) {eff=.02; err=.02;}
        if(cen==21) {eff=.03; err=.02;}

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

//  v2A1s->Scale(res) ;
//  v2A2s->Scale(res) ;
  
  v2A1s->SetMarkerStyle(sym[0]) ;
  v2A1s->SetMarkerColor(col[0]) ;
  v2A2s->SetMarkerStyle(sym[1]) ;
  v2A2s->SetMarkerColor(col[1]) ;

  TCanvas * c2 = new TCanvas("v2","v2") ;

  char bname[255];

  char cname[255];
  if(cen==0)sprintf(cname,"0-5%%");
  if(cen==1)sprintf(cname,"5-10%%");
  if(cen==2)sprintf(cname,"10-20%%");
  if(cen==3)sprintf(cname,"20-30%%");
  if(cen==4)sprintf(cname,"30-40%%");
  if(cen==5)sprintf(cname,"40-50%%");
  if(cen==10)sprintf(cname,"0-10%%");
  if(cen==11)sprintf(cname,"20-40%%");

  sprintf(bname,"v^{#gamma}_{2}{EP} for centrality %s, EP %s", cname, d);

  TH1D * box = new TH1D("box",bname,100,0.,20.) ;
  box->SetMinimum(-0.1) ;
  box->SetMaximum(0.4) ;
  box->SetXTitle("p_{t}^{#gamma} (GeV/c)") ;
  box->SetYTitle("v_{2}") ;
  box->Draw() ;
  v2A1s->Draw("same") ;
  v2A2s->Draw("same") ;
    
  TLegend * l = new TLegend(0.6,0.7,0.9,0.9) ;
  l->AddEntry(v2A1s,"Fit v2","p") ;
  l->AddEntry(v2A2s,"Fit v2+v4","p") ;

  l->Draw() ;

  for(Int_t i=1; i<=v2A1s->GetNbinsX() ;i++){
    Double_t mean=
      v2A1s->GetBinContent(i)/v2A1s->GetBinError(i)/v2A1s->GetBinError(i)
      +v2A2s->GetBinContent(i)/v2A2s->GetBinError(i)/v2A2s->GetBinError(i);
    Double_t weight=
       1./v2A1s->GetBinError(i)/v2A1s->GetBinError(i)
      +1./v2A2s->GetBinError(i)/v2A2s->GetBinError(i);

     Double_t statErr = TMath::Min(v2A1s->GetBinError(i),v2A2s->GetBinError(i));

    mean/=weight ;
    Double_t rms = (v2A1s->GetBinContent(i)-mean)*(v2A1s->GetBinContent(i)-mean)
/v2A1s->GetBinError(i)/v2A1s->GetBinError(i) + (v2A2s->GetBinContent(i)-mean)*(v2A2s->GetBinContent(i)-mean)/v2A2s->GetBinError(i)/v2A2s->GetBinError(i);
    rms/=weight;   

    v2Astat->SetBinContent(i,mean) ;
    v2Astat->SetBinError(i,statErr) ;
    v2Asys->SetBinContent(i,mean) ;
    v2Asys->SetBinError(i,TMath::Sqrt(rms)) ;

  }

  if(effcor){
        v2Astat->Add(effcorrect2);
        v2Asys->Add(effcorrect);
  }

  v2Asys->Scale(1/res) ;
  v2Astat->Scale(1/res) ;

  TCanvas * c5 = new TCanvas("v2_res","v2 res") ;

  v2Asys->SetFillStyle(1001) ;
  v2Asys->SetFillColor(kYellow) ;
 
  v2Astat->SetMarkerColor(kOrange+7) ;
  v2Astat->SetMarkerStyle(20) ;
  v2Astat->SetLineColor(kOrange+7) ;

  v2Asys->SetXTitle("p_{t} (GeV/c)") ;
  v2Asys->SetYTitle("v_{2}") ;

  box->Draw() ;
  v2Asys->Draw("sameE1") ;
  //  v2Fsys->Draw("E2same") ;
  v2Astat->Draw("same") ;
  //  v2Fstat->Draw("same") ;

  TFile fout("v2_photons.root","update");
  char nname[100];

  sprintf(nname,"v2sys_m1_%s_%s_%d_%d",w,PID,cen,side) ;
  v2Asys->SetName(nname) ;
  sprintf(nname,"v2stat_m1_%s_%s_%d_%d",w,PID,cen,side) ;
  v2Astat->SetName(nname) ;


/*  sprintf(nname,"v2sys_%s_%s_%d",kind,PID,cen) ;
  v2Asys->SetName(nname) ;
  sprintf(nname,"v2stat_%s_%s_%d",kind,PID,cen) ;
  v2Astat->SetName(nname) ;
*/

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

