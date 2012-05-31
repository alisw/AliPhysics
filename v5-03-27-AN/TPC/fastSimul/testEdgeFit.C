/*
  .L $ALICE_ROOT/TPC/fastSimul/testEdgeFit.C+
  testEdge(50000);
*/

#include "TVectorD.h"
#include "TMath.h"
#include "TLinearFitter.h"
#include "TRandom.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TPad.h"

#include "TTreeStream.h"
//#include "TVectorD.h"


void testEdge(Int_t ntest){

  TTreeSRedirector *cstream = new TTreeSRedirector("testEdge.root");
  Double_t yth=1;
  Double_t ymc0=0;
  Double_t kx=0;
  Double_t sy=0.25;
  Int_t npoints=100;
  TVectorD vecX;
  TVectorD vecY0;
  TVectorD vecYC;
  TVectorD vecYT;

  TVectorD vecU;
  TVectorD vecYTa;
  TVectorD vecYTb;
  TVectorD vecYEa;
  TVectorD vecYEb;
  //
  TVectorD vecYTM5;
  TVectorD vecYTk5;
  //
  TVectorD vecE;
  //
  //
  for (Int_t it=0; it<ntest;it++){
    //
    ymc0 = gRandom->Rndm()*1;
    kx = gRandom->Rndm()*1.;
    npoints = TMath::Nint(gRandom->Rndm()*50);
    vecX.ResizeTo(npoints);
    vecY0.ResizeTo(npoints);
    vecYC.ResizeTo(npoints);
    vecYT.ResizeTo(npoints);
    //
    vecYTa.ResizeTo(npoints);
    vecYTb.ResizeTo(npoints);
    vecYEa.ResizeTo(npoints);
    vecYEb.ResizeTo(npoints);
    vecU.ResizeTo(npoints);
    //
    TLinearFitter fitter(2,"pol1");
    for (Int_t ipoint=npoints-1;ipoint>0;ipoint--){
      Double_t x[10];
      //
      x[0]=(Double_t)(ipoint)*0.2;
      Double_t y0   = ymc0+Double_t(x[0])*kx; 
      Double_t yc   = y0+gRandom->Gaus()*sy;
      //
      if (yc<yth)   yc=yth;
      vecX[ipoint]   = x[0];;
      vecY0[ipoint]  = y0;
      vecYC[ipoint]  = yc;
      vecU[ipoint]   = 0;
      vecYT[ipoint]  = -10;
      vecYTa[ipoint] = -10;
      vecYTb[ipoint] = -10;
      fitter.AddPoint(x,yc,sy);
      if (fitter.GetNpoints()>4){
	fitter.Eval();
	fitter.GetErrors(vecE);
	vecYTa[ipoint]=fitter.GetParameter(0);
	vecYTb[ipoint]=fitter.GetParameter(1);
	vecYEa[ipoint]=vecE[0];
	vecYEb[ipoint]=vecE[1];
	vecU[ipoint]=fitter.GetNpoints();
	vecYT[ipoint]=fitter.GetParameter(0)+fitter.GetParameter(1)*x[0];
      }
    }
    (*cstream)<<"Dump"<<
      "n="<<npoints<<
      "ymc="<<ymc0<<
      "kx="<<kx<<
      "vX.="<<&vecX<<
      "vU.="<<&vecU<<
      "vY0.="<<&vecY0<<
      "vYC.="<<&vecYC<<
      "vYT.="<<&vecYT<<
      "vYTax.="<<&vecYTa<<
      "vYTb.="<<&vecYTb<<
      "vYEa.="<<&vecYEa<<
      "vYEb.="<<&vecYEb<<
      "\n";
  }
  delete cstream;
}


void MakePicDy(){
  TFile f("testEdge.root");
  TTree * tree = (TTree*)f.Get("Dump");
  TObjArray arrfdy0,arrfdyT;
  TH2F *phisdY0 = new TH2F("hisdY0","hisdY0",25,1,3,200,-3,3);
  TH2F *phisdYT = new TH2F("hisdYT","hisdYT",25,1,3,200,-3,3);
  tree->Draw("10*(vYT.fElements-vY0.fElements):vY0.fElements>>hisdY0","vU.fElements>20&&kx>0.1","colz");
  phisdY0->SetDirectory(0);
  tree->Draw("10*(vYT.fElements-vY0.fElements):vYT.fElements>>hisdYT","vU.fElements>20&&kx>0.1","colz");
  phisdYT->SetDirectory(0);

  phisdY0->FitSlicesY(0,0,-1,0,"QNR",&arrfdy0); 
  phisdYT->FitSlicesY(0,0,-1,0,"QNR",&arrfdyT);
  //
  TH1 *hisdy0 = (TH1*)arrfdy0.At(1);
  TH1 *hisdyT = (TH1*)arrfdyT.At(1);
  hisdy0->SetLineColor(1);
  hisdyT->SetLineColor(2);
  hisdy0->SetMarkerColor(1);
  hisdyT->SetMarkerColor(2);
  hisdy0->SetMarkerStyle(22);
  hisdyT->SetMarkerStyle(24);
  hisdyT->SetXTitle("Y (cm)");
  hisdyT->SetYTitle("Y_{t}-Y_{mc} (mm)");
  hisdyT->SetMinimum(-hisdyT->GetMaximum());
  hisdyT->Draw();
  hisdy0->Draw("same");
  TLegend *legend = new TLegend(0.25,0.60,0.85,0.85, "Track residuals close to edge (Edge= 1 cm, #sigma_{cl}=0.25 cm)");
  legend->AddEntry(hisdy0,"Y=Y_{mc}");
  legend->AddEntry(hisdyT,"Y=Y_{t}");
  legend->Draw();
  gPad->SaveAs("picEdge/dYedgeMC.eps");
  gPad->SaveAs("picEdge/dYedgeMC.gif");
}

void MakePickDy(){
  TFile f("testEdge.root");
  TTree * tree = (TTree*)f.Get("Dump");
  TObjArray arrfdky0,arrfdkyT;
  TH2F *phisdKY0 = new TH2F("hisdKY0","hisdKY0",25,1,3,200,-150,150);
  TH2F *phisdKYT = new TH2F("hisdKYT","hisdKYT",25,1,3,200,-150,150);
  tree->Draw("1000*(vYTb.fElements-kx):vY0.fElements>>hisdKY0","vU.fElements>20","colz");
  phisdKY0->SetDirectory(0);
  tree->Draw("1000*(vYTb.fElements-kx):vYT.fElements>>hisdKYT","vU.fElements>20","colz");
  phisdKYT->SetDirectory(0);

  phisdKY0->FitSlicesY(0,0,-1,0,"QNR",&arrfdky0); 
  phisdKYT->FitSlicesY(0,0,-1,0,"QNR",&arrfdkyT);
  //
  TH1 *hisdky0 = (TH1*)arrfdky0.At(1);
  TH1 *hisdkyT = (TH1*)arrfdkyT.At(1);
  hisdky0->SetLineColor(1);
  hisdkyT->SetLineColor(2);
  hisdky0->SetMarkerColor(1);
  hisdkyT->SetMarkerColor(2);
  hisdky0->SetMarkerStyle(22);
  hisdkyT->SetMarkerStyle(24);
  hisdkyT->SetMinimum(-20);
  hisdkyT->SetMaximum(20);
  hisdkyT->SetXTitle("Y (cm)");
  hisdkyT->SetYTitle("k_{yt}-k_{ymc} (mrad)");
  hisdkyT->Draw();
  hisdky0->Draw("same");
  TLegend *legend = new TLegend(0.25,0.60,0.85,0.85, "Track residuals close to edge (Edge= 1 cm, #sigma_{cl}=0.25 cm)");
  legend->AddEntry(hisdky0,"Y=Y_{mc}");
  legend->AddEntry(hisdkyT,"Y=Y_{t}");
  legend->Draw();
  gPad->SaveAs("picEdge/dKYedgeMC.eps");
  gPad->SaveAs("picEdge/dKYedgeMC.gif");
}


