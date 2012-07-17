#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TList.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TPad.h"
#endif

TGraphErrors *gr[20];
double mnX=1e6,mnY=1e6,mxX=-1e6,mxY=-1e6;
int nbTot = 0;
int nbc = 0;
TH1* hbase = 0;

void Compose(char* fL="resMult/res137366_eta_m0875_p0875_zv_m7_p7_V0_10bins_CutEta-0.9_0.9_Zv-7.0_7.0_bgInj_Shape_wdst_mcLB0_cutSig1.5_cutBg5.0.root",
	     char* fM="resMult/res137366_eta_m1875_m0625_zv_p6_p13_V0_10bins_CutEta-1.9_-0.6_Zv6.0_13.0_bgInj_Shape_wdst_mcLB0_cutSig1.5_cutBg5.0.root",
	     char* fR="resMult/res137366_eta_p0625_p1875_zv_m13_m6_V0_10bins_CutEta0.6_1.9_Zv-13.0_-6.0_bgInj_Shape_wdst_mcLB0_cutSig1.5_cutBg5.0.root")	     

{
  TObjArray *arr[3];
  TFile* fl = TFile::Open(fL);
  arr[0] = (TObjArray*)fl->Get("TObjArray");
  TFile* fm = TFile::Open(fM);
  arr[1] = (TObjArray*)fm->Get("TObjArray");
  TFile* fr = TFile::Open(fR);
  arr[2] = (TObjArray*)fr->Get("TObjArray");
  //
  for (int ia=0;ia<3;ia++) nbTot += ((TH1*)arr[ia]->At(0))->GetNbinsX();
  nbc = arr[0]->GetEntriesFast();
  //
  for (int ic=0;ic<nbc;ic++) {
    gr[ic] = new TGraphErrors(nbTot);
    gr[ic]->SetMarkerStyle(20+ic);
    gr[ic]->SetMarkerColor(kRed);    
    gr[ic]->SetLineColor(kRed);    
    int npg = 0;
    for (int ia=0;ia<3;ia++) {
      TH1* hh = (TH1*)arr[ia]->At(ic);
      int nbh = hh->GetNbinsX();
      for (int ib=1;ib<=nbh;ib++) {
	double vl  = hh->GetBinContent(ib);
	double vle = hh->GetBinError(ib);
	double x   = hh->GetBinCenter(ib);
	double xe  = hh->GetBinWidth(ib)/2;
	//
	if (mnX>x-xe) mnX = x-xe;
	if (mxX<x+xe) mxX = x+xe;	
	if (mnY>vl-vle) mnY = vl-vle;
	if (mxY<vl+vle) mxY = vl+vle;	
	//
	gr[ic]->SetPoint(npg, x,vl);
	gr[ic]->SetPointError(npg, xe,vle);
	npg++;
	//
      }
    }
  }
  //
  double del = mxX-mnX;
  mnX -= 0.05*del;
  mxX += 0.05*del;
  del = mxY-mnY;
  mnY -= 0.1*del;
  mxY += 0.1*del;
  hbase = new TH1F("hbase","",100,mnX,mxX);
  hbase->SetMinimum(mnY);
  hbase->SetMaximum(mxY);
  //
  hbase->Draw();
  for (int i=0;i<nbc;i++) {
    gr[i]->Draw("p");
  }
  gPad->SetGrid();
}
