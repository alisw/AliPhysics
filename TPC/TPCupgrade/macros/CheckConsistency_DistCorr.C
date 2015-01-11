#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeStream.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "AliTPCCorrection.h"
#include "AliTPCCorrectionLookupTable.h"
#include <AliToyMCEventGenerator.h>



void makeComparisonTree(TString filename, TString addToName)
{
  
  AliTPCCorrectionLookupTable *fTPCCorrection2=0x0;

  Bool_t doScaling=kTRUE;
  if (filename.Contains(":")) {
    TObjArray *arr=filename.Tokenize(":");
    TString s2(arr->At(1)->GetName());
    if (s2.Contains("-scale")) {
      doScaling=kFALSE;
      s2.ReplaceAll("-scale","");
    }
    TFile f2(s2);
    gROOT->cd();
    fTPCCorrection2=(AliTPCCorrectionLookupTable*)f2.Get("map");
    f2.Close();
    filename=arr->At(0)->GetName();
    delete arr;
  }
  TFile fn(filename.Data());
  gROOT->cd();
  AliTPCCorrectionLookupTable *fTPCCorrection=(AliTPCCorrectionLookupTable*)fn.Get("map");
  fn.Close();
  if (fTPCCorrection2 && doScaling) {
    Float_t dummy=0;
    fTPCCorrection2->SetCorrScaleFactor(AliToyMCEventGenerator::GetSCScalingFactor(fTPCCorrection, fTPCCorrection2,dummy));
    
  }
//   fTPCCorrection->BuildExactInverse();

//   TFile f("/tmp/corrTest.Root","recreate");
//   fTPCCorrection->Write("map");
//   f.Close();
  
  TString outFile=addToName;
  outFile.Append(".root");
  TTreeSRedirector *sred=new TTreeSRedirector(outFile.Data());
  
  Float_t dx[3]={0,0,0};
  
  for (Float_t iz=-245; iz<=245; iz+=10) {
    Short_t roc=(iz>=0)?0:18;
    for (Float_t ir=86; ir<250; ir+=10) {
      for (Float_t iphi=0; iphi<TMath::TwoPi(); iphi+=10*TMath::DegToRad()){
        Float_t x=ir*(Float_t)TMath::Cos(iphi);
        Float_t y=ir*(Float_t)TMath::Sin(iphi);
        Float_t x3[3]={x,y,iz};
        (*sred) << "t" <<
        "r="   << ir <<
        "phi=" << iphi <<
        "x=" << x <<
        "y=" << y <<
        "z=" << iz;
        
        //distortions
        fTPCCorrection->GetDistortion(x3,roc,dx);
        Float_t xd   = x+dx[0];
        Float_t yd   = y+dx[1];
        Float_t zd   = iz+dx[2];
        Float_t rd   = TMath::Sqrt(xd*xd+yd*yd);
        Float_t phid = TMath::ATan2(yd,xd);
        if (phid<0) phid+=TMath::TwoPi();
        (*sred) << "t" <<
        "xd="   << xd <<
        "yd="   << yd <<
        "zd="   << zd <<
        "rd="   << rd <<
        "phid=" << phid;
        
        // correct back distorted point
        Float_t xd3[3]={xd,yd,zd};

        if (fTPCCorrection2) {
          fTPCCorrection2->GetCorrection(xd3,roc,dx);
        } else {
          fTPCCorrection->GetCorrection(xd3,roc,dx);
        }
        Float_t xdc   = xd+dx[0];
        Float_t ydc   = yd+dx[1];
        Float_t zdc   = zd+dx[2];
        Float_t rdc   = TMath::Sqrt(xdc*xdc+ydc*ydc);
        Float_t phidc = TMath::ATan2(ydc,xdc);
        if (phidc<0) phidc+=TMath::TwoPi();
        
        (*sred)  << "t" <<
        "xdc="   << xdc <<
        "ydc="   << ydc <<
        "zdc="   << zdc <<
        "rdc="   << rdc <<
        "phidc=" << phidc;
        
        // write current point
        (*sred) << "t" <<
        "\n";
        
      }
    }
  }
  
  delete sred;
}

void makeAllComparisonTrees()
{
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps5_50kHz_precal.lookup.root","LUT_05");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps10_50kHz_precal.lookup.root","LUT_10");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2_eps20_50kHz_precal.lookup.root","LUT_20");
}

void makeAllComparisonTreesNew()
{
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2N2_eps5_50kHz_precal.lookup.root","LUT_05");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2N2_eps10_50kHz_precal.lookup.root","LUT_10");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2N2_eps20_50kHz_precal.lookup.root","LUT_20");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2N2_eps20_50kHz_precal.lookup.root","LUT_25");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2N2_eps20_50kHz_precal.lookup.root","LUT_30");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2N2_eps20_50kHz_precal.lookup.root","LUT_35");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/SC_NeCO2N2_eps20_50kHz_precal.lookup.root","LUT_40");
}

void makeAllComparisonTreesOld()
{
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/old/SC_NeCO2N2_eps5_50kHz_precal.lookup.root","LUT_05");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/old/SC_NeCO2N2_eps10_50kHz_precal.lookup.root","LUT_10");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/old/SC_NeCO2N2_eps20_50kHz_precal.lookup.root","LUT_20");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/old/SC_NeCO2N2_eps20_50kHz_precal.lookup.root","LUT_25");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/old/SC_NeCO2N2_eps20_50kHz_precal.lookup.root","LUT_30");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/old/SC_NeCO2N2_eps20_50kHz_precal.lookup.root","LUT_35");
  makeComparisonTree("$ALICE_ROOT/TPC/Calib/maps/old/SC_NeCO2N2_eps20_50kHz_precal.lookup.root","LUT_40");
}

TCanvas *GetCanvas(TString addToName);

void makeHistos(TString addToName) {
  TString filename; //("test_");
  filename.Append(addToName.Data());
  filename.Append(".root");
  TFile f(filename.Data());
  gROOT->cd();
  TTree *t=(TTree*)f.Get("t");
  gStyle->SetTitleX(0.18);
  gStyle->SetTitleW(1-.18-.1);

  t->SetMarkerStyle(20);
  t->SetMarkerSize(.8);

  TCanvas *c=0x0;
  c=GetCanvas(addToName+"_zRes");
  t->Draw("zdc-z:z:r","","colz");
  c->SaveAs(Form("%s_zRes.png",addToName.Data()));
  //
  c=GetCanvas(addToName+"_rRes");
  t->Draw("rdc-r:z:r","","colz");
  c->SaveAs(Form("%s_rRes.png",addToName.Data()));
  //
  c=GetCanvas(addToName+"_phiRes");
  t->SetAlias("phiFix","-((phidc-phi)>4)*2*TMath::Pi()+((phidc-phi)<-4)*2*TMath::Pi()");
  t->Draw("phidc-phi+phiFix:z:r","","colz");
  c->SaveAs(Form("%s_phiRes.png",addToName.Data()));
  //
  c=GetCanvas(addToName+"_rphiRes");
  t->Draw("(phidc*rdc)-(phi*r):z+(r-84)/(254-84)*18:r","abs(phidc-phi)<1","colz");
  c->SaveAs(Form("%s_rphiRes.png",addToName.Data()));

  TCanvas *c2=0x0;
  c2=GetCanvas(addToName+"_Res_1D");
  c2->Divide(2,2);
  
  c2->cd(1);
  t->Draw("zdc-z","","");
  //
  c2->cd(2);
  t->Draw("rdc-r","","");
  //
  c2->cd(3);
  t->SetAlias("phiFix","-((phidc-phi)>4)*2*TMath::Pi()+((phidc-phi)<-4)*2*TMath::Pi()");
  t->SetAlias("phiRes","phidc-phi+phiFix");
  t->Draw("phiRes","","");
  //
  c2->cd(4);
  t->Draw("(phidc*rdc)-(phi*r)","abs(phidc-phi)<1","");

  c2->SaveAs(Form("%s_Res_1D.png",addToName.Data()));
  
  f.Close();
}

void makeHistosDist(TString addToName) {
  TString filename; //("test_");
  filename.Append(addToName.Data());
  filename.Append(".root");
  TFile f(filename.Data());
  gROOT->cd();
  TTree *t=(TTree*)f.Get("t");
  gStyle->SetTitleX(0.18);
  gStyle->SetTitleW(1-.18-.1);
  
  t->SetMarkerStyle(20);
  t->SetMarkerSize(.8);
  
  TCanvas *c=0x0;
  c=GetCanvas(addToName+"_zResDist");
  t->Draw("zd-z:z:r","","colz");
  c->SaveAs(Form("%s_zResDist.png",addToName.Data()));
  //
  c=GetCanvas(addToName+"_rResDist");
  t->Draw("rd-r:z:r","","colz");
  c->SaveAs(Form("%s_rResDist.png",addToName.Data()));
  //
  c=GetCanvas(addToName+"_phiResDist");
  t->Draw("phid-phi:z:r","abs(phid-phi)<1","colz");
  c->SaveAs(Form("%s_phiResDist.png",addToName.Data()));
  //
  c=GetCanvas(addToName+"_rphiResDist");
  t->Draw("(phid*rd)-(phi*r):z+(r-84)/(254-84)*18:r","abs(phid-phi)<1","colz");
  c->SaveAs(Form("%s_rphiResDist.png",addToName.Data()));
  
  f.Close();
}


void makeAllHistos() {
  makeHistos("LUT_05");
  makeHistos("LUT_10");
  makeHistos("LUT_20");
  makeHistos("LUT_25");
  makeHistos("LUT_30");
  makeHistos("LUT_35");
  makeHistos("LUT_40");
  
}

void makeAllHistosDist() {
  makeHistosDist("LUT_05");
  makeHistosDist("LUT_10");
  makeHistosDist("LUT_20");
  makeHistosDist("LUT_25");
  makeHistosDist("LUT_30");
  makeHistosDist("LUT_35");
  makeHistosDist("LUT_40");
  
}

TCanvas *GetCanvas(TString addToName)
{
  TString cName(addToName);
  cName.Prepend("c_");
  TCanvas *c=(TCanvas*)gROOT->GetListOfCanvases()->FindObject(cName.Data());
  if (!c) c=new TCanvas(cName.Data(),addToName.Data());
  c->Clear();
  c->cd();
  return c;
}



