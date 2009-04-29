/*
  //Analysis of the output of AliTPCtaskPID.
  //3 6D histograms  - THnSparse created in the task:
  //TPC raw dEdx
  //TPC normalized dEdx (dEdx_rec/dNdx_mc)
  //TPC PID probabilities
  //
  //The values are binned in following variables:
  // Some of them are correlated - but THnSpase handle it  
  //                               ~ 14 MBy per object needed
  //


  .x ~/NimStyle.C
  .x ~/UliStyle.C
  .L $ALICE_ROOT/PWG1/Macros/tpcQA.C+
  Init();
  
  // 0 - chi2
  // 1 - number of clusters
  // 2 - number of findable clusters
  // 3 - number of clusters/ findable clusters  
  // 4 - pt          - at the entrance of the TPC
  // 5 - eta         - at the entrance of the TPC
  // 6 - phi         - at the entrance of the TPC
  GetProjection(fTPCqa,0,4,  0,10, -1,1, -3.14,3.14)  ->ProfileX()->Draw();
  //
  MakeReport();
  
 
*/
#include "TFile.h"
#include "THnSparse.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"

#include "TLegend.h"
#include "TCanvas.h"

#include "AliPID.h"

Int_t kmicolors[10]={1,2,3,4,6,7,8,9,10,11};
Int_t kmimarkers[10]={21,22,23,24,25,26,27,28,29,30};



THnSparse * fTPCqa     = 0;
TString selName;
Int_t version=0;

void Init(){
  TFile *f = new TFile("OutputQA.root");
  TObjArray *array= (TObjArray*)f->Get("tpcTaskQA");
  delete f;
  new TCanvas("dEdx study");
  fTPCqa = (THnSparse*)array->At(0);
}


TH2F* GetProjection(THnSparse*his, Int_t i0, Int_t i1,  Float_t p0, Float_t p1, Float_t eta0, Float_t eta1, Float_t phi0, Float_t phi1){

  his->GetAxis(4)->SetRangeUser(p0,p1);
  his->GetAxis(5)->SetRangeUser(eta0,eta1);
  his->GetAxis(6)->SetRangeUser(phi0,phi1);
  TH2F * res = (TH2F*) his->Projection(i0,i1);
  res->SetXTitle(his->GetAxis(i1)->GetTitle());
  res->SetYTitle(his->GetAxis(i0)->GetTitle());
  return res;
}

TH1F* GetProjection(THnSparse*his, Int_t i0,  Float_t p0, Float_t p1, Float_t eta0, Float_t eta1, Float_t phi0, Float_t phi1){
  
  his->GetAxis(4)->SetRangeUser(p0,p1);
  his->GetAxis(5)->SetRangeUser(eta0,eta1);
  his->GetAxis(6)->SetRangeUser(phi0,phi1);
  TH1F * res = (TH1F*) his->Projection(i0);
  res->SetXTitle(his->GetAxis(i0)->GetTitle());
  return res;
}


void DrawChi2(){
  //
  //
  //
  TCanvas *canvas= new TCanvas("Chi2","Chi2");
  canvas->Divide(2,2);
  //
  canvas->cd(1);
  TH1 *hischi2  =  GetProjection(fTPCqa,0,  0,10, -0.9,0.9, -3.14,3.14); 
  hischi2->SetXTitle("#chi^{2}/N_{cl}");
  hischi2->Draw();
  canvas->cd(2)->SetLogx(kTRUE);
  TH1 *hischi2Pt  =  GetProjection(fTPCqa,0,4,  0,10, -0.9,0.9, -3.14,3.14)->ProfileX(); 
  hischi2Pt->SetXTitle("p_{t} (GeV/c)");
  hischi2Pt->SetYTitle("#chi^{2}/N_{cl}");
  hischi2Pt->SetMinimum(0);
  hischi2Pt->SetMaximum(4);
  hischi2Pt->Draw();
  //
  canvas->cd(3)->SetLogx(kFALSE);
  TH1 *hischi2Eta  =  GetProjection(fTPCqa,0,5,  0,10, -0.9,0.9, -3.14,3.14)->ProfileX(); 
  hischi2Eta->SetXTitle("#eta");
  hischi2Eta->SetYTitle("#chi^{2}/N_{cl}"); 
  hischi2Eta->SetMinimum(0);
  hischi2Eta->SetMaximum(4);
  hischi2Eta->Draw();
  
  canvas->cd(4)->SetLogx(kFALSE);
  TH1 *hischi2Phi  =  GetProjection(fTPCqa,0,6,  0,10, -0.9,0.9, -3.14,3.14)->ProfileX(); 
  hischi2Phi->SetXTitle("#phi");
  hischi2Phi->SetYTitle("#chi^{2}/N_{cl}");  
  hischi2Phi->SetMinimum(0);
  hischi2Phi->SetMaximum(4);
  hischi2Phi->Draw();
  
  canvas->SaveAs("pic/chi2.eps");
  canvas->SaveAs("pic/chi2.gif");
}



void DrawNclRatio(){
  //
  //
  //
  TCanvas *canvas= new TCanvas("Ncl_Nclf","Ncl_Nclf");
  canvas->Divide(2,2);
  //
  canvas->cd(1);
  TH1 *hischi2  =  GetProjection(fTPCqa,3,  0,10, -0.9,0.9, -3.14,3.14); 
  hischi2->SetXTitle("N_{cl}/N_{f}");
  hischi2->Draw();
  canvas->cd(2)->SetLogx(kTRUE);
  TH1 *hischi2Pt  =  GetProjection(fTPCqa,3,4,  0,10, -0.9,0.9, -3.14,3.14)->ProfileX(); 
  hischi2Pt->SetXTitle("p_{t} (GeV/c)");
  hischi2Pt->SetYTitle("N_{cl}/N_{f}");
  hischi2Pt->SetMinimum(0.6);
  hischi2Pt->SetMaximum(1.1);
  hischi2Pt->Draw();
  //
  canvas->cd(3)->SetLogx(kFALSE);
  TH1 *hischi2Eta  =  GetProjection(fTPCqa,3,5,  0,10, -0.9,0.9, -3.14,3.14)->ProfileX(); 
  hischi2Eta->SetXTitle("#eta");
  hischi2Eta->SetYTitle("N_{cl}/N_{f}"); 
  hischi2Eta->SetMinimum(0.6);
  hischi2Eta->SetMaximum(1.1);
  hischi2Eta->Draw();
  
  canvas->cd(4)->SetLogx(kFALSE);
  TH1 *hischi2Phi  =  GetProjection(fTPCqa,3,6,  0,10, -0.9,0.9, -3.14,3.14)->ProfileX(); 
  hischi2Phi->SetXTitle("#phi");
  hischi2Phi->SetYTitle("N_{cl}/N_{f}");  
  hischi2Phi->SetMinimum(0.6);
  hischi2Phi->SetMaximum(1.1);
  hischi2Phi->Draw();
  
  canvas->SaveAs("pic/nclratio.eps");
  canvas->SaveAs("pic/nclratio.gif");
}


void DrawNcl(){
  //
  //
  //
  TCanvas *canvas= new TCanvas("Ncl","Ncl");
  canvas->Divide(2,2);
  //
  canvas->cd(1);
  TH1 *hischi2  =  GetProjection(fTPCqa,1,  0.5,10, -0.9,0.9, -3.14,3.14); 
  hischi2->SetXTitle("N_{cl}");
  hischi2->Draw();
  canvas->cd(2)->SetLogx(kTRUE);
  TH1 *hischi2Pt  =  GetProjection(fTPCqa,1,4,  0.5,10, -0.9,0.9, -3.14,3.14)->ProfileX(); 
  hischi2Pt->SetXTitle("p_{t} (GeV/c)");
  hischi2Pt->SetYTitle("N_{cl}");
  hischi2Pt->SetMinimum(50);
  hischi2Pt->SetMaximum(160);
  hischi2Pt->Draw();
  //
  canvas->cd(3)->SetLogx(kFALSE);
  TH1 *hischi2Eta  =  GetProjection(fTPCqa,1,5,  0.5,10, -0.9,0.9, -3.14,3.14)->ProfileX(); 
  hischi2Eta->SetXTitle("#eta");
  hischi2Eta->SetYTitle("N_{cl}"); 
  hischi2Eta->SetMinimum(50);
  hischi2Eta->SetMaximum(160);
  hischi2Eta->Draw();
  
  canvas->cd(4)->SetLogx(kFALSE);
  TH1 *hischi2Phi  =  GetProjection(fTPCqa,1,6,  0.5,10, -0.9,0.9, -3.14,3.14)->ProfileX(); 
  hischi2Phi->SetXTitle("#phi");
  hischi2Phi->SetYTitle("N_{cl}");  
  hischi2Phi->SetMinimum(50);
  hischi2Phi->SetMaximum(160);
  hischi2Phi->Draw();

  canvas->SaveAs("pic/ncl.eps");
  canvas->SaveAs("pic/ncl.gif");
}

void MakeReport(){
  DrawNcl();
  DrawNclRatio();
  DrawChi2();
}
