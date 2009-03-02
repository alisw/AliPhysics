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
  // 0 - particle specie as defined in the AliPID - negatives+5 <0,9>
  // 1 - momenta - at the entrance of the TPC
  // 2 - tan lambda- fP[3]
  // 3 - betagamma
  // 4 - npoints
  // 5 - measurement - dEdx, dEdx/BB resp.  PID probability
  // 6 - BB



  .x ~/NimStyle.C
  .x ~/UliStyle.C
  .L $ALICE_ROOT/PWG1/Macros/tpcPIDResol.C+
  Init();
  //
  // custom draw  GetProjection(THnSparse*his, Int_t i0, Int_t i1, Float_t type0, Float_t type1, Float_t p0, Float_t p1, Float_t eta0, Float_t eta1, Int_t ncl0, Int_t ncl1)
  GetProjection(fTPCsignalNorm,5,2,  -1,11,   5,10,  -1.4,1.4,   120,160)  ->ProfileX()->Draw();
  //

  //
  // predefined draw
  //                         0          1          2             4         6
  //1. SetRange -             part      p          eta           ncl       par 6
  SetRange(fTPCsignal,      0,10,    0.1,100,     -0.9,0.9,    120,160,  0,6);
  SetRange(fTPCsignalNorm,  0,10,    0.1,100,     -0.9,0.9,    120,160,  0,6);
  SetRange(fTPCr,           0,10,    0.1,100,     -0.9,0.9,    120,160,  0,10);
  
  //
  // Draw
  //
  DrawPIDpt();
  DrawSignalTypeBG();
  DrawSignalTypeMom();
  DrawNormSignalTypeBG();
  DrawNormSignalTypeMom();
  DrawNormSignalTypeBB();

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



THnSparse * fTPCsignal     = 0;
THnSparse * fTPCsignalNorm = 0;
THnSparse * fTPCr = 0;
Float_t     fNorm=50;
TString selName;
Int_t version=0;

void Init(){
  TFile *f = new TFile("OutputPID.root");
  TObjArray *array= (TObjArray*)f->Get("tpcTaskPID");
  delete f;
  fTPCsignal     = (THnSparse*)array->At(0);
  fTPCsignalNorm = (THnSparse*)array->At(1);
  fTPCr          = (THnSparse*)array->At(2);

  fTPCsignal->GetAxis(5)->SetTitle("dEdx_{rec}");
  fTPCsignal->GetAxis(6)->SetTitle("dNdx_{mc}");
  fTPCsignal->GetAxis(5)->SetName("dEdx_{rec}");
  fTPCsignal->GetAxis(6)->SetName("dNdx_{mc}");

  fTPCsignalNorm->GetAxis(5)->SetTitle("dEdx_{rec}/dNdx_{mc}");
  fTPCsignalNorm->GetAxis(6)->SetTitle("dNdx_{mc}");
  fTPCsignalNorm->GetAxis(5)->SetName("dEdx_{rec}/dNdx_{mc}");
  fTPCsignalNorm->GetAxis(6)->SetName("dNdx_{mc}");
  new TCanvas("dEdx study");
}

void ResetRange(THnSparse* his){
  //
  // reset range to full
  for (Int_t idim=0;idim< his->GetNdimensions(); idim++){
    his->GetAxis(idim)->SetRange(0, his->GetAxis(idim)->GetNbins());
  }    
}

void SetRange(THnSparse*his, Float_t type0, Float_t type1, Float_t p0, Float_t p1, Float_t eta0, Float_t eta1, Int_t ncl0, Int_t ncl1, Float_t p60, Float_t p61){
  //
  //
  //
  ResetRange(his);
  selName=Form("pic/p%f_%f_eta_%f_%f_ncl_%d_%d_p6_%f_%f",p0,p1,eta0,eta1,ncl0,ncl1,p60,p61);
  his->GetAxis(0)->SetRangeUser(type0,type1);
  his->GetAxis(1)->SetRangeUser(p0,p1);
  his->GetAxis(2)->SetRangeUser(eta0,eta1);
  his->GetAxis(4)->SetRangeUser(ncl0,ncl1);
  his->GetAxis(6)->SetRangeUser(p60,p61);
  version++;
}




TH2F* GetProjection(THnSparse*his, Int_t i0, Int_t i1, Float_t type0, Float_t type1, Float_t p0, Float_t p1, Float_t eta0, Float_t eta1, Int_t ncl0, Int_t ncl1, Float_t p60, Float_t p61){

  SetRange(his,  type0, type1,   p0, p1,  eta0,eta1 ,   ncl0, ncl1,   p60,p61);
  TH2F * res = (TH2F*) his->Projection(i0,i1);
  res->SetXTitle(his->GetAxis(i1)->GetTitle());
  res->SetYTitle(his->GetAxis(i0)->GetTitle());
  return res;
}



void DrawSignalTypeBG(){
  //
  //
  //
  //
  TObjArray array;
  TH2F *hisdEdxBG_part[10];
  TH1  *hisdEdxBG_partProfile[10];
  //  TH1  *hisdEdxBG_partFit[10];
  for (Int_t ipart=0;ipart<10;ipart++){  
    fTPCsignal->GetAxis(0)->SetRange(1+ipart, 1+ipart);
    hisdEdxBG_part[ipart] = (TH2F*)fTPCsignal->Projection(5,3);
    if (!hisdEdxBG_part[ipart]) continue;
    hisdEdxBG_partProfile[ipart]=(TH1F*)(hisdEdxBG_part[ipart]->ProfileX(Form("SignalTypeBG_%d_%dProf",ipart,version)));
    //    hisdEdxBG_partProfile[ipart]->SetName(;
    hisdEdxBG_partProfile[ipart]->SetMaximum(3);
    hisdEdxBG_partProfile[ipart]->SetXTitle("#beta#gamma");
    hisdEdxBG_partProfile[ipart]->SetYTitle("dEdx/MIP");
    hisdEdxBG_partProfile[ipart]->Scale(1/fNorm);
    hisdEdxBG_partProfile[ipart]->SetMarkerColor(kmicolors[ipart%5]);
    hisdEdxBG_partProfile[ipart]->SetMarkerStyle(22+ipart);
    delete hisdEdxBG_part[ipart];
  }
  if (gPad) gPad->SetLogx(kTRUE);
  hisdEdxBG_partProfile[0]->Draw("");
  TLegend * legend = new TLegend(.7,.7, .99, .99, "TPC dEdx");
  for (Int_t ipart=0;ipart<10;ipart++)  {
    hisdEdxBG_partProfile[ipart]->Draw("same");
    legend->AddEntry(hisdEdxBG_partProfile[ipart],AliPID::ParticleName(ipart%5),"p");
  }
  legend->Draw();
  gPad->SaveAs(selName+"_hisdEdxBG.gif");
  gPad->SaveAs(selName+"_hisdEdxBG.eps");
}

void DrawSignalTypeMom(){
  //
  //
  //
  TObjArray array;
  TH2F *hisdEdxMom_part[10];
  TH1  *hisdEdxMom_partProfile[10];
  //  TH1  *hisdEdxMom_partFit[10];
  for (Int_t ipart=0;ipart<10;ipart++){  
    fTPCsignal->GetAxis(0)->SetRange(1+ipart, 1+ipart);
    hisdEdxMom_part[ipart] = (TH2F*)fTPCsignal->Projection(5,1);
    if (!hisdEdxMom_part[ipart]) continue;
    hisdEdxMom_partProfile[ipart]=(TH1F*)(hisdEdxMom_part[ipart]->ProfileX(Form("SignalTypeMom_%d_%dProf",ipart,version)));
    
    hisdEdxMom_partProfile[ipart]->SetMaximum(5);
    hisdEdxMom_partProfile[ipart]->SetMinimum(0);
    hisdEdxMom_partProfile[ipart]->SetXTitle("p (GeV/c)");
    hisdEdxMom_partProfile[ipart]->SetYTitle("dEdx/MIP");
    hisdEdxMom_partProfile[ipart]->Scale(1/fNorm);
    hisdEdxMom_partProfile[ipart]->SetMarkerColor(kmicolors[ipart%5]);
    hisdEdxMom_partProfile[ipart]->SetMarkerStyle(22+ipart);
    delete hisdEdxMom_part[ipart];
  } 
  if (gPad) gPad->SetLogx(kTRUE);
  hisdEdxMom_partProfile[0]->Draw("");
  TLegend * legend = new TLegend(.7,.7, .99, .99, "TPC dEdx");
  for (Int_t ipart=0;ipart<10;ipart++)  {
    hisdEdxMom_partProfile[ipart]->Draw("same");
    legend->AddEntry(hisdEdxMom_partProfile[ipart],AliPID::ParticleName(ipart%5),"p");
  }
  legend->Draw();
  gPad->SaveAs(selName+"_hisdEdxMom.gif");
  gPad->SaveAs(selName+"_hisdEdxMom.eps");
}


void DrawNormSignalTypeBG(){
  //
  //
  //
  //
  TObjArray array;
  TH2F *hisdEdxNormBG_part[10];
  TH1  *hisdEdxNormBG_partProfile[10];
  //  TH1  *hisdEdxNormBG_partFit[10];
  for (Int_t ipart=0;ipart<10;ipart++){  
    fTPCsignalNorm->GetAxis(0)->SetRange(1+ipart, 1+ipart);
    hisdEdxNormBG_part[ipart] = (TH2F*)fTPCsignalNorm->Projection(5,3);
    if (!hisdEdxNormBG_part[ipart]) continue;
    hisdEdxNormBG_partProfile[ipart]=(TH1F*)(hisdEdxNormBG_part[ipart]->ProfileX(Form("NormSignalTypeBG_%d_%dProf",ipart,version)));
    //    hisdEdxNormBG_partProfile[ipart]->SetName(;
    hisdEdxNormBG_partProfile[ipart]->SetMaximum(1.5);
    hisdEdxNormBG_partProfile[ipart]->SetMinimum(0.5);
    hisdEdxNormBG_partProfile[ipart]->SetXTitle("#beta#gamma");
    hisdEdxNormBG_partProfile[ipart]->SetYTitle("dEdx/MIP");
    hisdEdxNormBG_partProfile[ipart]->Scale(1/fNorm);
    hisdEdxNormBG_partProfile[ipart]->SetMarkerColor(kmicolors[ipart%5]);
    hisdEdxNormBG_partProfile[ipart]->SetMarkerStyle(22+ipart);
    delete hisdEdxNormBG_part[ipart];
  }
  if (gPad) gPad->SetLogx(kTRUE);
  hisdEdxNormBG_partProfile[0]->Draw("");
  TLegend * legend = new TLegend(.7,.7, .99, .99, "TPC dEdx");
  for (Int_t ipart=0;ipart<10;ipart++)  {
    hisdEdxNormBG_partProfile[ipart]->Draw("same");
    legend->AddEntry(hisdEdxNormBG_partProfile[ipart],AliPID::ParticleName(ipart%5),"p");
  }
  legend->Draw();
  gPad->SaveAs(selName+"_hisdEdxNormBG.gif");
  gPad->SaveAs(selName+"_hisdEdxNormBG.eps");    
}


void DrawNormSignalTypeMom(){
  //
  //
  //
  //
  //
  TObjArray array;
  TH2F *hisdEdxNormMom_part[10];
  TH1  *hisdEdxNormMom_partProfile[10];
  //  TH1  *hisdEdxNormMom_partFit[10];
  for (Int_t ipart=0;ipart<10;ipart++){  
    fTPCsignalNorm->GetAxis(0)->SetRange(1+ipart, 1+ipart);
    hisdEdxNormMom_part[ipart] = (TH2F*)fTPCsignalNorm->Projection(5,1);
    if (!hisdEdxNormMom_part[ipart]) continue;
    hisdEdxNormMom_partProfile[ipart]=(TH1F*)(hisdEdxNormMom_part[ipart]->ProfileX(Form("NormSignalTypeMom_%d_%dProf",ipart,version)));
    hisdEdxNormMom_partProfile[ipart]->Scale(1/fNorm);
    hisdEdxNormMom_partProfile[ipart]->SetMaximum(1.5);
    hisdEdxNormMom_partProfile[ipart]->SetMinimum(0.5);
    hisdEdxNormMom_partProfile[ipart]->SetXTitle("p (GeV/c)");
    hisdEdxNormMom_partProfile[ipart]->SetYTitle("dEdx_{rec}/dNdx_{mc}");
    hisdEdxNormMom_partProfile[ipart]->SetMarkerColor(kmicolors[ipart%5]);
    hisdEdxNormMom_partProfile[ipart]->SetMarkerStyle(22+ipart);
    delete hisdEdxNormMom_part[ipart];
  }
  if (gPad) gPad->SetLogx(kTRUE);
  hisdEdxNormMom_partProfile[0]->Draw("");
  TLegend * legend = new TLegend(.7,.7, .99, .99, "TPC dEdx");
  for (Int_t ipart=0;ipart<10;ipart++)  {
    hisdEdxNormMom_partProfile[ipart]->Draw("same");
    legend->AddEntry(hisdEdxNormMom_partProfile[ipart],AliPID::ParticleName(ipart%5),"p");
  }
  legend->Draw();
  gPad->SaveAs(selName+"_hisdEdxNormMom.gif");
  gPad->SaveAs(selName+"_hisdEdxNormMom.eps");  
}



void DrawNormSignalTypeBB(){
  //
  //
  //
  //
  //
  TObjArray array;
  TH2F *hisdEdxNormBG_part[10];
  TH1  *hisdEdxNormBG_partProfile[10];
  //  TH1  *hisdEdxNormBG_partFit[10];
  Double_t itypes[10];
  Double_t ival[10];
  Double_t ierr[10];
  for (Int_t ipart=0;ipart<10;ipart++){  
    fTPCsignalNorm->GetAxis(0)->SetRange(1+ipart, 1+ipart);
    hisdEdxNormBG_part[ipart] = (TH2F*)fTPCsignalNorm->Projection(5,6);
    if (!hisdEdxNormBG_part[ipart]) continue;
    hisdEdxNormBG_partProfile[ipart]=(TH1F*)(hisdEdxNormBG_part[ipart]->ProfileX(Form("NormSignalTypeBB_%d_%dProf",ipart,version)));
    hisdEdxNormBG_partProfile[ipart]->Scale(1/fNorm);
    hisdEdxNormBG_partProfile[ipart]->SetMaximum(1.5);
    hisdEdxNormBG_partProfile[ipart]->SetMinimum(0.9);
    hisdEdxNormBG_partProfile[ipart]->SetXTitle("dNdx_{mc}");
    hisdEdxNormBG_partProfile[ipart]->SetYTitle("dEdx_{rec}/dNdx_{mc}");
    hisdEdxNormBG_partProfile[ipart]->SetMarkerColor(kmicolors[ipart%5]);
    hisdEdxNormBG_partProfile[ipart]->SetMarkerStyle(22+ipart);
    itypes[ipart]=ipart;
    itypes[ipart]=0;
    itypes[ipart]=0;
    if (ipart%5!=0) {
      hisdEdxNormBG_partProfile[ipart]->Fit("pol1","q","",1,2);
      ival[ipart]= hisdEdxNormBG_partProfile[ipart]->GetFunction("pol1")->GetParameter(1);
      ierr[ipart]= hisdEdxNormBG_partProfile[ipart]->GetFunction("pol1")->GetParError(1);
    }
    printf("%d\t\t%s%f\t%f\n",ipart,AliPID::ParticleName(ipart%5), ival[ipart], ierr[ipart]);
    delete hisdEdxNormBG_part[ipart];
  }
  if (gPad) gPad->SetLogx(kFALSE);
  hisdEdxNormBG_partProfile[0]->Draw("");
  TLegend * legend = new TLegend(.7,.7, .99, .99, "TPC dEdx");
  for (Int_t ipart=0;ipart<10;ipart++)  {
    hisdEdxNormBG_partProfile[ipart]->Draw("same");
    legend->AddEntry(hisdEdxNormBG_partProfile[ipart],AliPID::ParticleName(ipart%5),"p");
  }
  legend->Draw();
  gPad->SaveAs(selName+"_hisdEdxNormBG.gif");
  gPad->SaveAs(selName+"_hisdEdxNormBG.eps");  
}


void DrawPIDpt(){
  //
  //
  //
  //
  TObjArray array;
  TH2F *hisPIDMCRCSame_part[10];
  TH1  *hisPIDMCRCSame_partProfile[10];

  //  TH1  *hisPIDMCRCSame_partFit[10];
  for (Int_t ipart=0;ipart<5;ipart++){  
    fTPCr->GetAxis(0)->SetRange(1+ipart,1+ipart);
    fTPCr->GetAxis(6)->SetRange(1+ipart,1+ipart);
    hisPIDMCRCSame_part[ipart] = (TH2F*)fTPCr->Projection(5,1);
    if (!hisPIDMCRCSame_part[ipart]) continue;
    hisPIDMCRCSame_partProfile[ipart]=(hisPIDMCRCSame_part[ipart]->ProfileX(Form("PIDpt_%d_%dProf",ipart,version)));
    hisPIDMCRCSame_partProfile[ipart]->SetDirectory(0);
    hisPIDMCRCSame_partProfile[ipart]->SetMaximum(1.0);
    hisPIDMCRCSame_partProfile[ipart]->SetMinimum(0.0);
    hisPIDMCRCSame_partProfile[ipart]->SetXTitle("p (GeV/c)");
    hisPIDMCRCSame_partProfile[ipart]->SetYTitle("TPCr");
    hisPIDMCRCSame_partProfile[ipart]->SetMarkerColor(kmicolors[ipart%5]);
    hisPIDMCRCSame_partProfile[ipart]->SetMarkerStyle(22+ipart);
    delete hisPIDMCRCSame_part[ipart];
  }
  if (gPad) gPad->SetLogx(kTRUE);
  hisPIDMCRCSame_partProfile[0]->Draw("");
  TLegend * legend = new TLegend(.7,.7, .99, .99, "TPC dEdx");
  for (Int_t ipart=0;ipart<5;ipart++)  {
    hisPIDMCRCSame_partProfile[ipart]->Draw("same");
    legend->AddEntry(hisPIDMCRCSame_partProfile[ipart],AliPID::ParticleName(ipart%5),"p");
  }
  legend->Draw();
  gPad->SaveAs(selName+"pidPIDMCRCSame.gif");
  gPad->SaveAs(selName+"pidPIDMCRCSame.eps");    
}
