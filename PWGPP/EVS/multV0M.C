#include "periodLevelQA.C"
#include "map"
using namespace std;

void multV0M(){
  gStyle->SetOptStat(0);
  gStyle->SetLineScalePS(1.5);
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(1);
//  TFile* fin = new TFile("EventStat_temp.root");
//  TList* statsout = (TList*) fin->Get("cstatsout");
//  AliPhysicsSelection* ps = statsout ? (AliPhysicsSelection*) statsout->FindObject("AliPhysicsSelection") : 0;
//  AliTriggerAnalysis* taMB = ps->GetTriggerAnalysis(2);
//  AliTriggerAnalysis* taHM = ps->GetTriggerAnalysis(20);
//  TList* listMB = taMB->GetHistList();
//  TList* listHM = taHM->GetHistList();
//  TH1F* hV0MOnAccMB = (TH1F*) listMB->FindObject("fHistV0MOfAcc");
//  TH1F* hV0MOnAccHM = (TH1F*) listHM->FindObject("fHistV0MOfAcc");

//  TFile* fV0M_MB = new TFile("MB_V0M.root"); 
//  TFile* fV0M_HM = new TFile("HM_V0M.root");
//  TH1F* hV0MOnAccMB = (TH1F*) fV0M_MB->Get("hOnNorm_kINT7_V0M");
//  TH1F* hV0MOnAccHM = (TH1F*) fV0M_HM->Get("hOnNorm_kHighMultV0_V0M");
  TFile* fV0M_MB = new TFile("MB_OFO.root"); 
  TFile* fV0M_HM = new TFile("HM_OFO.root");
  TH1F* hV0MOnAccMB = (TH1F*) fV0M_MB->Get("hOfNorm_kINT7_TKL");
  TH1F* hV0MOnAccHM = (TH1F*) fV0M_HM->Get("hOfNorm_kHighMultV0_TKL");

  Float_t mb = hV0MOnAccMB->Integral(0,hV0MOnAccMB->GetNbinsX()+1);
  Float_t hm = hV0MOnAccHM->Integral(0,hV0MOnAccHM->GetNbinsX()+1);
  
  TCanvas* c = new TCanvas("c1","c1",1000,800);
  c->SetMargin(0.08,0.02,0.08,0.02);
  c->SetLogy();
  hV0MOnAccMB->SetTitle("");
  hV0MOnAccMB->GetYaxis()->SetTitle("Events");
  hV0MOnAccMB->Draw();
  hV0MOnAccHM->SetLineColor(kRed);
  hV0MOnAccHM->Draw("same");
  
  TLegend* l = new TLegend(0.70,0.8,0.97,0.97);
  l->AddEntry(hV0MOnAccMB,Form("MB: %.0f M",mb/1000000));
  l->AddEntry(hV0MOnAccHM,Form("HM: %.0f M",hm/1000000));
  l->Draw();
  gPad->Print("statTKL.png");
}
