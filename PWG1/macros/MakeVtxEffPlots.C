#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TSystem.h"
#include "TROOT.h"
#include <TString.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TFile.h>
#endif


void setDrawAtt(Int_t markerstyle,Int_t markercolor,Int_t markersize,Int_t linecolor,Int_t linewidth,TGraphAsymmErrors *h1)
{ 
  h1->SetMarkerStyle(markerstyle);
  h1->SetMarkerColor(markercolor);
  h1->SetMarkerSize(markersize);
  h1->SetLineColor(linecolor);
  h1->SetLineWidth(linewidth);
  }
void setDrawAtt(Int_t markerstyle,Int_t markercolor,Int_t markersize,Int_t linecolor,Int_t linewidth,TH1 *h1)
{ 
  h1->SetMarkerStyle(markerstyle);
  h1->SetMarkerColor(markercolor);
  h1->SetMarkerSize(markersize);
  h1->SetLineColor(linecolor);
  h1->SetLineWidth(linewidth);
  }

void MakeVtxEffPlots(TString foldname="",Bool_t savefile=kFALSE){
  
  //Make efficiency plots for SPD vertex, TRK Vertex and TPC Vertex
  // Author: 
  // L. Milano, milano@to.infn.it
  
  TFile *infil=new TFile(Form("%sAnalysisResults.root",foldname.Data()),"read");
  TDirectory *dirFile=(TDirectory*)infil->Get("Vertex_Performance");
  TList *cOutput = (TList*)dirFile->Get("cOutputVtxESD");
    
  TH1F *fhTriggeredTrklets=(TH1F*)cOutput->FindObject("fhTriggeredTrklets");
  TH1F *fhSPDZTrklets=(TH1F*)cOutput->FindObject("fhSPDZTrklets");
  TH1F *fhSPD3DTrklets=(TH1F*)cOutput->FindObject("fhSPD3DTrklets");
  TH1F *fhTRKTrklets=(TH1F*)cOutput->FindObject("fhTRKTrklets");
  TH1F *fhTRKcTrklets=(TH1F*)cOutput->FindObject("fhTRKcTrklets");
  TH1F *fhTRKncTrklets=(TH1F*)cOutput->FindObject("fhTRKncTrklets");
  TH1F *fhTPCTrklets=(TH1F*)cOutput->FindObject("fhTPCTrklets");
  TH1F *fhTPCcTrklets=(TH1F*)cOutput->FindObject("fhTPCcTrklets");
  TH1F *fhTPCncTrklets=(TH1F*)cOutput->FindObject("fhTPCncTrklets");
  TH1F *fhSPDZZreco=(TH1F*)cOutput->FindObject("fhSPDZZreco");
  TH1F *fhSPD3DZreco=(TH1F*)cOutput->FindObject("fhSPD3DZreco");
  
  TGraphAsymmErrors *fhSPDZEffTrklets=new TGraphAsymmErrors(fhSPDZTrklets,fhTriggeredTrklets,"w");
  fhSPDZEffTrklets->SetName("fhSPDZEffTrklets");
  fhSPDZEffTrklets->SetDrawOption("AP");
  TGraphAsymmErrors *fhSPD3DEffTrklets=new TGraphAsymmErrors(fhSPD3DTrklets,fhTriggeredTrklets,"w");
  fhSPD3DEffTrklets->SetName("fhSPD3DEffTrklets");
  TH1F * fhSPDOverallTrklets=(TH1F*)fhSPDZTrklets->Clone("fhSPDOverallTrklets");
  fhSPDOverallTrklets->Add(fhSPD3DTrklets);
  TGraphAsymmErrors *fhSPDOverallEffTrklets=new TGraphAsymmErrors(fhSPDOverallTrklets,fhTriggeredTrklets,"w");
  fhSPDOverallEffTrklets->SetName("fhSPDOverallEffTrklets");
  TGraphAsymmErrors *fhTRKEffTrklets=new TGraphAsymmErrors(fhTRKTrklets,fhTriggeredTrklets,"w");
  fhTRKEffTrklets->SetName("fhTRKEffTrklets");
  TGraphAsymmErrors *fhTRKcEffTrklets=new TGraphAsymmErrors(fhTRKcTrklets,fhTriggeredTrklets,"w");
  fhTRKcEffTrklets->SetName("fhTRKcEffTrklets");
  TGraphAsymmErrors *fhTRKncEffTrklets=new TGraphAsymmErrors(fhTRKncTrklets,fhTriggeredTrklets,"w");
  fhTRKncEffTrklets->SetName("fhTRKncEffTrklets");
  TGraphAsymmErrors *fhTPCEffTrklets=new TGraphAsymmErrors(fhTPCTrklets,fhTriggeredTrklets,"w");
  fhTPCEffTrklets->SetName("fhTPCEffTrklets");
  TGraphAsymmErrors *fhTPCcEffTrklets=new TGraphAsymmErrors(fhTPCcTrklets,fhTriggeredTrklets,"w");
  fhTPCcEffTrklets->SetName("fhTPCcEffTrklets");
  TGraphAsymmErrors *fhTPCncEffTrklets=new TGraphAsymmErrors(fhTPCncTrklets,fhTriggeredTrklets,"w");
  fhTPCncEffTrklets->SetName("fhTPCncEffTrklets");
  TH1F * fhSPDOverallZreco=(TH1F*)fhSPDZZreco->Clone("fhSPDOverallZreco");
  fhSPDOverallZreco->Add(fhSPD3DZreco);
  TGraphAsymmErrors *fhSPDEffZreco=new TGraphAsymmErrors(fhSPD3DZreco,fhSPDOverallZreco,"w");
  fhSPDEffZreco->SetName("fhSPDEffZreco");

  TH1F *fhEff = new TH1F("hEff","hEff",9,0.5,9.5);
  Int_t count=1;
  if(fhSPDZTrklets->GetEntries()!=0 && fhTriggeredTrklets->GetEntries()!=0){
    fhEff->Fill(count,fhSPDZTrklets->GetEntries()/fhTriggeredTrklets->GetEntries());
    fhEff->SetBinError(count,fhEff->GetBinContent(count)*TMath::Sqrt(1/fhSPDZTrklets->GetEntries()+1/fhTriggeredTrklets->GetEntries()));
  }
  fhEff->GetXaxis()->SetBinLabel(count,"SPDZ");
  
  count++;
  if(fhSPD3DTrklets->GetEntries()!=0 && fhTriggeredTrklets->GetEntries()!=0){
    fhEff->Fill(count,fhSPD3DTrklets->GetEntries()/fhTriggeredTrklets->GetEntries());
    fhEff->SetBinError(count,fhEff->GetBinContent(count)*TMath::Sqrt(1/fhSPD3DTrklets->GetEntries()+1/fhTriggeredTrklets->GetEntries()));
  }
  fhEff->GetXaxis()->SetBinLabel(count,"SPD3D");
  
  count++;
  if(fhSPDOverallTrklets->GetEntries()!=0 && fhTriggeredTrklets->GetEntries()!=0){
    fhEff->Fill(count,fhSPDOverallTrklets->GetEntries()/fhTriggeredTrklets->GetEntries());
    fhEff->SetBinError(count,fhEff->GetBinContent(count)*TMath::Sqrt(1/fhSPDOverallTrklets->GetEntries()+1/fhTriggeredTrklets->GetEntries()));
  }
  fhEff->GetXaxis()->SetBinLabel(count,"SPD Overall");
  
  count++;
  if(fhTRKTrklets->GetEntries()!=0 && fhTriggeredTrklets->GetEntries()!=0){
    fhEff->Fill(count,fhTRKTrklets->GetEntries()/fhTriggeredTrklets->GetEntries());
    fhEff->SetBinError(count,fhEff->GetBinContent(count)*TMath::Sqrt(1/fhTRKTrklets->GetEntries()+1/fhTriggeredTrklets->GetEntries()));
  }
  fhEff->GetXaxis()->SetBinLabel(count,"TRK");
  
  count++;
  if(fhTRKcTrklets->GetEntries()!=0 && fhTriggeredTrklets->GetEntries()!=0){
    fhEff->Fill(count,fhTRKcTrklets->GetEntries()/fhTriggeredTrklets->GetEntries());
    fhEff->SetBinError(count,fhEff->GetBinContent(count)*TMath::Sqrt(1/fhTRKcTrklets->GetEntries()+1/fhTriggeredTrklets->GetEntries()));
  }
  fhEff->GetXaxis()->SetBinLabel(count,"TRKc");
  
  count++;
  if(fhTRKncTrklets->GetEntries()!=0 && fhTriggeredTrklets->GetEntries()!=0){
    fhEff->Fill(count,fhTRKncTrklets->GetEntries()/fhTriggeredTrklets->GetEntries());
    fhEff->SetBinError(count,fhEff->GetBinContent(count)*TMath::Sqrt(1/fhTRKncTrklets->GetEntries()+1/fhTriggeredTrklets->GetEntries()));
  }
  fhEff->GetXaxis()->SetBinLabel(count,"TRKnc");
  
  count++;
  if(fhTPCTrklets->GetEntries()!=0 && fhTriggeredTrklets->GetEntries()!=0){
    fhEff->Fill(count,fhTPCTrklets->GetEntries()/fhTriggeredTrklets->GetEntries());
    fhEff->SetBinError(count,fhEff->GetBinContent(count)*TMath::Sqrt(1/fhTPCTrklets->GetEntries()+1/fhTriggeredTrklets->GetEntries()));
  }
  fhEff->GetXaxis()->SetBinLabel(count,"TPC");
  
  count++;
  if(fhTPCcTrklets->GetEntries()!=0 && fhTriggeredTrklets->GetEntries()!=0){
    fhEff->Fill(count,fhTPCcTrklets->GetEntries()/fhTriggeredTrklets->GetEntries());
    fhEff->SetBinError(count,fhEff->GetBinContent(count)*TMath::Sqrt(1/fhTPCcTrklets->GetEntries()+1/fhTriggeredTrklets->GetEntries()));
  }
  fhEff->GetXaxis()->SetBinLabel(count,"TPCc");
  
  count++;
  if(fhTPCncTrklets->GetEntries()!=0 && fhTriggeredTrklets->GetEntries()!=0){
    fhEff->Fill(count,fhTPCncTrklets->GetEntries()/fhTriggeredTrklets->GetEntries());
    fhEff->SetBinError(count,fhEff->GetBinContent(count)*TMath::Sqrt(1/fhTPCncTrklets->GetEntries()+1/fhTriggeredTrklets->GetEntries()));
  }
  fhEff->GetXaxis()->SetBinLabel(count,"TPCnc");
  
  TCanvas *canvSPDTrklets=new TCanvas("canvSPDTrklets","SPD Eff vs tracklet multiplicy");
  canvSPDTrklets->SetBottomMargin(0.14);
  canvSPDTrklets->SetTopMargin(0.08);
  canvSPDTrklets->SetLeftMargin(0.14);
  canvSPDTrklets->SetRightMargin(0.08);
  fhSPDZEffTrklets->SetMinimum(0);
  fhSPDZEffTrklets->SetMaximum(1.2);
  fhSPDZEffTrklets->GetXaxis()->SetTitle("tracklet multiplicity");
  fhSPDZEffTrklets->GetXaxis()->SetTitleSize(0.05);
  fhSPDZEffTrklets->GetYaxis()->SetTitle("efficiency");
  fhSPDZEffTrklets->GetYaxis()->SetTitleSize(0.05);
  setDrawAtt(22,4,2,4,2,fhSPDZEffTrklets);
  setDrawAtt(23,2,2,2,2,fhSPD3DEffTrklets);
  setDrawAtt(24,1,2,1,2,fhSPDOverallEffTrklets);
  fhSPDZEffTrklets->SetTitle("SPDZ efficiency");
  fhSPDZEffTrklets->Draw("ALP");
  fhSPD3DEffTrklets->Draw("PLSAME");
  fhSPDOverallEffTrklets->Draw("PLSAME"); 
  TLegend *legSPD=new TLegend(0.6,0.2,0.9,0.4);
  legSPD->AddEntry(fhSPD3DEffTrklets,"SPD3D","P");
  legSPD->AddEntry(fhSPDZEffTrklets,"SPDZ","P");
  legSPD->AddEntry(fhSPDOverallEffTrklets,"SPDZ + SPD3D","P");
  legSPD->Draw();
  
  TCanvas *canvTRKTrklets=new TCanvas("canvTRKTrklets","TRK Eff vs tracklet multiplicy");
  canvTRKTrklets->SetBottomMargin(0.14);
  canvTRKTrklets->SetTopMargin(0.08);
  canvTRKTrklets->SetLeftMargin(0.14);
  canvTRKTrklets->SetRightMargin(0.08);
  fhTRKEffTrklets->SetMinimum(0);
  fhTRKEffTrklets->SetMaximum(1.2);
  fhTRKEffTrklets->GetXaxis()->SetTitle("tracklet multiplicity");
  fhTRKEffTrklets->GetXaxis()->SetTitleSize(0.05);
  fhTRKEffTrklets->GetYaxis()->SetTitle("efficiency");
  fhTRKEffTrklets->GetYaxis()->SetTitleSize(0.05);
  setDrawAtt(22,4,2,4,2,fhTRKEffTrklets);
  setDrawAtt(23,2,2,2,2,fhTRKcEffTrklets);
  setDrawAtt(24,1,2,1,2,fhTRKncEffTrklets);
  fhTRKEffTrklets->SetTitle("TRK efficiency");
  fhTRKEffTrklets->Draw("ALP");
  fhTRKcEffTrklets->Draw("PLSAME");
  fhTRKncEffTrklets->Draw("PLSAME"); 
  TLegend *legTRK=new TLegend(0.6,0.2,0.9,0.4);
  legTRK->AddEntry(fhTRKEffTrklets,"TRK","P");  
  legTRK->AddEntry(fhTRKcEffTrklets,"TRKc","P");
  legTRK->AddEntry(fhTRKncEffTrklets,"TRKnc","P");
  legTRK->Draw();
  
  TCanvas *canvTPCTrklets=new TCanvas("canvTPCTrklets","TPC Eff vs tracklet multiplicy");
  canvTPCTrklets->SetBottomMargin(0.14);
  canvTPCTrklets->SetTopMargin(0.08);
  canvTPCTrklets->SetLeftMargin(0.14);
  canvTPCTrklets->SetRightMargin(0.08);
  fhTPCEffTrklets->SetMinimum(0);
  fhTPCEffTrklets->SetMaximum(1.2);
  fhTPCEffTrklets->GetXaxis()->SetTitle("tracklet multiplicity");
  fhTPCEffTrklets->GetXaxis()->SetTitleSize(0.05);
  fhTPCEffTrklets->GetYaxis()->SetTitle("efficiency");
  fhTPCEffTrklets->GetYaxis()->SetTitleSize(0.05);
  setDrawAtt(22,4,2,4,2,fhTPCEffTrklets);
  setDrawAtt(23,2,2,2,2,fhTPCcEffTrklets);
  setDrawAtt(24,1,2,1,2,fhTPCncEffTrklets);
  fhTPCEffTrklets->SetTitle("TPC efficiency");
  fhTPCEffTrklets->Draw("ALP");
  fhTPCcEffTrklets->Draw("PLSAME");
  fhTPCncEffTrklets->Draw("PLSAME"); 
  TLegend *legTPC=new TLegend(0.6,0.2,0.9,0.4);
  legTPC->AddEntry(fhTPCEffTrklets,"TPC","P");  
  legTPC->AddEntry(fhTPCcEffTrklets,"TPCc","P");
  legTPC->AddEntry(fhTPCncEffTrklets,"TPCnc","P");
  legTPC->Draw();
  
  TCanvas *canvZ=new TCanvas("canvZ","3D/reco vs Z");
  canvZ->SetBottomMargin(0.14);
  canvZ->SetTopMargin(0.08);
  canvZ->SetLeftMargin(0.14);
  canvZ->SetRightMargin(0.08);
  fhSPDEffZreco->SetMinimum(0);
  fhSPDEffZreco->SetMaximum(1.2);
  fhSPDEffZreco->GetXaxis()->SetTitle("Z - <Z> [cm]");
  fhSPDEffZreco->GetXaxis()->SetTitleSize(0.05);
  fhSPDEffZreco->GetYaxis()->SetTitle("3D rec / (3D+Z rec)");
  fhSPDEffZreco->GetYaxis()->SetTitleSize(0.05);
  setDrawAtt(22,4,2,4,2,fhSPDEffZreco);
  fhSPDEffZreco->SetTitle("SPD3D/(SPD3D+Z) vs Zreco");
  fhSPDEffZreco->Draw("AP");
  TLegend *legZ=new TLegend(0.6,0.2,0.9,0.4);
  legZ->AddEntry(fhSPDEffZreco,"SPD3D/(SPD3D+Z) vs Zreco","P");
  legZ->Draw();
  
  TCanvas *canvOverall=new TCanvas("canvOverall","Eff integrated over multiplicy");
  canvOverall->SetBottomMargin(0.14);
  canvOverall->SetTopMargin(0.08);
  canvOverall->SetLeftMargin(0.14);
  canvOverall->SetRightMargin(0.08);
  fhEff->SetMinimum(0);
  fhEff->SetMaximum(1.2);
  fhEff->GetXaxis()->SetTitleSize(0.05);
  fhEff->GetYaxis()->SetTitle("efficiency");
  fhEff->GetYaxis()->SetTitleSize(0.05);
  fhEff->SetTitle("integrated over multiplicity");
  fhEff->Draw("");

  if(savefile){
  TFile* fileEff = new TFile("VtxEff.root","recreate");
  fhSPDZEffTrklets->Write();
  fhSPD3DEffTrklets->Write();
  fhSPDOverallEffTrklets->Write();
  fhTRKEffTrklets->Write();
  fhTRKcEffTrklets->Write();
  fhTRKncEffTrklets->Write();
  fhSPDEffZreco->Write();
  fhEff->Write();
  fileEff->Close();
  delete fileEff;
  }
  
  
}



