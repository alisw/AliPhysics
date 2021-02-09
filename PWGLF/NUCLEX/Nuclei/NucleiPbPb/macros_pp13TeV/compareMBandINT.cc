#include "src/Common.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TLegend.h"

void compareMBandINT(){
    TFile file_MB(Form("%spp13TeV.mb.fullpT.INEL.FINAL-2019-05-30.root",kBaseOutputDir.data()));
    TH1F* hMB = (TH1F*)file_MB.Get("hsys_pp13_mb_proton_sum");
    hMB->SetDirectory(0);
    hMB->SetName("hMB");
    hMB->Scale(0.5);
    hMB->SetMarkerStyle(21);
    hMB->SetMarkerColor(kBlue);
    hMB->SetLineColor(kBlue);
    //
    TH1F* hMB_scaled = (TH1F*)hMB->Clone("hMB_scaled");
    hMB_scaled->SetTitle("#times #sqrt{2}");
    hMB_scaled->SetDirectory(0);
    hMB_scaled->Scale(TMath::Sqrt(2));
    hMB_scaled->SetMarkerStyle(21);
    hMB_scaled->SetMarkerColor(kGreen+3);
    hMB_scaled->SetLineColor(kGreen+3);
    //
    TH1F* hMB_scaled_bis = (TH1F*)hMB->Clone("hMB_scaled_bis");
    hMB_scaled->SetTitle("/ 0.7448");
    hMB_scaled_bis->SetDirectory(0);
    hMB_scaled_bis->Scale(1/0.7448);
    hMB_scaled_bis->SetMarkerStyle(21);
    hMB_scaled_bis->SetMarkerColor(kBlack);
    hMB_scaled_bis->SetLineColor(kBlack);
    //
    TFile file_Mult(Form("%sB2.root",kBaseOutputDir.data()));
    TGraphErrors* hMult = (TGraphErrors*)file_Mult.Get("deuterons/9/grProtSyst_M_9");
    hMult->SetName("hMult");
    hMult->SetMarkerStyle(20);
    hMult->SetMarkerColor(kRed);
    hMult->SetLineColor(kRed);
    //
    TCanvas* can = new TCanvas("cComp","cComp");
    can->SetLeftMargin(0.15);
    TH1* hFrame = can->DrawFrame(0,1e-5,4,0.1,";#it{p}_{T} (GeV/#it{c});#frac{1}{#it{N}_{ev}} #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}} (GeV/#it{c})^{-1}");
    hFrame->GetYaxis()->SetTitleOffset(1.3);
    hMult->Draw("PE");
    hMB->Draw("esame");
    hMB_scaled->Draw("esame");
    hMB_scaled_bis->Draw("esame");

    TLegend* leg = new TLegend(0.5,0.5,0.7,0.7);
    leg->AddEntry(hMult,"Integrated Multiplicity", "pe");
    leg->AddEntry(hMB,"MB", "pe");
    leg->AddEntry(hMB_scaled,"MB #times #sqrt{2}", "pe");
    leg->AddEntry(hMB_scaled_bis,"MB / 0.7448", "pe");
    leg->Draw();

    TFile fOutput(Form("%sInvSpectraMBvsInt.root",kBaseOutputDir.data()),"recreate");
    hMB->Write();
    hMB_scaled->Write();
    hMB_scaled_bis->Write();
    hMult->Write();
    can->Write();

    can->Print(Form("%simages/InvSpectraMBvsInt.pdf",kBaseOutputDir.data()));


}