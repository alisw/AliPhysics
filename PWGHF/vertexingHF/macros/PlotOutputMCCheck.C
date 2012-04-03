#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TFile.h>
#endif

/* $Id$ */ 

// Macro to plot the output of AliAnalysisTaskCheckHFMCProd
// Author: F. Prino, prino@to.infn.it


void PlotOutputMCCheck(){
  TFile *fil=new TFile("AnalysisResultsMerged.root");
  TDirectoryFile* df=(TDirectoryFile*)fil->Get("HFMCCheck");
  TList* l=(TList*)df->Get("clistHFMCCheck");
  l->ls();

  TH1F* hNEvents=(TH1F*)l->FindObject("hNEvents");
  Int_t nAnalEv=hNEvents->GetBinContent(1);
  printf("Number of events= %d\n",nAnalEv);


  TCanvas* cv=new TCanvas("cv","Vertex");
  cv->Divide(3,3);
  cv->cd(1);
  TH1F* hSPD3DvX=(TH1F*)l->FindObject("hSPD3DvX");
  hSPD3DvX->Draw();
  cv->cd(2);
  TH1F* hSPD3DvY=(TH1F*)l->FindObject("hSPD3DvY");
  hSPD3DvY->Draw();
  cv->cd(3);
  TH1F* hSPD3DvZ=(TH1F*)l->FindObject("hSPD3DvZ");
  hSPD3DvZ->Draw();
  cv->cd(4);
  TH1F* hSPDZvX=(TH1F*)l->FindObject("hSPDZvX");
  hSPDZvX->Draw();
  cv->cd(5);
  TH1F* hSPDZvY=(TH1F*)l->FindObject("hSPDZvY");
  hSPDZvY->Draw();
  cv->cd(6);
  TH1F* hSPDZvZ=(TH1F*)l->FindObject("hSPDZvZ");
  hSPDZvZ->Draw();
  cv->cd(7);
  TH1F* hTRKvX=(TH1F*)l->FindObject("hTRKvX");
  hTRKvX->Draw();
  cv->cd(8);
  TH1F* hTRKvY=(TH1F*)l->FindObject("hTRKvY");
  hTRKvY->Draw();
  cv->cd(9);
  TH1F* hTRKvZ=(TH1F*)l->FindObject("hTRKvZ");
  hTRKvZ->Draw();

  TCanvas* c1=new TCanvas("c1","Multipl");
  c1->Divide(3,1);
  c1->cd(1);
  TH1F* hTracklets=(TH1F*)l->FindObject("hTracklets");
  hTracklets->Draw();
  c1->cd(2);
  TH1F* hTracks=(TH1F*)l->FindObject("hTracks");
  hTracks->Draw();
  c1->cd(3);
  TH1F* hSelTracks=(TH1F*)l->FindObject("hSelTracks");
  hSelTracks->Draw();

  TH1F* hncharmed=(TH1F*)l->FindObject("hncharmed");
  TCanvas* cn=new TCanvas("cn","ncharm");
  hncharmed->Draw("box");
  hncharmed->GetXaxis()->SetTitle("dNch/dy");
  hncharmed->GetYaxis()->SetTitle("N Charm hadrons in golden channels");
  cn->Update();

  TH1F* hnbvsnc=(TH1F*)l->FindObject("hnbvsnc");
  TCanvas* cnhf=new TCanvas("cnhf","nb/c");
  hnbvsnc->Draw("box");
  hnbvsnc->GetXaxis()->SetTitle("Nc");
  hnbvsnc->GetYaxis()->SetTitle("Nb");
  cnhf->Update();

  TH2F* hyptd0prompt=(TH2F*)l->FindObject("hyptd0prompt");
  TH2F* hyptd0feeddown=(TH2F*)l->FindObject("hyptd0feeddown");
  TH2F* hyptD02=(TH2F*)l->FindObject("hyptD02");
  TH2F* hyptD04=(TH2F*)l->FindObject("hyptD04");

  TCanvas* cd0=new TCanvas("cd0","D0");
  cd0->Divide(2,2);
  cd0->cd(1);
  hyptd0prompt->Draw("colz");
  cd0->cd(2);
  hyptd0feeddown->Draw("colz");
  cd0->cd(3);
  hyptD02->Draw("colz");
  cd0->cd(4);
  hyptD04->Draw("colz");

  TH2F* hyptdplusprompt=(TH2F*)l->FindObject("hyptdplusprompt");
  TH2F* hyptdplusfeeddown=(TH2F*)l->FindObject("hyptdplusfeeddown");
  TH2F* hyptDplusnonreson=(TH2F*)l->FindObject("hyptDplusnonreson");
  TH2F* hyptDplusreson=(TH2F*)l->FindObject("hyptDplusreson");

  TCanvas* cdplus=new TCanvas("cdplus","Dplus");
  cdplus->Divide(2,2);
  cdplus->cd(1);
  hyptdplusprompt->Draw("colz");
  cdplus->cd(2);
  hyptdplusfeeddown->Draw("colz");
  cdplus->cd(3);
  hyptDplusnonreson->Draw("colz");
  cdplus->cd(4);
  hyptDplusreson->Draw("colz");

  TH2F* hyptdsprompt=(TH2F*)l->FindObject("hyptdsprompt");
  TH2F* hyptdsfeeddown=(TH2F*)l->FindObject("hyptdsfeedown");
  TH2F* hyptdsphi=(TH2F*)l->FindObject("hyptdsphi");
  TH2F* hyptdsK0st=(TH2F*)l->FindObject("hyptdsk0st");
  
  TCanvas* cds=new TCanvas("cds","Ds");
  cds->Divide(2,2);
  cds->cd(1);
  hyptdsprompt->Draw("colz");
  cds->cd(2);
  hyptdsfeeddown->Draw("colz");
  cds->cd(3);
  hyptdsphi->Draw("colz");
  cds->cd(4);
  hyptdsK0st->Draw("colz");

  TH2F* hyptdstarprompt=(TH2F*)l->FindObject("hyptdstarprompt");
  TH2F* hyptdstarfeedown=(TH2F*)l->FindObject("hyptdstarfeedown");
  TH2F* hyptlcprompt=(TH2F*)l->FindObject("hyptlcprompt");
  TH2F* hyptlcfeedown=(TH2F*)l->FindObject("hyptlcfeedown");

  TCanvas* cdstlc=new TCanvas("cdstls","Dstar LambdaC");
  cdstlc->Divide(2,2);
  cdstlc->cd(1);
  hyptdstarprompt->Draw("colz");
  cdstlc->cd(2);
  hyptdstarfeedown->Draw("colz");
  cdstlc->cd(3);
  hyptlcprompt->Draw("colz");
  cdstlc->cd(4);
  hyptlcfeedown->Draw("colz");

}
