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
  TFile *fil=new TFile("AnalysisResults.root");
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
  c1->Divide(2,2);
  c1->cd(1);
  gPad->SetLogy();
  TH1F* hPhysPrim=(TH1F*)l->FindObject("hPhysPrim");
  hPhysPrim->Draw();
  c1->cd(2);
  gPad->SetLogy();
  TH1F* hTracklets=(TH1F*)l->FindObject("hTracklets");
  hTracklets->Draw();
  c1->cd(3);
  gPad->SetLogy();
  TH1F* hTracks=(TH1F*)l->FindObject("hTracks");
  hTracks->Draw();
  c1->cd(4);
  gPad->SetLogy();
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
  hnbvsnc->Draw("colztext");
  hnbvsnc->GetXaxis()->SetTitle("Nc");
  hnbvsnc->GetYaxis()->SetTitle("Nb");
  cnhf->Update();
  TH2F* hyptD0all=(TH2F*)l->FindObject("hyptD0AllDecay");
  TH2F* hyptD0promptall=(TH2F*)l->FindObject("hyptD0promptAllDecay");
  TH2F* hyptD0feeddownall=(TH2F*)l->FindObject("hyptD0feeddownAllDecay");
  TH2F* hyptD0prompt=(TH2F*)l->FindObject("hyptD0prompt");
  TH2F* hyptD0feeddown=(TH2F*)l->FindObject("hyptD0feeddown");
  TH2F* hyptD02=(TH2F*)l->FindObject("hyptD02");
  TH2F* hyptD04=(TH2F*)l->FindObject("hyptD04");
  TH1D* hptD0all=hyptD0all->ProjectionX();
  TH1D* hptD0promptall=hyptD0promptall->ProjectionX();
  TH1D* hptD0feeddownall=hyptD0feeddownall->ProjectionX();
  TH2F* hyptDplusall=(TH2F*)l->FindObject("hyptDplusAllDecay");
  TH2F* hyptDpluspromptall=(TH2F*)l->FindObject("hyptDpluspromptAllDecay");
  TH2F* hyptDplusfeeddownall=(TH2F*)l->FindObject("hyptDplusfeeddownAllDecay");
  TH2F* hyptDplusprompt=(TH2F*)l->FindObject("hyptDplusprompt");
  TH2F* hyptDplusfeeddown=(TH2F*)l->FindObject("hyptDplusfeeddown");
  TH2F* hyptDplusnonreson=(TH2F*)l->FindObject("hyptDplusnonreson");
  TH2F* hyptDplusreson=(TH2F*)l->FindObject("hyptDplusreson");
  TH2F* hyptDsall=(TH2F*)l->FindObject("hyptDsAllDecay");
  TH2F* hyptDsprompt=(TH2F*)l->FindObject("hyptDsprompt");
  TH2F* hyptDsfeeddown=(TH2F*)l->FindObject("hyptDsfeedown");
  TH2F* hyptDsphi=(TH2F*)l->FindObject("hyptDsphi");
  TH2F* hyptDsK0st=(TH2F*)l->FindObject("hyptDsk0st");
  TH2F* hyptDstarall=(TH2F*)l->FindObject("hyptDstarAllDecay");
  TH2F* hyptDstarprompt=(TH2F*)l->FindObject("hyptDstarprompt");
  TH2F* hyptDstarfeedown=(TH2F*)l->FindObject("hyptDstarfeedown");
  TH2F* hyptLcprompt=(TH2F*)l->FindObject("hyptLcprompt");
  TH2F* hyptLcfeedown=(TH2F*)l->FindObject("hyptLcfeedown");
  TH2F* hyptLcall=(TH2F*)l->FindObject("hyptLcAllDecay");

  TH2F* hyptB0all=(TH2F*)l->FindObject("hyptB0AllDecay");
  TH2F* hyptBplusall=(TH2F*)l->FindObject("hyptBplusAllDecay");
  TH2F* hyptBsall=(TH2F*)l->FindObject("hyptBsAllDecay");
  TH2F* hyptBstarall=(TH2F*)l->FindObject("hyptBstarAllDecay");
  TH2F* hyptLball=(TH2F*)l->FindObject("hyptLbAllDecay");



  TH1F*	hD0fonll7=HistoFONLL7TeV();
  hD0fonll7->Scale(hptD0all->GetMaximum()/hD0fonll7->GetMaximum());  
  hD0fonll7->SetLineColor(kGreen+1);
  hD0fonll7->SetLineWidth(2);

  TH1F*	hD0fonll2=HistoFONLL2_76TeV();
  hD0fonll2->Scale(hptD0all->GetMaximum()/hD0fonll2->GetMaximum());  
  hD0fonll2->SetLineColor(kBlue+1);
  hD0fonll2->SetLineWidth(2);

  TH1F* hptD0pythia=HistoPYTHIA7();
  hptD0pythia->Scale(hptD0all->GetMaximum()/hptD0pythia->GetMaximum());
  hptD0pythia->SetLineColor(kRed+1);
  hptD0pythia->SetLineWidth(2);

  hptD0all->SetLineWidth(2);
  hptD0all->SetMarkerStyle(20);
  hptD0all->SetTitle("");
  hptD0all->GetYaxis()->SetTitle("D^{0} dN/dp_{T} (a.u.)");
  hptD0all->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hptD0all->SetStats(0);  
  hptD0all->GetYaxis()->SetTitleOffset(1.2);
  hptD0all->GetXaxis()->SetTitleOffset(1.2);

  TCanvas* cd0a=new TCanvas("cd0a","D0 spectra",700,700);
  cd0a->SetLogy();
  cd0a->SetLeftMargin(0.13);
  cd0a->SetRightMargin(0.07);
  hptD0all->Draw();
  hD0fonll7->Draw("lsame");
  hD0fonll2->Draw("lsame");
  hptD0pythia->Draw("lsame");
  TLegend* leg=new TLegend(0.45,0.6,0.89,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hptD0all,"MC production","LP");
  leg->AddEntry(hD0fonll7,"FONLL, #sqrt{s}=7 TeV","L");
  leg->AddEntry(hD0fonll2,"FONLL, #sqrt{s}=2.76 TeV","L");
  leg->AddEntry(hptD0pythia,"PYTHIA Perugia0, #sqrt{s}=7 TeV","L");
  leg->Draw();

  // Prompt and Feeddown
  TH1F* hOriginPrompt=(TH1F*)l->FindObject("hOriginPrompt");
  TH1F* hOriginFeeddown=(TH1F*)l->FindObject("hOriginFeeddown");

  hptD0promptall->SetLineColor(4);
  hptD0promptall->SetMarkerColor(4);
  hptD0promptall->SetMarkerStyle(26);
  hptD0feeddownall->SetLineColor(2);
  hptD0feeddownall->SetMarkerColor(2);
  hptD0feeddownall->SetMarkerStyle(23)
;
  TCanvas* cprf1=new TCanvas("cprf1","Prompt/Feeddown",700,700);
  gPad->SetLogy();
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.07);
  hptD0all->Draw();
  hptD0promptall->Draw("same");
  hptD0feeddownall->Draw("same");
  TLegend* leg2=new TLegend(0.4,0.5,0.89,0.85);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->AddEntry(hptD0all,"All D0","LP");
  leg2->AddEntry("",Form("Entries=%.0f  <pt>=%.2f GeV/c",hptD0all->Integral(),hptD0all->GetMean()),"");
  leg2->AddEntry(hptD0promptall,"Prompt D0","LP");
  leg2->AddEntry("",Form("Entries=%.0f  <pt>=%.2f GeV/c",hptD0promptall->Integral(),hptD0promptall->GetMean()),"");
  leg2->AddEntry(hptD0feeddownall,"Feeddown D0","LP");
  leg2->AddEntry("",Form("Entries=%.0f  <pt>=%.2f GeV/c",hptD0feeddownall->Integral(),hptD0feeddownall->GetMean()),"");
  leg2->Draw();
  hOriginPrompt->GetXaxis()->SetTitle("Distance of prompt D meson origin to vertex (cm)");
  hOriginFeeddown->GetXaxis()->SetTitle("Distance of feed-down D meson origin to vertex (cm)");

  TCanvas* cprf2=new TCanvas("cprf2","Origin",1000,600);
  cprf2->Divide(2,1);
  cprf2->cd(1);
  hOriginPrompt->Draw();
  cprf2->cd(2);
  hOriginFeeddown->Draw();

  // Hadrons
  Double_t nev=hNEvents->GetBinContent(1);
  Double_t nD0=hyptD0all->GetEntries();
  Double_t nDp=hyptDplusall->GetEntries();
  Double_t nDs=hyptDsall->GetEntries();
  Double_t nDst=hyptDstarall->GetEntries();
  Double_t nLc=hyptLcall->GetEntries();
  Double_t nB0=hyptB0all->GetEntries();
  Double_t nBp=hyptBplusall->GetEntries();
  Double_t nBs=hyptBsall->GetEntries();
  Double_t nBst=hyptBstarall->GetEntries();
  Double_t nLb=hyptLball->GetEntries();


  Double_t nccbar=0;
  Double_t nbbbar=0;
  for(Int_t iBinx=1; iBinx<=hnbvsnc->GetNbinsX(); iBinx++){
    for(Int_t iBiny=1; iBiny<=hnbvsnc->GetNbinsY(); iBiny++){
      Double_t bincentx=hnbvsnc->GetXaxis()->GetBinCenter(iBinx);
      Double_t bincenty=hnbvsnc->GetYaxis()->GetBinCenter(iBiny);
      Double_t bincont=hnbvsnc->GetBinContent(iBinx,iBiny);
      nccbar+=(bincentx*bincont);
      nbbbar+=(bincenty*bincont);
    }
  }


  printf("Events =%f\n",nev);
  printf("c+cbar =%f\n",nccbar);
  printf("D0 =%f\n",nD0);
  printf("D+ =%f\n",nDp);
  printf("D*+ =%f\n",nDst);
  printf("Ds =%f\n",nDs);
  printf("Lc =%f\n",nLc);
  printf("----------\n");
  printf("b+bbar =%f\n",nbbbar);
  printf("B0 =%f\n",nB0);
  printf("B+ =%f\n",nBp);
  printf("B*0 =%f\n",nBst);
  printf("Bs =%f\n",nBs);
  printf("Lb =%f\n",nLb);

  if(nccbar==0) nccbar=nD0+nDp+nDs+nLc;
  if(nbbbar==0) nbbbar=nB0+nBp+nBs+nLb;


  TH1F* hCharmHad=new TH1F("hCharmHad","",5,-0.5,4.5);
  hCharmHad->GetXaxis()->SetBinLabel(1,"D0");
  hCharmHad->SetBinContent(1,nD0/nccbar);
  hCharmHad->GetXaxis()->SetBinLabel(2,"D+");
  hCharmHad->SetBinContent(2,nDp/nccbar);
  hCharmHad->GetXaxis()->SetBinLabel(3,"Ds");
  hCharmHad->SetBinContent(3,nDs/nccbar);
  hCharmHad->GetXaxis()->SetBinLabel(4,"Lc");
  hCharmHad->SetBinContent(4,nLc/nccbar);
  hCharmHad->GetXaxis()->SetBinLabel(5,"D*+");
  hCharmHad->SetBinContent(5,nDst/nccbar);
  hCharmHad->SetMinimum(0);
  hCharmHad->SetStats(0);
  hCharmHad->GetYaxis()->SetTitle("N(species)/N(c+#bar{c})");
  hCharmHad->GetYaxis()->SetTitleOffset(1.3);

  TH1F* hBeautyHad=new TH1F("hBeautyHad","",5,-0.5,4.5);
  hBeautyHad->GetXaxis()->SetBinLabel(1,"B0");
  hBeautyHad->SetBinContent(1,nB0/nbbbar);
  hBeautyHad->GetXaxis()->SetBinLabel(2,"B+");
  hBeautyHad->SetBinContent(2,nBp/nbbbar);
  hBeautyHad->GetXaxis()->SetBinLabel(3,"Bs");
  hBeautyHad->SetBinContent(3,nBs/nbbbar);
  hBeautyHad->GetXaxis()->SetBinLabel(4,"Lb");
  hBeautyHad->SetBinContent(4,nLb/nbbbar);
  hBeautyHad->GetXaxis()->SetBinLabel(5,"B*0");
  hBeautyHad->SetBinContent(5,nBst/nbbbar);
  hBeautyHad->SetMinimum(0);
  hBeautyHad->SetStats(0);
  hBeautyHad->GetYaxis()->SetTitle("N(species)/N(b+#bar{b})");
  hBeautyHad->GetYaxis()->SetTitleOffset(1.3);



  TCanvas* chad=new TCanvas("chad","Hadrons",1200,600);
  chad->Divide(2,1);
  chad->cd(1);
  hCharmHad->Draw();
  chad->cd(2);
  hBeautyHad->Draw();
  


  TCanvas* cd0=new TCanvas("cd0","D0");
  cd0->Divide(2,2);
  cd0->cd(1);
  hyptD0prompt->Draw("colz");
  cd0->cd(2);
  hyptD0feeddown->Draw("colz");
  cd0->cd(3);
  hyptD02->Draw("colz");
  cd0->cd(4);
  hyptD04->Draw("colz");


  TCanvas* cdplus=new TCanvas("cdplus","Dplus");
  cdplus->Divide(2,2);
  cdplus->cd(1);
  hyptDplusprompt->Draw("colz");
  cdplus->cd(2);
  hyptDplusfeeddown->Draw("colz");
  cdplus->cd(3);
  hyptDplusnonreson->Draw("colz");
  cdplus->cd(4);
  hyptDplusreson->Draw("colz");

   
  TCanvas* cds=new TCanvas("cds","Ds");
  cds->Divide(2,2);
  cds->cd(1);
  hyptDsprompt->Draw("colz");
  cds->cd(2);
  hyptDsfeeddown->Draw("colz");
  cds->cd(3);
  hyptDsphi->Draw("colz");
  cds->cd(4);
  hyptDsK0st->Draw("colz");


  TCanvas* cdstlc=new TCanvas("cdstls","Dstar LambdaC");
  cdstlc->Divide(2,2);
  cdstlc->cd(1);
  hyptDstarprompt->Draw("colz");
  cdstlc->cd(2);
  hyptDstarfeedown->Draw("colz");
  cdstlc->cd(3);
  hyptLcprompt->Draw("colz");
  cdstlc->cd(4);
  hyptLcfeedown->Draw("colz");

}



TH1F* HistoFONLL7TeV(){

  TH1F* hFONLL7=new TH1F("hFONLLD07TeV","",61,0.,30.5);
  Float_t val[61]={
    1390542.31,3512269.33,4572233.65,4116353.65,3104057.40,
    2185147.21,1507632.40,1038687.03,721889.43,509311.55,
    365094.01,265684.43,196609.08,147415.64,112019.94,
    86170.85,66997.46,52651.54,41787.55,33458.64,
    27012.62,21981.48,18020.69,14873.73,12354.91,
    10324.90,8677.32,7331.56,6225.67,5311.74,
    4552.30,3918.10,3385.73,2936.79,2556.52,
    2233.25,1957.23,1720.64,1517.12,1341.44,
    1189.32,1057.17,942.05,841.47,753.33,
    675.89,607.68,547.45,494.14,446.87,
    404.83,367.39,333.96,304.08,277.30,
    253.25,231.62,212.13,194.55,178.66,
    164.27};
  for(Int_t i=0; i<61; i++) hFONLL7->SetBinContent(i+1,val[i]);
  return hFONLL7;
}

TH1F* HistoFONLL2_76TeV(){
  TH1F* hFONLL2=new TH1F("hFONLLD02_76TeV","",61,0.,30.5);
  Float_t val[61]={
    1154043.73,2596175.89,2995937.57,2442988.08,1701598.07,
    1122452.81,732896.42,482314.10,322062.75,218493.05,
    151653.77,107053.80,76999.08,56239.64,41687.89,
    31322.25,23818.86,18326.64,14252.10,11192.03,
    8869.29,7089.11,5711.50,4635.50,3788.31,
    3116.19,2578.84,2146.51,1796.18,1510.67,
    1276.72,1083.93,924.20,791.18,679.88,
    586.36,507.49,440.73,383.99,335.54,
    294.07,258.43,227.68,201.11,178.07,
    158.04,140.58,125.32,111.95,100.20,
    89.86,80.73,72.66,65.51,59.16,
    53.51,48.48,43.98,39.96,36.36,
    33.12};
  for(Int_t i=0; i<61; i++) hFONLL2->SetBinContent(i+1,val[i]);
  return hFONLL2;
}

TH1F* HistoPYTHIA7(){
  TH1F* hPYTHIA7=new TH1F("hPYTHIAD07TeV","",40,0.,20.);
  Float_t val[40]={
    826307.00,1753264.00,1877464.00,1634664.00,1278586.00,
    959137.00,713389.00,535745.00,410250.00,316284.00,
    249107.00,198235.00,158832.00,129936.00,106957.00,
    89098.00,73690.00,62752.00,53247.00,45004.00,
    38475.00,32691.00,28510.00,24516.00,21204.00,
    18276.00,15890.00,13702.00,12326.00,10612.00,
    9184.00,8028.00,7194.00,6384.00,5767.00,
    5102.00,4505.00,3939.00,3578.00,3288.00
  };
  for(Int_t i=0; i<40; i++) hPYTHIA7->SetBinContent(i+1,val[i]);
  return hPYTHIA7;
}
