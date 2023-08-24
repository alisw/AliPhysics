/** @file rDrawEmbeddedAndInputCellClusterTrack.C
* @brief Example to check output of embedding analysis
*
* @ingroup CaloTrackCorrMacros
*
* Get the cells, clusters and tracks energy distribution for the data, embedded MC, 
* and combined output of embedding. Compare them also to the simple sum of data and MC.
* For tracks it the sum and combined should be the same, 
* for clusters and cells it should be close but not exactly.
*
* The input are the histograms of the example file runCaloTrackCorrEmbeddingAnalysis.C
* 
* @author Gustavo Conesa Balbastre <gustavo.conesa.balbastre@cern.ch>, LPSC-Grenoble
* @date Mar 11, 2020
*/


void DrawEmbeddedAndInputCellClusterTrack()
{
  gStyle->SetOptTitle(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(000000);
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.1);
  
  TFile * file = TFile::Open("AnalysisResults.root");
  
  TList * lData = (TList*) file->Get("CTC_EMCAL_Trig_default_SPDPileUp");
  TList * lComb = (TList*) file->Get("CTC_EMCAL_Trig_default_SPDPileUp_EmbedMC");
  TList * lSimu = (TList*) file->Get("CTC_EMCAL_Trig_default_SPDPileUp_EmbedMCInput");
  
  //printf("list pointer: %p, %p %p\n",lData,lComb,lSimu);
  
  const Int_t nCases = 4;
  
  TH1F* hCe[nCases];
  TH1F* hCl[nCases];
  TH1F* hTr[nCases];
  
  hCe[0] = (TH1F*) lData->FindObject("QA_hAmplitude"); 
  hCl[0] = (TH1F*) lData->FindObject("QA_hE"); 
  hTr[0] = (TH1F*) lData->FindObject("AnaHadrons_hPt");
  
  //printf("histo0 pointer: %p, %p %p\n",hCe[0],hCl[0],hTr[0]);

  hCe[1] = (TH1F*) lSimu->FindObject("QA_hAmplitude"); 
  hCl[1] = (TH1F*) lSimu->FindObject("QA_hE"); 
  hTr[1] = (TH1F*) lSimu->FindObject("AnaHadrons_hPt");
  
  //printf("histo0 pointer: %p, %p %p\n",hCe[1],hCl[1],hTr[1]);
  
  hCe[2] = (TH1F*) lComb->FindObject("QA_hAmplitude"); 
  hCl[2] = (TH1F*) lComb->FindObject("QA_hE");
  hTr[2] = (TH1F*) lComb->FindObject("AnaHadrons_hPt");
  
  //printf("histo0 pointer: %p, %p %p\n",hCe[2],hCl[2],hTr[2]);
  
  hCe[3] = (TH1F*) hCe[0]->Clone("hCeSum");
  hCl[3] = (TH1F*) hCl[0]->Clone("hClSum");
  hTr[3] = (TH1F*) hTr[0]->Clone("hTrSum");
  
  hCe[3]->Add(hCe[1]);
  hCl[3]->Add(hCl[1]);
  hTr[3]->Add(hTr[1]);
  
  // Style  
  hCe[0]->SetLineColor(1);
  hCl[0]->SetLineColor(1);
  hTr[0]->SetLineColor(1);
  
  hCe[1]->SetLineColor(4);
  hCl[1]->SetLineColor(4);
  hTr[1]->SetLineColor(4);
  
  hCe[2]->SetLineColor(2);
  hCl[2]->SetLineColor(2);
  hTr[2]->SetLineColor(2);
  
  hCe[3]->SetLineColor(8);
  hCl[3]->SetLineColor(8);
  hTr[3]->SetLineColor(8);
 
  hCe[3]->SetLineWidth(3);
  hCl[3]->SetLineWidth(3);
  hTr[3]->SetLineWidth(3);
 
//  hCe[0]->SetMarkerColor(1);
//  hCl[0]->SetMarkerColor(1);
//  hTr[0]->SetMarkerColor(1);
//  
//  hCe[1]->SetMarkerColor(4);
//  hCl[1]->SetMarkerColor(4);
//  hTr[1]->SetMarkerColor(4);
//  
//  hCe[2]->SetMarkerColor(2);
//  hCl[2]->SetMarkerColor(2);
//  hTr[2]->SetMarkerColor(2);
//  
//  hCe[3]->SetMarkerColor(8);
//  hCl[3]->SetMarkerColor(8);
//  hTr[3]->SetMarkerColor(8);
  
//  hCe[0]->SetMarkerStyle(20);
//  hCl[0]->SetMarkerStyle(20);
//  hTr[0]->SetMarkerStyle(20);
//  
//  hCe[1]->SetMarkerStyle(20);
//  hCl[1]->SetMarkerStyle(20);
//  hTr[1]->SetMarkerStyle(20);
//  
//  hCe[2]->SetMarkerStyle(24);
//  hCl[2]->SetMarkerStyle(24);
//  hTr[2]->SetMarkerStyle(24);
//  
//  hCe[3]->SetMarkerStyle(24);
//  hCl[3]->SetMarkerStyle(24);
//  hTr[3]->SetMarkerStyle(24);
//  
//  hCe[0]->Sumw2();
//  hCl[0]->Sumw2();
//  hTr[0]->Sumw2();
//  
//  hCe[1]->Sumw2();
//  hCl[1]->Sumw2();
//  hTr[1]->Sumw2();
//  
//  hCe[2]->Sumw2();
//  hCl[2]->Sumw2();
//  hTr[2]->Sumw2();
//  
//  hCe[3]->Sumw2();
//  hCl[3]->Sumw2();
//  hTr[3]->Sumw2();
  
  TLegend l(0.5,0.65,0.9,0.85);
  l.SetTextSize(0.04);
  l.AddEntry(hCe[0],"Data","L");
  l.AddEntry(hCe[1],"MC","L");
  l.AddEntry(hCe[2],"Embedded Data+MC","L");
  l.AddEntry(hCe[3],"Direct Data+MC","L");
  
  // Plot
  TCanvas * c = new TCanvas("c","",3*1000,1*600);
  c->Divide(3,1);
  
  // cell
  c->cd(1);
  gPad->SetLogy();
  gPad->SetLogx();
  
  hCe[3]->SetTitle("cell");
  hCe[3]->SetAxisRange(1,100,"X");
  hCe[3]->Draw();
  
  hCe[1]->Draw("same");
  hCe[0]->Draw("same");
  hCe[2]->Draw("same");
  l.Draw();
  // cluster
  c->cd(2);
  gPad->SetLogy();
  gPad->SetLogx();

  hCl[3]->SetTitle("cluster");
  hCl[3]->SetAxisRange(1,100,"X");
  hCl[3]->Draw();
  
  hCl[1]->Draw("same");
  hCl[0]->Draw("same");
  hCl[2]->Draw("same");

  // tracl
  c->cd(3);
  gPad->SetLogy();
  gPad->SetLogx();

  hTr[3]->SetTitle("track");
  hTr[3]->SetAxisRange(1,100,"X");
  hTr[3]->Draw();
  
  hTr[1]->Draw("same");
  hTr[0]->Draw("same");
  hTr[2]->Draw("same");

  c->Print("CellClusterTrackSpectra.eps");
}
