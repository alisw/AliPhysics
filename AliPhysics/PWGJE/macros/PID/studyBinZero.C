#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLegend.h"

#include <iostream>

#include "THnSparseDefinitions.h"

const Int_t iMult   = 0;
const Int_t iPt     = 1;
const Int_t iEta    = 2;

const Bool_t drawPileUp = kFALSE;


//___________________________________________________________________
void SetupLegend(TLegend* leg)
{
  if (!leg)
    return;
  
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetTextSize(0.06);
  
}


//___________________________________________________________________
void SetCanvasMargins(TCanvas* c)
{
  if (!c)
    return;
  
  c->SetRightMargin(0.02);
  c->SetTopMargin(0.04);
  c->SetLeftMargin(0.13);
  c->SetBottomMargin(0.167);
}


//___________________________________________________________________
void normaliseHist(TH1* h, Double_t nEvt)
{
  if (h->GetSumw2N() <= 0)
    h->Sumw2();
  
  if (nEvt <= 0) {
    printf("Error: nEvt <= 0!\n");
    return;
  }
  
  h->Scale(1. / nEvt);
  
  for (Int_t i = 1; i <= h->GetNbinsX(); i++) {
    Double_t normFactor = h->GetBinWidth(i);
    h->SetBinContent(i, h->GetBinContent(i) / normFactor);
    h->SetBinError(i, h->GetBinError(i) / normFactor);
  }
}


//___________________________________________________________________
void setSparseErrors(THnSparse* h)
{
  if (h && !h->GetCalculateErrors()) {
    std::cout << "Re-calculating errors of " << h->GetName() << "..." << std::endl;
    h->Sumw2();
    Long64_t nBinsTHnSparse = h->GetNbins();
    Double_t binContent = 0;
    
    for (Long64_t bin = 0; bin < nBinsTHnSparse; bin++) {
      binContent = h->GetBinContent(bin);
      h->SetBinError(bin, TMath::Sqrt(binContent));
    }
  }
}


//___________________________________________________________________
void setupHist(TH1* h, Int_t colour, TString histTitle = "")
{
  h->SetLineWidth(2);
  h->SetMarkerColor(colour);
  h->SetLineColor(colour);
  
  if (histTitle != "")
    h->SetTitle(histTitle.Data());
  
  h->GetXaxis()->SetRangeUser(0.15, 50);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetXaxis()->SetNoExponent(kTRUE);
  h->GetXaxis()->SetTitleOffset(1.05);
  h->GetXaxis()->SetTitleSize(0.07);
  h->GetXaxis()->SetLabelSize(0.06);
  
  h->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N}/d#it{p}_{T}^{gen} (GeV/#it{c})^{-1}");
  h->GetYaxis()->SetTitleSize(0.07);
  h->GetYaxis()->SetTitleOffset(0.9);
  h->GetYaxis()->SetLabelSize(0.06);

  h->SetStats(kFALSE);
}


//___________________________________________________________________
Int_t studyBinZero(TString pathData, TString fileNameData, TString listName = "",
                   Double_t etaLow = -0.9, Double_t etaUp = 0.9,
                   Double_t lowerCentrality = -2, Double_t upperCentrality = -2)
{
  PrintSettingsAxisRangeForMultiplicityAxisForMB();
  
  TString pathNameData = Form("%s/%s", pathData.Data(), fileNameData.Data());
  
  if (listName == "") {
    listName = pathNameData;
    listName.Replace(0, listName.Last('/') + 1, "");
    listName.ReplaceAll(".root", "");
  }
  
  TObjArray* histList = 0x0;
  
  TFile* f = TFile::Open(pathNameData.Data());
  if (!f)  {
    std::cout << std::endl;
    std::cout << "Failed to open file \"" << pathNameData.Data() << "\"!" << std::endl;
    return -1;
  }
  
  histList = (TObjArray*)(f->Get(listName.Data()));
  if (!histList) {
    std::cout << std::endl;
    std::cout << "Failed to load list \"" << listName.Data() << "\"!" << std::endl;
    return -1;
  }
  
  // Extract the data histograms
  TH1* hNumEventsTriggerSel = dynamic_cast<TH1*>(histList->FindObject("fhEventsTriggerSel"));
  TH1* hNumEventsTriggerSelVtxCut = dynamic_cast<TH1*>(histList->FindObject("fhEventsTriggerSelVtxCut"));
  TH1* hNumEventsTriggerSelVtxCutZ= dynamic_cast<TH1*>(histList->FindObject("fhEventsProcessedNoPileUpRejection"));
  TH1* hNumEventsTriggerSelVtxCutZPileUpRej= dynamic_cast<TH1*>(histList->FindObject("fhEventsProcessed"));
  
  THnSparse* hDataTriggerSel = dynamic_cast<THnSparse*>(histList->FindObject("fChargedGenPrimariesTriggerSel"));
  THnSparse* hDataTriggerSelVtxCut = dynamic_cast<THnSparse*>(histList->FindObject("fChargedGenPrimariesTriggerSelVtxCut"));
  THnSparse* hDataTriggerSelVtxCutZ = dynamic_cast<THnSparse*>(histList->FindObject("fChargedGenPrimariesTriggerSelVtxCutZ"));
  THnSparse* hDataTriggerSelVtxCutZPileUpRej = dynamic_cast<THnSparse*>(histList->FindObject("fChargedGenPrimariesTriggerSelVtxCutZPileUpRej"));
  
  setSparseErrors(hDataTriggerSel);
  setSparseErrors(hDataTriggerSelVtxCut);
  setSparseErrors(hDataTriggerSelVtxCutZ);
  setSparseErrors(hDataTriggerSelVtxCutZPileUpRej);
  
  
  
  // Set multiplicity range, if desired. Note that mult and pT axes are the same for all histos
  Int_t lowerCentralityBinLimit = -1;
  Int_t upperCentralityBinLimit = -2;
  Bool_t restrictCentralityAxis = kFALSE;
  Double_t actualLowerCentrality = -1.;
  Double_t actualUpperCentrality = -1.;
  
  if (lowerCentrality >= -1 && upperCentrality >= -1) {
    // Add subtract a very small number to avoid problems with values right on the border between to bins
    lowerCentralityBinLimit = hNumEventsTriggerSel->GetXaxis()->FindFixBin(lowerCentrality + 0.001);
    upperCentralityBinLimit = hNumEventsTriggerSel->GetXaxis()->FindFixBin(upperCentrality - 0.001);
    
    // Check if the values look reasonable
    if (lowerCentralityBinLimit <= upperCentralityBinLimit && lowerCentralityBinLimit >= 1
        && upperCentralityBinLimit <= hNumEventsTriggerSel->GetXaxis()->GetNbins()) {
      actualLowerCentrality = hNumEventsTriggerSel->GetXaxis()->GetBinLowEdge(lowerCentralityBinLimit);
      actualUpperCentrality = hNumEventsTriggerSel->GetXaxis()->GetBinUpEdge(upperCentralityBinLimit);

      restrictCentralityAxis = kTRUE;
    }
    else {
      std::cout << std::endl;
      std::cout << "Requested centrality range out of limits or upper and lower limit are switched!" << std::endl;
      return -1;
    }
  }
  
  if (!restrictCentralityAxis)
    GetAxisRangeForMultiplicityAxisForMB(hNumEventsTriggerSel->GetXaxis(), lowerCentralityBinLimit, upperCentralityBinLimit);
  
  hNumEventsTriggerSel->GetXaxis()->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
  
  actualLowerCentrality = hNumEventsTriggerSel->GetXaxis()->GetBinLowEdge(hNumEventsTriggerSel->GetXaxis()->GetFirst());
  actualUpperCentrality = hNumEventsTriggerSel->GetXaxis()->GetBinUpEdge(hNumEventsTriggerSel->GetXaxis()->GetLast());
  
  std::cout << "centrality: ";
  if (restrictCentralityAxis) {
    std::cout << actualLowerCentrality << " - " << actualUpperCentrality << std::endl;
    std::cout << "WARNING: Does it really make sense to restrict the centrality? Usually, centrality estimators imply centrality >= 0!"
              << std::endl;
  }
  else
    std::cout << "MB (" << actualLowerCentrality << " - " << actualUpperCentrality << ")" << std::endl;
    
  if (restrictCentralityAxis) {
    hDataTriggerSel->GetAxis(iMult)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
    hDataTriggerSelVtxCut->GetAxis(iMult)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
    hDataTriggerSelVtxCutZ->GetAxis(iMult)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
    hDataTriggerSelVtxCutZPileUpRej->GetAxis(iMult)->SetRange(lowerCentralityBinLimit, upperCentralityBinLimit);
  }
  
  
  // Restrict eta range
  const Int_t lowerEtaBinLimit = hDataTriggerSel->GetAxis(iEta)->FindFixBin(etaLow + 0.0001);
  const Int_t upperEtaBinLimit = hDataTriggerSel->GetAxis(iEta)->FindFixBin(etaUp  - 0.0001);
  
  const Double_t actualLowerEta = hDataTriggerSel->GetAxis(iEta)->GetBinLowEdge(lowerEtaBinLimit);
  const Double_t actualUpperEta = hDataTriggerSel->GetAxis(iEta)->GetBinUpEdge(upperEtaBinLimit);
  
  std::cout << "Eta: ";
  std::cout << actualLowerEta << " - " << actualUpperEta << std::endl;
  
  
  
  Bool_t symmetricInEta = TMath::Abs(TMath::Abs(actualLowerEta) - TMath::Abs(actualUpperEta)) < 1e-6;
  TString etaRange = symmetricInEta ? Form("|#it{#eta}| < %.1f", TMath::Abs(actualLowerEta)) 
                                    : Form("%.1f < #it{#eta} < %.1f", actualLowerEta, actualUpperEta); 
  
  hDataTriggerSel->GetAxis(iEta)->SetRange(lowerEtaBinLimit, upperEtaBinLimit);
  hDataTriggerSelVtxCut->GetAxis(iEta)->SetRange(lowerEtaBinLimit, upperEtaBinLimit);
  hDataTriggerSelVtxCutZ->GetAxis(iEta)->SetRange(lowerEtaBinLimit, upperEtaBinLimit);
  hDataTriggerSelVtxCutZPileUpRej->GetAxis(iEta)->SetRange(lowerEtaBinLimit, upperEtaBinLimit);
  
  // Obtain number of events for each class
  // Note: Under- and overflow automatically taken into account for unrestricted axis by using -1 and -2 for limits
  const Double_t numEvtsTriggerSel = hNumEventsTriggerSel->Integral(lowerCentralityBinLimit, upperCentralityBinLimit);
  const Double_t numEvtsTriggerSelVtxCut = hNumEventsTriggerSelVtxCut->Integral(lowerCentralityBinLimit, upperCentralityBinLimit);
  const Double_t numEvtsTriggerSelVtxCutZ = hNumEventsTriggerSelVtxCutZ->Integral(lowerCentralityBinLimit, upperCentralityBinLimit);
  const Double_t numEvtsTriggerSelVtxCutZPileUpRej = hNumEventsTriggerSelVtxCutZPileUpRej->Integral(lowerCentralityBinLimit,
                                                                                                    upperCentralityBinLimit);
  
  // Project pT spectra
  TH1D* hSpectraTriggerSel = (TH1D*)hDataTriggerSel->Projection(iPt, "e");
  TH1D* hSpectraTriggerSelVtxCut = (TH1D*)hDataTriggerSelVtxCut->Projection(iPt, "e");
  TH1D* hSpectraTriggerSelVtxCutZ = (TH1D*)hDataTriggerSelVtxCutZ->Projection(iPt, "e");
  TH1D* hSpectraTriggerSelVtxCutZPileUpRej = (TH1D*)hDataTriggerSelVtxCutZPileUpRej->Projection(iPt, "e");
  
  // Normalise histos to 1/Nevt dN/dPt
  normaliseHist(hSpectraTriggerSel, numEvtsTriggerSel);
  normaliseHist(hSpectraTriggerSelVtxCut, numEvtsTriggerSelVtxCut);
  normaliseHist(hSpectraTriggerSelVtxCutZ, numEvtsTriggerSelVtxCutZ);
  normaliseHist(hSpectraTriggerSelVtxCutZPileUpRej, numEvtsTriggerSelVtxCutZPileUpRej);
  
  setupHist(hSpectraTriggerSel, kBlack, "Trigger selection");
  setupHist(hSpectraTriggerSelVtxCut, kRed, "& vertex cut");
  setupHist(hSpectraTriggerSelVtxCutZ, kBlue, "& vertex #it{z} cut");
  setupHist(hSpectraTriggerSelVtxCutZPileUpRej, kGreen, "& pile-up rejection");
  
  TCanvas* canvSpectra = new TCanvas("canvSpectra", "Spectra", 760, 420);
  canvSpectra->SetLogx();
  canvSpectra->SetLogy();
  canvSpectra->SetGrid(0, 0);
  SetCanvasMargins(canvSpectra);
  
  TLegend* leg = new TLegend(0.14, 0.26, 0.62, 0.55);
  leg->SetHeader(Form("MC pp #sqrt{s}=7 TeV, inclusive, %s", etaRange.Data()));
  leg->AddEntry(hSpectraTriggerSel, "", "l");
  leg->AddEntry(hSpectraTriggerSelVtxCut, "", "l");
  leg->AddEntry(hSpectraTriggerSelVtxCutZ, "", "l");
  if (drawPileUp)
    leg->AddEntry(hSpectraTriggerSelVtxCutZPileUpRej, "", "l");
  SetupLegend(leg);
  leg->SetMargin(0.15);
  
  hSpectraTriggerSel->Draw();
  hSpectraTriggerSelVtxCut->Draw("same");
  hSpectraTriggerSelVtxCutZ->Draw("same");
  if (drawPileUp)
    hSpectraTriggerSelVtxCutZPileUpRej->Draw("same");
  
  leg->Draw("same");
  
  // Ratios (take binomial errors, since real sub-samples)
  
  TH1D* hRatioVtxCut = new TH1D(*hSpectraTriggerSelVtxCut);
  hRatioVtxCut->SetName("hRatioVtxCut");
  setupHist(hRatioVtxCut, kBlack, "MB & vtx / MB");
  hRatioVtxCut->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N}/d#it{p}_{T}^{gen} ratio");
  hRatioVtxCut->Divide(hRatioVtxCut, hSpectraTriggerSel, 1., 1., "B");
  hRatioVtxCut->GetYaxis()->SetRangeUser(1.08, 1.14);
  
  // Mean ratio of spectra and of Nevt
  TF1* funcRatioVtxCut = new TF1("funcRatioVtxCut", "pol0",
                                 0.15, hRatioVtxCut->GetXaxis()->GetBinUpEdge(hSpectraTriggerSelVtxCut->FindLastBinAbove(0)));
  funcRatioVtxCut->SetLineWidth(2);
  funcRatioVtxCut->SetLineColor(kRed);
  hRatioVtxCut->Fit(funcRatioVtxCut, "N");
  const Double_t meanRatioVtxCut = funcRatioVtxCut->GetParameter(0);
  const Double_t ratioNevtVtxCut = numEvtsTriggerSelVtxCut > 0 ? numEvtsTriggerSel / numEvtsTriggerSelVtxCut : -1.;
  const Double_t doubleRatioMinusOne = meanRatioVtxCut > 0 ? (ratioNevtVtxCut / meanRatioVtxCut - 1.) : -999.;
  
  
  
  TH1D* hRatioVtxCutZ = new TH1D(*hSpectraTriggerSelVtxCutZ);
  hRatioVtxCutZ->SetName("hRatioVtxCutZ");
  setupHist(hRatioVtxCutZ, kBlack, "MB & vtx & #it{z} vtx / MB & vtx");
  hRatioVtxCutZ->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N}/d#it{p}_{T}^{gen} ratio");
  hRatioVtxCutZ->Divide(hRatioVtxCutZ, hSpectraTriggerSelVtxCut, 1., 1., "B");
  hRatioVtxCutZ->GetYaxis()->SetRangeUser(0.95, 1.05);
  
  // Mean ratio of spectra and of Nevt
  TF1* funcRatioVtxCutZ = new TF1("funcRatioVtxCutZ", "pol0",
                                  0.15, hRatioVtxCutZ->GetXaxis()->GetBinUpEdge(hSpectraTriggerSelVtxCutZ->FindLastBinAbove(0)));
  hRatioVtxCutZ->Fit(funcRatioVtxCutZ, "N");
  funcRatioVtxCutZ->SetLineWidth(2);
  funcRatioVtxCutZ->SetLineColor(kRed);
  const Double_t meanRatioVtxCutZ = funcRatioVtxCutZ->GetParameter(0);
  const Double_t ratioNevtVtxCutZ = numEvtsTriggerSelVtxCutZ > 0 ? numEvtsTriggerSelVtxCut / numEvtsTriggerSelVtxCutZ : -1.;
  
  TH1D* hRatioPileUp = new TH1D(*hSpectraTriggerSelVtxCutZPileUpRej);
  hRatioPileUp->SetName("hRatioPileUp");
  setupHist(hRatioPileUp, kBlack, "MB & vtx & #it{z} vtx & pile-up / MB & vtx & #it{z} vtx");
  hRatioPileUp->GetYaxis()->SetTitle("1/#it{N}_{evt} d#it{N}/d#it{p}_{T}^{gen} ratio");
  hRatioPileUp->Divide(hRatioPileUp, hSpectraTriggerSelVtxCutZ, 1., 1., "B");
  
  /*
  TF1* funcRatioPileUp = new TF1("funcRatioPileUp", "pol0",
                                 0.15, hRatioPileUp->GetXaxis()->GetBinUpEdge(hSpectraTriggerSelVtxCutZPileUpRej->FindLastBinAbove(0)));
  hRatioPileUp->Fit(funcRatioPileUp, "N");*/
  
  ClearTitleFromHistoInCanvas(canvSpectra);
  
  TCanvas* canvRatioVtx = new TCanvas("canvRatioVtx", "Ratio vertex", 760, 420);
  canvRatioVtx->SetLogx();
  canvRatioVtx->SetGrid(0, 1);
  SetCanvasMargins(canvRatioVtx);
  
  hRatioVtxCut->Draw();
  funcRatioVtxCut->Draw("same");
  
  TLegend* leg2 = new TLegend(0.22, 0.7, drawPileUp ? 0.87 : 0.72, 0.9);
  leg2->SetHeader(Form("MC pp #sqrt{s}=7 TeV, inclusive, %s", etaRange.Data()));
  leg2->AddEntry(hRatioVtxCut, "", "l");
  SetupLegend(leg2);
  leg2->SetMargin(0.1);
  
  leg2->Draw("same");
  
  ClearTitleFromHistoInCanvas(canvRatioVtx);
  
  TCanvas* canvRatioOther = new TCanvas("canvRatioOther", "Ratio others", 760, 420);
  canvRatioOther->SetLogx();
  canvRatioOther->SetGrid(0, 1);
  SetCanvasMargins(canvRatioOther);
  
  hRatioVtxCutZ->Draw();
  if (drawPileUp)
    hRatioPileUp->Draw("same");

  TLegend* leg3 = new TLegend(0.22, drawPileUp ? 0.63 : 0.7, drawPileUp ? 0.87 : 0.72, 0.9);
  leg3->SetHeader(Form("MC pp #sqrt{s}=7 TeV, inclusive, %s", etaRange.Data()));
  leg3->AddEntry(hRatioVtxCutZ, "", "l");
  if (drawPileUp)
    leg3->AddEntry(hRatioPileUp, "", "l");
  SetupLegend(leg3);
  leg3->SetMargin(0.1);
  
  leg3->Draw("same");
  
  ClearTitleFromHistoInCanvas(canvRatioOther);

  
  printf("meanRatioVtxCut %f <-> ratioNevtVtxCut %f => meanRatioVtxCutZ/ratioNevtVtxCutZ - 1 = %f\nmeanRatioVtxCutZ %f <-> ratioNevtVtxCutZ %f\n",
         meanRatioVtxCut, ratioNevtVtxCut, doubleRatioMinusOne, meanRatioVtxCutZ, ratioNevtVtxCutZ);
  
  
  TNamed* settings = new TNamed(
      Form("Settings: Data file \"%s\", lowerEta %.2f, uppEta %.2f, lowerCentrality %.3f, upperCentrality %.3f, doubleRatioMinusOne %.4f\n",
           pathNameData.Data(), etaLow, etaUp, lowerCentrality, upperCentrality, doubleRatioMinusOne), "");
  
  // Save results to file
  TString saveFileName = pathNameData;
  saveFileName = Form("%s_binZeroStudy.root", saveFileName.ReplaceAll(".root", "").Data());
  TFile *saveFile = TFile::Open(saveFileName.Data(), "RECREATE");
  saveFile->cd();
  hSpectraTriggerSel->Write();
  hSpectraTriggerSelVtxCut->Write();
  hSpectraTriggerSelVtxCutZ->Write();
  hSpectraTriggerSelVtxCutZPileUpRej->Write();
  
  hRatioVtxCut->Write();
  hRatioVtxCutZ->Write();
  hRatioPileUp->Write();
  
  funcRatioVtxCut->Write();
  funcRatioVtxCutZ->Write();
  
  canvSpectra->Write();
  canvRatioVtx->Write();
  canvRatioOther->Write();
  
  settings->Write();
  
  saveFile->Close();
  
  TString temp = saveFileName;
  canvRatioVtx->SaveAs(Form("%s", temp.ReplaceAll(".root", "_ratioVtx.pdf").Data()));
  temp = saveFileName;
  canvRatioOther->SaveAs(Form("%s", temp.ReplaceAll(".root", "_ratioOther.pdf").Data()));
  
  PrintSettingsAxisRangeForMultiplicityAxisForMB();
  
  return 0;
}