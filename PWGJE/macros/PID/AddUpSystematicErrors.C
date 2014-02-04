#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLegend.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TMath.h"

#include "AliPID.h"

#include "THnSparseDefinitions.h"

#include <iostream>
#include <iomanip>

const Int_t numSpecies = 5;

//________________________________________________________
TCanvas* DrawFractionHistos(TString canvName, TString canvTitle, Double_t pLow, Double_t pHigh, TGraphAsymmErrors** gr,  TH1F** histRef)
{
  TCanvas* canv = new TCanvas(canvName.Data(), canvTitle.Data(),100,10,1200,800);
  canv->SetGridx(1);
  canv->SetGridy(1);
  canv->SetLogx(1);
  
  for (Int_t i = 0; i < numSpecies; i++) {
    histRef[i]->GetYaxis()->SetRangeUser(0.0, 1.0);
    histRef[i]->GetYaxis()->SetTitle(canvTitle.Data());
    histRef[i]->GetXaxis()->SetRangeUser(pLow, pHigh);
    //histRef[i]->SetFillStyle(3004 + i);
    //histRef[i]->SetFillColor(kGray);
    histRef[i]->SetFillStyle(0);
    histRef[i]->SetFillColor(histRef[i]->GetMarkerColor());
    histRef[i]->SetLineColor(histRef[i]->GetMarkerColor());
  }
  histRef[2]->SetMarkerStyle(20);
  histRef[2]->Draw("e p");
  histRef[0]->SetMarkerStyle(21);
  histRef[0]->Draw("e p same");
  histRef[1]->SetMarkerStyle(22);
  histRef[1]->Draw("e p same");
  histRef[3]->SetMarkerStyle(29);
  histRef[3]->Draw("e p same");
  histRef[4]->SetMarkerStyle(30);
  histRef[4]->Draw("e p same");
  
  gr[0]->GetHistogram()->GetXaxis()->SetRangeUser(pLow, pHigh);
  gr[0]->GetHistogram()->GetYaxis()->SetRangeUser(0., 1.0);
  gr[0]->Draw("2 same");
  gr[1]->Draw("2 same");
  gr[2]->Draw("2 same");
  gr[3]->Draw("2 same");
  gr[4]->Draw("2 same");
  
  TLegend* legend = new TLegend(0.622126, 0.605932, 0.862069, 0.855932);    
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->AddEntry(histRef[2], "#pi", "flp");
  legend->AddEntry(histRef[0], "e", "flp");
  legend->AddEntry(histRef[1], "K", "flp");
  legend->AddEntry(histRef[3], "p", "flp");
  legend->AddEntry(histRef[4], "#mu", "flp");
  legend->Draw();
  
  ClearTitleFromHistoInCanvas(canv);
  
  return canv;
}


//________________________________________________________
TGraphAsymmErrors* loadGraph(const TString graphName, TFile* f)
{
  if (!f) {
    std::cout << "No file. Cannot load graph \"" << graphName.Data() << "\n!" << std::endl;
    return 0x0;
  }
  
  TGraphAsymmErrors* grTemp = dynamic_cast<TGraphAsymmErrors*>(f->Get(graphName.Data()));
  if (!grTemp) {
    std::cout << "Failed to load histo \"" << graphName.Data() << "\"!" << std::endl;
    return 0x0;
  } 
  
  return grTemp;
}


//________________________________________________________
TH1F* loadHisto(const TString histName, TFile* f)
{
  if (!f) {
    std::cout << "No file. Cannot load hist \"" << histName.Data() << "\n!" << std::endl;
    return 0x0;
  }
  
  TH1F* hTemp = dynamic_cast<TH1F*>(f->Get(histName.Data()));
  if (!hTemp) {
    std::cout << "Failed to load histo \"" << histName.Data() << "\"!" << std::endl;
    return 0x0;
  } 
  
  return hTemp;
}


//________________________________________________________
Int_t AddUpSystematicErrors(const TString path, const TString outFileTitle, const TString* fileNames, const Int_t numFiles,
                            const TString fileNameReference) 
{ 
  if (!fileNames || numFiles < 1)
    return -1;
  
  TFile* f[numFiles];
  TGraphAsymmErrors** grSysError[numSpecies];
  TString grNames[numSpecies] = {"systematicError_electron", "systematicError_kaon", "systematicError_pion", "systematicError_proton",
                                 "systematicError_muon" };
  
  const TString histNames[numSpecies] = {"hFractionElectrons", "hFractionKaons", "hFractionPions", "hFractionProtons", "hFractionMuons" };
  
  
  TGraphAsymmErrors** grSysErrorYields[numSpecies];
  TString grNamesYields[numSpecies] = {"systematicErrorYields_electron", "systematicErrorYields_kaon", "systematicErrorYields_pion",
                                       "systematicErrorYields_proton", "systematicErrorYields_muon" };
                                 
  const TString histNamesYields[numSpecies] = {"hYieldElectrons", "hYieldKaons", "hYieldPions", "hYieldProtons", "hYieldMuons" };
  
  for (Int_t i = 0; i < numSpecies; i++) {
    grSysError[i] = new TGraphAsymmErrors*[numFiles];
    grSysErrorYields[i] = new TGraphAsymmErrors*[numFiles];
  }
    
  for (Int_t iFile = 0; iFile < numFiles; iFile++) {
    f[iFile] = TFile::Open(fileNames[iFile].Data());
    if (!f[iFile])  {
      std::cout << "Failed to open file \"" << fileNames[iFile].Data() << "\"!" << std::endl;
      return -1;
    }
        
    // Extract the data graphs
    for (Int_t i = 0; i < numSpecies; i++) {
      grSysError[i][iFile] = loadGraph(grNames[i], f[iFile]);
      if (!grSysError[i][iFile])
        return -1;
      
      grSysError[i][iFile]->SetName(Form("%s_fileID%d", grSysError[i][iFile]->GetName(), iFile));
      
      
      grSysErrorYields[i][iFile] = loadGraph(grNamesYields[i], f[iFile]);
      if (!grSysErrorYields[i][iFile])
        return -1;
      
      grSysErrorYields[i][iFile]->SetName(Form("%s_fileID%d", grSysErrorYields[i][iFile]->GetName(), iFile));
    }
  }
  
  // Load data points with statistical errors
  TFile* fReferenceData = TFile::Open(fileNameReference.Data());
  if (!fReferenceData) {
    std::cout << "Failed to open file \"" << fileNameReference.Data() << "\"!" << std::endl;
    return -1;
  }
    
  TH1F* hReferenceFractions[numSpecies];
  TH1F* hReferenceYields[numSpecies];
  
  for (Int_t i = 0; i < numSpecies; i++) {
    hReferenceFractions[i] = loadHisto(histNames[i], fReferenceData);
    if (!hReferenceFractions[i])
      return -1;
    
    hReferenceYields[i] = loadHisto(histNamesYields[i], fReferenceData);
    if (!hReferenceYields[i])
      return -1;
  }
  
  const Int_t reference = 0;
  
  // The x,y coordinates should be those of the reference graph. Then, all corresponding sys. errors are added in quadrature.
  TGraphAsymmErrors* grTotSysError[numSpecies];
  TGraphAsymmErrors* grTotSysErrorYields[numSpecies];
  Double_t sysErrorTotalSquared = 0;
  Double_t temp = 0;
  
  for (Int_t i = 0; i < numSpecies; i++) {
    grTotSysError[i] = new TGraphAsymmErrors(*grSysError[i][reference]);
    grTotSysError[i]->SetName(grNames[i]);
    
    for (Int_t iPoint = 0; iPoint < grTotSysError[i]->GetN(); iPoint++) {
      sysErrorTotalSquared = 0;
      for (Int_t iFile = 0; iFile < numFiles; iFile++) {
        // Already averages high and low value -> Since they are by now the same, this is ok.
        temp = grSysError[i][iFile]->GetErrorY(iPoint); 
        if (temp > 0) {
          sysErrorTotalSquared += temp * temp;
        }
      }
      
      grTotSysError[i]->SetPointEYlow(iPoint, TMath::Sqrt(sysErrorTotalSquared));
      grTotSysError[i]->SetPointEYhigh(iPoint, TMath::Sqrt(sysErrorTotalSquared));
    }
    
    // Same for the yields
    grTotSysErrorYields[i] = new TGraphAsymmErrors(*grSysErrorYields[i][reference]);
    grTotSysErrorYields[i]->SetName(grNamesYields[i]);
    
    for (Int_t iPoint = 0; iPoint < grTotSysErrorYields[i]->GetN(); iPoint++) {
      sysErrorTotalSquared = 0;
      for (Int_t iFile = 0; iFile < numFiles; iFile++) {
        // Already averages high and low value -> Since they are by now the same, this is ok.
        temp = grSysErrorYields[i][iFile]->GetErrorY(iPoint); 
        if (temp > 0) {
          sysErrorTotalSquared += temp * temp;
        }
      }
      
      grTotSysErrorYields[i]->SetPointEYlow(iPoint, TMath::Sqrt(sysErrorTotalSquared));
      grTotSysErrorYields[i]->SetPointEYhigh(iPoint, TMath::Sqrt(sysErrorTotalSquared));
    }
  }
  
  const Double_t pLow = 0.15;
  const Double_t pHigh = 50.;
  TCanvas* cFractionsWithTotalSystematicError = DrawFractionHistos("cFractionsWithTotalSystematicError", "Particle fractions", pLow, pHigh,
                                                                   grTotSysError, hReferenceFractions);
  
    
  // Output file
  TFile* fSave = 0x0;
  TDatime daTime;
  TString saveFileName;
    
  saveFileName = Form("outputSystematicsTotal_%s__%04d_%02d_%02d.root", outFileTitle.Data(), daTime.GetYear(),
                      daTime.GetMonth(), daTime.GetDay());
    
  fSave = TFile::Open(Form("%s/%s", path.Data(), saveFileName.Data()), "recreate");
  if (!fSave) {
    std::cout << "Failed to open save file \"" << Form("%s/%s", path.Data(), saveFileName.Data()) << "\"!" << std::endl;
    return -1;
  }
  
  // Save final results
  fSave->cd();
  
  if (cFractionsWithTotalSystematicError)
      cFractionsWithTotalSystematicError->Write();
    
    for (Int_t i = 0; i < numSpecies; i++) {
      if (grTotSysError[i])
        grTotSysError[i]->Write();
      if (hReferenceFractions[i])
        hReferenceFractions[i]->Write();
      
      if (grTotSysErrorYields[i])
        grTotSysErrorYields[i]->Write();
      if (hReferenceYields[i])
        hReferenceYields[i]->Write();
  }
  
  // Save list of file names in output file
  TString listOfFileNames = "";
  for (Int_t i = 0; i < numFiles; i++) {
    listOfFileNames.Append(Form("%s%d: %s", i == 0 ? "" : ", ", i, fileNames[i].Data()));
  }
  
  TNamed* settings = new TNamed(Form("Used files for systematics: %s\n", listOfFileNames.Data()),
                                Form("Used files for systematics: %s\n", listOfFileNames.Data()));
  settings->Write();
  
  delete cFractionsWithTotalSystematicError;
    
  fSave->Close();
     
  return 0;
}
